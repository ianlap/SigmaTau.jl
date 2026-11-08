#!/usr/bin/env python3
"""Regenerate the allantools / scipy reference tables for the cross-check suite.

    pip install allantools==2024.6 scipy numpy
    python3 tools/validation/gen_allantools.py            # writes test/fixtures/validation/allantools/
    python3 tools/validation/gen_allantools.py --out DIR  # write somewhere else (CI)

Inputs are the committed phase records, so SigmaTau and allantools read the
same numbers:

  noise/alpha_*.txt   power-law records, alpha = +2 ... -4, N = 4096
  PHASE.DAT           Stable32's PHASE.DAT (the NBS 1000-point set as phase)
  stable32gen.DAT     8192-point Stable32-generated record

Outputs (CSV, values at %.17e so Float64 round-trips):

  devs.csv     dataset, estimator, m, tau, n, dev       every estimator allantools has
  edf.csv      family, alpha, d, m, N, edf              Greenhall / SP1065 EDF on a grid
  ci.csv       edf, confidence, lo_factor, hi_factor    chi-squared interval factors
  noiseid.csv  dataset, m, dmax, alpha_int, alpha, d, rho   lag-1 ACF noise ID
  welch.csv    dataset, kind, nperseg, f, psd           scipy.signal.welch S_y and S_x
  dynamic.csv  dataset, kind, window, start, m, dev     per-window OADEV / OHDEV
  gaps.csv     dataset, m, n, dev                       gap-tolerant ADEV (gradev)
"""
import argparse
import csv
import io
import math
import os
import sys
import contextlib

import numpy as np
import scipy.signal
import scipy.special
import allantools as at

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.normpath(os.path.join(HERE, "..", "..", "test", "fixtures", "validation"))

ALPHAS = [2, 1, 0, -1, -2, -3, -4]


def alpha_file(a):
    return "alpha_p%d.txt" % a if a > 0 else ("alpha_0.txt" if a == 0 else "alpha_m%d.txt" % -a)


def read_phase(path, skip_header_end=False):
    vals = []
    started = not skip_header_end
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if not started:
                if s.startswith("# Header End"):
                    started = True
                continue
            if s.startswith("#"):
                continue
            vals.append(float(s.split()[-1]))
    return np.array(vals)


def datasets():
    out = []
    for a in ALPHAS:
        out.append(("alpha_%+d" % a, read_phase(os.path.join(FIX, "noise", alpha_file(a)))))
    out.append(("phasedat", read_phase(os.path.join(FIX, "stable32_phasedat", "PHASE.DAT"))))
    out.append(("stable32gen", read_phase(os.path.join(FIX, "stable32gen.DAT"), skip_header_end=True)))
    return out


def octave(n_max):
    m, out = 1, []
    while m <= n_max:
        out.append(m)
        m *= 2
    return out


@contextlib.contextmanager
def quiet():
    # allantools prints from some EDF helpers; keep the CSVs clean.
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def fmt(v):
    return "nan" if (v is None or (isinstance(v, float) and math.isnan(v))) else "%.17e" % v


DEV_FNS = [
    ("adev", at.adev), ("oadev", at.oadev), ("mdev", at.mdev), ("tdev", at.tdev),
    ("hdev", at.hdev), ("ohdev", at.ohdev), ("totdev", at.totdev),
    ("mtotdev", at.mtotdev), ("htotdev", at.htotdev),
    ("tierms", at.tierms), ("mtie", at.mtie), ("pdev", at.pdev),
]
TOTAL_FAMILY = ("mtotdev", "htotdev")


def _devs_for(item):
    """Rows for one dataset (runs in a worker process)."""
    name, x = item
    N = len(x)
    out = []
    ms = sorted(set(octave(N // 2) + [3, 10, 100]))
    for est, fn in DEV_FNS:
        mlist = ms
        if est in TOTAL_FAMILY:
            # allantools' MTOT/HTOT are O(N m) Python loops: cap at m <= 256, and
            # skip the 8192-point record (its total family is already covered by
            # allantools_out/allantools_data_full.csv).
            if name == "stable32gen":
                continue
            mlist = [m for m in ms if m <= 256 and 3 * m <= N]
        if est == "mtie":
            mlist = [m for m in ms if m <= 256]
        taus, devs, errs, ns = fn(x, rate=1.0, data_type="phase", taus=[float(m) for m in mlist])
        for t, d, n in zip(taus, devs, ns):
            out.append([name, est, int(round(t)), fmt(float(t)), int(n), fmt(float(d))])
            if est == "mtotdev":
                # allantools.ttotdev is exactly mtotdev * tau / sqrt(3); derive it
                # instead of paying for a second O(N m) pass.
                out.append([name, "ttotdev", int(round(t)), fmt(float(t)), int(n),
                            fmt(float(d) * float(t) / math.sqrt(3.0))])
    # Theo1: even m only; allantools reports tau = m*tau0 (SP1065 uses 0.75 m).
    tm = [m for m in [2, 4, 8, 10, 16, 32, 64, 100, 128, 256] if m < N]
    taus, devs, errs, ns = at.theo1(x, rate=1.0, data_type="phase", taus=[float(m) for m in tm])
    for t, d, n in zip(taus, devs, ns):
        out.append([name, "theo1", int(round(t)), fmt(float(t)), int(n), fmt(float(d))])
    print("devs:", name, file=sys.stderr, flush=True)
    return out


def gen_devs(w):
    w.writerow(["dataset", "estimator", "m", "tau", "n", "dev"])
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor() as ex:
        for rows in ex.map(_devs_for, datasets()):
            for r in rows:
                w.writerow(r)


def gen_edf(w):
    w.writerow(["family", "alpha", "d", "m", "N", "edf"])
    grid_m = [1, 2, 3, 4, 8, 10, 16, 32, 64, 100, 128, 256, 512]
    for N in [1001, 4096, 8192]:
        for m in grid_m:
            for a in ALPHAS:
                for d in (2, 3):
                    if a + 2 * d <= 1:
                        continue
                    for fam, ov, mod in (("overlapping", True, False),
                                         ("modified", True, True),
                                         ("nonoverlapping", False, False)):
                        L = m / (1 if mod else m) + m * d
                        if N < L + m:
                            continue
                        try:
                            with quiet():
                                e = at.edf_greenhall(alpha=a, d=d, m=m, N=N,
                                                     overlapping=ov, modified=mod)
                        except (NotImplementedError, AssertionError):
                            e = float("nan")
                        w.writerow([fam, a, d, m, N, fmt(float(e))])
                if a in (0, -1, -2):
                    w.writerow(["totdev", a, 2, m, N, fmt(float(at.edf_totdev(N, m, a)))])
                if a in (2, 1, 0, -1, -2):
                    with quiet():
                        w.writerow(["mtotdev", a, 2, m, N, fmt(float(at.edf_mtotdev(N, m, a)))])
                    w.writerow(["simple", a, 2, m, N, fmt(float(at.edf_simple(N, m, a)))])


def gen_ci(w):
    w.writerow(["edf", "confidence", "lo_factor", "hi_factor"])
    for e in [1.0, 1.5, 2.0, 3.7, 5.0, 10.0, 33.3, 100.0, 1000.0, 12345.0]:
        for c in [0.683, at.ONE_SIGMA_CI if hasattr(at, "ONE_SIGMA_CI") else 0.6826894921370859, 0.90, 0.95, 0.99]:
            lo, hi = at.confidence_interval(dev=1.0, edf=e, ci=c)
            w.writerow([fmt(e), fmt(c), fmt(float(lo)), fmt(float(hi))])


def gen_noiseid(w):
    w.writerow(["dataset", "m", "dmax", "alpha_int", "alpha", "d", "rho"])
    for name, x in datasets():
        for m in octave(len(x) // 30):
            for dmax in (2, 3):
                try:
                    with quiet():
                        ai, af, d, rho = at.autocorr_noise_id(x, m, data_type="phase", dmin=0, dmax=dmax)
                except NotImplementedError:
                    continue
                w.writerow([name, m, dmax, int(ai), fmt(float(af)), int(d), fmt(float(rho))])


def gen_welch(w):
    w.writerow(["dataset", "kind", "nperseg", "f", "psd"])
    for name, x in datasets():
        y = np.diff(x)
        for nps in (256,):
            if nps > len(y):
                continue
            f, p = scipy.signal.welch(y, fs=1.0, window="hann", nperseg=nps, noverlap=nps // 2,
                                      detrend="constant", scaling="density", return_onesided=True)
            for fi, pi in zip(f, p):
                w.writerow([name, "Sy", nps, fmt(float(fi)), fmt(float(pi))])
            f, p = scipy.signal.welch(x, fs=1.0, window="hann", nperseg=nps, noverlap=nps // 2,
                                      detrend="constant", scaling="density", return_onesided=True)
            for fi, pi in zip(f, p):
                w.writerow([name, "Sx", nps, fmt(float(fi)), fmt(float(pi))])


def gen_dynamic(w):
    w.writerow(["dataset", "kind", "window", "start", "m", "dev"])
    for name, x in datasets():
        if not name.startswith("alpha_"):
            continue
        win, step = 1024, 512
        for s in range(0, len(x) - win + 1, step):
            seg = x[s:s + win]
            for kind, fn in (("dadev", at.oadev), ("dhdev", at.ohdev)):
                taus, devs, _, _ = fn(seg, rate=1.0, data_type="phase", taus=[1.0, 2.0, 4.0, 16.0, 64.0, 256.0])
                for t, d in zip(taus, devs):
                    w.writerow([name, kind, win, s + 1, int(round(t)), fmt(float(d))])


# Samples blanked for the gap-tolerant check, 1-based and inclusive; the same
# spans are in test/validation/CrossCheck.jl.
GAP_SPANS = [(101, 160), (777, 777), (1500, 1503), (3001, 3100)]


def gen_gaps(w):
    w.writerow(["dataset", "m", "n", "dev"])
    for name, x in datasets():
        if name not in ("alpha_+2", "alpha_+0", "alpha_-2"):
            continue
        x = x.copy()
        for lo, hi in GAP_SPANS:
            x[lo - 1:hi] = np.nan
        taus, devs, _, ns = at.gradev(x, rate=1.0, data_type="phase", taus=[float(m) for m in octave(len(x) // 4)])
        for t, d, n in zip(taus, devs, ns):
            w.writerow([name, int(round(t)), int(n), fmt(float(d))])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(FIX, "allantools"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for fname, gen in (("edf.csv", gen_edf), ("ci.csv", gen_ci), ("noiseid.csv", gen_noiseid),
                       ("welch.csv", gen_welch), ("dynamic.csv", gen_dynamic), ("gaps.csv", gen_gaps),
                       ("devs.csv", gen_devs)):
        with open(os.path.join(args.out, fname), "w", newline="") as f:
            gen(csv.writer(f, lineterminator="\n"))
        print("wrote", fname, file=sys.stderr)
    with open(os.path.join(args.out, "VERSION"), "w") as f:
        f.write("allantools %s\nnumpy %s\nscipy %s\n" % (at.__version__, np.__version__,
                                                         __import__("scipy").__version__))


if __name__ == "__main__":
    main()

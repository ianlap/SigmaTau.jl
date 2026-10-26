#!/usr/bin/env python3
"""Time allantools on the record written by bench_sigmatau.jl and print the
side-by-side table.

    python bench/bench_allantools.py [N] [--report-only]
    python bench/bench_allantools.py --input <path> [--tau0 S] [--tag NAME] [--report-only]

Synthetic (default): reads bench/results/input_<N>.txt and
bench/results/sigmatau_<N>.json. Each kernel is timed as the minimum of 5
calls (one call for mtotdev and htotdev).

--input: reads the last whitespace-separated column of <path> as phase at
--tau0 (default 1 s), checks it against bench/results/input_<tag>.txt, and
times each kernel with a single call. <tag> defaults to "measured", matching
bench_sigmatau.jl.

Writes bench/results/allantools_<key>.json and prints a Markdown table (and the
end-to-end default-call table when the SigmaTau JSON carries one). With
--report-only it skips the timing and prints the tables from the saved JSON.

Times are printed to three significant figures and each speed-up is the ratio
of the two printed times: whole numbers from 10× up, one decimal below.
"""
import argparse
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
import allantools

ap = argparse.ArgumentParser()
ap.add_argument("N", nargs="?", type=int, default=25_000)
ap.add_argument("--input", type=Path)
ap.add_argument("--tau0", type=float, default=1.0)
ap.add_argument("--tag", default="measured")
ap.add_argument("--report-only", action="store_true")
opt = ap.parse_args()

OUT = Path(__file__).resolve().parent / "results"
key = opt.tag if opt.input else str(opt.N)
sig = json.loads((OUT / f"sigmatau_{key}.json").read_text())


def read_last_column(path):
    vals = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in "#%":
                continue
            vals.append(float(s.split()[-1]))
    return np.array(vals)


if opt.input:
    x = read_last_column(opt.input)
    ref = np.loadtxt(OUT / f"input_{key}.txt")
    assert np.array_equal(x, ref), "samples differ from those SigmaTau timed"
else:
    x = np.loadtxt(OUT / f"input_{key}.txt")
assert len(x) == sig["N"], f"N = {len(x)} but SigmaTau timed N = {sig['N']}"
assert math.isclose(opt.tau0, sig.get("tau0", 1.0)), "tau0 differs from SigmaTau's"
taus = np.array(sig["m_values"], dtype=float) * opt.tau0

KERNELS = {
    "adev": allantools.oadev,     # SigmaTau's adev is the overlapping estimator
    "mdev": allantools.mdev,
    "hdev": allantools.ohdev,
    "tdev": allantools.tdev,
    "totdev": allantools.totdev,
    "mtotdev": allantools.mtotdev,
    "htotdev": allantools.htotdev,
}


def best_of(fn, reps):
    best = float("inf")
    for _ in range(reps):
        t0 = time.perf_counter()
        fn(x, rate=1.0 / opt.tau0, data_type="phase", taus=taus)
        best = min(best, time.perf_counter() - t0)
    return best


def fmt(t):
    """Three significant figures in µs, ms or s; returns (text, value shown in s)."""
    for unit, scale in (("µs", 1e-6), ("ms", 1e-3), ("s", 1.0)):
        v = t / scale
        if v < 1000 or unit == "s":
            break
    dec = max(0, 2 - math.floor(math.log10(v)))
    v = round(v, dec)
    if v >= 10 ** (3 - dec) and dec > 0:   # rounding crossed a decade, e.g. 9.996 -> 10.0
        dec -= 1
    return f"{v:,.{dec}f} {unit}", v * scale


def speedup(a, s):
    r = a / s
    return f"{r:,.0f}×" if r >= 10 else f"{r:.1f}×"


if opt.report_only:
    saved = json.loads((OUT / f"allantools_{key}.json").read_text())
    for field in ("N", "tau0", "m_values"):
        assert saved.get(field) == sig.get(field), \
            f"saved allantools run has {field} = {saved.get(field)}, SigmaTau run has {sig.get(field)}; re-time"
    res = saved["seconds"]
else:
    print(f"allantools {allantools.__version__}, numpy {np.__version__}, "
          f"Python {platform.python_version()}, N = {len(x)}, tau0 = {opt.tau0} s")
    res, reps_used = {}, {}
    for name, fn in KERNELS.items():
        reps = 1 if opt.input or name in ("mtotdev", "htotdev") else 5
        t = best_of(fn, reps)
        res[name], reps_used[name] = t, reps
        print(f"  {name:<8s} {1e3 * t:12.3f} ms (min of {reps})", flush=True)
    (OUT / f"allantools_{key}.json").write_text(json.dumps(
        {"library": f"allantools {allantools.__version__}", "numpy": np.__version__,
         "python": platform.python_version(), "N": len(x), "tau0": opt.tau0,
         "m_values": sig["m_values"], "reps": reps_used, "seconds": res}, indent=2) + "\n")

print()
print(f"| Kernel | SigmaTau ({sig['threads']} threads) | allantools | Speed-up |")
print("|---|---:|---:|---:|")
for name in KERNELS:
    (s_txt, s), (a_txt, a) = fmt(sig["seconds"][name]), fmt(res[name])
    print(f"| {name} | {s_txt} | {a_txt} | {speedup(a, s)} |")

if "seconds_default_call" in sig:
    print()
    print(f"| Kernel | SigmaTau kernel-only | SigmaTau default call | allantools | "
          f"Speed-up (allantools ÷ default call) |")
    print("|---|---:|---:|---:|---:|")
    for name in KERNELS:
        k_txt, _ = fmt(sig["seconds"][name])
        (d_txt, d), (a_txt, a) = fmt(sig["seconds_default_call"][name]), fmt(res[name])
        print(f"| {name} | {k_txt} | {d_txt} | {a_txt} | {speedup(a, d)} |")

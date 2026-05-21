#!/usr/bin/env python3
"""Empirical scaling fits for allantools, mirroring `scaling.jl`.

Same shape: T(N) ≈ a · N^b, fit in log-log space across an N grid.
Uses `timeit.Timer.repeat` (closest analog to Julia's `@belapsed` — best
of repeated measurements) so the numbers line up with the Julia side.

What's measured: the bare kernel call (e.g. `at.oadev(x, …)`). allantools
does NOT auto-compute confidence intervals or bias correction — those are
separate functions the user would call manually — so this measures less
work than the SigmaTau side's full-API timing. The renderer notes the
asymmetry.

Run from repo root:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
      .bench-venv/bin/python benchmarks/bench/scaling_allantools.py

Output: benchmarks/bench/scaling_allantools.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import timeit
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import allantools as at  # type: ignore[import-not-found]


# Per-kernel mapping. Names follow the Julia side; allantools function
# objects are looked up at module load. `None` entries (mhdev on older
# allantools, htdev / mhtotdev never) are skipped silently.
KERNELS = [
    ("adev",     at.oadev),
    ("mdev",     at.mdev),
    ("tdev",     at.tdev),
    ("hdev",     at.ohdev),
    ("mhdev",    getattr(at, "mhdev", None)),
    ("htdev",    None),                          # not implemented in allantools
    ("totdev",   at.totdev),
    ("mtotdev",  at.mtotdev),
    ("htotdev",  at.htotdev),
    ("mhtotdev", None),                          # not implemented in allantools
]

# N grid for the display table. Mirrors `scaling.jl`. The very-slow
# kernels (mtotdev, htotdev) need a lower cap on the allantools side
# because per-call time runs to minutes well before SigmaTau's 65 k cap.
NS = [1 << k for k in range(10, 19, 2)]    # 1 024, 4 096, 16 k, 65 k, 262 k
VERY_SLOW_MAX = 4096                       # mtotdev, htotdev
SLOW_MAX = 65_536                          # totdev (allantools' totdev is mild)
VERY_SLOW_KERNELS = {"mtotdev", "htotdev"}
SLOW_KERNELS = {"totdev"}


def _human_n(n: int) -> str:
    return f"{n >> 10}k" if n >= 1024 else str(n)


def _m_values(n: int, kernel: str) -> np.ndarray:
    """Octave-spaced m grid matching `_default_m_values` on the Julia side."""
    if kernel in ("adev", "totdev"):
        m_max = (n - 1) // 2
    elif kernel in ("mdev", "tdev", "mtotdev"):
        m_max = n // 3
    elif kernel in ("hdev", "htotdev"):
        m_max = (n - 1) // 3
    elif kernel in ("mhdev",):
        m_max = n // 4
    else:
        m_max = (n - 1) // 2
    if m_max < 1:
        return np.array([1], dtype=int)
    cap = int(math.floor(math.log2(m_max)))
    return np.array([1 << k for k in range(0, cap + 1)], dtype=int)


def _fit_loglog(ns, ts):
    lx = np.log(np.asarray(ns, dtype=float))
    ly = np.log(np.asarray(ts, dtype=float))
    n = len(lx)
    mx, my = lx.mean(), ly.mean()
    Sxx = np.sum((lx - mx) ** 2)
    Sxy = np.sum((lx - mx) * (ly - my))
    b = Sxy / Sxx
    log_a = my - b * mx
    a = math.exp(log_a)
    pred = log_a + b * lx
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - my) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return a, b, r2


def _belapsed(fn, x, taus, rate, budget_s, max_evals=100):
    """Closest analog of Julia's @belapsed: best per-call time, picking
    evals-per-sample so each sample takes ~50 ms, then `repeat(samples)`
    so total work fits the budget."""
    # Pilot run to size evals.
    t0 = timeit.default_timer()
    fn(x, rate=rate, data_type="phase", taus=taus)
    single = max(timeit.default_timer() - t0, 1e-9)
    evals = max(1, min(max_evals, int(0.05 / single)))
    samples = max(1, int(budget_s / (single * evals)))
    timer = timeit.Timer(
        lambda fn=fn: fn(x, rate=rate, data_type="phase", taus=taus)
    )
    times = timer.repeat(repeat=samples, number=evals)
    return min(times) / evals


def fit_scaling(seconds: float, predict_at: int, out_path: Path) -> dict:
    print(f"Single-thread (OMP/OpenBLAS/MKL = 1). "
          f"allantools = {getattr(at, '__version__', '?')}. "
          f"Per-measurement budget = {seconds:.1f}s.\n")
    # Header
    hdr_n = "  ".join(f"N={_human_n(n):>5s}" for n in NS)
    print(f"{'kernel':<9} {hdr_n}  |  T ≈ a · N^b           R²    | predict T(N={_human_n(predict_at)})")
    print("─" * (10 + len(hdr_n) + 50))

    fits = []
    for name, fn in KERNELS:
        row = f"{name:<9} "
        if fn is None:
            row += "  ".join(["         —"] * len(NS))
            row += "  |  (not in allantools)"
            print(row)
            fits.append({"kernel": name, "implemented": False})
            continue

        if name in VERY_SLOW_KERNELS:
            cap = VERY_SLOW_MAX
        elif name in SLOW_KERNELS:
            cap = SLOW_MAX
        else:
            cap = NS[-1]

        ns_measured = []
        ts_measured = []
        cells = []
        for n in NS:
            if n > cap:
                cells.append(f"         —")
                continue
            rng = np.random.default_rng(0)
            x = rng.standard_normal(n)
            taus = _m_values(n, name).astype(float)
            t = _belapsed(fn, x, taus, rate=1.0, budget_s=seconds)
            ns_measured.append(n)
            ts_measured.append(t)
            cells.append(f"{t:8.2e}s")
        row += "  ".join(cells)

        if len(ns_measured) >= 2:
            a, b, r2 = _fit_loglog(ns_measured, ts_measured)
            t_pred = a * float(predict_at) ** b
            row += f"  |  {a:.3e} · N^{b:.3f}  {r2:5.3f} | {t_pred:9.3f} s"
            fits.append({
                "kernel": name, "implemented": True,
                "a": a, "b": b, "r2": r2,
                "Ns": ns_measured, "times": ts_measured,
            })
        else:
            row += "  |  (too few points to fit)"
            fits.append({"kernel": name, "implemented": True,
                         "Ns": ns_measured, "times": ts_measured})
        print(row)

    payload = {
        "meta": {
            "library": "allantools",
            "python_version": sys.version.split()[0],
            "allantools_version": getattr(at, "__version__", "?"),
            "threads": int(os.environ.get("OMP_NUM_THREADS", "1")),
            "kwargs": "bare kernel; CI and bias correction NOT auto-computed",
            "predict_at": predict_at,
            "seconds_budget": seconds,
        },
        "fits": fits,
    }
    out_path.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {out_path}")
    return payload


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--seconds", type=float, default=1.0,
                   help="Per-measurement budget in seconds.")
    p.add_argument("--predict-at", type=int, default=1_000_000,
                   help="N value to report a predicted runtime for.")
    p.add_argument("--out", type=Path,
                   default=Path(__file__).parent / "scaling_allantools.json")
    args = p.parse_args()
    fit_scaling(args.seconds, args.predict_at, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

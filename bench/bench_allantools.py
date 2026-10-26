#!/usr/bin/env python3
"""Time allantools on the record written by bench_sigmatau.jl and print the
side-by-side table.

    python bench/bench_allantools.py [N] [--report-only]

Reads bench/results/input_<N>.txt and bench/results/sigmatau_<N>.json, writes
bench/results/allantools_<N>.json, prints a Markdown table. With
--report-only it skips the timing and prints the table from the saved JSON.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import allantools

args = [a for a in sys.argv[1:] if not a.startswith("--")]
N = int(args[0]) if args else 25_000
REPORT_ONLY = "--report-only" in sys.argv
OUT = Path(__file__).resolve().parent / "results"

x = np.loadtxt(OUT / f"input_{N}.txt")
sig = json.loads((OUT / f"sigmatau_{N}.json").read_text())
taus = np.array(sig["m_values"], dtype=float)  # tau0 = 1 s, so tau = m

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
        fn(x, rate=1.0, data_type="phase", taus=taus)
        best = min(best, time.perf_counter() - t0)
    return best


if REPORT_ONLY:
    res = json.loads((OUT / f"allantools_{N}.json").read_text())["seconds"]
else:
    print(f"allantools {allantools.__version__}, numpy {np.__version__}, N = {N}")
    res = {}
    for name, fn in KERNELS.items():
        reps = 1 if name in ("mtotdev", "htotdev") else 5
        t = best_of(fn, reps)
        res[name] = t
        print(f"  {name:<8s} {1e3 * t:12.3f} ms", flush=True)
    (OUT / f"allantools_{N}.json").write_text(json.dumps(
        {"library": f"allantools {allantools.__version__}", "numpy": np.__version__,
         "N": N, "m_values": sig["m_values"], "seconds": res}, indent=2))

print()
print(f"| Kernel | SigmaTau ({sig['threads']} threads) | allantools | Speedup |")
print("|---|---:|---:|---:|")
for name in KERNELS:
    s, a = sig["seconds"][name], res[name]
    print(f"| {name} | {1e3 * s:.3g} ms | {1e3 * a:.3g} ms | {a / s:,.0f}× |")

#!/usr/bin/env python3
"""Minimal `@btime`-style bench for allantools: oadev, mdev, ohdev.

Mirrors the Julia script `btime_sigmatau.jl` so the numbers line up
directly. Synthetic standard-normal phase data; no file IO. Pinned to a
single BLAS thread for a fair comparison.

Python has no exact `@btime` equivalent, but `timeit.Timer.repeat`
(minimum of repeated runs) gives the same "best steady-state time per
call" reading. We also report the per-call median allocation footprint
via `tracemalloc` peak so the output shape matches BenchmarkTools'
"time / bytes / allocs".

Usage:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
      python3 benchmarks/bench/btime_allantools.py            # N = 2**15

    python3 benchmarks/bench/btime_allantools.py 131072       # bigger record
"""
from __future__ import annotations

import os
import sys
import timeit
import tracemalloc

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import allantools as at  # type: ignore[import-not-found]

KERNELS = [
    ("oadev", at.oadev),
    ("mdev",  at.mdev),
    ("ohdev", at.ohdev),
]


def _m_grid(N: int) -> np.ndarray:
    cap = int(np.floor(np.log2(N / 3)))
    return np.array([1 << k for k in range(0, cap + 1)], dtype=int)


def _fmt_time(s: float) -> str:
    if s < 1e-6:
        return f"{s*1e9:7.2f} ns"
    if s < 1e-3:
        return f"{s*1e6:7.2f} μs"
    if s < 1.0:
        return f"{s*1e3:7.2f} ms"
    return f"{s:7.3f}  s"


def _fmt_bytes(b: int) -> str:
    if b < 1024:
        return f"{b:6d}  B"
    if b < 1024**2:
        return f"{b/1024:6.1f} KiB"
    if b < 1024**3:
        return f"{b/1024**2:6.1f} MiB"
    return f"{b/1024**3:6.1f} GiB"


def btime_run(N: int = 1 << 15, seed: int = 0,
              budget_s: float = 5.0, max_evals: int = 100) -> None:
    rng = np.random.default_rng(seed)
    x = rng.standard_normal(N)
    ms = _m_grid(N)
    taus = ms.astype(float)  # rate = 1.0, so taus = m
    print(f"N = {N}, {len(ms)} m values ({ms[0]} … {ms[-1]}), "
          f"OMP_NUM_THREADS = {os.environ.get('OMP_NUM_THREADS','?')}, "
          f"allantools = {getattr(at, '__version__', '?')}")
    print("-" * 54)

    for name, fn in KERNELS:
        # One untimed warmup to pay any one-off imports / caches.
        fn(x, rate=1.0, data_type="phase", taus=taus)

        # Pilot run to pick an evals-per-sample so each sample is ~50 ms.
        t0 = timeit.default_timer()
        fn(x, rate=1.0, data_type="phase", taus=taus)
        single = max(timeit.default_timer() - t0, 1e-9)
        evals = max(1, min(max_evals, int(0.05 / single)))
        samples = max(1, int(budget_s / (single * evals)))

        timer = timeit.Timer(
            lambda fn=fn: fn(x, rate=1.0, data_type="phase", taus=taus)
        )
        times = timer.repeat(repeat=samples, number=evals)
        best_per_call = min(times) / evals

        tracemalloc.start()
        fn(x, rate=1.0, data_type="phase", taus=taus)
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        print(f"{name:<6s} {_fmt_time(best_per_call)}  "
              f"(best of {samples}×{evals})   peak alloc ≈ {_fmt_bytes(peak)}")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 1 << 15
    btime_run(N=n)

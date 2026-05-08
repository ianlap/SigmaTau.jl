#!/usr/bin/env python3
"""Wall-clock benchmark for allantools on a real clock recording.

Mirrors `bench_sigmatau.jl`: same per-octave m grid, same kernels, same
file. Pinned to single-thread (OMP/MKL/OpenBLAS) so the comparison is
fair against SigmaTau (BLAS=1) and Stable32 (single-threaded).

Run from repo root:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
      python3 validation/bench/bench_allantools.py \\
        reference/clock_data/6krb25apr.txt

Optional `--kernels adev,mtotdev` to time a subset.

Note: allantools.tdev / ldev derive from mdev / mhdev. We time the
public allantools entry points, which on long records may include some
redundant differencing — that's faithful to "what a user would call".
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

# ANSI tag stripped from the env vars when imported from Julia harness.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import allantools as at  # type: ignore[import-not-found]

# (name, function, allantools-supports-modified-hadamard?)
KERNELS = [
    ("adev",     at.oadev,                 True),
    ("mdev",     at.mdev,                  True),
    ("tdev",     at.tdev,                  True),
    ("hdev",     at.ohdev,                 True),
    ("mhdev",    getattr(at, "mhdev", None), False),
    ("ldev",     None,                       False),
    ("totdev",   at.totdev,                True),
    ("mtotdev",  at.mtotdev,               True),
    ("htotdev",  at.htotdev,               True),
    ("mhtotdev", None,                       False),
]


def load_phase_2col(path: Path) -> tuple[np.ndarray, int, float]:
    """Loader matching `_loader.jl`. Returns (phase, N, tau0_seconds)."""
    arr = np.loadtxt(path, dtype=np.float64)
    if arr.ndim != 2 or arr.shape[1] != 2:
        raise ValueError(f"{path}: expected 2-column whitespace-separated, got shape {arr.shape}")
    t_mjd = arr[:, 0]
    x = arr[:, 1]
    tau0 = float(np.median(np.diff(t_mjd))) * 86400.0
    return x, x.size, tau0


def bench_m_values(N: int) -> np.ndarray:
    cap = int(np.floor(np.log2(N / 3)))
    return np.array([1 << k for k in range(0, cap + 1)], dtype=int)


def warmup() -> None:
    """Touch every kernel once on a tiny array — pays imports/JIT/cache costs
    out of band so the timed runs measure steady-state work."""
    print("Warming up allantools (2048-sample dummy run) …", flush=True)
    rng = np.random.default_rng(0)
    x = rng.standard_normal(2048)
    taus = np.asarray([1.0, 2.0, 4.0, 8.0])
    for name, fn, _ in KERNELS:
        if fn is None:
            continue
        t0 = time.monotonic()
        try:
            fn(x, rate=1.0, data_type="phase", taus=taus)
        except Exception as e:  # noqa: BLE001
            print(f"  warm {name:<9s} FAILED ({type(e).__name__}: {e})", flush=True)
            continue
        print(f"  warm {name:<9s} {time.monotonic() - t0:.3f}s", flush=True)


def bench_one(path: Path, kernels_filter: set[str] | None) -> list[tuple[str, float]]:
    print(f"\n=== {path.name} ===", flush=True)
    x, N, tau0 = load_phase_2col(path)
    m_values = bench_m_values(N)
    taus = (m_values * tau0).astype(float)
    print(f"  N = {N}, τ₀ = {tau0:.6g} s, {len(m_values)} m values "
          f"(1 … {m_values[-1]})", flush=True)
    print(f"  rate = {1.0/tau0:.6g}", flush=True)
    print("  " + "-" * 46, flush=True)

    results: list[tuple[str, float]] = []
    for name, fn, _ in KERNELS:
        if kernels_filter is not None and name not in kernels_filter:
            continue
        if fn is None:
            print(f"  {name:<9s}  (not implemented in allantools)", flush=True)
            continue
        t0 = time.monotonic()
        fn(x, rate=1.0/tau0, data_type="phase", taus=taus)
        elapsed = time.monotonic() - t0
        results.append((name, elapsed))
        print(f"  {name:<9s}  {elapsed:>10.3f} s", flush=True)

    print("  " + "-" * 46, flush=True)
    print(f"  total: {sum(t for _, t in results):.3f} s", flush=True)
    return results


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("path", type=Path)
    p.add_argument("--kernels", type=str, default=None,
                   help="Comma-separated subset, e.g. 'adev,mtotdev'.")
    args = p.parse_args()

    kernels_filter = (
        {k.strip() for k in args.kernels.split(",")} if args.kernels else None
    )

    warmup()
    bench_one(args.path, kernels_filter)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

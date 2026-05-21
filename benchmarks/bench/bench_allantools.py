#!/usr/bin/env python3
"""Wall-clock benchmark for allantools on a real clock recording.

Mirrors `bench_sigmatau.jl`: same per-octave m grid, same kernels, same
file. Pinned to single-thread (OMP/MKL/OpenBLAS) so the comparison is
fair against SigmaTau (BLAS=1) and Stable32 (single-threaded).

Run from repo root:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
      python3 benchmarks/bench/bench_allantools.py \\
        reference/clock_data/6krb25apr.txt

Optional `--kernels adev,mtotdev` to time a subset.

Note: allantools.tdev / ldev derive from mdev / mhdev. We time the
public allantools entry points, which on long records may include some
redundant differencing — that's faithful to "what a user would call".
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import threading
import time
from pathlib import Path

# ANSI tag stripped from the env vars when imported from Julia harness.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import psutil  # type: ignore[import-not-found]
import allantools as at  # type: ignore[import-not-found]


class RSSPeakSampler:
    """Background sampler tracking peak RSS for the current process.

    Run a kernel inside `with sampler:` and read `sampler.peak_delta` after.
    Polls every 25 ms; ample resolution for our 1-300 s kernels.
    """

    def __init__(self, interval_s: float = 0.025):
        self._proc = psutil.Process()
        self._interval = interval_s
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self.baseline = 0
        self.peak = 0

    def __enter__(self):
        self.baseline = self._proc.memory_info().rss
        self.peak = self.baseline
        self._stop.clear()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, *_):
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
        # Final read in case the kernel finished between samples.
        self.peak = max(self.peak, self._proc.memory_info().rss)

    def _run(self):
        while not self._stop.is_set():
            try:
                rss = self._proc.memory_info().rss
                if rss > self.peak:
                    self.peak = rss
            except psutil.Error:
                pass
            self._stop.wait(self._interval)

    @property
    def peak_delta(self) -> int:
        return self.peak - self.baseline

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


def bench_one(path: Path, kernels_filter: set[str] | None,
              m_max: int | None = None,
              per_m: bool = False) -> tuple[list[dict], dict]:
    print(f"\n=== {path.name} ===", flush=True)
    x, N, tau0 = load_phase_2col(path)
    m_values = bench_m_values(N)
    if m_max is not None:
        m_values = m_values[m_values <= m_max]
    taus = (m_values * tau0).astype(float)
    print(f"  N = {N}, tau0 = {tau0:.6g} s, {len(m_values)} m values "
          f"({m_values[0]} ... {m_values[-1]})", flush=True)
    print(f"  rate = {1.0/tau0:.6g}", flush=True)
    print(f"  threads (OMP/MKL/OpenBLAS) = "
          f"{os.environ.get('OMP_NUM_THREADS','?')}", flush=True)
    print(f"  mode = {'per-m sweep' if per_m else 'one-shot (all taus per call)'}",
          flush=True)
    print("  " + "-" * 60, flush=True)

    results: list[dict] = []
    for name, fn, _ in KERNELS:
        if kernels_filter is not None and name not in kernels_filter:
            continue
        if fn is None:
            print(f"  {name:<9s}  (not implemented in allantools)", flush=True)
            results.append({"kernel": name, "time_s": None,
                            "rss_delta_bytes": None, "implemented": False})
            continue
        sampler = RSSPeakSampler()
        kernel_t0 = time.monotonic()
        if per_m:
            per_m_times: list[tuple[int, float]] = []
            print(f"  {name:<9s}  starting per-m sweep ...", flush=True)
            with sampler:
                for m, tau in zip(m_values.tolist(), taus.tolist()):
                    t0 = time.monotonic()
                    fn(x, rate=1.0/tau0, data_type="phase",
                       taus=np.array([tau], dtype=float))
                    dt = time.monotonic() - t0
                    per_m_times.append((m, dt))
                    cumulative = time.monotonic() - kernel_t0
                    print(f"    m={m:<7d} t={dt:>10.3f} s   "
                          f"(cum {cumulative:>10.3f} s)", flush=True)
            extra = {"per_m": [{"m": m, "t_s": t} for m, t in per_m_times]}
        else:
            print(f"  {name:<9s}  running ...", flush=True)
            with sampler:
                fn(x, rate=1.0/tau0, data_type="phase", taus=taus)
            extra = {}
        elapsed = time.monotonic() - kernel_t0
        results.append({"kernel": name, "time_s": elapsed,
                        "rss_delta_bytes": sampler.peak_delta,
                        "implemented": True, **extra})
        print(f"  {name:<9s}  TOTAL {elapsed:>10.3f} s   "
              f"peak deltaRSS={sampler.peak_delta/2**20:>8.1f} MiB", flush=True)

    print("  " + "-" * 60, flush=True)
    total = sum(r["time_s"] for r in results if r["time_s"] is not None)
    print(f"  total: {total:.3f} s", flush=True)
    meta = {
        "file": path.name,
        "N": int(N),
        "tau0_s": float(tau0),
        "n_m": int(len(m_values)),
        "m_min": int(m_values[0]),
        "m_max": int(m_values[-1]),
        "threads": int(os.environ.get("OMP_NUM_THREADS", "1")),
        "python_version": sys.version.split()[0],
        "allantools_version": getattr(at, "__version__", "?"),
    }
    return results, meta


def bench_synth(synth_dir: Path, kernels_filter: set[str] | None,
                tau0: float = 1.0,
                m_max: int | None = None) -> tuple[list[dict], dict]:
    """Statistical bench: run all KERNELS on each .txt realization file,
    record per-realization wall + peak ΔRSS."""
    files = sorted(synth_dir.glob("*.txt"))
    if not files:
        raise FileNotFoundError(f"no .txt files in {synth_dir}")
    # probe N
    x0 = np.loadtxt(files[0], dtype=np.float64)
    N = x0.size
    m_values = bench_m_values(N)
    if m_max is not None:
        m_values = m_values[m_values <= m_max]
    taus = (m_values * tau0).astype(float)

    print(f"\n=== synthetic bench ===", flush=True)
    print(f"  dir = {synth_dir}", flush=True)
    print(f"  reals = {len(files)}, N = {N}, tau0 = {tau0:.6g} s", flush=True)
    print(f"  {len(m_values)} m values ({m_values[0]} ... {m_values[-1]})", flush=True)
    print(f"  threads (OMP/MKL/OpenBLAS) = "
          f"{os.environ.get('OMP_NUM_THREADS','?')}", flush=True)
    print("  " + "-" * 60, flush=True)

    per_kernel: dict[str, list[dict]] = {}
    for name, fn, _ in KERNELS:
        if kernels_filter is not None and name not in kernels_filter:
            continue
        per_kernel[name] = []

    for i, path in enumerate(files):
        x = np.loadtxt(path, dtype=np.float64)
        for name, fn, _ in KERNELS:
            if kernels_filter is not None and name not in kernels_filter:
                continue
            if fn is None:
                # Skip silently for synth mode; will appear as empty per_realization
                continue
            sampler = RSSPeakSampler()
            t0 = time.monotonic()
            with sampler:
                fn(x, rate=1.0/tau0, data_type="phase", taus=taus)
            dt = time.monotonic() - t0
            per_kernel.setdefault(name, []).append({
                "realization": i,
                "time_s": dt,
                "rss_delta_bytes": sampler.peak_delta,
            })
        print(f"  realization {i+1:3d} / {len(files)} done", flush=True)

    print("  " + "-" * 60, flush=True)
    for name in per_kernel:
        ts = [r["time_s"] for r in per_kernel[name]]
        if not ts:
            continue
        mean = sum(ts) / len(ts)
        var = sum((t - mean) ** 2 for t in ts) / max(len(ts) - 1, 1)
        std = var ** 0.5
        med = sorted(ts)[len(ts) // 2]
        print(f"  {name:<9s}  mean={mean:>9.4f} s   "
              f"median={med:>9.4f} s   std={std:>9.4f} s", flush=True)

    results = [{"kernel": name,
                "per_realization": per_kernel.get(name, []),
                "implemented": next((not (fn is None)
                                     for n, fn, _ in KERNELS if n == name), True)}
               for name in [k for k, _, _ in KERNELS]
               if kernels_filter is None or name in kernels_filter]
    meta = {
        "mode": "synth",
        "synth_dir": str(synth_dir),
        "n_reals": len(files),
        "N": int(N),
        "tau0_s": float(tau0),
        "m_values": m_values.tolist(),
        "n_m": int(len(m_values)),
        "m_min": int(m_values[0]),
        "m_max": int(m_values[-1]),
        "threads": int(os.environ.get("OMP_NUM_THREADS", "1")),
        "python_version": sys.version.split()[0],
        "allantools_version": getattr(at, "__version__", "?"),
    }
    return results, meta


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("path", type=Path, nargs="?",
                   help="Single record .txt file (legacy mode). Omit when using --synth.")
    p.add_argument("--synth", type=Path, default=None,
                   help="Directory of realization .txt files (statistical bench mode).")
    p.add_argument("--kernels", type=str, default=None,
                   help="Comma-separated subset, e.g. 'adev,mtotdev'.")
    p.add_argument("--out", type=Path,
                   default=Path(__file__).parent / "results_allantools.json",
                   help="Where to write JSON results.")
    p.add_argument("--m-max", type=int, default=None,
                   help="Cap m grid at this value (e.g. 512).")
    p.add_argument("--per-m", action="store_true",
                   help="Diagnostic mode: time each m separately (single-record only).")
    args = p.parse_args()

    kernels_filter = (
        {k.strip() for k in args.kernels.split(",")} if args.kernels else None
    )

    warmup()
    if args.synth is not None:
        results, meta = bench_synth(args.synth, kernels_filter, m_max=args.m_max)
    else:
        if args.path is None:
            p.error("either <path> or --synth <dir> is required")
        results, meta = bench_one(args.path, kernels_filter,
                                  m_max=args.m_max, per_m=args.per_m)
        meta["mode"] = "per-m" if args.per_m else "one-shot"
    if args.m_max is not None:
        meta["m_max_cap"] = args.m_max
    args.out.write_text(json.dumps({"meta": meta, "results": results}, indent=2))
    print(f"wrote {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

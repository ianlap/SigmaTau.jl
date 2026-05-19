#!/usr/bin/env python3
"""Time allantools and numba kernels on the shared phase datasets.

Outputs one CSV row per (dataset, kernel, implementation) into
results/python_results.csv, plus a small JSON sidecar with the actual
deviation values for cross-checking against the Julia run.
"""
import os, sys, time, json, csv
import numpy as np

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)

import allantools  # noqa: E402
from numba_kernels import mtotdev_numba, htotdev_numba  # noqa: E402

DATA = os.path.join(HERE, "data")
RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

TAU0 = 1.0
# allantools mtotdev/htotdev cost grows ~ O(N*M) with M ∝ N, so N=16000
# would take ~8 min/kernel. Run allantools only at the smaller sizes; numba
# and Julia run at all three.
NS_ALL = (1000, 4000, 16000)
NS_ALLANTOOLS = (1000, 4000)


def octave_ms(N, max_m_frac=1.0 / 3.0):
    m = 1
    out = []
    while m <= int(N * max_m_frac):
        out.append(m)
        m *= 2
    return np.array(out, dtype=np.int64)


def time_call(fn, repeats=3):
    best = float("inf")
    result = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = fn()
        dt = time.perf_counter() - t0
        if dt < best:
            best = dt
    return best, result


def allantools_mtotdev(x, taus_s, tau0):
    t, devs, _, _ = allantools.mtotdev(
        x, rate=1.0 / tau0, data_type="phase", taus=taus_s
    )
    return devs


def allantools_htotdev(x, taus_s, tau0):
    t, devs, _, _ = allantools.htotdev(
        x, rate=1.0 / tau0, data_type="phase", taus=taus_s
    )
    return devs


def main():
    rows = []
    devs_dump = {}

    # Numba warmup on tiny input to avoid charging JIT to the smallest dataset
    print("# warming up numba kernels...", flush=True)
    _warm_x = np.cumsum(np.random.default_rng(0).standard_normal(64))
    _ = mtotdev_numba(_warm_x, np.array([1, 2], dtype=np.int64), TAU0)
    _ = htotdev_numba(_warm_x, np.array([1, 2], dtype=np.int64), TAU0)
    print("# warmup done", flush=True)

    for N in NS_ALL:
        path = os.path.join(DATA, f"phase_N{N}.txt")
        x = np.loadtxt(path)
        m_values = octave_ms(N)
        taus_s = m_values.astype(np.float64) * TAU0
        print(f"\n## N={N}  m_values={list(m_values)}", flush=True)

        for kernel, at_fn, nb_fn in (
            ("mtotdev", allantools_mtotdev, mtotdev_numba),
            ("htotdev", allantools_htotdev, htotdev_numba),
        ):
            # allantools (skipped at largest N to keep runtime bounded)
            if N in NS_ALLANTOOLS:
                label = f"allantools.{kernel}"
                print(f"  timing {label} ...", flush=True)
                try:
                    t_at, devs_at = time_call(lambda: at_fn(x, taus_s, TAU0), repeats=2)
                except Exception as e:
                    print(f"   FAILED: {e}", flush=True)
                    t_at, devs_at = float("nan"), np.full(len(m_values), np.nan)
                print(f"   {t_at:.4f}s", flush=True)
                rows.append(dict(N=N, kernel=kernel, impl="allantools", seconds=t_at))
                devs_dump[f"{kernel}.allantools.N{N}"] = list(devs_at)
            else:
                print(f"  allantools.{kernel} skipped at N={N} (would take ~8 min)", flush=True)

            # numba (single-thread, post-warmup)
            label = f"numba.{kernel}"
            print(f"  timing {label} ...", flush=True)
            t_nb, devs_nb = time_call(
                lambda: nb_fn(x, m_values, TAU0), repeats=3
            )
            print(f"   {t_nb:.4f}s", flush=True)
            rows.append(dict(N=N, kernel=kernel, impl="numba", seconds=t_nb))
            devs_dump[f"{kernel}.numba.N{N}"] = list(devs_nb)

    out_csv = os.path.join(RESULTS, "python_results.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["N", "kernel", "impl", "seconds"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out_csv}")

    out_json = os.path.join(RESULTS, "python_devs.json")
    with open(out_json, "w") as f:
        json.dump(devs_dump, f, indent=2)
    print(f"wrote {out_json}")


if __name__ == "__main__":
    main()

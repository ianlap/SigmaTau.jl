#!/usr/bin/env python3
"""Stats renderer for the synthetic statistical bench.

Reads results_sigmatau_synth.json + results_allantools_synth.json,
produces a per-kernel stats table and a per-realization CSV export
(suitable for histogram/violin plotting).

Usage:
    python render_stats.py results_sigmatau_synth.json \
                           results_allantools_synth.json \
                           --csv-out per_realization.csv \
                           --kernels adev,mdev,hdev,tdev,totdev,mtotdev,htotdev
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean, median, pstdev, quantiles


def stats(values):
    if not values:
        return None
    n = len(values)
    m = mean(values)
    med = median(values)
    sd = pstdev(values) if n > 1 else 0.0
    if n >= 5:
        q = quantiles(values, n=20, method="inclusive")  # ventiles: 5th, 95th
        p5, p95 = q[0], q[18]
    else:
        sv = sorted(values)
        p5, p95 = sv[0], sv[-1]
    return {"n": n, "mean": m, "median": med, "std": sd,
            "p5": p5, "p95": p95, "min": min(values), "max": max(values)}


def fmt(x, unit="s", digits=4):
    if x is None:
        return "n/a"
    if unit == "MiB":
        return f"{x / 2**20:.2f}"
    return f"{x:.{digits}f}"


def kernel_rows(kernel, sigmatau_results, allantools_results):
    s_per = next((r["per_realization"] for r in sigmatau_results
                  if r["kernel"] == kernel), [])
    a_per = next((r["per_realization"] for r in allantools_results
                  if r["kernel"] == kernel), [])
    s_t = stats([r["time_s"] for r in s_per])
    a_t = stats([r["time_s"] for r in a_per])
    s_b = stats([r["bytes"] for r in s_per])
    a_b = stats([r.get("rss_delta_bytes", 0) for r in a_per])
    return s_t, a_t, s_b, a_b, s_per, a_per


def main():
    p = argparse.ArgumentParser()
    p.add_argument("sigmatau_json", type=Path)
    p.add_argument("allantools_json", type=Path)
    p.add_argument("--kernels", type=str,
                   default="adev,mdev,hdev,tdev,totdev,mtotdev,htotdev")
    p.add_argument("--csv-out", type=Path,
                   default=Path(__file__).parent / "per_realization.csv",
                   help="Long-format CSV: kernel,library,realization,time_s,mem_bytes")
    args = p.parse_args()

    sj = json.loads(args.sigmatau_json.read_text())
    aj = json.loads(args.allantools_json.read_text())
    kernels = [k.strip() for k in args.kernels.split(",")]

    print()
    print(f"=== Statistical bench (synthetic, {sj['meta']['n_reals']} realizations of "
          f"N={sj['meta']['N']}) ===")
    print(f"SigmaTau threads: {sj['meta']['threads']}   |   "
          f"allantools threads: {aj['meta']['threads']}")
    print(f"Julia {sj['meta']['julia_version']}   |   "
          f"Python {aj['meta']['python_version']}, "
          f"allantools {aj['meta']['allantools_version']}")
    print(f"m grid: {sj['meta']['m_min']}..{sj['meta']['m_max']} "
          f"({sj['meta']['n_m']} values)")
    print()

    # ---- time stats table ----
    print("Time per kernel call (seconds, across realizations)")
    hdr = (f"  {'kernel':<9} | "
           f"{'SigmaTau mean':>14} {'± std':>10} {'median':>10} {'p5':>10} {'p95':>10} | "
           f"{'allantools mean':>16} {'± std':>10} {'median':>10} {'p5':>10} {'p95':>10} | "
           f"{'speedup (mean)':>14}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for k in kernels:
        s_t, a_t, _, _, _, _ = kernel_rows(k, sj["results"], aj["results"])
        if s_t is None and a_t is None:
            continue
        speed = (a_t["mean"] / s_t["mean"]) if (s_t and a_t and s_t["mean"] > 0) else None
        speed_s = f"{speed:>14.1f}x" if speed is not None else f"{'-':>14}"
        st_cols = ([fmt(s_t[k2]) for k2 in ("mean", "std", "median", "p5", "p95")]
                   if s_t else ["n/a"]*5)
        at_cols = ([fmt(a_t[k2]) for k2 in ("mean", "std", "median", "p5", "p95")]
                   if a_t else ["n/a"]*5)
        print(f"  {k:<9} | "
              f"{st_cols[0]:>14} {st_cols[1]:>10} {st_cols[2]:>10} "
              f"{st_cols[3]:>10} {st_cols[4]:>10} | "
              f"{at_cols[0]:>16} {at_cols[1]:>10} {at_cols[2]:>10} "
              f"{at_cols[3]:>10} {at_cols[4]:>10} | "
              f"{speed_s}")
    print()

    # ---- memory stats table ----
    print("Memory per kernel call (MiB; SigmaTau = cumulative @timed .bytes; "
          "allantools = peak ΔRSS via psutil polling)")
    hdr2 = (f"  {'kernel':<9} | "
            f"{'SigmaTau mean':>14} {'± std':>10} {'max':>10} | "
            f"{'allantools mean':>16} {'± std':>10} {'max':>10}")
    print(hdr2)
    print("  " + "-" * (len(hdr2) - 2))
    for k in kernels:
        _, _, s_b, a_b, _, _ = kernel_rows(k, sj["results"], aj["results"])
        if s_b is None and a_b is None:
            continue
        sb_cols = ([fmt(s_b[k2], unit="MiB") for k2 in ("mean", "std", "max")]
                   if s_b else ["n/a"]*3)
        ab_cols = ([fmt(a_b[k2], unit="MiB") for k2 in ("mean", "std", "max")]
                   if a_b else ["n/a"]*3)
        print(f"  {k:<9} | "
              f"{sb_cols[0]:>14} {sb_cols[1]:>10} {sb_cols[2]:>10} | "
              f"{ab_cols[0]:>16} {ab_cols[1]:>10} {ab_cols[2]:>10}")
    print()

    # ---- CSV export (long format for plotting) ----
    with args.csv_out.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kernel", "library", "realization", "time_s", "mem_bytes"])
        for k in kernels:
            _, _, _, _, s_per, a_per = kernel_rows(k, sj["results"], aj["results"])
            for r in s_per:
                w.writerow([k, "sigmatau", r["realization"], r["time_s"], r["bytes"]])
            for r in a_per:
                w.writerow([k, "allantools", r["realization"], r["time_s"],
                            r.get("rss_delta_bytes", 0)])
    print(f"wrote {args.csv_out} (long-format per-realization data for plotting)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Merge bench_sigmatau.jl + bench_allantools.py JSON outputs into one table.

Usage:
    python render_table.py validation/bench/results_sigmatau.json \
                           validation/bench/results_allantools.json
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

KERNEL_ORDER = ["adev", "mdev", "tdev", "hdev", "mhdev", "ldev",
                "totdev", "mtotdev", "htotdev", "mhtotdev"]


def parse_kernel_filter(argv):
    """If --kernels k1,k2,... is in argv, return (filtered_argv, filter_list)."""
    for i, a in enumerate(argv):
        if a == "--kernels" and i + 1 < len(argv):
            ks = [k.strip() for k in argv[i + 1].split(",")]
            return argv[:i] + argv[i + 2:], ks
        if a.startswith("--kernels="):
            ks = [k.strip() for k in a.split("=", 1)[1].split(",")]
            return argv[:i] + argv[i + 1:], ks
    return argv, None


def fmt_time(t):
    return "n/a" if t is None else f"{t:.3f}"


def fmt_mem_mib(b):
    return "n/a" if b is None else f"{b/2**20:.1f} MiB"


def main():
    argv, kernel_filter = parse_kernel_filter(sys.argv[1:])
    sj = json.loads(Path(argv[0]).read_text())
    aj = json.loads(Path(argv[1]).read_text())
    order = kernel_filter if kernel_filter else KERNEL_ORDER

    sm = {r["kernel"]: r for r in sj["results"]}
    am = {r["kernel"]: r for r in aj["results"]}

    meta_s = sj["meta"]
    meta_a = aj["meta"]

    print()
    print(f"=== {meta_s['file']} (N={meta_s['N']}, "
          f"tau0={meta_s['tau0_s']:.6g} s, "
          f"{meta_s['n_m']} m values {meta_s['m_min']}..{meta_s['m_max']}) ===")
    print(f"SigmaTau threads: {meta_s['threads']}   |   "
          f"allantools threads: {meta_a['threads']}")
    print(f"Julia {meta_s['julia_version']}   |   "
          f"Python {meta_a['python_version']}, "
          f"allantools {meta_a['allantools_version']}")
    print("Memory metrics: SigmaTau = Julia @timed .bytes (cumulative alloc);")
    print("                allantools = peak deltaRSS via psutil polling.")
    print("                Not directly comparable; both indicate memory pressure.")
    print()

    header = f"  {'kernel':<9} {'SigmaTau (s)':>13} {'allantools (s)':>15} " \
             f"{'speedup':>8} {'SigmaTau alloc':>17} {'allantools dRSS':>17}"
    print(header)
    print("  " + "-" * (len(header) - 2))

    total_s = 0.0
    total_a = 0.0
    for k in order:
        srow = sm.get(k)
        arow = am.get(k)
        st = srow["time_s"] if srow else None
        sb = srow["bytes"] if srow else None
        at_t = arow["time_s"] if arow and arow.get("implemented") else None
        ar = arow["rss_delta_bytes"] if arow and arow.get("implemented") else None

        if st is not None:
            total_s += st
        if at_t is not None:
            total_a += at_t

        if st is not None and at_t is not None and st > 0:
            speed = f"{at_t/st:.1f}x"
        else:
            speed = "-"

        st_s = fmt_time(st)
        if at_t is not None:
            at_s = fmt_time(at_t)
        elif arow is not None and arow.get("dnf"):
            at_s = "DNF"
        elif arow is None or not arow.get("implemented"):
            at_s = "(n/a)"
        else:
            at_s = "n/a"
        sb_s = fmt_mem_mib(sb)
        if ar is not None:
            ar_s = fmt_mem_mib(ar)
        elif arow is not None and arow.get("dnf"):
            ar_s = "DNF"
        elif arow is None or not arow.get("implemented"):
            ar_s = "-"
        else:
            ar_s = "n/a"

        print(f"  {k:<9} {st_s:>13} {at_s:>15} {speed:>8} {sb_s:>17} {ar_s:>17}")

    print("  " + "-" * (len(header) - 2))
    if total_s > 0 and total_a > 0:
        print(f"  {'total wall':<9} {total_s:>13.3f} {total_a:>15.3f} "
              f"{total_a/total_s:>7.1f}x")
    else:
        print(f"  {'total wall':<9} {total_s:>13.3f} {total_a:>15.3f}")


if __name__ == "__main__":
    main()

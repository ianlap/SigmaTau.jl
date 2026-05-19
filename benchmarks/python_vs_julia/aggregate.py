#!/usr/bin/env python3
"""Combine Python + Julia CSV outputs into a single table and cross-check
that all three implementations agree on the deviation values.
"""
import csv, json, os, sys, math

HERE = os.path.dirname(__file__)
R = os.path.join(HERE, "results")


def load_csv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    py = load_csv(os.path.join(R, "python_results.csv"))
    jsingle = load_csv(os.path.join(R, "julia_results_single.csv"))
    jthreaded = load_csv(
        os.path.join(R, "julia_results_threaded4.csv")
    )

    rows = py + jsingle + jthreaded
    if not rows:
        print("no results found", file=sys.stderr)
        sys.exit(1)

    # canonicalise dtypes
    for r in rows:
        r["N"] = int(r["N"])
        r["seconds"] = float(r["seconds"])

    by_key = {}
    for r in rows:
        key = (r["N"], r["kernel"], r["impl"])
        by_key[key] = r["seconds"]

    impls = ["allantools", "numba", "sigmatau_single", "sigmatau_threaded4"]
    ns = sorted({r["N"] for r in rows})
    kernels = ["mtotdev", "htotdev"]

    print("\n# Wall-clock times (best of N repeats, seconds)\n")
    header = f"{'N':>6}  {'kernel':<8}  " + "  ".join(f"{i:>20}" for i in impls)
    print(header)
    print("-" * len(header))
    for N in ns:
        for k in kernels:
            cells = []
            for i in impls:
                v = by_key.get((N, k, i))
                cells.append("--" if v is None else f"{v:.4f}s")
            print(f"{N:>6}  {k:<8}  " + "  ".join(f"{c:>20}" for c in cells))

    print("\n# Speedup vs allantools (single-threaded)\n")
    header = f"{'N':>6}  {'kernel':<8}  " + "  ".join(f"{i:>20}" for i in impls)
    print(header)
    print("-" * len(header))
    for N in ns:
        for k in kernels:
            base = by_key.get((N, k, "allantools"))
            cells = []
            for i in impls:
                v = by_key.get((N, k, i))
                if v is None or base is None:
                    cells.append("--")
                else:
                    cells.append(f"{base / v:.1f}x")
            print(f"{N:>6}  {k:<8}  " + "  ".join(f"{c:>20}" for c in cells))

    # cross-check deviations agree
    print("\n# Deviation cross-check (max relative diff per (kernel, N))\n")
    py_devs = {}
    p = os.path.join(R, "python_devs.json")
    if os.path.exists(p):
        with open(p) as f:
            py_devs = json.load(f)

    jl_devs = {}
    for tag in ("single", "threaded4"):
        p = os.path.join(R, f"julia_devs_{tag}.txt")
        if os.path.exists(p):
            with open(p) as f:
                for ln in f:
                    if "=" not in ln:
                        continue
                    name, vals = ln.strip().split("=", 1)
                    jl_devs[(tag, name)] = [float(v) for v in vals.split(",")]

    def relmax(a, b):
        if a is None or b is None or len(a) != len(b):
            return float("nan")
        worst = 0.0
        for x, y in zip(a, b):
            if math.isnan(x) or math.isnan(y) or y == 0.0:
                continue
            r = abs(x - y) / abs(y)
            if r > worst:
                worst = r
        return worst

    for k in kernels:
        for N in ns:
            jl = jl_devs.get(("single", f"{k}.N{N}"))
            nb = py_devs.get(f"{k}.numba.N{N}")
            at = py_devs.get(f"{k}.allantools.N{N}")
            r_nb_jl = relmax(nb, jl)
            r_at_jl = relmax(at, jl) if at else float("nan")
            r_at_nb = relmax(at, nb) if at else float("nan")
            print(
                f"  N={N:>5} {k:<8}  numba vs julia: {r_nb_jl:.2e}  "
                f"allantools vs julia: {r_at_jl:.2e}  "
                f"allantools vs numba: {r_at_nb:.2e}"
            )


if __name__ == "__main__":
    main()

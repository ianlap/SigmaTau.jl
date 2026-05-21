#!/usr/bin/env python3
"""Side-by-side scaling comparison: SigmaTau (full API, with CI/bias) vs
allantools (bare kernel). Reads the two JSONs produced by `scaling.jl` and
`scaling_allantools.py`.

Caveat printed at the top: the two libraries do different work in their
"default user call". SigmaTau auto-computes confidence intervals and bias
correction; allantools doesn't. The comparison is "wall-time a user pays
for one call", not "kernel-vs-kernel".

Usage:

    python benchmarks/bench/render_scaling.py
    python benchmarks/bench/render_scaling.py --predict-at 10000000
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def _fmt_t(s: float | None) -> str:
    if s is None:
        return "       —"
    if s < 1e-3:
        return f"{s*1e6:6.1f} μs"
    if s < 1.0:
        return f"{s*1e3:6.2f} ms"
    if s < 60:
        return f"{s:6.2f}  s"
    if s < 3600:
        return f"{s/60:6.2f}  m"
    return f"{s/3600:6.2f}  h"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--sigmatau", type=Path,
                   default=Path(__file__).parent / "scaling_sigmatau.json")
    p.add_argument("--allantools", type=Path,
                   default=Path(__file__).parent / "scaling_allantools.json")
    p.add_argument("--predict-at", type=int, default=None,
                   help="Override the predicted-N from the JSON metadata.")
    args = p.parse_args()

    sj = json.loads(args.sigmatau.read_text())
    aj = json.loads(args.allantools.read_text())

    predict_at = args.predict_at or sj["meta"].get("predict_at", 1_000_000)

    print("Scaling comparison — SigmaTau (full API) vs allantools (bare kernel)")
    print("─" * 78)
    print(f"  SigmaTau:   {sj['meta']['kwargs']}")
    print(f"              Julia {sj['meta']['julia_version']}, threads = {sj['meta']['threads']}")
    print(f"  allantools: {aj['meta']['kwargs']}")
    print(f"              Python {aj['meta']['python_version']}, "
          f"allantools {aj['meta']['allantools_version']}, "
          f"threads = {aj['meta']['threads']}")
    print()
    print("Note: SigmaTau measurement includes auto-computed confidence intervals,")
    print("noise identification, and bias correction. allantools doesn't auto-compute")
    print("those — the user calls `at.ci.*` separately. Comparison is user-facing")
    print("wall-time per call, not pure-kernel work.")
    print()

    # Build a kernel -> (a, b, r2) map for each library
    s_fits = {f["kernel"]: f for f in sj["fits"]}
    a_fits = {f["kernel"]: f for f in aj["fits"]}

    kernels = [f["kernel"] for f in sj["fits"]]   # SigmaTau ordering

    print(f"{'kernel':<9} │ {'ST fit':<26}    {'AT fit':<26} │ "
          f"{'ST T(N=' + str(predict_at) + ')':>14}   {'AT T(N=' + str(predict_at) + ')':>14}   speedup")
    print("─" * 130)

    for k in kernels:
        sf = s_fits.get(k, {})
        af = a_fits.get(k, {})

        st_str = f"{sf['a']:.2e} · N^{sf['b']:.3f} (R²={sf['r2']:.3f})" if "a" in sf else "—"
        if af.get("implemented") is False:
            at_str = "(not in allantools)"
        elif "a" in af:
            at_str = f"{af['a']:.2e} · N^{af['b']:.3f} (R²={af['r2']:.3f})"
        else:
            at_str = "(too few points)"

        st_t = sf["a"] * predict_at ** sf["b"] if "a" in sf else None
        at_t = af["a"] * predict_at ** af["b"] if "a" in af else None

        if st_t is not None and at_t is not None and st_t > 0:
            speedup = f"{at_t / st_t:>7.1f}×"
        else:
            speedup = "—"

        print(f"{k:<9} │ {st_str:<30}  {at_str:<30} │ "
              f"{_fmt_t(st_t):>14}   {_fmt_t(at_t):>14}   {speedup:>9}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Generate a deterministic set of synthetic phase-noise realization files
for the statistical bench. Each file is a single-column whitespace-separated
.txt of N standard-normal values; seed = realization index.

Usage:
    python gen_synth.py --N 25000 --reals 30 --out-dir validation/bench/synth
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--N", type=int, default=25000)
    p.add_argument("--reals", type=int, default=30)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--seed-offset", type=int, default=0,
                   help="Add this to the realization index for the seed.")
    args = p.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    width = max(4, len(str(args.reals - 1)))
    for i in range(args.reals):
        rng = np.random.default_rng(args.seed_offset + i)
        x = rng.standard_normal(args.N).astype(np.float64)
        path = args.out_dir / f"realization_{i:0{width}d}.txt"
        np.savetxt(path, x, fmt="%.17g")
    print(f"wrote {args.reals} files of N={args.N} to {args.out_dir}")


if __name__ == "__main__":
    raise SystemExit(main())

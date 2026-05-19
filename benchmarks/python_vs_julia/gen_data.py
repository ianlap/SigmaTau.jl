#!/usr/bin/env python3
"""Generate synthetic phase data for cross-language benchmarks.

A mix of white phase noise and a small flicker-like contribution. Same seed,
written to plain text so Julia and Python can both load it without depending
on each other's serialisation.
"""
import numpy as np
import sys, os

OUT = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(OUT, exist_ok=True)

rng = np.random.default_rng(20260519)

for N in (1000, 4000, 16000):
    # phase = cumulative white frequency + white phase, scaled to seconds
    y_wfm = rng.standard_normal(N) * 1e-12
    x_wpm = rng.standard_normal(N) * 1e-13
    x = np.cumsum(y_wfm) + x_wpm
    path = os.path.join(OUT, f"phase_N{N}.txt")
    np.savetxt(path, x, fmt="%.17g")
    print(f"wrote {path}  N={N}  range=({x.min():.3e}, {x.max():.3e})")

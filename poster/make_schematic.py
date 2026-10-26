"""F6: how MHTOTDEV is built from one 4m-sample subsequence (src/kernels.jl, _mhtotdev_greenhall)."""
import os, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_figures import CRIMSON, BLUE, AMBER, GRAY, INK, INK2, GRID, save, rcParams

def main():
    rng = np.random.default_rng(7)
    m = 3; Lp = 4*m
    # a short phase record with a slope and noise
    x = np.cumsum(rng.normal(0, 1.0, 40)) + 0.6*np.arange(40)
    n0 = 14; seg = x[n0:n0+Lp].copy()
    half = Lp//2
    slope = (seg[half:].mean() - seg[:half].mean()) / (Lp/2)
    det = seg - slope*np.arange(Lp)
    ext = np.concatenate((det[::-1], det, det[::-1]))

    fig, axes = plt.subplots(3, 1, figsize=(9.0, 7.2), gridspec_kw=dict(height_ratios=[1.1, 1, 1], hspace=0.5))
    # (a) the record and one window of length 4m
    ax = axes[0]
    ax.plot(np.arange(40), x, color=GRAY, lw=2)
    ax.plot(np.arange(n0, n0+Lp), seg, color=CRIMSON, lw=3)
    ax.axvspan(n0-0.5, n0+Lp-0.5, color=CRIMSON, alpha=0.08, lw=0)
    ax.text(n0+Lp/2-0.5, max(seg)+6, 'one subsequence, length $L = 4m$\n(slides over all $N - 4m + 1$ positions)', ha='center', va='bottom', fontsize=17, color=CRIMSON)
    ax.set_ylim(min(x)-2, max(x)+18); ax.set_xlim(-1, 40)
    ax.set_yticks([]); ax.set_xticks([]); ax.grid(False)
    ax.set_title('a  phase record $x_k$', loc='left', fontsize=19, color=INK2)
    # (b) half-mean slope removal then 3x time-reverse extension
    ax = axes[1]
    j = np.arange(3*Lp)
    ax.plot(j[:Lp], ext[:Lp], color=BLUE, lw=2.5, marker='o', ms=5)
    ax.plot(j[Lp:2*Lp], ext[Lp:2*Lp], color=CRIMSON, lw=3, marker='o', ms=5)
    ax.plot(j[2*Lp:], ext[2*Lp:], color=BLUE, lw=2.5, marker='o', ms=5)
    for b in (Lp-0.5, 2*Lp-0.5): ax.axvline(b, color=AXIS if False else GRAY, lw=1.2, ls=':')
    ax.text(Lp/2-0.5, max(ext)+1.0, 'reversed copy', ha='center', fontsize=16, color=BLUE)
    ax.text(1.5*Lp-0.5, max(ext)+1.0, 'window, slope removed', ha='center', fontsize=16, color=CRIMSON)
    ax.text(2.5*Lp-0.5, max(ext)+1.0, 'reversed copy', ha='center', fontsize=16, color=BLUE)
    ax.set_ylim(min(ext)-1.5, max(ext)+4.5); ax.set_xlim(-1, 3*Lp); ax.set_yticks([]); ax.set_xticks([]); ax.grid(False)
    ax.set_title('b  half-mean slope removal, then 3x time-reverse extension (length 3L)', loc='left', fontsize=19, color=INK2)
    # (c) third differences at lag m, averaged over m: the "modified Hadamard" step
    ax = axes[2]
    L3 = 3*Lp - 3*m
    d3 = np.array([ext[k] - 3*ext[k+m] + 3*ext[k+2*m] - ext[k+3*m] for k in range(L3)])
    ax.bar(np.arange(L3), d3, color=GRAY, width=0.7, alpha=0.5)
    k0 = 6
    ax.bar(np.arange(k0, k0+m), d3[k0:k0+m], color=AMBER, width=0.7)
    ax.axhline(0, color=INK2, lw=1); ax.set_xlim(-1, L3); ax.set_yticks([]); ax.grid(False)
    ax.set_xlabel(f'third differences $\\Delta^3_m x^*_j$ along the extension (gray); one block of $m$ = {m} averaged (amber)', fontsize=16)
    ax.set_xticks([])
    ax.set_title(r'c  modified Hadamard: $\Delta^3$ at lag $m$, averaged over $m$, squared, $\div\,6m^2\tau^2$', loc='left', fontsize=19, color=INK2)
    save(fig, 'F6_mhtot_construction')

if __name__ == '__main__':
    main()

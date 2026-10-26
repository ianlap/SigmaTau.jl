"""Poster figures for SigmaTau.jl. All data come from the repository:
  data/validation.json   <- docs/src/validation.md (SigmaTau, Stable32, allantools per tau)
  tools/artifacts/mhtotdev_mc_full.json  <- Monte Carlo fits
  bench/README.md        <- recorded timings (typed in below, verbatim)
  test/fixtures/validation/stable32gen.DAT  <- the 8192-point fixture (numpy port, stab_np.py)
Run: python3 make_figures.py   (writes figures/*.png and figures/*.pdf)
"""
import json, sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import stab_np as S

# ---- palette (validated with the dataviz validator, adjacent pairs, white surface) ----
CRIMSON = '#9E1B32'          # Alabama crimson: brand + emphasis series
BLUE    = '#2a78d6'
AMBER   = '#eda100'
GREEN   = '#008300'
VIOLET  = '#6250d6'
GRAY    = '#898781'          # context / de-emphasis
INK     = '#0b0b0b'; INK2 = '#52514e'; GRID = '#e1e0d9'; AXIS = '#c3c2b7'
SERIES  = [CRIMSON, BLUE, AMBER, GREEN, VIOLET]
MARKERS = ['o', 's', '^', 'D', 'v']

rcParams.update({
    'font.family': 'Carlito',          # Calibri-metric; the poster body font
    'font.size': 20, 'axes.titlesize': 22, 'axes.labelsize': 22,
    'xtick.labelsize': 18, 'ytick.labelsize': 18, 'legend.fontsize': 18,
    'axes.edgecolor': AXIS, 'axes.linewidth': 1.2, 'axes.labelcolor': INK,
    'xtick.color': INK2, 'ytick.color': INK2, 'text.color': INK,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 1.0,
    'axes.spines.top': False, 'axes.spines.right': False,
    'lines.linewidth': 2.5, 'lines.markersize': 9,
    'legend.frameon': False, 'figure.dpi': 100, 'savefig.dpi': 300,
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
})
OUT = os.path.join(HERE, 'figures'); os.makedirs(OUT, exist_ok=True)
def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(OUT, f'{name}.{ext}'), bbox_inches='tight', pad_inches=0.08)
    plt.close(fig); print('wrote', name)

V = json.load(open(os.path.join(HERE, 'data', 'validation.json')))
MC = json.load(open(os.path.join(ROOT, 'tools', 'artifacts', 'mhtotdev_mc_full.json')))

# =====================================================================
# F1  Deviation suite on the shipped validation fixture, with 68.3 % chi^2 intervals
# =====================================================================
def fig_suite():
    fig, ax = plt.subplots(figsize=(9.0, 6.6))
    order = [('adev', 'ADEV (overlapping)'), ('mdev', 'MDEV'), ('hdev', 'HDEV'), ('totdev', 'TOTDEV'), ('htotdev', 'HTOTDEV')]
    for k, (key, label) in enumerate(order):
        rows = V[key]
        tau = np.array([r['tau'] for r in rows]); dev = np.array([r['sigmatau'] for r in rows])
        lo = np.array([r['lo'] for r in rows]); hi = np.array([r['hi'] for r in rows])
        ax.errorbar(tau, dev, yerr=[dev - lo, hi - dev], fmt=MARKERS[k] + '-', color=SERIES[k],
                    capsize=4, capthick=1.5, elinewidth=1.5, label=label, zorder=5 - k * 0.1)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(r'averaging time $\tau$ (s)'); ax.set_ylabel(r'$\sigma_y(\tau)$')
    ax.set_xlim(0.7, 6000); ax.set_ylim(1.2e-3, 1.6)
    # noise regimes identified per tau (Stable32 column alpha; SigmaTau agrees except tau = 256 s)
    ax.text(3.2, 0.62, 'white PM\n(slope -1)', color=INK2, fontsize=18, ha='left', va='bottom')
    ax.text(1300, 0.013, 'flicker FM /\nrandom-walk FM', color=INK2, fontsize=18, ha='center', va='bottom')
    ax.legend(loc='lower left', ncol=1, handlelength=2.2)
    save(fig, 'F1_suite')

# =====================================================================
# F2  Validation: relative difference to Stable32, per tau, per estimator
# =====================================================================
def fig_validation():
    fig, ax = plt.subplots(figsize=(9.6, 5.6))
    names = [('adev', 'ADEV'), ('mdev', 'MDEV'), ('tdev', 'TDEV'), ('hdev', 'HDEV'), ('totdev', 'TOTDEV'), ('mtotdev', 'MTOTDEV'), ('htotdev', 'HTOTDEV')]
    rng = np.random.default_rng(1)
    for i, (key, label) in enumerate(names):
        rows = V[key]
        col = 'sigmatau_raw' if key == 'mtotdev' else 'sigmatau'
        for r in rows:
            rel = abs(r[col] - r['stable32']) / r['stable32']
            x = i + rng.uniform(-0.18, 0.18)
            if r['noise_flip']:
                ax.plot(x, max(rel, 2e-7), marker='o', mfc='white', mec=GRAY, mew=2, ms=10, ls='none', zorder=4)
            else:
                ax.plot(x, max(rel, 2e-7), marker='o', color=CRIMSON, ms=9, ls='none', alpha=0.85, zorder=5)
    ax.axhline(5e-5, color=INK2, lw=1.8, ls='--')
    ax.text(6.45, 6.5e-5, 'Stable32 prints 5 significant figures: floor 5 x 10$^{-5}$', ha='right', va='bottom', fontsize=16, color=INK2)
    ax.set_yscale('log'); ax.set_ylim(1.2e-7, 1.5e-2)
    ax.set_xticks(range(len(names))); ax.set_xticklabels([n for _, n in names], fontsize=16)
    ax.text(5, 2.5e-7, 'raw kernel', ha='center', fontsize=14, color=INK2)
    ax.set_ylabel(r'$|\sigma_{\rm SigmaTau} - \sigma_{\rm Stable32}|\,/\,\sigma_{\rm Stable32}$')
    ax.plot([], [], marker='o', color=CRIMSON, ls='none', ms=9, label='one point per $\\tau$ (1 s to 4096 s)')
    ax.plot([], [], marker='o', mfc='white', mec=GRAY, mew=2, ms=10, ls='none', label=r'$\tau$ = 256 s: noise-ID convention differs')
    ax.legend(loc='upper left', fontsize=16)
    save(fig, 'F2_validation')

# =====================================================================
# F3  Speed: bench/README.md recorded results (N = 25 000 synthetic record)
# =====================================================================
BENCH_12 = [('adev', 15e-6, 393e-6, 20), ('mdev', 56e-6, 1.13e-3, 20), ('hdev', 25e-6, 341e-6, 14), ('tdev', 51e-6, 991e-6, 19),
            ('totdev', 254e-6, 684e-6, 2.7), ('mtotdev', 30.3e-3, 107.3, 3543), ('htotdev', 31.5e-3, 114.3, 3632)]
BENCH_4  = {'adev': 11, 'mdev': 15, 'hdev': 11, 'tdev': 9, 'totdev': 74, 'mtotdev': 1641, 'htotdev': 1860}
def fig_speed():
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    names = [b[0] for b in BENCH_12]; sp12 = [b[3] for b in BENCH_12]; sp4 = [BENCH_4[n] for n in names]
    y = np.arange(len(names))[::-1]
    ax.barh(y, sp12, height=0.58, color=CRIMSON, label='12-core workstation', zorder=3)
    ax.plot(sp4, y, marker='D', color=GRAY, ls='none', ms=10, label='4-core cloud container, re-run', zorder=4)
    for yi, s in zip(y, sp12):
        txt = f'{s:,.0f}x' if s >= 10 else f'{s:g}x'
        ax.text(s * 1.25, yi, txt, va='center', ha='left', fontsize=19, color=INK)
    ax.set_xscale('log'); ax.set_xlim(1, 6e4)
    ax.set_yticks(y); ax.set_yticklabels([n.upper() for n in names])
    ax.set_xlabel('speed-up over allantools (same samples, same $\\tau$ grid, same statistic)')
    ax.grid(axis='y', visible=False)
    ax.legend(loc='upper right', bbox_to_anchor=(1.0, 1.0), fontsize=15, borderaxespad=0.2)
    save(fig, 'F3_speed')

# =====================================================================
# F4/F5  MHTOTDEV Monte Carlo calibration: bias factor and EDF models
# =====================================================================
NOISE = [('2', 'white PM'), ('1', 'flicker PM'), ('0', 'white FM'), ('-1', 'flicker FM'), ('-2', 'random-walk FM')]
def fig_mhtot_bias():
    fig, ax = plt.subplots(figsize=(5.0, 4.3))
    r = np.linspace(0.0, 0.25, 200)   # tau/T; the kernel's reach limit is tau = T/4
    for k, (a, label) in enumerate(NOISE):
        b = MC['bias_fit'][a]; B = b['b0'] + b['b1'] * r
        ax.plot(r, B, color=SERIES[k], label=f'{label} (α = {a})')
        ax.plot(r[::40], B[::40], marker=MARKERS[k], color=SERIES[k], ls='none', ms=8)
        if a in ('-1', '-2'):
            ax.text(0.258, {'-1':1.20,'-2':B[-1]}[a], f'{label}', color=SERIES[k], va='center', ha='left', fontsize=14)
    ax.axhline(1.0, color=INK2, lw=1.5, ls='--'); ax.text(0.004, 0.80, 'B = 1: unbiased', color=INK2, fontsize=14)
    ax.text(0.258, 1.0, 'PM noises,\nwhite FM', color=INK2, va='center', ha='left', fontsize=14)
    ax.set_xlim(0, 0.40); ax.set_ylim(0.65, 2.3); ax.set_xticks([0, 0.1, 0.2, 0.25]); ax.set_xticklabels(['0', '0.1', '0.2', '1/4'])
    ax.set_xlabel(r'$\tau\,/\,T$'); ax.set_ylabel(r'bias factor $B(\alpha, \tau/T)$')
    save(fig, 'F4_mhtot_bias')

def fig_mhtot_edf():
    fig, ax = plt.subplots(figsize=(5.0, 4.3))
    Ttau = np.logspace(np.log10(4), np.log10(2048), 300)
    for k, (a, label) in enumerate(NOISE):
        e = MC['edf_fit'][a]   # the shipped model is the Mod-Totvar form for every alpha (src/edf.jl _coeff_mhtot)
        edf = e['modtot']['b'] * Ttau - e['modtot']['c']
        ok = edf > 0
        ax.plot(Ttau[ok], edf[ok], color=SERIES[k], label=label)
        ax.plot(Ttau[ok][::50], edf[ok][::50], marker=MARKERS[k], color=SERIES[k], ls='none', ms=8)
    ax.axvspan(4, 6, color=GRAY, alpha=0.18, lw=0); ax.text(6.6, 2600, 'edge zone\n$T/\\tau$ < 5–6', fontsize=14, color=INK2, va='top', ha='left')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(3.5, 2500); ax.set_ylim(0.8, 4000)
    ax.set_xlabel(r'$T/\tau$'); ax.set_ylabel('EDF  ν')
    ax.legend(loc='lower right', fontsize=14, handlelength=1.6, labelspacing=0.3)
    save(fig, 'F5_mhtot_edf')

if __name__ == '__main__':
    which = sys.argv[1:] or ['suite', 'validation', 'speed', 'mhtot_bias', 'mhtot_edf']
    for w in which:
        globals()['fig_' + w]()

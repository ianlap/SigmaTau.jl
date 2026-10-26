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
    ax.text(3.2, 0.62, 'white / flicker PM\n(slope -1)', color=INK2, fontsize=18, ha='left', va='bottom')
    ax.text(1300, 0.013, 'flicker FM', color=INK2, fontsize=18, ha='center', va='bottom')
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


# =====================================================================
# F1 (package output)  Suite on the fixture from poster/data/suite_fixture.csv (julia_figdata.jl)
# =====================================================================
import csv
def _suite_rows():
    f = os.path.join(HERE, 'data', 'suite_fixture.csv')
    return list(csv.DictReader(open(f))) if os.path.exists(f) else None
NOISE_LABEL = {'WHPM': 'white PM', 'FLPM': 'flicker PM', 'WHFM': 'white FM', 'FLFM': 'flicker FM', 'RWFM': 'random-walk FM'}
def fig_suite_pkg():
    rows = _suite_rows()
    if rows is None:
        print('no suite_fixture.csv; run julia_figdata.jl'); return
    fig, ax = plt.subplots(figsize=(9.0, 6.4))
    order = [('adev', 'ADEV'), ('mdev', 'MDEV'), ('hdev', 'HDEV'), ('totdev', 'TOTDEV'), ('mhtotdev', 'MHTOTDEV (new)')]
    for k, (key, label) in enumerate(order):
        rr = [r for r in rows if r['dev'] == key]
        tau = np.array([float(r['tau']) for r in rr]); dev = np.array([float(r['sigma']) for r in rr])
        lo = np.array([float(r['lo']) for r in rr]); hi = np.array([float(r['hi']) for r in rr])
        ax.errorbar(tau, dev, yerr=[dev - lo, hi - dev], fmt=MARKERS[k] + '-', color=SERIES[k], capsize=4, capthick=1.5, elinewidth=1.5, label=label, zorder=5 - k * 0.1)
    # noise type identified per tau (ADEV row), as a strip of labels above the curves
    rr = [r for r in rows if r['dev'] == 'adev']
    taus = [float(r['tau']) for r in rr]; nts = [r['noise_type'] for r in rr]
    runs_ = []; start = 0
    for i in range(1, len(nts) + 1):
        if i == len(nts) or nts[i] != nts[start]:
            runs_.append((taus[start], taus[i - 1], nts[start])); start = i
    ytop = 2.6
    for (a, b, nt) in runs_:
        ax.text(np.sqrt(a * b), ytop, nt, ha='center', va='bottom', fontsize=15, color=INK2, rotation=0)
        ax.plot([a / 1.3, b * 1.3], [ytop * 0.93, ytop * 0.93], color=GRID, lw=2)
    ax.text(0.75, ytop * 1.9, 'noise type identified per τ (SP1065 codes):', ha='left', va='bottom', fontsize=15, color=INK2)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.7, 6000); ax.set_ylim(1.2e-3, 9)
    ax.set_xlabel(r'averaging time $\tau$ (s)'); ax.set_ylabel(r'$\sigma_y(\tau)$')
    ax.legend(loc='lower left', ncol=1, handlelength=2.2, fontsize=17)
    save(fig, 'F1_suite')

# =====================================================================
# F7  Dynamic Allan deviation map from poster/data/dynamic_map.csv (julia_figdata.jl)
# =====================================================================
from matplotlib.colors import LinearSegmentedColormap, LogNorm
CRIMSON_RAMP = LinearSegmentedColormap.from_list('crimson_seq', ['#FBEEF0', '#E9A3AE', '#C9566C', '#9E1B32', '#5A0F1D'])
def fig_dynamic():
    """Classic offset waterfall (SP1065 dynamic-stability look): one sigma_y(tau) curve per analysis
    window, each shifted up and right by its position in time; nearer windows occlude farther ones."""
    from matplotlib.patches import Polygon
    f = os.path.join(HERE, 'data', 'dynamic_map.csv')
    if not os.path.exists(f):
        print('no dynamic_map.csv; run julia_figdata.jl'); return
    m = list(csv.reader(open(f)))
    taus = np.array([float(h.split('=')[1]) for h in m[0][1:]])
    t = np.array([float(r[0]) for r in m[1:]]); Z = np.array([[float(v) for v in r[1:]] for r in m[1:]])
    xs = np.log10(taus); lz = np.log10(Z); n = len(t)
    dx, dy = 0.017, 0.021                      # per-window shift in log10(tau) and log10(sigma) units
    floor = np.nanmin(lz) - 0.12
    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    for k in range(n - 1, -1, -1):             # back to front
        x = xs + k * dx; y = lz[k] + k * dy
        verts = [(x[0], floor + k * dy)] + list(zip(x, y)) + [(x[-1], floor + k * dy)]
        face = '#FBEEF0' if t[k] < 4096 else '#E9A3AE'
        ax.add_patch(Polygon(verts, closed=True, facecolor=face, edgecolor=CRIMSON, linewidth=1.3, zorder=10 + (n - k)))
    # tau axis along the front curve
    for v, lab in [(1, '1'), (10, '10'), (100, '100')]:
        ax.plot([np.log10(v)] * 2, [floor - 0.06, floor], color=INK2, lw=1.2, zorder=100)
        ax.text(np.log10(v), floor - 0.09, lab, ha='center', va='top', fontsize=15, color=INK2)
    ax.text(xs.mean(), floor - 0.30, r'averaging time $\tau$ (s)', ha='center', va='top', fontsize=17, color=INK)
    # time axis along the offset direction (right-hand ends of the floors)
    for tt in (1000, 3000, 5000, 7000):
        k = int(np.argmin(np.abs(t - tt)))
        ax.plot([xs[-1] + k * dx, xs[-1] + k * dx + 0.06], [floor + k * dy] * 2, color=INK2, lw=1.2, zorder=100)
        ax.text(xs[-1] + k * dx + 0.09, floor + k * dy, f'{tt}', ha='left', va='center', fontsize=15, color=INK2)
    kmid = n // 2
    ax.text(xs[-1] + kmid * dx + 0.62, floor + kmid * dy - 0.16, 'window centre $t$ (s)', rotation=51, ha='center', va='center', fontsize=17, color=INK)
    # sigma axis at the front left
    for v in (-13, -12):
        ax.plot([xs[0] - 0.07, xs[0]], [v, v], color=INK2, lw=1.2, zorder=100)
        ax.text(xs[0] - 0.1, v, r'$10^{%d}$' % v, ha='right', va='center', fontsize=15, color=INK2)
    ax.text(xs[0] - 0.88, (floor + lz[0].max()) / 2, r'$\sigma_y(t,\tau)$', rotation=90, ha='center', va='center', fontsize=17, color=INK)
    ax.text(xs[0] - 0.62, lz.max() + (n - 1) * dy + 0.05, 'noise level changes at $t$ = 4096 s\n(lighter: before, darker: after)', fontsize=14, color=INK2, ha='left', va='bottom')
    ax.set_xlim(xs[0] - 1.02, xs[-1] + (n - 1) * dx + 0.95); ax.set_ylim(floor - 0.45, lz.max() + (n - 1) * dy + 0.42)
    ax.axis('off')
    save(fig, 'F7_dynamic')

if __name__ == '__main__':
    which = sys.argv[1:] or ['suite_pkg', 'validation', 'speed', 'mhtot_bias', 'mhtot_edf', 'dynamic']
    for w in which:
        globals()['fig_' + w]()

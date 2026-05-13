#!/usr/bin/env python3
"""Compare SigmaTau deviation outputs against Stable32 fixture.

Reads:
- sigmatau_deviations.csv  (SigmaTau output: dev, ci_lower, ci_upper, alpha)
- all_deviations.csv       (Stable32 fixture: sigma, sigma_min, sigma_max, alpha)

Produces:
- relative_diffs.csv  — relative differences for sigma, sigma_min, sigma_max
- plots/overlay_<deviation>.png — log-log overlay with CI ribbons
- plots/overlay_all.png — 2x4 grid of all 8 deviations
- summary text on stdout

Run from this directory:
    /home/ian/SigmaTau.jl/.bench-venv/bin/python compare.py
"""

import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
PLOTS_DIR = HERE / "plots"
PLOTS_DIR.mkdir(exist_ok=True)

DEVIATIONS = ["adev", "mdev", "tdev", "hdev",
              "totdev", "mtotdev", "ttotdev", "htotdev"]


def load_sigmatau(path):
    rows = {d: {} for d in DEVIATIONS}
    with open(path) as f:
        for r in csv.DictReader(f):
            d = r["deviation"]
            af = int(r["AF"])
            rows.setdefault(d, {})[af] = {
                "tau":       float(r["tau"]),
                "dev":       float(r["dev"]),
                "ci_lower":  float(r["ci_lower"]),
                "ci_upper":  float(r["ci_upper"]),
                "alpha":     int(r["alpha"]),
            }
    return rows


def load_stable32(path):
    rows = {d: {} for d in DEVIATIONS}
    with open(path) as f:
        for r in csv.DictReader(f):
            d = r["deviation"]
            af = int(r["AF"])
            rows.setdefault(d, {})[af] = {
                "tau":       float(r["tau"]),
                "dev":       float(r["sigma"]),
                "ci_lower":  float(r["sigma_min"]),
                "ci_upper":  float(r["sigma_max"]),
                "alpha":     int(r["alpha"]),
            }
    return rows


def write_diffs(st, s32, path):
    """Per-row relative diff = (sigmatau - stable32) / stable32, for dev,
    ci_lower, ci_upper. Stable32 reports to ~5 sig figs (rtol ~1e-4)."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "deviation", "AF", "tau",
            "s32_sigma", "st_sigma", "rel_diff_sigma",
            "s32_sigma_min", "st_sigma_min", "rel_diff_sigma_min",
            "s32_sigma_max", "st_sigma_max", "rel_diff_sigma_max",
            "s32_alpha", "st_alpha",
        ])
        for d in DEVIATIONS:
            common = sorted(set(st[d]) & set(s32[d]))
            for af in common:
                a = st[d][af]
                b = s32[d][af]
                w.writerow([
                    d, af, b["tau"],
                    b["dev"],      a["dev"],      (a["dev"]      - b["dev"])      / b["dev"]      if b["dev"]      else "",
                    b["ci_lower"], a["ci_lower"], (a["ci_lower"] - b["ci_lower"]) / b["ci_lower"] if b["ci_lower"] else "",
                    b["ci_upper"], a["ci_upper"], (a["ci_upper"] - b["ci_upper"]) / b["ci_upper"] if b["ci_upper"] else "",
                    b["alpha"],    a["alpha"],
                ])


def diff_stats(st, s32):
    """Summarize max|rel diff| per deviation for σ, σ_min, σ_max."""
    summary = {}
    for d in DEVIATIONS:
        common = sorted(set(st[d]) & set(s32[d]))
        rec = {"n": len(common)}
        for field, key in [("sigma", "dev"),
                           ("sigma_min", "ci_lower"),
                           ("sigma_max", "ci_upper")]:
            diffs = []
            for af in common:
                num = st[d][af][key] - s32[d][af][key]
                den = s32[d][af][key]
                if den != 0.0:
                    diffs.append(num / den)
            if diffs:
                a = np.abs(diffs)
                rec[field + "_max_abs"]  = float(np.max(a))
                rec[field + "_med_abs"]  = float(np.median(a))
        summary[d] = rec
    return summary


def plot_one(ax, d, st, s32, show_legend=False):
    common = sorted(set(st[d]) & set(s32[d]))
    if not common:
        ax.set_title(f"{d.upper()} (no data)")
        return
    tau = np.array([s32[d][af]["tau"]      for af in common])
    s32_sig = np.array([s32[d][af]["dev"]      for af in common])
    s32_lo  = np.array([s32[d][af]["ci_lower"] for af in common])
    s32_hi  = np.array([s32[d][af]["ci_upper"] for af in common])
    st_sig  = np.array([st[d][af]["dev"]       for af in common])
    st_lo   = np.array([st[d][af]["ci_lower"]  for af in common])
    st_hi   = np.array([st[d][af]["ci_upper"]  for af in common])

    # Stable32: filled band + line
    ax.fill_between(tau, s32_lo, s32_hi, color="#1f77b4", alpha=0.18,
                    label="Stable32 CI" if show_legend else None)
    ax.plot(tau, s32_sig, "-o", color="#1f77b4", markersize=4, lw=1.5,
            label="Stable32" if show_legend else None)
    # SigmaTau: dashed line + edge-only band
    ax.plot(tau, st_sig, "--s", color="#d62728", markersize=4, lw=1.5,
            label="SigmaTau" if show_legend else None)
    ax.plot(tau, st_lo, ":",  color="#d62728", lw=0.9)
    ax.plot(tau, st_hi, ":",  color="#d62728", lw=0.9,
            label="SigmaTau CI" if show_legend else None)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\tau$ (s)")
    ax.set_ylabel(r"$\sigma$")
    ax.set_title(d.upper())
    ax.grid(True, which="both", alpha=0.3)


def main():
    st  = load_sigmatau(HERE / "sigmatau_deviations.csv")
    s32 = load_stable32(HERE / "all_deviations.csv")

    write_diffs(st, s32, HERE / "relative_diffs.csv")
    summary = diff_stats(st, s32)

    # Per-deviation single-panel plots
    for d in DEVIATIONS:
        fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=140)
        plot_one(ax, d, st, s32, show_legend=True)
        ax.legend(fontsize=8, loc="best")
        fig.tight_layout()
        fig.savefig(PLOTS_DIR / f"overlay_{d}.png")
        plt.close(fig)

    # 2x4 grid
    fig, axes = plt.subplots(2, 4, figsize=(18, 9), dpi=130)
    for ax, d in zip(axes.flat, DEVIATIONS):
        plot_one(ax, d, st, s32, show_legend=False)
    # one shared legend
    handles = [
        plt.Line2D([0], [0], color="#1f77b4", marker="o", linestyle="-",  label="Stable32"),
        plt.Line2D([0], [0], color="#1f77b4", linestyle="-", alpha=0.4,  lw=8, label="Stable32 CI"),
        plt.Line2D([0], [0], color="#d62728", marker="s", linestyle="--", label="SigmaTau"),
        plt.Line2D([0], [0], color="#d62728", linestyle=":",              label="SigmaTau CI"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=10,
               bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("SigmaTau vs Stable32 — composite-noise fixture (N=32768)",
                 y=1.05, fontsize=13)
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "overlay_all.png", bbox_inches="tight")
    plt.close(fig)

    # Stdout summary
    print(f"\nRelative differences vs Stable32 (rtol target ≈ 1e-4 for 5 sig figs)\n")
    hdr = f"{'deviation':<10s} {'n':>3s}   {'σ max|Δ|':>10s} {'σ med|Δ|':>10s}   {'σmin max|Δ|':>12s} {'σmax max|Δ|':>12s}"
    print(hdr)
    print("-" * len(hdr))
    for d in DEVIATIONS:
        r = summary[d]
        print(f"{d:<10s} {r['n']:>3d}   "
              f"{r.get('sigma_max_abs', float('nan')):>10.2e} "
              f"{r.get('sigma_med_abs', float('nan')):>10.2e}   "
              f"{r.get('sigma_min_max_abs', float('nan')):>12.2e} "
              f"{r.get('sigma_max_max_abs', float('nan')):>12.2e}")

    print(f"\nFull table → {HERE / 'relative_diffs.csv'}")
    print(f"Plots      → {PLOTS_DIR}/  (overlay_*.png, overlay_all.png)")


if __name__ == "__main__":
    main()

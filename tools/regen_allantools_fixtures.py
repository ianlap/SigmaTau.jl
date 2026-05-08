#!/usr/bin/env python3
"""Regenerate allantools cross-validation fixture for SigmaTau.jl.

Reads the Stable32 phase fixture (`reference/validation/stable32gen.DAT`)
and the Stable32 row schema (`reference/validation/stable32out/
stable32_data_full.csv`), runs the matching allantools deviation for
each (Type, AF) pair, and writes the results as
`reference/validation/allantools_out/allantools_data_full.csv` with the
same columns Stable32 reports.

Run once when allantools is updated; the output is checked into the
repo so the Julia test does not depend on a runtime Python invocation.

Usage: python3 tools/regen_allantools_fixtures.py
"""
from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import allantools as at  # type: ignore[import-not-found]

REPO = Path(__file__).resolve().parents[1]
PHASE_FILE = REPO / "reference" / "validation" / "stable32gen.DAT"
STABLE32_CSV = REPO / "reference" / "validation" / "stable32out" / "stable32_data_full.csv"
OUT_DIR = REPO / "reference" / "validation" / "allantools_out"
OUT_CSV = OUT_DIR / "allantools_data_full.csv"

TAU0 = 1.0
HEADER_LINES = 10  # stable32gen.DAT has a 10-line header

# Stable32 "Type" column → allantools function. Functions not implemented
# by allantools (or not relevant to SigmaTau's exported surface) are
# skipped via the sentinel `None`.
TYPE_TO_AT = {
    "Overlapping Allan": at.oadev,
    "Modified Allan":    at.mdev,
    "Overlapping Hadamard": at.ohdev,
    "Time":              at.tdev,
    "Total":             at.totdev,
    "Hadamard Total":    at.htotdev,
    "Modified Total":    at.mtotdev,
    # not exercised against allantools (either not in our API or
    # not in allantools): "Allan", "Hadamard", "ThêoH", "Time Total"
}


def load_phase() -> np.ndarray:
    with open(PHASE_FILE) as f:
        lines = [ln.strip() for ln in f.readlines()]
    return np.asarray([float(s) for s in lines[HEADER_LINES:] if s], dtype=float)


def load_stable32_rows() -> list[dict[str, str]]:
    with open(STABLE32_CSV, newline="") as f:
        return list(csv.DictReader(f))


def _emit(msg: str) -> None:
    """Print a progress line and flush so it shows up immediately when
    piped to a file or captured by a wrapper."""
    print(msg, flush=True)


def main() -> int:
    t_start = time.monotonic()
    x = load_phase()
    if x.size != 8192:
        print(f"warning: expected 8192 phase samples, got {x.size}", file=sys.stderr)

    rows = load_stable32_rows()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Pre-filter: how many rows we'll actually run vs skip. Tells the
    # user up-front how long this will take.
    runnable = sum(1 for r in rows if r["Type"] in TYPE_TO_AT)
    skipped_count: dict[str, int] = {}
    for r in rows:
        if r["Type"] not in TYPE_TO_AT:
            skipped_count[r["Type"]] = skipped_count.get(r["Type"], 0) + 1

    _emit(f"allantools fixture regen — N={x.size}, {runnable} rows to compute, "
          f"{sum(skipped_count.values())} skipped")
    if skipped_count:
        for k, n in sorted(skipped_count.items()):
            _emit(f"  (skipping {n} rows of type {k!r})")
    _emit("")

    written = 0
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Type", "AF", "Tau", "N", "Sigma"])
        for row in rows:
            kind = row["Type"]
            fn = TYPE_TO_AT.get(kind)
            if fn is None:
                continue
            m = int(row["AF"])
            tau = m * TAU0
            written += 1
            t0 = time.monotonic()
            (_taus, devs, _errs, ns) = fn(x, rate=1.0/TAU0, data_type="phase",
                                          taus=np.asarray([tau]))
            elapsed = time.monotonic() - t0

            if devs.size == 0 or not np.isfinite(devs[0]):
                # allantools refuses some short-tau / long-m combinations
                # for total-family kernels — record as NaN and let the
                # Julia side decide whether to assert.
                w.writerow([kind, m, f"{tau:.6e}", int(ns[0]) if ns.size else 0, "nan"])
                _emit(f"[{written:>3}/{runnable}]  {kind:<22} m={m:<6} → nan ({elapsed:.2f}s)")
            else:
                w.writerow([kind, m, f"{tau:.6e}", int(ns[0]), f"{devs[0]:.6e}"])
                _emit(f"[{written:>3}/{runnable}]  {kind:<22} m={m:<6} → "
                      f"σ={devs[0]:.4e}  ({elapsed:.2f}s)")
            f.flush()  # so the CSV is recoverable if the run is interrupted

    total = time.monotonic() - t_start
    _emit("")
    _emit(f"wrote {written} rows to {OUT_CSV.relative_to(REPO)} in {total:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

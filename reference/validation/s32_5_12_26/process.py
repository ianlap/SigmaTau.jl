#!/usr/bin/env python3
"""Parse the Stable32 dumps in this directory into clean CSVs.

Inputs:
- output_notfiltered.txt : multi-deviation stability run (8 types x tau ladder)
- s32_sinlgetau          : single-tau (tau=1000, AF=1000) detail run with EDF/CI

Outputs:
- all_deviations.csv  : long format, one row per (deviation, AF, tau)
- singletau_details.csv : one row per deviation at tau=1000

Run from this directory:  python3 process.py
"""

import csv
import re
from pathlib import Path

HERE = Path(__file__).parent

# Stable32 "Sigma Type" string -> canonical short name
SIGMA_TYPE_MAP = {
    "Overlapping Allan":    "adev",
    "Modified Allan":       "mdev",
    "Time":                 "tdev",
    "Overlapping Hadamard": "hdev",
    "Total":                "totdev",
    "Modified Total":       "mtotdev",
    "Time Total":           "ttotdev",
    "Hadamard Total":       "htotdev",
}

VARIANCE_TYPE_MAP = {
    "Overlap Allan": "adev",
    "Mod Allan":     "mdev",
    "Time":          "tdev",
    "Overlap Had":   "hdev",
    "Total":         "totdev",
    "Mod Total":     "mtotdev",
    "Time Total":    "ttotdev",
    "Had Total":     "htotdev",
}

# Manual overrides for fields that the Stable32 sigma-dialog dump misreports
# or omits. Keyed by deviation; values overlay parsed fields.
# - hdev: Stable32 wrote "Hadamard Dev=0.000000e+00" but a fresh sigma-app
#   read returned 5.905421e-02. Use the corrected value.
SINGLETAU_OVERRIDES = {
    "hdev": {"dev_value": 5.905421e-02},
}


def parse_stability_dump(path: Path):
    """Yield dicts {deviation, AF, tau, n, alpha, sigma_min, sigma, sigma_max}."""
    rows = []
    current = None
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("Sigma Type:"):
            kind = s.split(":", 1)[1].strip()
            current = SIGMA_TYPE_MAP.get(kind)
            if current is None:
                raise ValueError(f"Unknown Sigma Type: {kind!r}")
            continue
        # Data rows: leading integer AF then 6 more numbers
        # Robust parse: a row with 7 whitespace-separated tokens, first 4 ints/floats
        toks = s.split()
        if len(toks) != 7:
            continue
        try:
            af = int(toks[0])
            tau = float(toks[1])
            n = int(toks[2])
            alpha = int(toks[3])
            smin = float(toks[4])
            sig = float(toks[5])
            smax = float(toks[6])
        except ValueError:
            continue
        if current is None:
            continue
        rows.append({
            "deviation": current,
            "AF":        af,
            "tau":       tau,
            "n":         n,
            "alpha":     alpha,
            "sigma_min": smin,
            "sigma":     sig,
            "sigma_max": smax,
        })
    return rows


# Single-tau parser: blocks start with "SIGMA FOR FILE:" and contain "key=value" lines.
SINGLETAU_BLOCK_RE = re.compile(r"SIGMA FOR FILE:", re.IGNORECASE)


def parse_singletau(path: Path):
    text = path.read_text()
    # Split into blocks at each "SIGMA FOR FILE:" header. Skip the leading prefix.
    parts = SINGLETAU_BLOCK_RE.split(text)[1:]
    rows = []
    for part in parts:
        fields = {}
        for raw in part.splitlines():
            line = raw.strip()
            if not line:
                continue
            # Match "Key=value" or "Key: value"
            m = re.match(r"([A-Za-z #][A-Za-z0-9 #/_]*?)\s*[=:]\s*(.+)$", line)
            if not m:
                continue
            key = m.group(1).strip()
            val = m.group(2).strip()
            fields[key] = val
        vt = fields.get("Variance Type")
        if vt is None:
            continue
        deviation = VARIANCE_TYPE_MAP.get(vt)
        if deviation is None:
            # buffer-leak garbage block, skip
            continue

        def fget(*keys, cast=float, default=None):
            for k in keys:
                if k in fields:
                    try:
                        return cast(fields[k])
                    except ValueError:
                        return default
            return default

        # The per-type deviation value lives under different field names.
        dev_value = fget(
            "Mod Sigma", "Time Sigma", "Hadamard Dev", "Total Dev",
            "ModTotdev", "Mod Totdev",
        )
        # For overlap-Allan, the reported "Sigma" *is* the deviation.
        if dev_value is None and deviation == "adev":
            dev_value = fget("Sigma")
        # Time Total / Had Total dumps in this file lack an explicit dev field.
        # We pull "Max R" / "Min R" but leave dev_value blank for those.

        record = {
            "deviation":     deviation,
            "AF":            fget("Avg Factor", cast=int),
            "tau":           fget("Tau"),
            "sigma_input":   fget("Sigma"),
            "dev_value":     dev_value,
            "n":             fget("# Analysis Points", cast=int),
            "edf":           fget("Chi Square DF"),
            "sigma_min":     fget("Min ADEV", "Min MDEV", "Min TDEV",
                                  "Min TOTDEV", "Min Mod Totdev", "Min R"),
            "sigma_max":     fget("Max ADEV", "Max MDEV", "Max TDEV",
                                  "Max TOTDEV", "Max Mod Totdev", "Max R"),
            "B1_ratio":      fget("B1 Ratio"),
            "Rn_ratio":      fget("Rn Ratio"),
            "BW_factor":     fget("BW Factor"),
            "noise_type":    fields.get("Noise Type"),
            "alpha":         fget("Alpha", cast=int),
            "mu":            fget("Mu", cast=int),
            "confidence":    fget("Confidence Factor"),
        }
        record.update(SINGLETAU_OVERRIDES.get(deviation, {}))
        rows.append(record)
    return rows


def write_csv(path: Path, rows, fieldnames):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    stab_rows = parse_stability_dump(HERE / "output_notfiltered.txt")
    st_rows = parse_singletau(HERE / "s32_sinlgetau")

    write_csv(
        HERE / "all_deviations.csv",
        stab_rows,
        ["deviation", "AF", "tau", "n", "alpha",
         "sigma_min", "sigma", "sigma_max"],
    )
    write_csv(
        HERE / "singletau_details.csv",
        st_rows,
        ["deviation", "AF", "tau", "sigma_input", "dev_value",
         "n", "edf", "sigma_min", "sigma_max",
         "B1_ratio", "Rn_ratio", "BW_factor",
         "noise_type", "alpha", "mu", "confidence"],
    )

    print(f"all_deviations.csv: {len(stab_rows)} rows")
    print(f"singletau_details.csv: {len(st_rows)} rows")


if __name__ == "__main__":
    main()

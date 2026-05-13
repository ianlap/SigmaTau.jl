# s32_5_12_26 — Stable32 fixture for composite-noise validation

Stable32 output dumps for the N=32768 composite-noise reference dataset
(`s32_5_12_freq.DAT`), captured on 2026-05-12.

## Files

| File | Source | Description |
|------|--------|-------------|
| `output_notfiltered.txt` | Stable32 stability run | Raw multi-deviation dump, 8 deviation types × τ ladder |
| `s32_sinlgetau`          | Stable32 sigma function | Single-τ (AF=1000, τ=1000) EDF/CI details |
| `process.py`             | parser | Idempotent script to regenerate the CSVs below |
| `all_deviations.csv`     | derived | Long format: one row per (deviation, AF, τ) |
| `singletau_details.csv`  | derived | One row per deviation at τ=1000 with EDF / CI / noise-ID |
| `compute_sigmatau.jl`    | driver | Runs SigmaTau on the phase record and writes `sigmatau_deviations.csv` |
| `compare.py`             | comparator | Joins SigmaTau and Stable32 CSVs, writes `relative_diffs.csv`, renders overlay plots into `plots/` |

Regenerate the CSVs with `python3 process.py` (Stable32 dumps → CSV) and
`julia --project=. compute_sigmatau.jl` (SigmaTau values), then
`python3 compare.py` (relative diffs + plots).

The `plots/` directory is gitignored — `compare.py` writes overlay PNGs
there for inspection but they're regenerable and not part of the
fixture.

## Coverage

Eight deviations are present (matches Stable32's complete deviation set):
`adev`, `mdev`, `tdev`, `hdev`, `totdev`, `mtotdev`, `ttotdev`, `htotdev`.

**Not in this fixture (Stable32 does not implement):** `mhdev`, `htdev`, `mhtotdev` — these are SigmaTau-only and must be validated by other means (legacy MATLAB parity / internal consistency).

## τ ladder

- ADEV / MDEV / TDEV / HDEV: AF ∈ {1, 2, 4, 10, 20, 40, 100, 200, 400, 1000, 2000, 4000}
- Total family (TOTDEV / MTOTDEV / TTOTDEV / HTOTDEV): extends one octave further with AF=10000

## Composite-noise design (input recap)

The data was generated with σ(τ=1) values chosen to give each noise type its own decade of τ:

| Noise | Target σ(τ=1) | Expected dominant region (τ) |
|-------|---------------|------------------------------|
| WPM/FPM | 3e-10 | 1 – 10 |
| WFM   | 1e-10 | 10 – 100 |
| FFM   | 1e-11 | 100 – 1000 |
| RWFM  | 3e-13 | 1000 – 10000 |

Stable32's B1-based noise ID in `all_deviations.csv` (column `alpha`) tracks this almost exactly: FPM/WPM at τ ≤ 4, WFM at τ=10–40, FFM at τ=100–400, RWFM at τ=1000.

## Known Stable32 quirks (handle on the SigmaTau side)

1. **HDEV single-τ block reports `Hadamard Dev=0.000000e+00`** — Stable32 UI bug in the single-tau dialog. A fresh sigma-app read returned the correct `5.905421e-02` (matches `all_deviations.csv` `hdev` row AF=1000). `process.py` overrides this value in `SINGLETAU_OVERRIDES`; if you re-dump the file, drop the override or update it.
2. **TTOTDEV / HTOTDEV single-τ blocks omit the explicit deviation value field** (only `Max R`/`Min R` were captured). The parser leaves `dev_value` blank for those two rows; the true values are in `all_deviations.csv` (ttotdev@1000 = 31.223, htotdev@1000 = 0.061872).
3. **Noise-ID flips back to WFM (α=0) at τ ≥ 2000 for several deviations** despite RWFM being the truly dominant component. Expected: B1 / Greenhall–Riley noise ID becomes unreliable near the long-τ edge where EDF drops below ~10. Not a fixture error, just a Stable32 behavior to be aware of when cross-checking SigmaTau's noise-ID output.
4. **`sigma_input` in `singletau_details.csv` is the σ value Stable32 used to compute the χ² CI**, not necessarily the deviation value itself. For ADEV it equals the deviation; for other types it's a separate user-supplied number. EDF (`edf` column) is the figure to cross-check against SigmaTau's `_*_edf` routines.

## Cross-validation tolerance

Stable32 prints to 5 significant figures → use `rtol ≈ 1e-4` when comparing SigmaTau output against `all_deviations.csv`.

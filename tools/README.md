# tools/

`mc_mhtotdev.jl` is the Monte Carlo harness that measured the MHTOTDEV bias
factor `B(α, τ/T)` and EDF coefficients shipped in `src/edf.jl`. MHTOTDEV has
no published coefficient table, so this is the provenance for those numbers;
the method is described in `docs/src/mhtotdev.md`.

```bash
julia --project=. --threads=auto tools/mc_mhtotdev.jl          # full sweep (hours)
julia --project=. --threads=auto tools/mc_mhtotdev.jl quick    # laptop check (minutes)
```

`artifacts/mhtotdev_mc_full.json` is the full-sweep fit summary (seed, grid,
R = 3000 realizations per cell, per-α coefficients and R²) that the shipped
tables were taken from.

## validation/

The reference cross-check suite behind `test/validation/`:

```bash
julia --project=. tools/validation/make_noise_fixtures.jl   # power-law records (once; committed)
pip install allantools==2024.6 scipy numpy
python3 tools/validation/gen_allantools.py                   # allantools / scipy tables (~5 min)
julia --project=. tools/validation/report.jl                 # REPORT.md: agreement and discrepancies
```

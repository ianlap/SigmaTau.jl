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

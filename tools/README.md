# tools/

Development and provenance tools. Nothing here is part of the public API; it is
tracked so the shipped coefficient tables and test fixtures can be reproduced.

## `mc_mhtotdev.jl`

Monte Carlo harness that measures the MHTOTDEV bias factor `B(α, τ/T)` and
equivalent degrees of freedom. MHTOTDEV is novel to SigmaTau, so no published
coefficient table exists; the values shipped in `src/edf.jl` (`_coeff_mhtot`
and `bias_correction(:mhtot, …)`) were fitted by this harness. Methodology is
documented in `docs/src/theory/mhtotdev_bias_edf.md`.

```bash
julia --project=. --threads=auto tools/mc_mhtotdev.jl            # full sweep (hours)
julia --project=. --threads=auto tools/mc_mhtotdev.jl quick      # laptop check (minutes)
```

Outputs a per-cell CSV and a JSON summary (metadata plus fitted coefficients)
under `artifacts/`.

## `artifacts/mhtotdev_mc_full.json`

Tracked provenance for the shipped MHTOTDEV coefficients: the authoritative
full-sweep fit summary (git d8a7ca7, master seed 20260531, R = 3000
realizations per cell, N = 1024…32768, all fits R² ≥ 0.998). Per-cell CSVs and
quick-run outputs are gitignored and stay local.

## `regen_allantools_fixtures.py`

Regenerates the allantools cross-validation fixtures under
`test/fixtures/validation/` from a Python environment with allantools
installed. Run only when intentionally refreshing the reference values; the
fixtures are otherwise frozen so test failures mean SigmaTau changed, not the
reference.

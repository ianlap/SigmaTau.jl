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

## `nnls_noise_id.jl`

Monte Carlo study of weighted-NNLS power-law noise identification across MDEV,
MHDEV, HDEV and MHTOTDEV. A composite record carrying WPM + WFM + RWFM + RRFM
at known amplitudes is fitted against an *empirical* basis — the ensemble-mean
σ̂²(τ) of synthesized unit noises, one basis per estimator — and the study
reports amplitude recovery, conditioning (κ₂), across-trial spread, a
parametric χ²(ν) bootstrap that isolates the EDF contribution, and per-noise-type
misidentification rates — each type scored over the τ where the analytic
crossover design (with family-correct response exponents) makes it dominant by a
clear margin, so an estimator cannot score well on a type that owns a single
uncontested octave. The exact solver it exercises is the shipped
`SigmaTau._nnls_exact`; a local Lawson–Hanson implementation cross-checks every
solve.

```bash
julia --project=. --threads=auto tools/nnls_noise_id.jl            # full study (hours)
julia --project=. --threads=auto tools/nnls_noise_id.jl quick      # laptop check (~2 min)
```

Hard gates (nonzero exit on failure): the pure-noise MHDEV slopes must be
(−3/2, −1/2, +1/2, +3/2) ± 0.15 for α = (+2, 0, −2, −4); the two solvers must
agree; and MHDEV / HDEV / MHTOTDEV must each recover a single-type record with
≥ 90 % of the fitted variance on the correct column. MDEV's failure on RRFM is
expected and documented, not gated: MVAR does not converge for α ≤ −3, so its
α = −4 column grows with the record length instead of following τ³.

Writes basis / trials / summary CSVs and a JSON (metadata, gate outcomes,
per-estimator metrics) under `artifacts/`; `paper/scripts/make_figures.jl`
reads them for Fig 4.

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

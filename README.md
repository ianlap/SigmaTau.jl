# SigmaTau.jl

[![CI](https://github.com/ianlap/SigmaTau.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ianlap/SigmaTau.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs — stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ianlap.github.io/SigmaTau.jl/stable/)
[![Docs — dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ianlap.github.io/SigmaTau.jl/dev/)
[![Julia ≥ 1.11](https://img.shields.io/badge/julia-%E2%89%A5%201.11-9558B2.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

SigmaTau.jl is a Julia package for time-and-frequency stability analysis. It
implements the Allan, Modified Allan, Hadamard, and Total deviation families,
MTIE, parabolic deviation, lag-1 ACF noise identification, Greenhall-Riley
EDF / chi-squared confidence intervals, and calibrated power-law noise
generation in one flat namespace.

The stability estimators are cross-checked against Stable32 and allantools.
Coming from Stable32? See the
[migration guide](https://ianlap.github.io/SigmaTau.jl/stable/migration_from_stable32/).

## Scope

- Deviation kernels are checked against Stable32, allantools, and a MATLAB
  reference where external references exist.
- Against allantools 2024.06 on a real N≈407 000 record, the modified-total
  kernels run **~3,800–4,200× faster**
  (`mtotdev` 0.47 s vs ~30 min; `htotdev` 0.46 s vs ~32 min), and the full
  seven-kernel sweep finishes ~4,000× sooner. The cheap kernels
  (`adev`/`mdev`/`hdev`/`tdev`) run 13–40× faster. See
  [Performance](https://ianlap.github.io/SigmaTau.jl/stable/performance/).
  *(Benchmarked against allantools; Stable32 is parity-verified for numerics,
  not speed-tested.)*
- `using SigmaTau` brings the shared types (`PhaseData`, `FrequencyData`,
  `StabilityResult`) and every public deviation, noise-ID, EDF/CI, MTIE, PDEV,
  and IO function into scope.

## Install

```julia
pkg> add https://github.com/ianlap/SigmaTau.jl
```

Or, working from a clone of this repo:

```julia
pkg> activate .
pkg> instantiate
```

## Quickstart

```julia
using SigmaTau

# Phase residuals (seconds), sampled every τ₀ = 1 s.
x = randn(10_000)
data = PhaseData(x, 1.0)

# One deviation, octave-spaced τ grid (Stable32's default spacing).
result = adev(data)

result.tau          # τ values (s)
result.dev          # σ_y(τ)
result.noise_type   # :WHPM / :FLPM / :WHFM / :FLFM / :RWFM, per τ
result.ci           # χ² (or Gaussian fallback) confidence bounds — one
                    # (lo, hi) tuple per τ: result.ci[1].lo / result.ci[1].hi;
                    # ci_lower(result) / ci_upper(result) give plain vectors
result.edf          # equivalent degrees of freedom (empty when ci=false)
result.neff         # number of analysis windows per τ (always populated)
```

Compute a whole suite in one call and index it by deviation:

```julia
suite = stability(data; devs=(:adev, :mdev, :hdev, :mhdev))
suite[:adev].dev    # the ADEV curve
keys(suite)         # which deviations are present, in order
```

`nch` separates per-clock noise from pairwise comparisons (the
three-cornered hat, generalized to N clocks) — see the three-cornered-hat
tutorial in the docs:

```julia
clock_a, clock_b, clock_c = nch(r_ab, r_bc, r_ca)
```

Read a Stable32-style data file and choose a τ grid explicitly or by spacing.
A single-column `.DAT` needs `time_col=0` (no time column) and an explicit
`tau0`:

```julia
p = read_phase("clock.DAT"; time_col=0, value_col=1, tau0=1.0)
adev(p, [1, 2, 4, 8, 16, 32, 64])   # explicit averaging factors
adev(p, Decade)                     # decade-spaced grid
```

Every public function also accepts `FrequencyData`:

```julia
y = randn(10_000) .* 1e-9
adev(FrequencyData(y, 1.0), [1, 2, 4])
```

### Available deviations

```
adev  mdev  tdev
hdev  mhdev  htdev
totdev  mtotdev  ttotdev  htotdev  mhtotdev
mtie  tierms  pdev
theo1  theoh
dadev  dhdev
```

### Plotting

`SigmaTau` ships a `RecipesBase` package extension. Loading any
`Plots`-compatible backend brings in a log-log τ–σ recipe with optional error
bars from the result's CI bounds, plus overlays for a whole suite:

```julia
using Plots, SigmaTau
plot(adev(data))                    # single deviation with CI error bars
plot(stability(data))              # overlay the default suite
```

### Streaming

`StreamingStability` computes `adev`/`mdev`/`hdev`/`mhdev` in real time as
phase samples arrive — O(1) per sample per τ, bounded memory (Dobrogowski &
Kasznia 2007):

```julia
acc = StreamingStability(:adev, 1.0, [1, 10, 100])
push!(acc, x_next)        # feed samples as they are measured
snapshot(acc)              # running StabilityResult, matches the batch kernel
```

## Documentation

Full docs — tutorials, theory, the migration guide, performance numbers, and
the API reference — are at
**[ianlap.github.io/SigmaTau.jl](https://ianlap.github.io/SigmaTau.jl/stable/)**.

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Tests are organised under `test/` (`types/`, `stab/`, `io/`,
`umbrella_smoke.jl`); the command above runs all sub-suites under a single
top-level `test/runtests.jl`.

## Contributing

Bug reports, new deviations, validation fixtures, and benchmarks are welcome.
See [`CONTRIBUTING.md`](CONTRIBUTING.md). The main references for the underlying
math are NIST SP1065 (Riley & Howe), Greenhall & Riley 2003, and IEEE 1139-2022.

## Citing

If you use SigmaTau.jl in published work, cite it as:

```bibtex
@software{lapinski-sigmatau,
  author  = {Lapinski, Ian},
  title   = {{SigmaTau.jl}: time-and-frequency stability analysis in {Julia}},
  url     = {https://github.com/ianlap/SigmaTau.jl},
  version = {0.5.0},
  year    = {2026},
}
```

Machine-readable citation metadata is in [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).

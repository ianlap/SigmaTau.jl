# SigmaTau.jl

[![CI](https://github.com/ianlap/SigmaTau.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ianlap/SigmaTau.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs — stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ianlap.github.io/SigmaTau.jl/stable/)
[![Docs — dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ianlap.github.io/SigmaTau.jl/dev/)
[![Julia ≥ 1.11](https://img.shields.io/badge/julia-%E2%89%A5%201.11-9558B2.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Frequency-stability analysis in Julia: Allan / Modified Allan /
Hadamard / Total deviations, lag-1 ACF noise identification,
Greenhall–Riley EDF / χ² confidence intervals, MTIE, parabolic
deviation, and a calibrated power-law noise generator.

`SigmaTau.jl` ships a single flat module covering the full stability
surface. Shared types (`PhaseData`, `FrequencyData`, `StabilityResult`)
plus every deviation, noise-ID, EDF/CI, MTIE, PDEV, and IO function
land in one namespace under `using SigmaTau`.

For clock state-space estimation (Kalman filter, PID steering, holdover
projection) see the sister package
[ClockEnsemble.jl](https://github.com/ianlap/ClockEnsemble.jl) (formerly
the `SigmaTau.Est` submodule, split out at v0.3.0).

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

# Phase residuals (in seconds), sampled every τ₀ = 1 s.
x = randn(10_000)
data = PhaseData(x, 1.0)

# Overlapping Allan deviation across a power-of-two τ grid.
result = adev(data, [1, 2, 4, 8, 16, 32, 64])

result.tau          # τ values (s)
result.dev          # σ_y(τ)
result.noise_type   # :WHPM / :FLPM / :WHFM / :FLFM / :RWFM
result.ci_lower     # χ² (or Gaussian fallback) confidence bounds
result.ci_upper
result.edf          # equivalent degrees of freedom (empty when calc_ci=false)
```

`adev` and the other public functions also accept `FrequencyData`:

```julia
y = randn(10_000) .* 1e-9
adev(FrequencyData(y, 1.0), [1, 2, 4])
```

### Available deviations

```
adev  mdev  tdev
hdev  mhdev  htdev
totdev  mtotdev  ttotdev  htotdev  mhtotdev
mtie  pdev
```

### Plotting

`SigmaTau` ships a `RecipesBase` package extension. Loading any
`Plots`-compatible backend brings in a log-log τ–σ recipe with optional
error bars from the result's CI bounds:

```julia
using Plots, SigmaTau
plot(adev(data, [1, 2, 4, 8, 16, 32, 64]))
```

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Tests are organised under `test/` (`types/`, `stab/`, `io/`,
`umbrella_smoke.jl`); the command above runs all sub-suites under a
single top-level `test/runtests.jl`.

## Status

See [`project_overview.md`](project_overview.md) for the per-component
status matrix and [`TODO.md`](TODO.md) for outstanding work. Notable
references for the underlying math: NIST SP1065 (Riley & Howe), Greenhall &
Riley 2003, IEEE 1139-2022.

## License

MIT. See [`LICENSE`](LICENSE).

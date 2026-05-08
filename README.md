# SigmaTau.jl

Frequency-stability analysis and clock-ensemble estimation in Julia.

`SigmaTau` is an umbrella package that re-exports three focused subpackages:

| Subpackage | Purpose |
|---|---|
| [`SigmaTauBase`](lib/SigmaTauBase) | Core types (`PhaseData`, `FrequencyData`, `StabilityResult`) shared across the rest of the suite. |
| [`SigmaTauStability`](lib/SigmaTauStability) | Stability deviations (Allan / Hadamard / Total families), noise identification, equivalent degrees of freedom, χ² confidence intervals. |
| [`SigmaTauEnsemble`](lib/SigmaTauEnsemble) | Clock state-space models and Kalman filtering for atomic-clock ensembles. AD-friendly out-of-place math via `StaticArrays`. |

## Install

```julia
pkg> add https://github.com/ianlap/SigmaTau.jl
```

Or, working from a clone of this repo:

```julia
pkg> activate .
pkg> instantiate
```

The repo is structured as a Julia 1.11 workspace; the root `Project.toml`
pulls each subpackage from its local `lib/…` path via `[sources]`.

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
totdev  mtotdev  htotdev  mhtotdev
```

### Plotting

`SigmaTau` ships a `RecipesBase` package extension. Loading any
`Plots`-compatible backend brings in a log-log τ–σ recipe with optional
error bars from the result's CI bounds:

```julia
using Plots, SigmaTau
plot(adev(data, [1, 2, 4, 8, 16, 32, 64]))
```

## Clock ensemble (Kalman filter)

```julia
using SigmaTau

model = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
est   = StandardKalmanFilter([0.0, 0.0, 0.0], 1e-12 * I(3))

for z in measurements
    predict!(est, model, model.tau)
    update!(est, model, z)
end
```

Pass `legacy_compat=true` to the constructor to reproduce the MATLAB-era
`safe_sqrt` diagonal clamping bit-for-bit.

## Running tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Each subpackage has its own `test/runtests.jl` and can be tested in
isolation, e.g. `julia --project=lib/SigmaTauStability -e 'using Pkg; Pkg.test()'`.

## Status

See [`project_overview.md`](project_overview.md) for the per-component
status matrix and [`TODO.md`](TODO.md) for outstanding work. Notable
references for the underlying math: NIST SP1065 (Riley & Howe), Greenhall &
Riley 2003, IEEE 1139-2022.

## License

MIT. See [`LICENSE`](LICENSE).

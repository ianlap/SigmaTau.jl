# Records and results

A record goes in and a result comes out. Every estimator takes a
[`PhaseData`](@ref) or a [`FrequencyData`](@ref) record, which carries the
sample interval with the samples, and returns a [`StabilityResult`](@ref).
[Conventions](../explanation/conventions.md) defines the quantities and
their units.

## Records

```@docs
AbstractTimingData
PhaseData
FrequencyData
```

## Results

| Type | Returned by |
|:-----|:------------|
| [`StabilityResult`](@ref) | Every estimator, [`deviation`](@ref) and [`snapshot`](@ref) |
| [`StabilitySuite`](@ref) | [`stability`](@ref) |
| [`SpectralResult`](@ref) | [`psd_y`](@ref), [`psd_x`](@ref) and [`ssb_phase_noise`](@ref) |
| [`DynamicStabilityResult`](@ref) | [`dadev`](@ref) and [`dhdev`](@ref) |

The other result types are documented with the functions that return them:
[`DriftFit`](@ref), [`PowerLawFit`](@ref), [`NoiseFit`](@ref),
[`CrossStabilityResult`](@ref), [`CrossSpectralResult`](@ref),
[`KalmanResult`](@ref), [`HoldoverResult`](@ref) and
[`TimescaleResult`](@ref).

```@docs
StabilityResult
ci_lower
ci_upper
StabilitySuite
SpectralResult
DynamicStabilityResult
```

## Results as tables

With [Tables.jl](https://tables.juliadata.org/stable/) loaded, a
`StabilityResult` and a `StabilitySuite` are row tables, with one row for
each deviation and averaging time:

| Column | Type | Content |
|:-------|:-----|:--------|
| `deviation_type` | `Symbol` | The estimator, such as `:adev` |
| `tau` | `Float64` | The averaging time, in seconds |
| `dev` | `Float64` | The deviation |
| `noise_type` | `Symbol` or `missing` | The noise type |
| `ci_lower`, `ci_upper` | `Float64` or `missing` | The confidence limits |
| `edf` | `Float64` or `missing` | The equivalent degrees of freedom |
| `neff` | `Int` | The number of terms averaged |

A column reads `missing` where the result holds no value, as the limits
and the EDF do after a call with `ci=false`.

```julia
using SigmaTau, Tables

rows = Tables.rows(adev(data))
cols = Tables.columntable(stability(data; devs=(:adev, :mdev)))
```

`Tables` is a weak dependency. Its methods load with it, and
`using SigmaTau` alone loads no table package.
[Plotting and tables](../guide/plotting.md#Tables) hands a result to
`DataFrames` and `CSV`.

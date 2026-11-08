# Types

Records go in, results come out. Every estimator takes a
[`PhaseData`](@ref) or [`FrequencyData`](@ref) and returns a
[`StabilityResult`](@ref); [`stability`](@ref) returns a
[`StabilitySuite`](@ref), the spectra a [`SpectralResult`](@ref) and the
dynamic deviations a [`DynamicStabilityResult`](@ref).

## Records

```@docs
AbstractTimingData
PhaseData
FrequencyData
```

## Results

```@docs
StabilityResult
ci_lower
ci_upper
StabilitySuite
SpectralResult
DynamicStabilityResult
```

## Tables.jl

`StabilityResult` and `StabilitySuite` implement the
[Tables.jl](https://tables.juliadata.org/stable/) row-table interface when
`Tables` is loaded:

```julia
using SigmaTau, Tables

rows = Tables.rows(adev(data; ci=true))
cols = Tables.columntable(stability(data; devs=(:adev, :mdev)))
```

The extension is optional. `Tables` is a weak dependency, so loading
`SigmaTau` alone does not load the wider tabular-data stack. Empty confidence
fields from `ci=false` are exposed as `missing` table values.

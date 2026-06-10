# Shared types

Core data types used across the package and re-exported from
`SigmaTau`.

## Types

```@docs
AbstractTimingData
PhaseData
FrequencyData
StabilityResult
StabilitySuite
SpectralResult
DynamicStabilityResult
```

## Accessors

```@docs
ci_lower
ci_upper
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

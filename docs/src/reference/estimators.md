# Estimators

Every estimator takes a [`PhaseData`](@ref) or [`FrequencyData`](@ref) record
and has three call forms:

```julia
adev(data)                    # octave τ grid up to the estimator's limit
adev(data, Decade)            # a TauMode grid; see Suites, τ grids and streaming
adev(data, [1, 10, 100])      # explicit averaging factors m, τ = m·τ₀
```

They return a [`StabilityResult`](@ref). The common keywords are

- `ci=true`: identify the noise type and compute the EDF and confidence
  interval at each τ. `ci=false` returns the deviation and `neff` only, and
  is faster.
- `confidence=DEFAULT_CONFIDENCE` (0.683): the interval's confidence level.
- `correct_bias=true` (total family and `theo1` only): divide out the
  estimator's known bias; see [Noise, EDF and confidence](../explanation/noise_and_confidence.md#Bias-correction).

[Choosing an estimator](../guide/choosing.md) compares them.

## Allan and Hadamard

```@docs
adev
mdev
tdev
hdev
mhdev
htdev
pdev
```

## Total

```@docs
totdev
mtotdev
ttotdev
htotdev
mhtotdev
```

## Thêo1 and ThêoH

```@docs
theo1
theoh
```

## Time interval error

```@docs
mtie
tierms
```

## Dynamic (time-resolved)

`dadev` and `dhdev` slide a window across the record and return σ(t, τ) as a
[`DynamicStabilityResult`](@ref). [Plotting](../guide/plotting.md#Dynamic-deviations)
shows how to draw the map.

```@docs
dadev
dhdev
```

## N-cornered hat

```@docs
nch
```

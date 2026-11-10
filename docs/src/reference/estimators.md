# Estimators

Every estimator takes a [`PhaseData`](@ref) or [`FrequencyData`](@ref) record
and has three call forms:

```julia
adev(data)                    # octave τ grid up to the estimator's limit
adev(data, Decade)            # a TauMode grid; see Suites, τ grids and streaming
adev(data, [1, 10, 100])      # explicit averaging factors m, τ = m·τ₀ (a range works too)
```

They return a [`StabilityResult`](@ref). The common keywords are

- `ci=true`: identify the noise type and compute the EDF and confidence
  interval at each τ. `ci=false` returns the deviation and `neff` only, and
  is faster (with `correct_bias=true` the total estimators still identify
  the noise for the bias correction).
- `confidence=DEFAULT_CONFIDENCE` (0.683): the interval's confidence level.
- `correct_bias=true` (total family and `theo1` only): divide out the
  estimator's known bias. `:published` keeps to the published tables,
  `:stable32` does what Stable32 does for that estimator, and `false`
  returns the raw statistic, as allantools does; see
  [Noise, EDF and confidence](../explanation/noise_and_confidence.md#Matching-a-reference).

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

## Telecom masks

The ITU-T recommendations specify a clock's or an interface's wander as
limits on MTIE and TDEV over the observation interval τ. `mask_check`
compares an [`mtie`](@ref) or [`tdev`](@ref) result with one; `sigmatau
analyze --mask` and the explorer's Stability page draw the same limits.
Measure the way the recommendations do: phase against a reference better
than the mask, sampled at 1/30 s or faster through a 10 Hz low-pass filter.

```@docs
TELECOM_MASKS
telecom_mask
mask_check
```

## Dynamic (time-resolved)

`dadev` and `dhdev` slide a window across the record and return σ(t, τ) as a
[`DynamicStabilityResult`](@ref). [Plotting](../guide/plotting.md#Dynamic-deviations)
shows how to draw the map.

```@docs
dadev
dhdev
```

## Clock ensembles

```@docs
ncornered_hat
groslambert
```

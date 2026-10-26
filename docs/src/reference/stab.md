# Stability functions

## Allan family deviations

```@docs
adev
mdev
hdev
mhdev
tdev
htdev
```

## Total family deviations

```@docs
totdev
mtotdev
ttotdev
htotdev
mhtotdev
```

## MTIE, TIE rms, and PDEV

```@docs
mtie
tierms
pdev
```

## Theo1 and ThêoH

```@docs
theo1
theoh
```

## Dynamic (time-resolved) deviations

`dadev` and `dhdev` slide a window across the record and return a σ(t, τ)
map as a [`DynamicStabilityResult`](@ref).

```@docs
dadev
dhdev
```

## Streaming accumulators

[`StreamingStability`](@ref) computes `adev`/`mdev`/`hdev`/`mhdev` as phase
samples arrive, at O(1) cost per sample per averaging factor; at every
sample count the streamed estimate equals the batch call on the samples
seen so far.

```julia
acc = StreamingStability(:adev, 1.0, [1, 10, 100])
push!(acc, x_next)        # one phase sample, O(length(m_values))
append!(acc, xs)          # a chunk, in order
snapshot(acc)              # StabilityResult over everything seen so far
nsamples(acc)             # samples streamed so far
```

```@docs
StreamingStability
Base.push!(::StreamingStability, ::Real)
Base.append!(::StreamingStability, ::AbstractVector{<:Real})
snapshot
nsamples
```

## Tau grids and suite API

```@docs
TauMode
tau_values
DEFAULT_DEVIATIONS
stability
```

The exported [`TauMode`](@ref) instances are `AllTaus`, `Octave`,
`HalfOctave`, `QuarterOctave`, `Decade`, and `HalfDecade`.

## N-cornered hat

```@docs
nch
```

## Noise identification

```@docs
identify_noise
noise_gen
```

## Spectral densities

```@docs
Sy
Sx
L
```

## EDF, bias, and confidence intervals

```@docs
calculate_edf
bias_correction
confidence_intervals
DEFAULT_CONFIDENCE
```

## I/O and preprocessing

```@docs
read_phase
read_frequency
detrend
fillgaps
find_outliers
remove_outliers
save
save_result
load_result
save_suite
load_suite
```

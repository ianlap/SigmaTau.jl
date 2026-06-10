# SigmaTau — Stability analysis

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

DADEV and DHDEV slide an analysis window across the record and return a
2-D σ(t, τ) map as a [`DynamicStabilityResult`](@ref) rather than a
`StabilityResult` curve. They are therefore not available through
[`stability`](@ref).

```@docs
dadev
dhdev
```

## Streaming accumulators

[`StreamingStability`](@ref) computes `adev`/`mdev`/`hdev`/`mhdev` in real
time as phase samples arrive, using the Dobrogowski–Kasznia running-sum
scheme [dobrogowski-2007-realtime-adev](@cite): O(1) work per new sample
per averaging factor, memory bounded by a ring buffer of at most
`4·maximum(m_values) + 1` samples. At every sample count the streamed
estimate matches the batch call on the samples seen so far. The
total-family and Thêo estimators cannot stream — their whole-record
extension/sampling schemes have no O(1) per-sample update.

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

## Advanced / research kernels

These are not exported, but remain reachable as `SigmaTau._adev_core` and
friends for benchmarking, parity validation against reference implementations,
and research use. The underscore prefix signals that the calling convention is
not stability-guaranteed across minor versions; prefer the wrappers above for
application code.

```@docs
SigmaTau._adev_core
SigmaTau._mdev_core
SigmaTau._tdev_core
SigmaTau._hdev_core
SigmaTau._mhdev_core
SigmaTau._totdev_core
SigmaTau._mtotdev_core
SigmaTau._htotdev_core
SigmaTau._mhtotdev_core
SigmaTau._mtie_core
SigmaTau._tierms_core
SigmaTau._pdev_core
SigmaTau._theo1_core
SigmaTau._theobr_ratio
SigmaTau._kernel_m_max
SigmaTau._default_m_values
```

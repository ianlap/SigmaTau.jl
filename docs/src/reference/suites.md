# Suites, τ grids and streaming

## Several estimators at once

```@docs
stability
DEFAULT_DEVIATIONS
```

## Estimators by name

[`estimators`](@ref) describes every estimator (name, label, family, axis,
whether it takes `correct_bias`, a one-line description), and
[`deviation`](@ref) runs one chosen by its name. The explorer and the
command-line tool use both, so they list the same estimators as the library.

```@docs
estimators
deviation
```

## τ grids

```@docs
TauMode
tau_mode
tau_values
```

The exported [`TauMode`](@ref) instances are `AllTaus`, `Octave`,
`HalfOctave`, `QuarterOctave`, `Decade` and `HalfDecade`.

## Streaming

[`StreamingStability`](@ref) computes `adev`, `mdev`, `hdev` or `mhdev` as
phase samples arrive, at O(1) cost per sample per averaging factor. At every
sample count the streamed estimate equals the batch call on the samples seen
so far.

```julia
acc = StreamingStability(:adev, 1.0, [1, 10, 100])
push!(acc, x_next)        # one phase sample
append!(acc, xs)          # a chunk, in order
snapshot(acc)             # StabilityResult over everything seen so far
nsamples(acc)             # samples streamed so far
```

```@docs
StreamingStability
Base.push!(::StreamingStability, ::Real)
Base.append!(::StreamingStability, ::AbstractVector{<:Real})
snapshot
nsamples
```

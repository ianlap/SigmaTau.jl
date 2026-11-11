# Suites, grids and streaming

## Several estimators at once

[`stability`](@ref) runs several estimators on one record with one set of
keywords and returns a [`StabilitySuite`](@ref), which is indexed by the
name of the estimator.

```@docs
stability
DEFAULT_DEVIATIONS
```

## Estimators by name

[`deviation`](@ref) runs the estimator of a given name, for a choice that
is made at run time. [`estimators`](@ref) lists the names with the label,
the family and the unit of each. The command-line program and the explorer
offer the same list.

```@docs
estimators
deviation
```

## τ grids

An estimator takes its averaging factors as a vector, as a range or as a
grid. The six grids are the instances of [`TauMode`](@ref): `AllTaus`,
`Octave`, `HalfOctave`, `QuarterOctave`, `Decade` and `HalfDecade`.
[Averaging times](../guide/choosing.md#Averaging-times) says which to
choose.

```@docs
TauMode
tau_mode
tau_values
```

## Streaming

[`StreamingStability`](@ref) keeps running sums for ADEV, MDEV, HDEV or
MHDEV and updates them with each phase sample. The deviation is available
while the measurement runs, and the record need not be kept.

```julia
acc = StreamingStability(:adev, 1.0, [1, 10, 100])
push!(acc, x_next)        # one phase sample
append!(acc, xs)          # several, in order
snapshot(acc)             # the StabilityResult of the samples so far
nsamples(acc)             # their number
```

```@docs
StreamingStability
Base.push!(::StreamingStability, ::Real)
Base.append!(::StreamingStability, ::AbstractVector{<:Real})
snapshot
nsamples
```

## Live measurement

[`LiveStability`](@ref) holds one accumulator for each deviation and takes
the readings as an instrument gives them: phase, fractional frequency, or
the frequency of a counter in hertz. [`read_reading`](@ref) returns the
next reading of a source, and `NaN` once its stream has ended:

```julia
using Sockets

counter = InstrumentSource(connect("192.168.1.20", 5025); query="READ?")
live = LiveStability([:adev, :mdev], 1.0, [1, 2, 4, 8, 16];
                     quantity=:hz, nominal=10e6)
while true
    reading = read_reading(counter)
    isnan(reading) && break
    push!(live, reading)
end
snapshot(live)            # one StabilityResult for each deviation
```

[`SimulatedSource`](@ref) stands in for the instrument while a setup is
tried out. The command line follows a measurement with
[`sigmatau live`](cli.md), and the explorer on its
[Live](../explorer/live.md) page.

```@docs
LiveStability
SimulatedSource
InstrumentSource
read_reading
```

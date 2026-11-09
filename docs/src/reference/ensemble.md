# Clock ensembles

A timescale starts from several clocks read against one reference at the
same instants: a counter with a channel per clock, or a switch that steps
through them. [`ClockComparisons`](@ref) holds those readings, and
[`simulate_ensemble`](@ref) makes a set with known true states, each clock
from its own [`ClockModel`](@ref), to check the algorithms against.

```julia
ms = [ClockModel(; tau0=1.0, q1=1e-24), ClockModel(; tau0=1.0, q1=4e-24, q2=1e-30)]
cc, truth = simulate_ensemble(ms, 86400)
adev(cc[2, 1])       # clock 2 − clock 1, free of the reference
```

```@docs
ClockComparisons
simulate_ensemble
```

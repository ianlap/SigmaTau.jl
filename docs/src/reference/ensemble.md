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

## Timescales

[`timescale`](@ref) averages the clocks' predicted phases into one scale.
`:ewfa` gives every clock the same weight, the baseline any weighting has to
beat; `:at1` learns each clock's prediction-error variance and weights by
its inverse, which on clocks of steady noise converges to the optimal
inverse-variance weights. [`timescale_phase`](@ref) returns the scale
against the reference, for the usual stability analysis.

```julia
r = timescale(cc; algorithm=:at1, tau_y=3600.0)
r.weights[end, :]                     # each clock's weight now
adev(timescale_phase(r))              # the scale's stability
```

```@docs
timescale
TimescaleResult
timescale_phase
```

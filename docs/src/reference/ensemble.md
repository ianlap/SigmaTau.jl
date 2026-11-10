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
beat. `:at1` is the NIST algorithm as Weiss & Weissert describe it
[weiss-weissert-1991-at2](@cite): it learns each clock's prediction-error
variance and weights by its inverse, which on clocks of steady noise
converges to the optimal inverse-variance weights. The equation numbers in
the docstring of [`timescale`](@ref) are the paper's. Three choices are
the paper's and can be swapped: the bias term `K` that makes up for each
clock being part of the scale (`bias=:paper`, or the exact
`bias=:exact`), no cap on the weights (`wmax`), and no outlier test
(`outlier`). Each clock's frequency constant comes from its own τ_min,
which `tau_min` sets and which is otherwise read from the record.
The paper leaves its starting values open, so these are this package's
choices: the τ_min estimate, the starting variance τ²σ_y²(τ) taken from the
record's first part, each clock's starting frequency (its fitted slope), and
a variance memory of 20 days (`memory`) capped at a quarter of the record.
[`timescale_phase`](@ref) returns the scale against the reference, for the
usual stability analysis.

```julia
r = timescale(cc; algorithm=:at1)
r.weights[end, :]                     # each clock's weight now
adev(timescale_phase(r))              # the scale's stability
```

AT1's weights are optimal for the interval between readings, since that is
what each clock's prediction error measures; at much longer averaging times,
with clocks whose long-term behaviour differs from their short-term, the scale
can be worse than the best clock.

`:greenhall` is Greenhall's Kalman timescale. One Kalman filter tracks every
clock's phase and frequency from the clock models; after each update the
phase part of its covariance is reduced about a weighted mean of the
corrected clocks, which keeps the frequency estimates and makes each step
of the scale as small as it can be. Unlike AT1 it needs a model per clock,
and it uses each clock's white PM level `R`. At `R = 0` it is the reduced Kalman timescale of the 2003 paper; at `R > 0` the filter is IEMx (2006). With reading noise the corrected clocks no longer agree and the implicit ensemble mean cannot be observed, so the scale is read at one clock, `pivot` (default 1): its reading less its estimated offset, which carries only that clock's reading noise. The realized scale carries that clock's reading noise as white PM, so at `R > 0` it is not quieter than a weighted readout (with three equal clocks, √3 noisier at short τ). Pivots agree exactly at `R = 0`. While the pivot has no reading the first clock read stands in, and a pivot that starts late steps the scale at its first reading.

```julia
timescale(cc; algorithm=:greenhall, models=ms)
```

`:jones_tryon` is the same filter without the covariance reduction, the natural Jones–Tryon scale [jones-tryon-1983-estimating-time](@cite): a clock's reading less its estimated phase. It has the same frequency estimates and, with noiseless readings, the same scale at every clock, but it is noisier at short τ, which the reduction removes. JT83's error detection, step corrections and drift state are not included.

```@docs
timescale
TimescaleResult
timescale_phase
```

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
clock_models
```

## Timescales

[`timescale`](@ref) averages the clocks' predicted phases into one scale.
`:ewfa` gives every clock the same weight (Kettering's equal-weights frequency
averaging [kettering-2023-clock-ensemble](@cite)), the baseline any weighting
has to beat. `:at1` is the NIST algorithm as Weiss & Weissert describe it
[weiss-weissert-1991-at2](@cite): it learns each clock's prediction-error
variance and weights by its inverse, which on clocks of steady noise
converges to the optimal inverse-variance weights. The equation numbers in
the docstring of [`timescale`](@ref) are the paper's. Three choices are
the paper's and can be swapped: the bias term `K` that makes up for each
clock being part of the scale (`bias=:paper`, or the exact
`bias=:exact`), no cap on the weights (`wmax`), and no outlier test
(`outlier`). Each clock's frequency constant comes from its own τ\_min,
which `tau_min` sets and which is otherwise read from the record.
The paper leaves its starting values open, so these are this package's
choices: the τ\_min estimate, the starting variance τ²σ\_y²(τ) taken from the
record's first part, each clock's starting frequency (its fitted slope), and
a variance memory of 20 days (`memory`) capped at a quarter of the record.
`:at2` keeps AT1's weights and replaces each clock's frequency filter by a
Kalman update that carries a frequency variance, which adapts as the clock's
prediction error does; in steady state it is AT1's filter, so the scales are
close. `steps=:best` or `:real` adds the paper's frequency-step detection (eqs 19–22):
a clock whose frequency jumps is left out for τ\_min, while its new frequency
is learned, and the scale rerun from the step; the threshold's units are
corrected from the printed equation (see the docstring), and off by default.
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

A clock without a reading keeps its prediction and has no part in the
average. When it comes back it joins as a new clock does: its phase is taken
from the scale at that sample and it weighs in from the next, so the scale
does not step at the return. The papers give no rule for a returning clock.

The clocks need not be read at every sample of the record. The interval of
the equations is the time since the last sample with a reading
([Weiss & Weissert (1991), §2.8](@cite weiss-weissert-1991-at2)), so clocks
read together at every tenth sample give the scale of those readings alone.
They do have to be read together: a clock read less often than the others
joins anew at each of its readings and never has a weight. For clocks on
different schedules use `:greenhall`, whose filter takes each reading when it
comes.

`:greenhall` is Greenhall's Kalman timescale. One Kalman filter tracks every
clock's phase and frequency from the clock models; after each update the
phase part of its covariance is reduced about a weighted mean of the
corrected clocks, which keeps the frequency estimates and makes each step
of the scale as small as it can be. Unlike AT1 it needs a model per clock,
and it uses each clock's reading noise `R`. At `R = 0` it is the reduced Kalman timescale of the 2003 paper; at `R > 0` the filter is IEMx (2006). With reading noise the corrected clocks no longer agree and the implicit ensemble mean cannot be observed, so the scale is read at one clock, `pivot` (default 1): its reading less its estimated offset, which carries only that clock's reading noise. The realized scale carries that clock's reading noise as white PM, so at `R > 0` it is not quieter than a weighted readout (with three equal clocks, √3 noisier at short τ). Pivots agree exactly at `R = 0`. While the pivot has no reading the first clock read stands in. A clock without a reading has no weight in the mean and runs on its prediction ([Greenhall (2006), §10](@cite greenhall-2006-measurement-noise)), so the scale does not step when it returns or when a clock joins late.

```julia
timescale(cc; algorithm=:greenhall, models=ms)
```

In simulation the scale's own stability is that of the implicit ensemble mean `Σ aᵢ (xᵢ − x̂ᵢ)`, the weighted mean of the corrected clocks: each clock's true phase less the filter's estimate of it, `r.xhat` ([`implicit_mean`](@ref)). `adev(timescale_phase(r))` adds the pivot clock's reading noise to it at short τ; at `R = 0` the readout is the mean with its sign reversed.

`:jones_tryon` is the same filter without the covariance reduction, the natural Jones–Tryon scale [jones-tryon-1983-estimating-time](@cite): a clock's reading less its estimated phase. It has the same frequency estimates and, with noiseless readings, the same scale at every clock, but it is slightly noisier at short τ (about 4% at τ = 1 on four white-FM clocks over ten seeds; single records up to 40% higher), which the reduction removes. Its reported weights are `NaN`. `R` is each clock's reading noise against the reference, so differences share the reference's noise; JT83 takes independent errors in each difference. JT83's error detection, step corrections and drift state are not included. Its scale steps when a clock comes back from a gap or joins late; `:greenhall`'s does not.

`:stein` is Stein's Kalman timescale [stein-2003-demystified](@cite), the KAS-2 algorithm [levine-2012-statistical-modeling, erickson-2011-symmetricom-time-scale](@cite). Every clock is differenced against the pivot and each difference has a Kalman filter; the pivot's own shocks follow from the constraint that the weighted shocks of all clocks sum to zero, with separate weights for the phase, frequency and aging, inversely proportional to the variance of the matching noise. `wmax` caps the weights as in KAS-2. A clock back from a gap has no part in the pivot's shock at its first reading, so the scale does not step at its return. The pivot should be read throughout: while it has no reading every clock runs on its prediction and the scale steps when the pivot returns.

```julia
timescale(cc; algorithm=:stein, models=ms, wmax=0.5)
```

```@docs
timescale
TimescaleResult
timescale_phase
implicit_mean
```

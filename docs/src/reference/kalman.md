# Kalman filter

[`kalman_filter`](@ref) tracks one clock's phase, frequency and, with a
drift state, drift through a phase record, on a [`ClockModel`](@ref)
from [`clock_model`](@ref). [`kalman_smooth`](@ref) refines every estimate
with the samples after it, [`holdover_error`](@ref) says how far the clock
can wander once the reference is gone, and [`innovation_check`](@ref) says
whether the model fits the record at all. [`simulate_clock`](@ref) makes a
record with a known truth to try all of this on.

```julia
m  = clock_model(powerlaw_fit(p; alphas=[2, 0, -2]))
kf = kalman_filter(p, m)
innovation_check(kf).ok                 # white, unit-variance innovations?
holdover_error(kf, [3600.0, 86400.0])   # 1σ time error after 1 h and 1 day
ks = kalman_smooth(kf)
```

Gaps (`NaN` samples) are predicted across without an update. Flicker noise
has no state in the default model, so on a record dominated by it the
innovation check fails. Either give the model flicker states, Markov
processes placed where the flicker matters
([The clock model](../explanation/clock_model.md)), or let
[`tune_q`](@ref) retune the white and random-walk FM levels for the best
phase prediction over the horizons you care about, starting from the fit.

```julia
mf = clock_model(powerlaw_fit(p; alphas=[2, 1, 0, -1, -2]); flicker=true)
kalman_filter(p, mf)                    # flicker FM and PM as states
t = tune_q(p, m; horizons=[10, 100, 1000])
kalman_filter(p, t.model)
```

```@docs
kalman_filter
KalmanResult
kalman_smooth
holdover_error
innovation_check
simulate_clock
prediction_error
tune_q
```

# Kalman filter

[`kalman_filter`](@ref) estimates the phase and the frequency of one clock
from its phase record, and the frequency drift where the model has a drift
state. It runs on a [`ClockModel`](@ref), whose states and noise levels
[The clock model](../explanation/clock_model.md) derives. The same filter
is the predictor of the [holdover test](holdover.md) and the estimator of
the [steering loop](steering.md), and the Kalman
[time scales](ensemble.md) run one over several clocks.

## From a record to a filter

A fit of the noise gives the model, and the model gives the filter.

```@example kalman
using SigmaTau, Random

# white FM of 1e-12 at 1 s on a flicker floor of 1e-13, read with 1 ps of white PM
h = Dict(2 => 8π^2 * 1e-24, 0 => 2e-24, -1 => 1e-26 / 2log(2))
p = noise_gen(PhaseData, 16_384, 1.0; h, rng = Xoshiro(1))
fit = powerlaw_fit(p)
```

```@example kalman
m = clock_model(fit)
```

The model keeps the flicker FM that the fit found, as Markov states over
the averaging times at which it matters ([Flicker states](@ref)).

```@example kalman
kf = kalman_filter(p, m)
kf.x[end, 1:2]                  # phase (s) and frequency at the last sample
```

Each row of `kf.x` is the state at one sample, and `kf.P[:, :, k]` is its
covariance at sample ``k``. A `NaN` sample is a gap, which the filter
predicts across without an update. The filter needs no starting values: it
takes them from a line through the first readings.

[`kalman_smooth`](@ref) runs back over the result, so that every estimate
uses the readings after it as well as those before. That suits a record
analyzed afterwards. The filter alone gives what a clock running in real
time can know.

```@example kalman
ks = kalman_smooth(kf)
# 1σ of the phase at sample 8000, smoothed over filtered
sqrt(ks.P[1, 1, 8000] / kf.P[1, 1, 8000])
```

## Testing the model

The record tests the model. Under the right model the innovations, each
reading less its prediction, are uncorrelated and have the variance that
the filter expects, and [`innovation_check`](@ref) tests both.

```@example kalman
innovation_check(kf)
```

`nis` is the mean square of the innovations over their expected variance,
and `acf` is their lag-1 autocorrelation. Here both are inside their
intervals. With `flicker = false` the model of the same fit leaves the
flicker FM out and fails:

```@example kalman
m0 = clock_model(fit; flicker = false)
innovation_check(kalman_filter(p, m0))
```

An `nis` above its interval says that the model has too little noise. On a
long record the intervals are narrow, and the test also fails models that
predict well ([What the states buy](@ref)). The holdover test asks the
more direct question, whether the prediction errors are as large as the
filter claims. [`holdover`](@ref) with a model predicts from many points of
the record and returns the RMS error beside the claim:

```@example kalman
tau = [10, 100, 1000]
a, b = holdover(p, m; tau), holdover(p, m0; tau)
[a.rms ./ a.claimed  b.rms ./ b.claimed]    # with the states, and without
```

With the states the errors are those claimed. Without them the filter
claims too little, and more so the longer the prediction.

## Prediction

[`holdover_error`](@ref) carries the last estimate forward and returns the
standard deviation of the time error that the model gives after each
holdover time.

```@example kalman
holdover_error(kf, [60.0, 3600.0])    # after a minute and after an hour, in seconds
```

[`prediction_error`](@ref) measures the same thing inside the record: the
RMS difference between the readings and the predictions made a number of
samples earlier.

## Tuning in place of states

[`tune_q`](@ref) is a remedy that adds no states. It keeps a model of
white and random-walk FM and searches their two levels for the smallest
prediction errors at chosen horizons, starting from a fit without flicker
terms.

```@example kalman
m3 = clock_model(powerlaw_fit(p; alphas = [2, 0, -2]))
t = tune_q(p, m3)
c, d = holdover(p, t.model; tau), holdover(p, m3; tau)
# RMS error over that with the states: tuned, and as fitted
[c.rms ./ a.rms  d.rms ./ a.rms]
```

Tuning trades the short predictions for the long ones. The tuned levels
are no longer those of the noise. `tune_q` sets their size so that the
error claimed one sample ahead is the error made (`t.scale` is the factor
it applied); farther ahead the covariance of the tuned filter, and with
it [`holdover_error`](@ref), is not to be relied on. In simulation the
predictions of tuned models were 5 to 29 % worse than those of models
with flicker states, and the errors they claimed from 10 to 1000 s were
1.2 to 2.8 times too small ([What the states buy](@ref)). `t.edge` says
whether the best pair lies on the edge of the grid, where a larger
`decades` is called for.

## Simulated clocks

[`simulate_clock`](@ref) draws a record from a model and returns the true
states behind it, so that an estimate can be compared with the truth. The
record follows the model exactly, Markov states included, which
[`noise_gen`](@ref) records of power-law noise do only approximately.

```@example kalman
truth_model = ClockModel(; tau0 = 1.0, R = 1e-22, q1 = 1e-24, q2 = 1e-30)
q, X = simulate_clock(truth_model, 4096; rng = Xoshiro(2))
kq = kalman_filter(q, truth_model)
err = kq.x[1001:end, 1] .- X[1001:end, 1]
# RMS phase error over the filter's 1σ
sqrt(sum(abs2, err) / length(err)) / sqrt(kq.P[1, 1, end])
```

```@docs
kalman_filter
KalmanResult
kalman_smooth
innovation_check
prediction_error
holdover_error
tune_q
simulate_clock
```

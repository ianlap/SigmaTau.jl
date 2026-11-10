# Holdover

When a clock loses its reference it runs on what it learned while it had
one. [`holdover`](@ref) measures how well that works on a phase record: it
learns the frequency, and with `drift = true` the drift, from a training
stretch, lets the clock run from the last reading of the stretch, and
compares with what the record did next,

```math
e(\tau) = x(t+\tau) - x(t) - \hat y\,\tau - \tfrac12 \hat D\,\tau^2 ,
```

the time interval error x(t + τ) − x(t) of ITU-T G.810 [ITU-G810](@cite)
less the learned frequency and drift. Every sample the record allows is a
start, the sliding test of Delporte et al.
[delporte-2000-extrapolation](@cite), so the statistics at each holdover
time come from the whole record: the RMS, the quantile at a chosen
confidence, the worst case and the mean.

```julia
h = holdover(p; train=86400, tau=[60, 3600, 14400], limit=1e-6)
h.rms              # RMS time error at each h.tau
h.quantile         # the error 95 % of the starts stay within
h.time_to_limit    # holdover time at which that reaches 1 µs
```

No clock model enters. The learned frequency is the slope of a
least-squares line through the stretch, or of a parabola at its last
sample, and what limits it depends on the noise. For white FM `q1` and
random-walk FM `q2` (see [`ClockModel`](@ref)) and a line over T seconds,

```math
\sigma_e^2(\tau) = q_1\tau\left(1 + \frac{6\tau}{5T}\right)
                 + q_2\left(\frac{\tau^3}{3} + \frac{13\,T\tau^2}{35}\right).
```

In each bracket the first term is the clock's own noise over τ and the
second the error of the learned frequency: a longer stretch averages white
FM down but lets the frequency walk away, so there is a best `train`,
`√(42 q1 / 13 q2)` for this pair, and trying a few shows it on a real
record. A mean error that grows as τ² under a line is drift the line
cannot follow.

With `anchor = false` the fit's own end value replaces the last reading,
which is the extrapolated fit Delporte et al. analyse; their closed-form
variances for white, flicker and random-walk FM are what the tests check
the function against. [`holdover_error`](@ref) gives the error a
[`ClockModel`](@ref) predicts from a Kalman filter's last estimate.

```@docs
holdover
HoldoverResult
```

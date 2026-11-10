# Holdover

When a clock loses its reference it runs on what it learned while it had
one. [`holdover`](@ref) measures how well that works on a phase record: it
learns the frequency, and with `drift = true` the drift, from a training
stretch, lets the clock run from the last reading of the stretch, and
compares with what the record did next,

```math
e(\tau) = x(t+\tau) - x(t) - \hat y\,\tau - \tfrac12 \hat D\,\tau^2 ,
```

the time interval error x(t + τ) − x(t) of ITU-T G.810, clause 4.5.14
[ITU-G810](@cite), less the learned frequency and drift. Every sample the record allows is a
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
which is the extrapolated fit Delporte et al. analyze.

## Predicted from a noise model

The same call on a noise model predicts what the record should show: a
[`PowerLawFit`](@ref), a [`ClockModel`](@ref), or a `Dict` of h_α such as a
data sheet's levels.

```julia
fit = powerlaw_fit(p)
m = holdover(p;   train=86400, limit=1e-6)      # measured
q = holdover(fit; train=86400, limit=1e-6)      # predicted
q.sigma            # against m.rms
sqrt.(q.budget)    # what each noise type adds, one column per q.alphas
```

The error is a weighted sum of phase samples, so its variance follows
exactly from the spectrum of the sampled record, the one
[`powerlaw_sigma`](@ref) uses. Summing the weights once turns them into
weights on the frequency samples, twice into weights on the frequency
steps; that leaves white FM and random-walk FM as plain sums of squares
and the flicker noises as a quadratic form with the kernel
1/(1 − 4k²), the covariance of a sequence whose spectrum is sin(πfτ₀).
The flicker noises are therefore predicted on the same footing as the
rest, with no state model standing in for them.

For long stretches the prediction gives the closed forms of Delporte et
al. (their eqs 11 to 13 for a parabola and 17 to 19 for a line, on white,
flicker and random-walk FM, all with `anchor = false`); the tests hold it
to them within 0.4 % at 2048 samples in the stretch, and hold the measured
test to the prediction on all five noise types. The first term of each
bracket above, q₁τ and q₂τ³/3, is τσ_y(τ), the error of the best possible
predictor on that noise [allan-1987-prediction](@cite): what remains is
the price of learning the frequency from a finite stretch.

The prediction holds noise only. A drift D that a line does not follow
adds the mean ½ D τ (τ + T). [`holdover_error`](@ref) gives the error a
[`ClockModel`](@ref) predicts from a Kalman filter's last estimate
instead of a fitted line.

## Holdover and the time deviation

The time variance is the mean square of an error of the same kind. With x̄
the phase averaged over τ, the definitions of Mod σ_y(τ) and σ_x(τ)
[ieee1139-1999-definitions](@cite) (eqs A.22 and A.23) give

```math
\sigma_x^2(\tau) = \tfrac16 \left\langle \left( \bar x(t+2\tau) - 2\,\bar x(t+\tau) + \bar x(t) \right)^2 \right\rangle ,
```

and the bracket is the next τ-average less the line through the two before
it. So √6 [`tdev`](@ref) is the RMS error of a holdover test that learns
the frequency from two τ-averages and predicts the third. The third
difference does the same with three averages and a parabola, drift
learned: its mean square is 6τ² Mod σ²_H(τ), which makes that error
√20 [`htdev`](@ref).

`grow = k` in place of `train` gives the measured test the matching
stretch, k times each holdover time, so that a longer holdover is learned
over a longer stretch:

```julia
g = holdover(p; grow=2)                         # a line over 2τ, then τ ahead
g.rms ./ (sqrt(6) .* tdev(p, Int.(g.tau ./ p.tau0); ci=false).dev)
```

The two differ in the readings: TDEV averages them over τ, the test
compares one reading with a least-squares line. At τ = τ₀ an average is a
reading and the two are equal. Above it the ratio depends on the noise
type alone, from the exact variance of the last section:

| Noise | line over 2τ, over √6 TDEV | parabola over 3τ, over √20 HTDEV |
|:---|:---|:---|
| white FM | √(8/5) = 1.26 | 1.22 |
| flicker FM | 1.31 | 1.31 |
| random-walk FM | √(452/231) = 1.40 | 1.42 |
| flicker PM | 1.3 to 1.8, rising with τ | 1.1 to 1.4 |
| white PM | rising as √τ | rising as √τ |

On FM noise the test therefore runs parallel to the scaled time deviation,
a quarter to two fifths above it, and √6 TDEV read as a holdover error
understates it by that factor. On white PM the averages smooth what a
single reading does not, and the two part. The ratios are this package's
calculation; the tests hold the two white FM and random-walk FM values to
their closed forms and the equality at τ₀ to rounding. With a fixed
`train` there is no such relation, because the stretch no longer scales
with τ.

```@docs
holdover
HoldoverResult
HoldoverPrediction
```

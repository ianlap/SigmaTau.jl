# Holdover

A clock that loses its reference runs on what it learned while it had one.
[`holdover`](@ref) measures how well that works on a phase record. It
learns the frequency from a training stretch, and the drift as well with
`drift = true`, lets the clock run from the last reading of the stretch,
at time ``t``, and compares with what the record did next:

```math
e(τ) = x(t + τ) - x(t) - \hat{y}\,τ - \tfrac{1}{2} \hat{D}\,τ^2 .
```

That is the time interval error ``x(t + τ) - x(t)`` of
[ITU-T G.810, clause 4.5.14](@cite ITU-G810), less the learned frequency
and drift. Every sample that the record allows is a start, as in the
sliding test of [Delporte et al. (2000)](@cite delporte-2000-extrapolation),
so the statistics at each holdover time come from the whole record.

```@example holdover
using SigmaTau, Random

# white FM of 1e-11 at 1 s that gives way to random-walk FM near 400 s
p = noise_gen(PhaseData, 65_536, 1.0; h = Dict(0 => 2e-22, -2 => 1e-28),
              rng = Xoshiro(1))
h = holdover(p; train = 512, tau = [16, 64, 256, 1024], limit = 1e-9)
```

The columns are the holdover time, the RMS error, the error that 95 % of
the starts stay within, the largest error and the number of starts. With
`limit`, the line above the table gives the holdover time at which the
95 % column reaches that error.

## The training stretch

No clock model enters the test. The learned frequency is the slope of a
least-squares line through the stretch, and what limits it depends on the
noise. For white FM of level ``q_1`` and random-walk FM of level ``q_2``
([Noise levels](@ref)) and a line over ``T`` seconds, the variance of the
error is

```math
σ_e^2(τ) = q_1 τ \left(1 + \frac{6τ}{5T}\right)
         + q_2 \left(\frac{τ^3}{3} + \frac{13\,T τ^2}{35}\right) .
```

In each bracket the first term is the noise of the clock over ``τ`` and
the second the error of the learned frequency. The first terms,
``q_1 τ`` and ``q_2 τ^3/3``, are ``τ^2 σ_y^2(τ)``, the mean square of the
optimum prediction on each noise
([Allan (1987), Table III](@cite allan-1987-prediction)). The second are
the price of learning the frequency from a finite stretch: a longer one
averages the white FM down and lets the frequency walk away, so there is a
best `train`, ``\sqrt{42 q_1 / 13 q_2}``, about 400 s for this clock.

```@example holdover
q1, q2 = 1e-22, 2π^2 * 1e-28              # h₀/2 and 2π²h₋₂
σe(τ, T) = sqrt(q1 * τ * (1 + 6τ / 5T) + q2 * (τ^3 / 3 + 13T * τ^2 / 35))
T = [64, 128, 256, 512, 1024, 4096]
# the stretch, the RMS error measured 256 s on, and the form
[T  [holdover(p; train, tau = [256]).rms[1] for train in T]  σe.(256, T)]
```

The form is that of a line through a continuous stretch. With 64 samples
in the stretch it holds within 0.1 % on white FM and 2 % on random-walk
FM, and with 1024 within 0.2 % on both. One record scatters about it: over
200 records like this one the measured RMS was within 1 % of the form in
the median at every stretch, and four records in five were within 3 % of
it with 64 samples in the stretch and within 17 % with 4096.

A parabola, `drift = true`, also learns a frequency drift. Under a line a
drift ``D`` leaves a mean error of ``\tfrac{1}{2} D τ (τ + T)``, which
shows in the `mean` of the result. The parabola pays for its third
coefficient in noise. Its RMS error over that of a line, on a clock
without drift and with 256 samples in the stretch:

| Noise | ``τ = T/16`` | ``T/4`` | ``T`` | ``4T`` |
|:---|---:|---:|---:|---:|
| white FM | 1.13 | 1.51 | 2.96 | 8.65 |
| flicker FM | 1.11 | 1.31 | 2.07 | 4.81 |
| random-walk FM | 0.82 | 0.92 | 1.27 | 2.18 |

On random-walk FM the parabola predicts better than the line while the
holdover is shorter than two fifths of the stretch.

With `anchor = false` the value of the fit at ``t`` takes the place of the
last reading, which is the extrapolated fit that Delporte et al. analyze.
It averages the white PM of the readings down, and on a clock with
random-walk FM it starts from the misfit of the fit.

## Predicted from a noise model

The same call on a noise model gives the error that a record of that noise
should show. The model is a [`PowerLawFit`](@ref), a [`NoiseFit`](@ref), a
[`ClockModel`](@ref), or a `Dict` of the levels ``h_α``, such as those of
a data sheet.

```@example holdover
fit = powerlaw_fit(p)
pred = holdover(fit; train = 512, tau = [16, 64, 256, 1024], limit = 1e-9)
```

The columns after the 95 % error are the shares of the noise types, as
standard deviations. Measured over predicted, for this record:

```@example holdover
h.rms ./ pred.sigma
```

The error is a weighted sum of phase samples, so its variance follows
without approximation from the spectrum of the sampled record, the one
that [`powerlaw_sigma`](@ref) integrates. The flicker noises are predicted
on the same footing as the others, and no state model stands in for them.
The tables of this page are such predictions. For a long stretch and
`anchor = false` the prediction gives the closed forms of
[Delporte et al. (2000), Eqs. (11)–(13) and (17)–(19)](@cite delporte-2000-extrapolation)
for a parabola and for a line on white, flicker and random-walk FM, within
0.2 % at 4096 samples in the stretch.

The prediction holds noise only. A drift that a line does not learn adds
its mean error.

## With a Kalman filter as the predictor

A record and a [`ClockModel`](@ref) together run the test with the
[Kalman filter](kalman.md) in place of the fit. The filter passes once
over the record, and from every sample past the first `train` seconds its
estimate is carried ``τ`` ahead with no further readings and compared
with the reading there.

```@example holdover
m = clock_model(fit; flicker = true)
k = holdover(p, m; tau = [16, 64, 256, 1024])
```

`claimed` is the standard deviation that the filter itself gives the
error: [`holdover_error`](@ref) with the noise of the reading it is
compared with. The two columns side by side test the model where it
matters for holdover. An `rms` well above `claimed` says that the model
holds too little noise at that ``τ``.

A model made from a fit should keep the flicker FM that the fit finds,
`clock_model(fit; flicker = true)`. Without it that noise is left out of
the model, and in simulation the filter of a clock on its flicker floor
claimed 1.4 to 8.7 times less than the error it made
([What the states buy](@ref)). `sigmatau holdover --kalman` and the
explorer fit the model to the first training stretch alone, so that the
levels are not learned from the readings that the predictions are
compared with, and keep the flicker.

With the right model the filter's prediction is the conditional mean of
the phase, which no predictor betters in mean square
([Särkkä and Svensson (2023), §2.4 and Theorem 6.6](@cite sarkka-2023-bayesian-filtering)).
A least-squares line comes close at its best stretch and falls behind away
from it. For a clock of white and random-walk FM, the standard deviation
of the line's error over that of the filter's, at holdover times from a
tenth of the best stretch to six times it:

| Stretch, over the best | 1/4 | 1/2 | 1 | 2 | 4 | 16 |
|:---|:---|:---|:---|:---|:---|:---|
| line over filter | 1.12 to 1.35 | 1.04 to 1.14 | 1.02 to 1.07 | 1.04 to 1.14 | 1.12 to 1.35 | 1.53 to 2.33 |

The measured RMS errors of a simulated record of that clock gave the same
ratios within 5 %.

!!! note "At the longest holdover times"
    The ratio of `rms` to `claimed` is least certain at the longest
    holdover times. The RMS there rests on few independent trials, and
    the claim on what one stretch shows of the slow noise. In simulation,
    on records of 8192 samples with the first quarter as the stretch, the
    ratio at 2048 s was, in four records out of five:

    | Clock | With the true model | With the model fitted to the stretch |
    |:---|:---|:---|
    | white FM | 0.56 to 1.56 | 0.33 to 1.29 |
    | white FM, random-walk FM from 390 s | 0.55 to 1.37 | 0.60 to 8.6 |

    At ``τ_0`` it was 0.98 to 1.04 with the fitted model. On white FM the
    fit found some flicker or random-walk FM in 61 of the 100 records, and
    in the median the claim was twice that of the true model. On the
    second clock the fit found less than 1 % of the random-walk FM in half
    of the records, and in the median the claim was half that of the true
    model. A ratio well above 1 at long ``τ`` says that the stretch was
    too short to show the slow noise.

## Holdover and the time deviation

The time variance is the mean square of an error of the same kind. With
``\bar{x}`` the phase averaged over ``τ``, the definitions of
``\mathrm{Mod}\,σ_y(τ)`` and ``σ_x(τ)``
([IEEE Std 1139-2022, Eqs. (A.26) and (A.29)](@cite ieee1139-2022-definitions))
give

```math
σ_x^2(τ) = \tfrac{1}{6} \left\langle \left( \bar{x}(t + 2τ) - 2\,\bar{x}(t + τ)
           + \bar{x}(t) \right)^2 \right\rangle ,
```

and the bracket is the next ``τ``-average less the line through the two
before it. So ``\sqrt{6}`` [`tdev`](@ref) is the RMS error of a holdover
test that learns the frequency from two ``τ``-averages and predicts the
third. The third difference does the same with three averages and a
parabola, drift learned: its mean square is
``6τ^2\,\mathrm{Mod}\,σ_H^2(τ)``, which makes that error ``\sqrt{20}``
[`htdev`](@ref). [The Hadamard time deviation](@ref) has the steps.

`grow = k` in place of `train` gives the measured test the matching
stretch, ``k`` times each holdover time, so that a longer holdover is
learned over a longer stretch:

```@example holdover
g = holdover(p; grow = 2)                 # a line over 2τ, then τ ahead
mτ = Int.(g.tau ./ p.tau0)
[g.tau  g.rms ./ (sqrt(6) .* tdev(p, mτ; ci = false).dev)]
```

The two differ in the readings: TDEV averages them over ``τ``, and the
test compares one reading with a least-squares line. At ``τ = τ_0`` an
average is a reading and the two are equal. Above it the ratio depends on
the noise type. On FM noise it settles within a few octaves, and on PM
noise it keeps rising. At long ``τ``:

| Noise | Line over ``2τ``, over ``\sqrt{6}`` TDEV | Parabola over ``3τ``, over ``\sqrt{20}`` HTDEV |
|:---|:---|:---|
| white FM | ``\sqrt{8/5}`` = 1.26 | 1.22 |
| flicker FM | 1.31 | 1.31 |
| random-walk FM | ``\sqrt{452/231}`` = 1.40 | 1.42 |
| flicker PM | 1.3 at ``4τ_0``, 1.9 at ``4096τ_0`` | 1.1 at ``4τ_0``, 1.4 at ``4096τ_0`` |
| white PM | ``\sqrt{τ/3τ_0}`` | ``\sqrt{τ/10τ_0}`` |

The record above goes from the white FM value toward that of random-walk
FM as that noise takes over. On FM noise the test runs parallel to the
scaled time deviation, 1.2 to 1.4 times it, so a scaled time deviation
read as a holdover error understates the error of a least-squares fit by
that factor. On white PM the averages smooth what a single reading does
not, and the two part. With a fixed `train` there is no such relation,
because the stretch no longer scales with ``τ``.

```@docs
holdover
HoldoverResult
HoldoverPrediction
```

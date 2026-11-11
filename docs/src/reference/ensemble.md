# Clock ensembles

A time scale is a clock on paper: at each instant a weighted average of
several clocks, kept as the offset of each clock from the average. It can
be more stable than its best clock, and it goes on when a clock stops.
[`timescale`](@ref) forms one, by any of six algorithms, from clocks that
are read against a common reference.

## Comparisons

The clocks are read against one reference at the same instants, by a
counter with a channel for each clock or by a switch that steps through
them. [`ClockComparisons`](@ref) holds the readings, one column for each
clock. [`simulate_ensemble`](@ref) draws such a set from a
[`ClockModel`](@ref) for each clock and returns the true states with it,
so that a scale can be held against true time.

```@example ensemble
using SigmaTau, Random

# white FM of 1, 2 and 4 × 10⁻¹² at 1 s, each clock with random-walk FM
ms = [ClockModel(; tau0 = 1.0, q1 = q, q2 = 1e-6q) for q in (1e-24, 4e-24, 1.6e-23)]
cc, truth = simulate_ensemble(ms, 20_000; reference = 0, rng = Xoshiro(1))
cc
```

`cc[i]` is the record of clock ``i`` against the reference, and `cc[i, j]`
that of clock ``i`` against clock ``j``, in which the reference cancels.
With `reference = 0` the reference is perfect, which a simulation can have
and a laboratory cannot.

```@docs
ClockComparisons
simulate_ensemble
clock_models
```

## Time scales

[`timescale`](@ref) returns the phase and the frequency of each clock
against the scale, and its weight, at every sample. Three algorithms
average the predictions of the clocks and learn what they need from the
comparisons:

- `:ewfa` gives every clock the same weight
  ([Kettering (2023), §3.1](@cite kettering-2023-clock-ensemble)).
- `:at1`, the NIST algorithm, weighs each clock by the inverse of its mean
  squared prediction error
  ([Weiss and Weissert (1991), Eqs. (1)–(10)](@cite weiss-weissert-1991-at2)).
- `:at2` adds a variance to each frequency estimate of AT1, and with it a
  test for frequency steps (Eqs. (11)–(22) of the same paper).

Three run a Kalman filter on a [`ClockModel`](@ref) of each clock:

- `:greenhall` reduces the covariance of the filter after each update,
  which gives the scale the weights that make its steps smallest
  ([Greenhall (2003)](@cite greenhall-2003-jones-tryon);
  [Greenhall (2006)](@cite greenhall-2006-measurement-noise)).
- `:jones_tryon` is the same filter without the reduction
  ([Jones and Tryon (1983)](@cite jones-tryon-1983-estimating-time)).
- `:stein` filters the difference of each clock from a pivot clock, with
  one set of weights each for phase, frequency and aging
  ([Stein (2003)](@cite stein-2003-demystified)).

The docstring of [`timescale`](@ref) has the equations of each.

```@example ensemble
r = timescale(cc; algorithm = :at1)
```

AT1 learns its weights from the record. The reduced Kalman scale takes
them from the models, here the true ones, as the inverse of the white FM
levels:

```@example ensemble
g = timescale(cc; algorithm = :greenhall, models = ms)
round.(g.weights[end, :]; digits = 3)
```

[`timescale_phase`](@ref) returns the reference against the scale as a
phase record. Against the perfect reference of this simulation that is the
error of the scale, and its Allan deviation is that of the scale itself:

```@example ensemble
using Plots

m = 2 .^ (0:10)
σ(algorithm; kw...) = adev(timescale_phase(timescale(cc; algorithm, kw...)), m)
plot(adev(cc[1], m); label = "clock 1, the best", color = :gray)
plot!(σ(:ewfa); label = ":ewfa")
plot!(σ(:at1); label = ":at1")
plot!(σ(:greenhall; models = ms); label = ":greenhall", legend = :bottomleft)
```

Equal weights let the noisiest clock set the level, and that scale is less
stable than its best clock. The two weighted scales lie on one another
here, below every clock.

### The scales compared

For independent clocks, no average with fixed weights has a smaller Allan
variance at a given ``τ`` than ``1/\sum_i σ_i^{-2}(τ)``, that of the
average weighted by the inverse Allan variances of the clocks at that
``τ``. The tables give the Allan deviation of each scale over that bound,
in simulation, on clocks without reading noise and with their true models.
Each scale is set against that average of the same records, which takes
the scatter of the records out of the ratio. The rows of `:ewfa`, whose
weights are equal, and of the best clock follow from the models.

With four clocks whose white FM and random-walk FM levels both rise by a
factor of 4 from clock to clock, one set of weights is the best at every
``τ`` (400 records of 32 768 samples; the last entry of `:at1` is good to
0.03):

| Scale | ``τ_0`` | ``16τ_0`` | ``256τ_0`` | ``4096τ_0`` |
|:---|---:|---:|---:|---:|
| `:greenhall`, `:stein`, `:jones_tryon` | 1.00 | 1.00 | 1.00 | 1.00 |
| `:at1` | 1.00 | 1.00 | 1.01 | 1.07 |
| `:ewfa` | 2.66 | 2.66 | 2.66 | 2.66 |
| best clock | 1.15 | 1.15 | 1.15 | 1.15 |

[Greenhall (2006), §7](@cite greenhall-2006-measurement-noise) takes three
clocks of white and random-walk FM with ``(q_1, q_2)`` = (4, 0.37),
(25, 0.031) and (4, 0.37) at ``τ_0 = 1``. The second is the noisiest at
short ``τ`` and the quietest at long ``τ``, so the best weights change with
``τ``. Without reading noise, as in
[Greenhall (2007), §3](@cite greenhall-2007-measurement-noise-corrections)
(400 records of 8192 samples; `:at1` and `:at2` from ``64τ_0`` are good
to 0.05):

| Scale | ``τ_0`` | ``4τ_0`` | ``16τ_0`` | ``64τ_0`` | ``256τ_0`` | ``1024τ_0`` |
|:---|---:|---:|---:|---:|---:|---:|
| `:greenhall` | 1.01 | 1.07 | 1.08 | 1.08 | 1.01 | 1.00 |
| `:stein` | 1.03 | 1.17 | 1.18 | 1.09 | 1.01 | 1.00 |
| `:jones_tryon` | 3.10 | 2.63 | 1.38 | 1.01 | 1.00 | 1.00 |
| `:at1` | 1.00 | 1.01 | 1.08 | 1.47 | 1.52 | 1.48 |
| `:at2` | 1.00 | 1.00 | 1.09 | 1.51 | 1.56 | 1.52 |
| `:ewfa` | 1.39 | 1.24 | 1.01 | 1.51 | 1.77 | 1.79 |
| best clock | 1.47 | 1.50 | 1.60 | 1.12 | 1.08 | 1.08 |

The reduced Kalman scale stays within 10 % of the bound, and the
unreduced one is three times above it at ``τ_0`` and reaches it from
``64τ_0``. The 2007 paper reports both: within 1.09 times the bound for
the first, and for the second at the bound from ``64τ_0`` and behind the
first below ``32τ_0``. `:at1` and `:at2` take their weights from the
prediction errors over one interval. They are at the bound at short ``τ``
and half above it at long ``τ``, where the best clock alone is 30 to 45 %
more stable than the scale.

The unreduced filter is the more exposed to the long-term levels of the
models. Its scale was within 1 % of `:greenhall` on the four clocks of one
shape, 5 to 12 % less stable on four clocks of white FM alone, and 5.7
times less stable at ``τ_0`` on four clocks ordered one way in white FM
and the other way in random-walk FM.

### Fitted models

The Kalman scales are as good as their models. [`clock_models`](@ref) fits
them from the comparisons of three or more clocks: it separates each clock
from the others by the Groslambert covariance ([`groslambert`](@ref)) and
fits a power law to its curve.

```@example ensemble
fitted = clock_models(cc)
# made with q1 = 1e-24, 4e-24, 1.6e-23 and q2 = 1e-6 q1
[[m.q1 for m in fitted] [m.q2 for m in fitted]]
```

The white FM levels come out within 5 %. The random-walk FM of these
clocks shows only at the longest ``τ`` of 20 000 samples, and its levels
do not: the first is 50 times too large here, and over 16 such records the
fitted ``q_2`` ran from zero to 32 times the truth. `:greenhall` and
`:stein` on the fitted models were then within 2 % of `:greenhall` on the
true models up to 100 s and a quarter less stable at 1000 s, where `:at1`,
which takes no model, was 1 % behind. `:jones_tryon` on the fitted models
was 1.9 times less stable at every ``τ``. With 200 000 samples the fitted
``q_2`` lay between a quarter of the truth and 2.4 times it, `:greenhall`
and `:stein` were within 5 % of the true models up to 10 000 s,
`:jones_tryon` 5 to 10 % behind, and `:at1` 16 % behind at 10 000 s. A
record has to be long enough to show the slow noise of each clock before a
Kalman scale gains on AT1.

### Reading noise and the pivot

With noise on the readings the clocks, each corrected by the filter's
estimate of its phase, no longer agree, and the scale is read at one of
them, `pivot`: the reading of that clock less the estimate of its phase.
The readout carries the white PM of that reading. In simulation, with
three equal clocks of white FM ``10^{-12}`` at 1 s read with 10 ps of
white PM, the Allan deviation of the readout was 29 times that of the
implicit ensemble mean at ``τ_0`` and 1.2 times at ``1024τ_0``.
[`implicit_mean`](@ref) returns that mean for a simulated ensemble, whose
true phases are known. No measurement gives it
([Greenhall (2006), §4](@cite greenhall-2006-measurement-noise)). The
average of the readouts at the three clocks, with the weights of the
scale, was 1.67 times quieter than one readout at ``τ_0`` and 1.15 times
at ``1024τ_0``.

Without reading noise every pivot gives the same `:greenhall` and
`:jones_tryon` scale. `:stein` runs the filters of one pivot only and
depends on it: by less than 0.1 % in Allan deviation on the four clocks of
one shape, and by up to 3.7 % on Greenhall's three.

### Two clocks

Two clocks give one difference, which cannot say which of them is the
better. In simulation with white FM levels of 1 and 4, two clocks read
against each other kept equal weights under `:at1`. Read against a
quieter third reference, they took their weights from the start of the
record, where each clock is seen alone, and learned nothing after it: the
weights held under `bias = :exact` and relaxed to equal within the
`memory` under `bias = :paper`.

## Frequency-step detection

With `steps = :best` or `:real`, `:at2` tests each clock for a step in
frequency. A clock that fails is taken out of the scale for its
``τ_{\min}``, the scale is run again from the step, and the clock goes on
being estimated while it is out
([Weiss and Weissert (1991), §4](@cite weiss-weissert-1991-at2)). The
docstring of [`timescale`](@ref) gives the test and where it departs from
the paper.

In simulation, a step was put on a clock with white FM of
``1.4 × 10^{-11}`` at the sample interval, one of four clocks (20 records
of 8000 samples, threshold 5):

| Step | Caught | Samples until caught: least, median, most |
|:---|:---|:---|
| ``1 × 10^{-10}`` | 20 of 20 | 1, 1, 2 |
| ``5 × 10^{-11}`` | 20 of 20 | 3, 6.5, 11 |
| ``2 × 10^{-11}`` | 18 of 20 | 24, 41.5, 82 |
| ``1 × 10^{-11}`` | 6 of 20 | 63, 159, 264 |

A step several times the white FM of the clock is caught at once, and one
of its size after tens of samples or not at all. On the same ensemble
without a step, a threshold of 4 took a clock out 15 times in 24 records
and raised the Allan deviation of the scale at ``τ_0`` by 1 % in the
median and 4 % at most. Thresholds of 5 and 6 took none out. The paper
tests at 4.

What the test buys shows on five equal clocks of white and random-walk
FM, with a ``τ_{\min}`` of 40 samples and a step on one clock (20
records). Steps of 3 and 10 times the white FM of the clock were caught,
and the clock came back after one ``τ_{\min}`` with 0.97 to 0.98 of the
step in its frequency estimate. Larger steps took two: the first reading
after a large step raises the prediction error of the clock, the filter
then trusts its readings less, and at its first return the clock had
learned 0.45 of a step of 30 times its white FM and 0.05 of a step of 100
times. It came back with a weight of 0.001 or less, was caught again, and
held 0.98 of the step at its second return. At a return the scale stepped
by 0.4 to 1.3 times its usual step. Two hundred samples after the step,
the scale had lost 44 to 56 ns to it, against 830 to 1680 ns with the
same steps left in.

## Clocks that drop out and return

A `NaN` in the comparisons is a missing reading. A clock without a reading
keeps its prediction and takes no part in the scale at that sample, and no
scale steps when a clock drops out. Over the gap the prediction of the
clock gathers an error, which every scale but one keeps out when the
readings return:

- In `:ewfa`, `:at1`, `:at2` and `:stein` the clock joins again as a new
  clock does. Its phase is taken from the scale at the sample of its
  return, and it weighs in from the next. The sources give no procedure
  for a clock that returns, and this one is the package's own.
- `:greenhall` gives a clock no weight in its mean while it is unread, the
  procedure of
  [Greenhall (2006), §10](@cite greenhall-2006-measurement-noise).
- `:jones_tryon` has no such rule. Its filter leaves a missing difference
  out of the measurement
  ([Jones and Tryon (1983), §2](@cite jones-tryon-1983-estimating-time)),
  and the scale takes a share of the error when the clock returns.

In simulation, three clocks with white FM of 1, 2 and ``4 × 10^{-12}`` at
1 s ran at frequency offsets of 0, ``1 × 10^{-12}`` and ``-2 × 10^{-12}``,
and one of them went unread for ``g`` samples. The table gives the step of
the scale at the return over its RMS step, in the median of 500 records,
which is good to 10 %. For an ordinary step that ratio is 0.67.

| Scale, and the clock unread | ``g = 10`` | ``g = 100`` | ``g = 2000`` |
|:---|---:|---:|---:|
| all but `:jones_tryon`, the noisiest clock | 0.7 to 0.8 | 0.7 to 0.8 | 0.7 to 0.8 |
| all but `:jones_tryon`, the best clock | 1.0 to 1.4 | 1.1 to 1.3 | 1.0 to 1.3 |
| `:jones_tryon`, the noisiest clock | 0.8 | 1.5 | 8.4 |
| `:jones_tryon`, the best clock | 3.1 | 12 | 61 |
| `:stein`, its pivot | 2.4 | 7.5 | 38 |

At the sample of a return the scale is carried by the other clocks, and
its step is theirs: larger than usual when the best clock is the one that
returns, and no larger for a longer gap. `:jones_tryon` steps with the
length of the gap, and by more for a better clock.

A clock whose first reading comes after the start is measured against the
scale. In `:ewfa`, `:at1` and `:at2` its starting frequency is its slope
over its first tenth of the record less the slope of the scale over the
same samples. In `:greenhall` and `:jones_tryon` it begins as the clock of
most weight read with it plus a line through its first 16 readings against
that clock, with a wide covariance. The reference of the readings cancels
in both, and the scale takes an ordinary step at the join. The sources
give no procedure
([Greenhall (2003), §5](@cite greenhall-2003-jones-tryon) names the
insertion of a new clock as a practical difficulty of a Kalman scale), and
these are the package's own. In simulation, with three clocks whose white
FM variances were as 1 : 4 : 16 and one of them first read at sample 1001
of 8192, the rate of the `:ewfa`, `:greenhall` and `:jones_tryon` scales
changed at the join by what the noise of the clocks gives. That of `:at1`
and `:at2` changed by 1.3 times as much when the best clock was the late
one: the error of its frequency, measured over 819 samples, times its
weight of 0.75. `:stein` starts a difference from 16 readings and passes
the later corrections of its filter to the pivot, and its rate changed by
9 times as much when the best clock was the late one.

`:stein` measures every clock against its pivot. While the pivot is unread
no difference is measured, every clock runs on its prediction, and the
scale steps when the pivot returns. A `:stein` scale needs a pivot that is
read throughout.

`:ewfa` keeps its phase through a gap and not its rate. It carries each
clock on its last frequency step alone, so when a clock drops out the
frequency of the scale shifts by the noise of that one step, and stays
there until the clock returns. With the noisiest clock unread for 2000
samples, the mean frequency of `:ewfa` over the gap differed from that of
the 2000 samples before it by 37 times as much as in the same records
without the gap. For the other scales the factor was 1.0 to 1.2.

The clocks need not be read at every sample of the record. In `:ewfa`,
`:at1` and `:at2` the interval of the equations is the time since the last
sample with a reading
([Weiss and Weissert (1991), §2.8](@cite weiss-weissert-1991-at2)), so
clocks read together at every tenth sample give the scale of those readings
alone. They do have to be read together: a clock read less often than the
others joins anew at each of its readings and never has a weight. For
clocks on different schedules, form the scale on the samples they share:
the Kalman scales take each reading when it comes, but then lean on the
clock read most often.

```@docs
timescale
TimescaleResult
timescale_phase
implicit_mean
```

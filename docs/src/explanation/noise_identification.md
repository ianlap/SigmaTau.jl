# Noise identification

The confidence interval of a deviation depends on which power-law noise
dominates at that averaging time, and so do the bias corrections of the
total and Thêo estimators ([SP1065 §5.5.1](@cite riley-2008-sp1065)).
Unless the noise type is stated, every estimator finds it from the record
with [`identify_noise`](@ref), one averaging time at a time. This page
describes the tests, where the package departs from the published
procedure, and how far a label can be trusted.

## What depends on the label

The label sets the equivalent degrees of freedom (EDF) of the estimate.
Here is one Allan deviation, at 16 averaged points of a 1024-point record,
under each label in turn:

```@example noiseid
using SigmaTau, Random

p = noise_gen(PhaseData, 1024, 1.0; sigma1 = Dict(:WHFM => 1e-12), rng = Xoshiro(1))
for n in (:WHPM, :FLPM, :WHFM, :FLFM, :RWFM)
    r = adev(p, [64]; noise = n)
    lo, hi = round.((ci_lower(r)[1], ci_upper(r)[1]) ./ r.dev[1]; digits = 2)
    println(rpad(n, 6), "EDF ", lpad(round(r.edf[1]; digits = 1), 5),
            "   interval ", lo, " to ", hi, " of the estimate")
end
```

Taking one FM noise for another moves the interval a little. Taking FM for
PM makes it several times too narrow: the overlapping terms of the estimate
are nearly uncorrelated under white PM, so its EDF follows the record
length and not the number of averaged points. The rules below are weighted
against that mistake.

## The lag-1 autocorrelation test

The main test is that of
[Riley and Greenhall (2004)](@cite riley-2004-lag1-acf), also given in
[SP1065 §5.5.5 and §5.5.6](@cite riley-2008-sp1065). For averaging factor
``m``:

1. The phase is decimated by ``m``, which leaves the phase record at sample
   interval ``mτ_0``, and a least-squares quadratic is removed from it, as
   allantools does.
2. The lag-1 autocorrelation ``r_1`` of the series gives
   ``δ = r_1/(1 + r_1)``.
3. If ``δ ≥ 0.25`` the series is differenced and step 2 repeated, up to
   `dmax` times: twice for the Allan family, three times for the Hadamard
   family.
4. After ``d`` differences the exponent of the noise is
   ``α = 2 - 2(δ + d)``.

The estimate is a real number. On the record as sampled (``m = 1``) it is
rounded to the nearest integer, as the paper does. The test is reliable
from about 30 points ([SP1065 §5.5.6](@cite riley-2008-sp1065)), and it is
used while ``N/m ≥ 30``.

Samples farther than five standard deviations from the mean of the record
are left out of the test.

## Decimated records

The test is derived for discrete-time noise with the spectrum
``(2 \sin πf)^{-2δ}``. A record decimated by ``m ≥ 2`` behaves as a
continuous process read at intervals, and on three noise types the
difference is large enough to change the label. The package allows for
each.

**Flicker FM and random-walk FM.** The second differences of decimated
phase have lag-1 correlations of ``(9 \ln 3 - 16 \ln 2)/(8 \ln 2) ≈ -0.217``
and 1/4 on these two, where the discrete-time model has −1/3 and 0. The
estimate therefore settles on ``α = -1.45`` and −2.4, not −1 and −2. After
a third difference it reads −3.0, −3.93 and −4.57 on random-walk,
flicker-walk and random-run FM. All five values follow from the generalized
autocovariances of
[Greenhall and Riley (2003)](@cite greenhall-2003-edf-stability). For
``m ≥ 2`` the limits between the FM labels are placed halfway between the
readings, at −0.72, −1.92, −3.46 and −4.25, where rounding puts them at the
half-integers. Rounding at −1.5 cuts through the middle of the flicker FM
readings: in simulation (4096 points) flicker FM keeps its label in 71 % of
records at ``m = 8`` and 48 % at ``m = 64`` with rounding, and in 100 % and
87 % with these limits.

**Flicker PM.** Decimation whitens it. On records of 16 384 points the
estimate averages 1.0 at ``m = 1``, 1.5 at ``m = 16`` and 1.8 at
``m = 256``, and the test alone calls flicker PM white in 45 % of records
at ``m = 16`` and in 94 % at ``m = 256``. From ``m = 4`` a white PM result
is therefore checked with ``R(n)``, the ratio of the modified to the plain
Allan variance ([SP1065 §5.5.2 and §5.10](@cite riley-2008-sp1065)), which
reads every sample of the record. ``R(n)`` is ``1/m`` on white PM and falls
only as ``1/\ln m`` on flicker PM. A result nearer the flicker PM value is
relabeled. Flicker PM then keeps its label in 99 % of those records or
more, and white PM loses its own in 0.3 % at most.

**PM or FM.** With few decimated points the test takes FM for PM, the
costly mistake: in simulation (8192 points) white FM in 8 % of records at
64 points and in 22 % at 32. ``R(n)`` is 1/2, 27/40 and 33/40 on white,
flicker and random-walk FM, against 0.26 on flicker PM at ``m = 16`` and
less beyond. From ``m = 16`` a PM result stands only if ``R(n)`` lies
nearer its flicker PM value than 1/2. One that fails is reported as white
FM. With this rule, and with the conditions of the next section below 30
points, none of the 1500 simulated records of white FM was labeled PM at
any averaging factor from 64 to 2048.

## Fewer than 30 averaged points

Below 30 points the lag-1 estimate is too noisy, and the package turns to
the ``B_1`` ratio, as [Stable32](@cite hamilton-2008-stable32-manual) does.
``B_1`` is the standard variance of the ``m``-averaged frequency values
divided by their Allan variance. Its expected value depends on the number
of values ``N`` and on the slope ``μ`` of the Allan variance, ``σ_y^2(τ) ∝
τ^μ`` ([SP1065 §5.5.2 and §5.7](@cite riley-2008-sp1065)):

```math
B_1(N, μ) = \frac{N\,(1 - N^{μ})}{2\,(N - 1)\,(1 - 2^{μ})} ,
```

which is ``N/2`` on random-walk FM, ``N \ln N / (2 (N-1) \ln 2)`` on
flicker FM, 1 on white FM and ``2(N+1)/(3N)`` on PM. The measured ratio is
compared with the geometric means of neighboring values, and a PM result is
split into white and flicker PM with ``R(n)``. The test stops at
random-walk FM.

Three choices here are the package's own.

- **The length ``N``.** The expected ratios are evaluated at the number of
  averaged values the ratio was computed from, which is the ``N`` of the
  definition. `identify_noise(...; b1_length=:record)` evaluates them at
  the length of the whole record, which follows Stable32's labels more
  closely. At that length the limit for random-walk FM lies far above what
  a few values can reach. In simulation (1024 points, 6 to 24 averaged
  values) the record length finds random-walk FM in at most 2 % of records
  and flicker FM in 23 to 46 %, where the default finds them in 59 to 74 %
  and 34 to 78 %. White FM pays for it: 71 to 96 % against 86 to 100 %.
- **A floor of 6 values.** With 2 to 4 values the ratio is close to a coin
  toss, so below 6 averaged values the result is the one found at
  ``⌊N/6⌋``, the longest averaging factor that still has 6. The value 6
  comes from simulation, not from the literature.
- **Two conditions for PM.** In a sum of power-law noises the dominant
  type only gets redder toward longer ``τ``. A result from fewer than 30
  points is therefore PM only if the result at ``⌊N/30⌋``, the longest
  averaging factor that has 30 points, is PM, and ``R(n)`` at ``m``
  agrees. With 6 points or more ``R(n)`` must lie below the midpoint of
  ``1/m`` and 1/2, the value it takes where white PM and white FM hold
  equal shares of the Allan variance, and it then decides between white
  and flicker PM. With fewer it must lie on the white PM side of the limit
  between the two PM types, and it is read at ``⌊N/4⌋`` where ``m`` is
  larger. Below ``m = 16``, which a record shorter than 480 samples
  reaches with fewer than 30 points, a PM result of the ``B_1`` test takes
  the place of ``R(n)``. Every other result is the FM type of the ``B_1``
  ratio, with white FM for a PM result of it.

The conditions guard both ways. FM noise is not labeled PM on the strength
of a few points, and white PM is not labeled FM: in simulation (1024
points) the ``B_1`` ratio alone takes white PM for FM in one record in six
at 16 averaged points and in one in three at 8, which puts the FM bias
factor of the total estimators on PM noise, and with the conditions white
PM keeps its label down to 2 points. The cost falls on flicker PM with
fewer than 6 averaged points, which was labeled white FM in four records
in five or more (78 % of the records of 1024 points, 86 % of those of
8192) and flicker PM in none. Its interval is then wide, and the total,
modified total and Hadamard total variances read 19, 10 and 25 % high at 4
points.

With `b1_length=:record` neither the floor nor the conditions apply: the
result of the ``B_1`` test stands, subject to the ``R(n)`` check of the
previous section. `identify_noise(...; b1=false)` leaves the ratio out,
so that an FM result from fewer than 30 points is the type found at
``⌊N/30⌋``.

The lag-1 test removes a quadratic from each decimated series, but the
standard variance of the ``B_1`` ratio is taken from the record as given. A
linear frequency drift that dominates the Allan deviation at the longest
averaging times is therefore labeled random-walk FM there by the Allan
family. The deviation itself rises with the drift, so remove the drift
first with [`detrend`](@ref) if those points matter.

## The Hadamard family

The Hadamard estimators converge for two more noise types, flicker-walk FM
(``α = -3``) and random-run FM (``α = -4``), so their labels are found with
up to three differences ([SP1065 §5.5.6](@cite riley-2008-sp1065)). They
are also blind to linear frequency drift, and by default their labels are
too: the noise type is identified on the record less its least-squares
quadratic.

On that record the ``B_1`` ratio is not used. Its limits are those of a
record with the drift left in: in simulation, at 16 averaged values,
removing the quadratic lowered the mean ratio from 7.8 to 3.2 on
random-walk FM and from 2.2 to 1.7 on flicker FM. The test also stops at
random-walk FM, two types short of what the family converges for. With
fewer than 30 averaged points an FM label is therefore the one found at
``⌊N/30⌋``, the longest averaging factor that has 30, and a PM label is
subject to the two conditions above. A change from one FM type to another
beyond that averaging time is not seen. State it with `noise`.

`remove_drift=false` on [`hdev`](@ref), [`mhdev`](@ref), [`htdev`](@ref)
and [`htotdev`](@ref) identifies the noise on the record as given, with the
``B_1`` ratio, which is what Stable32 does. The deviation is the same
either way. In simulation, on white FM with a drift that dominates the
Allan deviation at the longest averaging times, the default labeled the
points with 16, 8 and 4 averaged values white FM in 98 % of records, and
`remove_drift=false` labeled them random-walk FM in every record. On
flicker-walk and random-run FM it is the other way round: with
`remove_drift=false` neither was found in any record below 30 points,
where the default kept the labels it found at 30 (tables below).

## How far to trust a label

Share of simulated records of pure noise (8192 points, 1000 records per
type) in which the estimator reports the true type, by the number of
averaged points ``N/m``. Sampling error is about 1.5 points.

| Allan family | 128 | 64 | 32 | 16 | 8 | 4 |
|:---|---:|---:|---:|---:|---:|---:|
| white PM | 100 | 100 | 100 | 100 | 100 | 100 |
| flicker PM | 100 | 99 | 98 | 98 | 96 | 0 |
| white FM | 100 | 100 | 99 | 89 | 73 | 66 |
| flicker FM | 98 | 86 | 51 | 64 | 38 | 32 |
| random-walk FM | 100 | 98 | 91 | 67 | 59 | 60 |

| Hadamard family | 128 | 64 | 32 | 16 | 8 | 4 |
|:---|---:|---:|---:|---:|---:|---:|
| white PM | 100 | 100 | 100 | 100 | 100 | 100 |
| flicker PM | 100 | 99 | 98 | 98 | 96 | 0 |
| white FM | 100 | 100 | 99 | 100 | 100 | 100 |
| flicker FM | 98 | 86 | 51 | 51 | 51 | 51 |
| random-walk FM | 98 | 92 | 82 | 82 | 82 | 82 |
| flicker-walk FM | 98 | 90 | 72 | 72 | 72 | 72 |
| random-run FM | 100 | 98 | 88 | 88 | 88 | 88 |

From 64 points up the label is right in at least 86 % of records. Flicker
FM is the hard case throughout: at 32 points the lag-1 test finds it in
half of the records. Below 30 points the Hadamard family repeats the FM
label found at 30, an FM label of the Allan family is weak evidence of the
type, and flicker PM is not found with fewer than 6 points.

The intervals hold up better than the labels, because the wrong label is
usually a neighbor. In the same simulation the 68.3 % interval of the
Allan, modified Allan, Hadamard and modified Hadamard deviations, with the
label found automatically, holds the true value in 64 to 76 % of records on
every FM noise at every ``N/m`` in the tables. The low end is the Allan
deviation on random-walk FM at 8 and 4 points. On PM noise the figure is
66 to 73 % at 32 points or more and 66 to 78 % below, with one exception:
at 4 points the intervals of the Allan and Hadamard deviations of flicker
PM, which is labeled white FM there, hold the true value in 96 % of
records. They are on the safe side and wide: the upper limit is 2.2 and
2.6 times the lower, where the true label gives 1.25 and 1.3.

A record whose noise type changes is harder, because near the crossing no
label is right. On simulated records where white or flicker PM gives way
to FM between 1/64 and 1/8 of the record, the 68.3 % interval of the
Allan, Hadamard and total deviations held the true value in 64 to 96 % of
records at 16 to 4 averaged points, and that of the modified total and
Hadamard total deviations in 58 to 86 %.

Three things a label is not:

- It is not a component of the record. Where two noise types cross, the
  test reads the type between them, as the example below shows.
- It is not evidence of noise. A periodic or wandering term moves the
  labels near its time scale: in simulation a sinusoid on white FM was read
  as flicker FM and random-walk FM at averaging times a sixteenth and an
  eighth of its period.
- It is not shared between the families. The Allan and Hadamard families
  identify the noise separately and can differ at the same ``τ``.

```@example noiseid
p = noise_gen(PhaseData, 16384, 1.0; sigma1 = Dict(:WHFM => 1e-12, :RWFM => 1e-14),
              rng = Xoshiro(1))
r = adev(p)
foreach((τ, n) -> println(lpad(Int(τ), 5), " s  ", n), r.tau, r.noise_type)
```

The record holds white FM and random-walk FM that cross near
``τ = 100`` s, and no flicker FM.

## Stating the noise type

When the noise type is known, state it. Every estimator that reports an
interval takes `noise`:

```@example noiseid
r = adev(p, [64, 128, 256]; noise = :WHFM)
println(join(r.noise_type, " "))
```

| `noise` | Label at each ``τ`` |
|:---|:---|
| `:auto` | identified from the record (the default) |
| a label, such as `:WHFM` | that type everywhere |
| a vector of labels | one per averaging factor |
| a [`PowerLawFit`](@ref) | the strongest component of the fitted model at that ``τ`` |

A fitted model gives labels that change once, at the crossing, with no
flicker FM between:

```@example noiseid
f = powerlaw_fit(p; alphas = [0, -2])
println(join(adev(p; noise = f).noise_type, " "))
```

The command line takes `--noise auto`, `--noise fit` or a label, and the
explorer has the same choice under **Noise type**. A type redder than an
estimator converges for counts as the reddest it does: `:RRFM` given to
[`adev`](@ref) is taken as `:RWFM`.

[Coming from Stable32](../guide/from_stable32.md) lists the points at which
these labels differ from Stable32's, with the option that reproduces each
where there is one.

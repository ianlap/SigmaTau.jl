# Choosing an estimator

The estimator follows from the question asked of the record and from what
the record contains. The table gives the short answer and the sections
below it the reasons, after the advice of
[SP1065 §5.2 and §10.13](@cite riley-2008-sp1065). With no question yet,
run [`stability`](@ref), which computes `adev`, `mdev`, `hdev` and `mhdev`,
and compare the curves.

| Question | Use | Reason |
|:---|:---|:---|
| A first look, or a figure to hold against a specification | [`adev`](@ref) | The overlapping Allan deviation is the general-purpose measure and the most widely used. |
| Which phase noise limits the short averaging times | [`mdev`](@ref) | It separates white PM from flicker PM, which the Allan deviation does not. |
| The stability of a clock that drifts | [`hdev`](@ref), [`mhdev`](@ref) | The third difference does not respond to linear frequency drift. |
| The long averaging times, with better confidence | [`totdev`](@ref), [`mtotdev`](@ref), [`htotdev`](@ref), [`mhtotdev`](@ref) | A total estimator extends the data by reflection and gains degrees of freedom. |
| Averaging times beyond half the record | [`theo1`](@ref), [`theoh`](@ref) | Thêo1 reaches three quarters of the record. Its bias removal holds on FM noise of one type. |
| The time error of a clock or a link, in seconds | [`tdev`](@ref), [`htdev`](@ref), [`ttotdev`](@ref) | The time deviations are the modified deviations times ``τ``, within a constant. |
| The peak time error, or a telecom mask | [`mtie`](@ref), [`tierms`](@ref), [`mask_check`](@ref) | MTIE is the largest peak-to-peak time error within any interval of a given length. |
| The same, from the shortest record | [`pdev`](@ref) | The parabolic deviation has the two slopes of `mdev` and is the more efficient of the two. |
| Whether the stability changes during the record | [`dadev`](@ref), [`dhdev`](@ref) | The dynamic deviations follow ``σ_y(τ)`` over a sliding window. |
| Each clock, from comparisons among three or more | [`ncornered_hat`](@ref), [`groslambert`](@ref) | They separate the clocks that the pairwise deviations combine. |
| One device, measured through two channels | [`cross_adev`](@ref) | Noise in only one channel averages out of the covariance ([Two-channel measurements](cross.md)). |
| A measurement still in progress | [`StreamingStability`](@ref), [`LiveStability`](@ref) | They update the deviation with each sample. |

[`estimators`](@ref) returns the list as data, and [`deviation`](@ref) runs
an estimator chosen by name.

## Drift

A linear frequency drift ``D`` adds ``D^2τ^2/2`` to the Allan variance
and hides the noise at long averaging times
([Power-law noise](../explanation/noise.md)). There are two ways around it.

The Hadamard family (`hdev`, `mhdev`, `htdev`, `htotdev`, `mhtotdev`)
takes a third difference of the phase, to which a linear frequency drift
contributes nothing. It also converges on flicker-walk and random-run FM,
where the Allan family does not. The price is confidence: at the same
averaging time a Hadamard deviation has fewer degrees of freedom than the
Allan deviation ([SP1065 §5.2](@cite riley-2008-sp1065)).

The other way is to estimate the drift, remove it, and use the Allan
family on what remains. SP1065 advises it as the usual course: the
stability plot shows the noise, and the drift is analyzed separately
([§5.19 and §11.2](@cite riley-2008-sp1065)).

```julia
d = drift_fit(p)                          # the drift and its model
adev(detrend(p; method = :quadratic))     # the noise that remains
```

[Gaps, outliers and drift](preparing.md) covers both functions.

## Phase noise

White and flicker PM both give the Allan deviation a slope near ``τ^{-1}``.
The modified deviations average the phase over ``τ`` before differencing
and fall as ``τ^{-3/2}`` on white PM and ``τ^{-1}`` on flicker PM. Use
[`mdev`](@ref), or [`mhdev`](@ref) when the record also drifts, wherever a
phase noise is to be identified or its level read off.

[`pdev`](@ref) has the same two slopes. It estimates each frequency by
linear regression on the phase, the best estimate in white phase noise,
which makes it the more efficient of the two and lets it identify a noise
type from a shorter record
([IEEE Std 1139-2022, A.3.4](@cite ieee1139-2022-definitions)).

## Long averaging times

At the longest averaging times an estimator rests on a few terms. The
equivalent degrees of freedom (EDF) say how few. Here they are for three
estimators at ``τ ≈ T/32``, ``T/8`` and ``T/4`` on white FM, and for Thêo1
before its bias is removed:

```@example choosing
using SigmaTau, Random
p = noise_gen(PhaseData, 4096, 1.0; sigma1 = Dict(0 => 1e-12), rng = Xoshiro(1))
edf(est, m; kw...) = round.(est(p, m; noise = :WHFM, kw...).edf; digits = 1)

(adev = edf(adev, [128, 512, 1024]), totdev = edf(totdev, [128, 512, 1024]),
 theo1 = edf(theo1, [170, 682, 1366]),
 theo1_uncorrected = edf(theo1, [170, 682, 1366]; correct_bias = false))
```

The total estimators extend the record by reflection at both ends, which
adds half again to the EDF of the Allan deviation at a quarter of the
record. Each keeps the properties of the deviation it extends: `totdev`
goes with `adev`, `mtotdev` with `mdev`, `htotdev` with `hdev` and
`mhtotdev` with `mhdev`. They are biased by an amount that depends on the
noise type, which the package identifies and corrects for
([Bias correction](../explanation/bias.md)). State the noise type with
`noise=` when it is known, because the label is least certain where these
estimators are most useful.

Thêo1 has about three times the EDF of the Allan deviation and reaches
three quarters of the record, but it is biased. [`theo1`](@ref) removes
the bias by a ratio measured on the record, which holds on FM noise of a
single type and adds a scatter of its own, about 15 % in deviation. The
EDF it reports, and its interval, include that scatter
([Confidence intervals](../explanation/confidence.md)): in the table the
corrected estimate has fewer degrees of freedom than the Allan deviation
at ``T/32`` and more from ``T/8``. It reads wrong where PM
noise dominates and where the noise type changes within the record
([Bias correction](../explanation/bias.md#Thêo1) has the sizes). At short
averaging times that scatter leaves it less precise than the Allan
deviation, so for a whole curve use [`theoh`](@ref), which joins the Allan
deviation at short averaging times to the corrected Thêo1 at long ones.

## Time error

[`tdev`](@ref), [`htdev`](@ref) and [`ttotdev`](@ref) are the modified
deviations as a time, in seconds: the usual measure of a time-transfer
link or a network clock. [`mtie`](@ref) gives the largest peak-to-peak
time error within any interval of length ``τ``, and [`tierms`](@ref) the
rms change of the phase over ``τ``, which includes any frequency offset.
MTIE is a peak, not an average: a single outlier sets it, and it has no
confidence interval. The ITU-T masks for network clocks limit MTIE and
TDEV, and [`mask_check`](@ref) compares a result with one of them.

## Averaging times

Octave spacing, the default, gives points evenly spaced on a logarithmic
axis, and `Decade` gives round values. `AllTaus` computes every averaging
factor. It shows a periodic disturbance as a ripple in the curve
([SP1065 §5.23](@cite riley-2008-sp1065)), and costs a call for each
factor.

```julia
adev(p, AllTaus)
```

## Counter readings

What a deviation of frequency readings estimates depends on how the
counter averaged each reading
([IEEE Std 1139-2022, A.3.1.2](@cite ieee1139-2022-definitions)). A Π
counter weights the interval uniformly, and the two-sample deviation of
its readings is the Allan deviation. A Λ counter weights it with a
triangle and an Ω counter with a parabola. The two-sample deviation of
their readings is a modified Allan deviation and a parabolic deviation,
whatever the function that computed it is called. Phase samples avoid the
question. The standard also prefers them where there is dead time between
the readings, which this package does not correct for.

## What each estimator needs

| Estimators | Gaps | Cost for each ``τ`` |
|:---|:---|:---|
| `adev`, `mdev`, `tdev`, `hdev`, `mhdev`, `htdev` | skipped | ``O(N)`` |
| `mtie`, `tierms` | skipped | ``O(N)`` |
| `pdev` | must be filled | ``O(N)`` |
| `totdev` | must be filled | ``O(N)`` |
| `mtotdev`, `ttotdev`, `htotdev`, `mhtotdev` | must be filled | ``O(Nm)``, on all threads |
| `theo1`, `theoh` | must be filled | ``O(Nm)``, on all threads |
| `dadev`, `dhdev` | must be filled | one deviation for each window |

An estimator that skips gaps leaves out the terms that touch a missing
sample. The others need [`fillgaps`](@ref) first, which changes the noise
inside each gap, so prefer an estimator that skips them when the record
has many ([Gaps, outliers and drift](preparing.md)). Start Julia with
`--threads=auto` for the last rows of the table: on a long record the
modified total and Thêo1 estimators are the only ones that take
noticeable time.

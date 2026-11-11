# Coming from Stable32

Stable32 and this package compute the same statistics from the same
definitions, those of [SP1065](@cite riley-2008-sp1065), whose author wrote
Stable32. With Stable32's noise types and bias convention stated, the
deviations agree to the five digits that Stable32 prints
([Validation](../explanation/validation.md)). A default call can differ,
because this package identifies the noise and corrects bias differently in
places. This page maps the functions of Stable32, as its
[manual](@cite hamilton-2008-stable32-manual) describes them, onto this
package, lists each difference, and shows how to reproduce a Stable32
table.

## Reading a file

[`read_data`](@ref) reads the `.phd` and `.frd` files of Stable32 and text
files with its header:

```julia
using SigmaTau
p = read_data("clock.phd")      # phase, with τ₀ from the Tau: line
f = read_data("clock.frd")      # fractional frequency
```

Without a `Tau:` line the sample interval comes from the time tags. MJD
time tags give the record its date ([`start_time`](@ref)), and `scaling=`
does the work of the Multiplier. The `Addend:`, `Multiplier:`,
`Tag Scale:` and `Tag Offset:` lines of a header are not applied.

!!! note "Zeros are gaps in a Stable32 file"
    Stable32 marks a missing sample with a zero: every zero of a frequency
    record, and every zero of a phase record but its first and last sample
    ([Stable32 manual, "Gaps & Outliers"](@cite hamilton-2008-stable32-manual)).
    `read_data` reads a `.phd` or `.frd` file the same way, and gives each
    such zero as `NaN`, which marks a missing sample in this package.
    Where the zeros are readings, as quantized data can have them, keep
    them:

    ```julia
    f = read_data("counts.frd"; keep_zeros = true)
    ```

    A Stable32 file under another extension is read as plain columns, in
    which a zero is a sample. `format = :stable32` reads it as above.

## Functions

| Stable32 | SigmaTau |
|:---|:---|
| Open, Read Archive, multicolumn data | [`read_data`](@ref), [`read_phase`](@ref), [`read_frequency`](@ref) |
| Regularize | the readers, which place a `NaN` where the time tags skip a sample |
| Convert | `FrequencyData(p)`, `PhaseData(f)` |
| Normalize | `detrend(p; method = :mean)` |
| Average | for a phase record, `PhaseData(p.x[1:m:end], m * p.tau0)`. Every estimator takes the averaging factor itself |
| Part | `PhaseData(p.x[a:b], p.tau0)` |
| Scale | `scaling=` on the readers, or arithmetic on the samples |
| Fill | [`fillgaps`](@ref), by another method ([below](@ref stable32-differences)) |
| Check | [`find_outliers`](@ref), by the same test, with five scaled median absolute deviations as the default. [`remove_outliers`](@ref) puts a filled value where Stable32 leaves a gap |
| Statistics | [`summarize`](@ref) |
| Drift | [`drift_fit`](@ref) for the fit, [`detrend`](@ref) for the removal ([Drift types](@ref stable32-drift)) |
| Sigma | an estimator at one averaging factor, `adev(p, [16])` |
| Run | an estimator, or [`stability`](@ref) for several |
| Dynamic Stability | [`dadev`](@ref), [`dhdev`](@ref) |
| Power | [`psd_y`](@ref), [`psd_x`](@ref), [`ssb_phase_noise`](@ref) |
| Domain | [`powerlaw_sigma`](@ref), [`powerlaw_psd`](@ref) |
| Noise | [`noise_gen`](@ref) |
| 3-Cornered Hat | [`ncornered_hat`](@ref), for three clocks or more, and [`groslambert`](@ref) |
| Plot | the plot recipes ([Plotting and tables](plotting.md)) |

The [explorer](../explorer/index.md) does these steps in a window, and the
[command line](cli.md) in a shell.

Not in this package: the Filter, Histogram and Autocorrelation functions,
the detection of frequency jumps, the dead-time corrections and the scaling
of a result for one of two identical units. For the last, divide `r.dev` by
``\sqrt{2}``.

## [Drift types](@id stable32-drift)

| Drift Type | SigmaTau |
|:---|:---|
| Quadratic Fit, phase data | `drift_fit(p; method = :phase_fit)` |
| Avg of 2nd Diffs, phase data | `drift_fit(p; method = :second_difference)` |
| 3-Point Fit, phase data | `drift_fit(p; method = :bisection)` |
| Linear Fit, frequency data | `drift_fit(f)` |
| Bisection Fit, frequency data | `drift_fit(f; method = :bisection)` |
| Log Fit, Diffusion Fit | `drift_fit(f; model = :log)`, `model = :diffusion` |
| Linear, Avg of 1st Diffs, Endpoints: the frequency offset of phase data | `detrend(p; method = :linear)` removes the first and `method = :endpoint` the other two, which are one number, the mean of `FrequencyData(p).y` |
| Autoregression | not in this package |
| Remove Drift | `detrend(p, fit)` with the result of [`drift_fit`](@ref), or `detrend(p; method = :quadratic)` for the quadratic fit |

Every method takes a record of either kind. [`drift_fit`](@ref) also has
a quadratic model of the frequency and the four-point estimator of
[Greenhall 1997](@cite greenhall-1997-drift-estimator). Its slopes are per
second and its drift rate per day, where the Drift function gives a slope
per sample interval. On the 1000-point test suite of
[SP1065 §12.4](@cite riley-2008-sp1065) the linear, bisection and
first-difference slopes agree with the published values to their seven
digits.

For a record with an even number of phase samples, which has no middle
sample, `:bisection` takes the mean frequencies of the two halves of the
frequency record and leaves its middle sample out. The manual does not say
what the 3-Point Fit does there.

## Variance types

| Variance type of Sigma and Run | SigmaTau |
|:---|:---|
| Normal (Allan) | `adev(p; overlapping = false)` |
| Overlapping (Allan) | [`adev`](@ref) |
| Modified | [`mdev`](@ref) |
| Time | [`tdev`](@ref) |
| Hadamard | `hdev(p; overlapping = false)` |
| Overlapping Hadamard | [`hdev`](@ref) |
| Total | [`totdev`](@ref) |
| Modified Total | [`mtotdev`](@ref) |
| Time Total | [`ttotdev`](@ref) |
| Hadamard Total | [`htotdev`](@ref) |
| TIE rms, MTIE | [`tierms`](@ref), [`mtie`](@ref) |
| Thêo1, ThêoBR | [`theo1`](@ref) |
| ThêoH | [`theoh`](@ref) |

The modified Hadamard deviations ([`mhdev`](@ref), [`htdev`](@ref),
[`mhtotdev`](@ref)) and the parabolic deviation ([`pdev`](@ref)) have no
counterpart in Stable32.

| Control of Run | SigmaTau |
|:---|:---|
| Alpha, Auto or a value | `noise = :auto`, or a noise type such as `:WHFM` |
| Conf Factor | `confidence`, 0.683 by default |
| No Drift | `detrend(p; method = :quadratic)` before the call, or `:linear` for a frequency record |
| Octave | `Octave`, the default |
| Decade, which is 1, 2, 4, 10, 20, 40, … | the factors themselves, `adev(p, [1, 2, 4, 10, 20, 40])`. `Decade` here is 1, 10, 100, … |
| All Tau | `AllTaus` |
| Start AF | the factors themselves |

## [Where the results differ](@id stable32-differences)

### Noise type

The noise type at each averaging time sets the EDF, and with it the
interval, and for the total deviations the bias factor. Both programs use
the lag-1 autocorrelation while at least 30 points remain and the ``B_1``
ratio below that. They differ in these steps, for the reasons that
[Noise identification](../explanation/noise_identification.md) gives:

| Step | Stable32 | SigmaTau |
|:---|:---|:---|
| Before the lag-1 test | not described in its manual | a least-squares quadratic is removed from each decimated series |
| A PM result at ``τ ≥ 4τ_0`` | stands | white PM is checked against ``R(n)``, and from ``16τ_0`` any PM result stands only if ``R(n)`` agrees |
| Limits between the FM types at ``τ ≥ 2τ_0`` | the estimate is rounded to the nearest ``α`` | halfway between the values the estimate takes on decimated phase |
| ``B_1`` ratio, fewer than 30 averaged points | labels as with ``B_1(N, μ)`` at the length of the record | ``B_1(N, μ)`` at the number of averaged points, and below 6 points the label found at 6 |
| PM with fewer than 30 averaged points | the result of the ``B_1`` ratio and ``R(n)`` | where the result at the longest ``τ`` that has 30 points is PM and ``R(n)`` at the ``τ`` itself agrees, and not otherwise |
| Hadamard family | on the reference record, the noise type of the other overlapping estimators | its own, found with up to three differences on the record less its frequency drift, and without the ``B_1`` ratio: below 30 averaged points an FM type is the one found at the longest ``τ`` that has 30 |

`identify_noise(x, m; remove_drift = false, b1_length = :record)` undoes the
first, fourth and fifth steps, and its labels go to any estimator as
`noise`. `remove_drift = false` on a Hadamard estimator undoes the last. No
option undoes the second and third.

The labels differ in two places. One is the longest averaging times, where
few averaged points remain. The other is flicker PM beyond about
``16τ_0``: on the flicker PM record from Stable32's own noise generator
that is among the test fixtures, Stable32 labels 32 s, 64 s and 128 s
white PM and 256 s white FM, where this package reads flicker PM. Its
interval for the Allan deviation at 32 s is half as wide as the one here.

### Bias

| Estimator | Stable32 | SigmaTau by default |
|:---|:---|:---|
| Total, Hadamard total, on FM noise | corrected | the same factor, except that the Hadamard total takes the exact factor below ``16τ_0`` |
| Total, Hadamard total, on PM noise | not corrected | corrected ([Total estimators on PM noise](../explanation/total_pm.md)) |
| Modified total, time total | printed without correction | corrected |
| Thêo1 | the factor of the noise type set, or ThêoBR with Auto | ThêoBR |
| ThêoBR, ThêoH | the ratio comes from at most 3000 samples of the record | from the whole record |

`correct_bias = :stable32` on an estimator, on [`deviation`](@ref) or on
[`stability`](@ref) applies Stable32's convention to the total deviations
and to Thêo1. The command line takes `--bias stable32`, and the explorer
**As Stable32** under **Bias correction**.

For [`theo1`](@ref) the setting gives the ThêoBR of a Stable32 run with
automatic noise identification. It agrees with Stable32 to the five digits
that Stable32 prints on the two records compared
([Validation](../explanation/validation.md)). The default averages the
same ratio over the same range, up to half the record, from every sample,
where Stable32 takes 3000 at most. On the record of 1001 samples the two
agree within 0.1 %, and on the record of 8192 the default reads 1.8 %
lower. `correct_bias = :published` averages the ratio up to a tenth of the
record only and can differ by far more
([Bias correction](../explanation/bias.md#Thêo1)).
`correct_bias = false` gives Stable32's Thêo1 with white FM set. The
factor that Stable32 applies for another noise type set by hand has no
counterpart here.

[`theoh`](@ref) takes the default ratio. For the long averaging times of a
Stable32 ThêoH run on more than 3000 samples, call `theo1` with
`correct_bias = :stable32`.

### Gaps

Stable32 skips gaps in the Allan, Hadamard, modified total, time total and
Hadamard total deviations, Thêo1, ThêoH and TIE rms, and fills them itself
for the modified Allan, time and total deviations and MTIE. This package
never fills on its own. The Allan and Hadamard families, MTIE and TIE rms skip
gaps, and the other estimators stop with an error until the record is
filled ([Gaps, outliers and drift](preparing.md)). Its [`fillgaps`](@ref)
reflects the neighboring samples into a gap, where the Fill function
interpolates linearly, so filled records differ.

### Smaller points

- The longest averaging time of a default grid is each estimator's own
  limit and can differ from Stable32's by a step: on a record of 8192
  points Stable32's total deviation ends at 4096 s and `totdev(p)` at
  2048 s. Give the factors to get a row.
- TIE rms, MTIE and Thêo1 take a frequency record here, and integrate it.
  Stable32 computes them for phase data only.
- The intervals are two-sided. Stable32 also offers single-sided limits.

## Reproducing a Stable32 table

State Stable32's noise type at each averaging time, from the Alpha column
of its table, and its bias convention. The reference record of the
[Validation](../explanation/validation.md) page comes with Stable32's
tables. A default call for its total deviation at 1024 s:

```@example s32
using SigmaTau
file = joinpath(pkgdir(SigmaTau), "test", "fixtures", "validation",
                "stable32gen.DAT")
p = read_phase(file)
totdev(p, [1024])
```

Stable32 prints 2.8694e-03 for this row, with the limits 2.4276e-03 and
3.6939e-03 and ``α = 0``, white FM:

```@example s32
totdev(p, [1024]; noise = :WHFM, correct_bias = :stable32)
```

The deviation is Stable32's to the digits it prints, and the limits are
within 0.03 % of its own. `noise` also takes one type for each averaging
factor, so a whole table is one call. For the Allan and Hadamard families
the noise type is all there is to state.

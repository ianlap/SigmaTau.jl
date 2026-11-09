# Choosing an estimator

Start with the question you are asking of the clock, then pick the estimator
that answers it with the least bias and the narrowest interval. When in
doubt, run the default suite, `stability(p)` (ADEV, MDEV, HDEV and MHDEV), and
compare the curves.

## Quick guide

| Your situation | Use | Why |
|:---|:---|:---|
| A first look, or a number for a spec sheet | `adev` | The standard measure σ\_y(τ); every datasheet quotes it. |
| Short τ dominated by phase noise (white or flicker PM) | `mdev` | ADEV falls as τ⁻¹ for both white and flicker PM; MDEV separates them (τ⁻³ᐟ² vs τ⁻¹). |
| The clock drifts (rubidium, quartz, uncorrected hydrogen maser) | `hdev`, `mhdev` | The third difference removes linear frequency drift, which turns ADEV into a τ¹ ramp. |
| Long τ, where few averages remain and the intervals are wide | `totdev`, `mtotdev`, `htotdev`, `mhtotdev` | Total estimators extend the record by reflection, so the last τ values carry more degrees of freedom. |
| The longest τ the record can support | `theo1` (ThêoBR), `theoh` | Thêo1 reaches 75 % of the record span; ThêoH joins it to ADEV at short τ. Slow on long records. |
| Time error of a clock or a timing link, in seconds | `tdev`, `htdev`, `ttotdev` | σ\_x(τ): the time deviation scaled from MDEV, MHDEV or MTOTDEV. |
| Telecom timing masks (ITU-T G.810 family) | `mtie`, `tierms` | Maximum and RMS time interval error over each window. |
| Data from an averaging (Λ-type) frequency counter | `pdev` | The parabolic deviation is the matched estimator for averaged frequency readings. |
| Does the stability change during the record? | `dadev`, `dhdev` | σ(t, τ) over sliding windows; see [Plotting](plotting.md#Dynamic-deviations). |
| Several clocks compared in pairs, and you want each clock alone | `ncornered_hat` | The N-cornered hat separates the pairwise variances, assuming independent clocks. |
| Live data, updated sample by sample | `StreamingStability`, `LiveStability` | ADEV, MDEV, HDEV and MHDEV at O(1) cost per sample; equals the batch call at every step. `LiveStability` takes instrument readings (phase, frequency or Hz) from an `InstrumentSource`. |

## The estimators

`estimators()` returns the same list as data (name, label, family, axis and a
one-line description), and `deviation(:mdev, p, m)` runs any of them by name.

### σ\_y(τ): fractional-frequency stability

| Function | Name | Note |
|:---|:---|:---|
| [`adev`](@ref) | Allan deviation (overlapping) | The standard. Cannot tell white PM from flicker PM; diverges under drift. |
| [`mdev`](@ref) | Modified Allan deviation | Averages phase over each τ, which separates white and flicker PM. |
| [`hdev`](@ref) | Hadamard deviation (overlapping) | Insensitive to linear frequency drift. |
| [`mhdev`](@ref) | Modified Hadamard deviation | Drift-insensitive and separates white and flicker PM. |
| [`pdev`](@ref) | Parabolic deviation | Least-squares frequency estimate per window; matched to averaging counters. |
| [`totdev`](@ref) | Total Allan deviation | ADEV with a reflected record: more EDF at long τ. |
| [`mtotdev`](@ref) | Total modified Allan deviation | MDEV with per-window drift removal and reflection. |
| [`htotdev`](@ref) | Total Hadamard deviation | HDEV with reflection: drift-insensitive with good long-τ confidence. |
| [`mhtotdev`](@ref) | Total modified Hadamard deviation | MHDEV with reflection. New in SigmaTau; see [MHTOTDEV bias and EDF](../explanation/mhtotdev.md). |
| [`theo1`](@ref) | Thêo1, bias-corrected (ThêoBR) by default | Reaches τ up to 75 % of the span. O(N·m) per τ, so seconds on long records. |
| [`theoh`](@ref) | ThêoH | ADEV at short τ joined to ThêoBR at long τ: one curve over the full range. |

### σ\_x(τ): time stability, in seconds

| Function | Name | Note |
|:---|:---|:---|
| [`tdev`](@ref) | Time deviation | τ/√3 · MDEV. The usual measure for time transfer and network timing. |
| [`htdev`](@ref) | Hadamard time deviation | τ/√(10/3) · MHDEV: TDEV with drift removed. New in SigmaTau. |
| [`ttotdev`](@ref) | Total time deviation | τ/√3 · MTOTDEV. |
| [`mtie`](@ref) | Maximum time interval error | Largest peak-to-peak time error in any window of length τ. A deterministic envelope, so it has no confidence interval. |
| [`tierms`](@ref) | RMS time interval error | RMS of the time error over windows of length τ. |

### Time-resolved and multi-clock

| Function | Returns |
|:---|:---|
| [`dadev`](@ref), [`dhdev`](@ref) | a [`DynamicStabilityResult`](@ref): ADEV or HDEV in each window of a sliding window |
| [`ncornered_hat`](@ref) | one [`StabilityResult`](@ref) per clock from the pairwise results |
| [`StreamingStability`](@ref) | an accumulator; [`snapshot`](@ref) gives the current [`StabilityResult`](@ref) |

## Reading the slope

On a log-log σ\_y(τ) plot each power-law noise has its own slope, which is how
you read the noise mix off a curve (and how SigmaTau's noise identification
works; see [Noise, EDF and confidence](../explanation/noise_and_confidence.md)):

| Noise | ADEV slope | MDEV slope |
|:---|:---|:---|
| White PM | τ⁻¹ | τ⁻³ᐟ² |
| Flicker PM | τ⁻¹ | τ⁻¹ |
| White FM | τ⁻¹ᐟ² | τ⁻¹ᐟ² |
| Flicker FM | τ⁰ (flat) | τ⁰ |
| Random-walk FM | τ¹ᐟ² | τ¹ᐟ² |
| Linear frequency drift | τ¹ | τ¹ |

HDEV and MHDEV follow the ADEV and MDEV slopes for these noises, but a
linear drift does not appear in them at all.

## Practical notes

- **Drift.** If ADEV turns up as τ¹ at long τ, either use the Hadamard family
  or remove the drift first with [`detrend`](@ref)`(p; method=:quadratic)`
  (on phase, the quadratic term is the frequency drift).
- **Bias correction.** The total estimators and Thêo1 are bias-corrected by
  default. Pass `correct_bias=false` for the raw statistic, which is what
  Stable32 reports for MTOT and Thêo1.
- **Speed.** The Allan and Hadamard estimators are O(N) per τ. The total
  family is heavier but multithreaded (start Julia with `--threads=auto`).
  Thêo1 and ThêoH are O(N²) overall, so they are off by default in the
  explorer.
- **Gaps.** `adev`, `mdev`, `tdev`, `hdev`, `mhdev`, `htdev`, `tierms` and
  `mtie` skip `NaN` gaps, as allantools' `gradev` does for ADEV: on phase, a
  term counts when every sample it reads was measured; on frequency, when no
  missing sample falls inside its span. `neff` is the number of terms used
  and the EDF follows from it. The total estimators, Thêo1, PDEV and the
  dynamic deviations still need a gap-free record: fill it with
  [`fillgaps`](@ref) first. Filling changes the noise inside the gap, so
  prefer skipping where the estimator allows it.
- **Glitches.** Remove phase steps or frequency spikes with
  [`remove_outliers`](@ref) before the analysis; see
  [Reading and saving data](../tutorials/02_reading_your_data.md).

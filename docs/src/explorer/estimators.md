# Estimators and noise types

This page is a short reference for the names the explorer shows. The
Background pages give the definitions and their sources, and
[Choosing an estimator](../guide/choosing.md) says which estimator answers
which question.

## Estimators

Every estimator works on the phase of the record, at averaging times
``τ = mτ_0`` from the τ grid of the [side panel](panel.md#Options). A
frequency record is integrated to phase first.

The frequency estimators give ``σ_y(τ)``, without unit:

| Estimator | Name | What it gives |
|:----------|:-----|:--------------|
| ADEV | Overlapping Allan deviation, [`adev`](@ref) | The standard measure of frequency stability |
| MDEV | Modified Allan deviation, [`mdev`](@ref) | Separates white PM from flicker PM, which ADEV does not |
| HDEV | Overlapping Hadamard deviation, [`hdev`](@ref) | Does not respond to linear frequency drift |
| MHDEV | Modified Hadamard deviation, [`mhdev`](@ref) | Separates the two phase noises and does not respond to drift |
| PDEV | Parabolic deviation, [`pdev`](@ref) | The deviation of the readings of an Ω counter. It separates the two phase noises as MDEV does |
| TOTDEV | Total deviation, [`totdev`](@ref) | ADEV with more degrees of freedom at long averaging times |
| MTOTDEV | Modified total deviation, [`mtotdev`](@ref) | The total form of MDEV |
| HTOTDEV | Hadamard total deviation, [`htotdev`](@ref) | The total form of HDEV |
| MHTOTDEV | Modified Hadamard total deviation, [`mhtotdev`](@ref) | The total form of MHDEV |
| THEO1 | Thêo1, [`theo1`](@ref) | An estimate of ADEV out to three quarters of the record |
| THEOH | ThêoH, [`theoh`](@ref) | ADEV at short averaging times joined to Thêo1, with its bias removed, at long ones |

The time estimators give ``σ_x(τ)``, in seconds:

| Estimator | Name | What it gives |
|:----------|:-----|:--------------|
| TDEV | Time deviation, [`tdev`](@ref) | MDEV times ``τ/\sqrt{3}`` |
| HTDEV | Hadamard time deviation, [`htdev`](@ref) | MHDEV times ``τ/\sqrt{10/3}``, a statistic defined in this package |
| TTOTDEV | Time total deviation, [`ttotdev`](@ref) | MTOTDEV times ``τ/\sqrt{3}`` |
| MTIE | Maximum time interval error, [`mtie`](@ref) | The largest peak-to-peak time error within any interval of length ``τ`` |
| TIE rms | Rms time interval error, [`tierms`](@ref) | The rms change of the time error over ``τ`` |

ADEV, MDEV, HDEV, MHDEV, TDEV, HTDEV, MTIE and TIE rms leave out the terms
that touch a missing sample. The other estimators need a record without
gaps.

MTOTDEV, HTOTDEV, MHTOTDEV, TTOTDEV and the two Thêo estimators take
seconds on a record of ``10^5`` samples. The others take a fraction of a
second.

## Noise types

Clock noise is modeled as a sum of power laws,
``S_y(f) = \sum_α h_α f^α``. The label at an averaging time names the term
that dominates there. PM stands for phase modulation and FM for frequency
modulation:

| Label | In exports | Noise | ``α`` | ADEV | MDEV |
|:------|:-----------|:------|------:|:-----|:-----|
| WPM | `WHPM` | White PM | 2 | ``τ^{-1}`` | ``τ^{-3/2}`` |
| FPM | `FLPM` | Flicker PM | 1 | ``τ^{-1}``, nearly | ``τ^{-1}`` |
| WFM | `WHFM` | White FM | 0 | ``τ^{-1/2}`` | ``τ^{-1/2}`` |
| FFM | `FLFM` | Flicker FM | −1 | ``τ^{0}`` | ``τ^{0}`` |
| RWFM | `RWFM` | Random-walk FM | −2 | ``τ^{1/2}`` | ``τ^{1/2}`` |
| FWFM | `FWFM` | Flicker-walk FM | −3 | none | none |
| RRFM | `RRFM` | Random-run FM | −4 | none | none |

The last two columns give the slope of the deviation on logarithmic axes.
HDEV follows ADEV, MHDEV and PDEV follow MDEV, and the total and Thêo
estimators follow the deviation they extend. The Allan variances do not converge for
the last two noise types. The Hadamard family does, with the slopes ``τ^1``
and ``τ^{3/2}``, and these two labels appear for it alone.
[Power-law noise](../explanation/noise.md) has the full table.

The label comes from the
[noise identification](../explanation/noise_identification.md) of the run,
unless the option **Noise type** states it. It sets the equivalent degrees
of freedom, and with them the confidence interval, and it selects the bias
factor of the total and Thêo estimators.

## Confidence intervals

The interval of a deviation is that of a ``χ^2`` variable with the
equivalent degrees of freedom (EDF) of the estimate, at the level of the
option **Confidence intervals**. The EDF depends on the estimator, the
noise type, the averaging time and the length of the record, and it falls
toward long averaging times, where the intervals widen. The level of
68.3 % is that of one standard deviation of a normal variable.
[Confidence intervals](../explanation/confidence.md) gives the formulas and
says how well the intervals hold.

MTIE and TIE rms have no interval.

## Bias correction

The total estimators and Thêo1 reuse the record to gain degrees of freedom,
and each then differs from the deviation it stands for by a factor that
depends on the noise type. The option **Bias correction** chooses the
factors:

| Choice | Factors |
|:-------|:--------|
| **Exact** | Those of this package, for every noise type and averaging factor |
| **Published tables** | Those of SP1065 and of the papers that define the estimators |
| **As Stable32** | Those that Stable32 applies to each estimator |
| **None** | No correction, as allantools returns |

[Bias correction](../explanation/bias.md) lists the factors and their
sources, and shows what each choice does to each estimator.

# Noise, EDF and confidence

A deviation computed from a finite record is an estimate, and how far it can
be trusted depends on the kind of noise in the clock. This page explains what
SigmaTau does between the raw deviation and the interval it reports: τ grids,
noise identification, equivalent degrees of freedom, χ² intervals and bias
correction.

## Phase, frequency and τ

A clock record is either **phase** (time error) x(t) in seconds or
**fractional frequency** y(t), dimensionless, sampled every τ₀ seconds. They
are two views of one signal: y is the slope of x,
``y_k = (x_{k+1} - x_k)/\tau_0``. SigmaTau's estimators work on phase, so a
[`FrequencyData`](@ref) record is integrated first.

A deviation is evaluated at averaging times τ = m·τ₀, where the integer m is
the **averaging factor**. You choose the factors explicitly
(`adev(p, [1, 10, 100])`) or through a named grid, a [`TauMode`](@ref):
`Octave` (powers of two, the default), `HalfOctave`, `QuarterOctave`,
`Decade`, `HalfDecade` (steps of √10) or `AllTaus`. A grid stops at the
largest m the estimator supports on the record length:
[`tau_values`](@ref) shows the factors it gives.

## Power-law noise

Clock noise is modelled as a sum of power laws in the frequency spectrum,
``S_y(f) = \sum_\alpha h_\alpha f^\alpha``. Five values of α cover almost every oscillator:

| Label | Noise | α in S\_y(f) ∝ f^α | σ\_y²(τ) ∝ τ^μ | Typical source |
|:---|:---|:---|:---|:---|
| `:WHPM` | White phase modulation | 2 | μ = −2 | Measurement noise, short τ |
| `:FLPM` | Flicker phase modulation | 1 | μ = −2 | Amplifiers and electronics |
| `:WHFM` | White frequency modulation | 0 | μ = −1 | The floor of passive atomic standards |
| `:FLFM` | Flicker frequency modulation | −1 | μ = 0 | The flicker floor |
| `:RWFM` | Random-walk frequency modulation | −2 | μ = 1 | Environment, long-term wander |

The slope of σ\_y(τ) on a log-log plot is μ/2, which is why the noise mix can
be read off a stability plot; [Choosing an estimator](../guide/choosing.md#Reading-the-slope)
tabulates the slopes. [`noise_gen`](@ref) synthesizes records with any mix of
these noises, which is how the examples on this site are built.

### Synthesized noise

[`noise_gen`](@ref) draws each component with the Kasdin–Walter discrete
fractional filter: white Gaussian noise convolved with the impulse response
of ``(1 - z^{-1})^{\alpha/2}``. It also covers the two steeper noises that the
Hadamard-class estimators converge on, flicker-walk FM (`:FWFM`, α = −3) and
random-run FM (`:RRFM`, α = −4). A component requested by `h_α` has the
spectrum

```math
S_y(f) = h_\alpha \left(\frac{\sin \pi f \tau_0}{\pi \tau_0}\right)^{\alpha},
```

which is ``h_\alpha f^\alpha`` below about a tenth of the sampling rate and
rolls off as ``\sin`` toward Nyquist. For α ≤ −1 the process is not
stationary; its d-th difference is, with spectrum
``(2 \sin \pi f \tau_0)^{2d} S_y(f)``. The test suite generates 40
records per α from +2 to −4, differences where needed, and checks the
averaged Welch spectrum against this formula bin by bin across the band: the
band-average ratio is 1 to within 1 % and no bin is off by more than 12 %,
the scatter expected of the estimate.

## Noise identification

With `ci=true` (the default) every estimator except MTIE and TIE rms identifies the dominant noise
type at each τ ([`identify_noise`](@ref)), because the interval depends on it.
SigmaTau uses the lag-1 autocorrelation of the phase decimated by m, after a
quadratic detrend of each decimated series, which is reliable while at least
30 averaged samples remain. Beyond that it falls back to the B1 and R(n)
ratios (SP1065 §5.5.2 [riley-2008-sp1065](@cite)). B1 is the ratio of the
standard variance of the averaged frequency to its Allan variance, and it is
compared with the B1(N, μ) of each noise type at the number of averaged
points N, the length the ratio is defined for. With 8 to 24 points no test
separates neighbouring noise types better than about two times in three, so
labels at the longest τ are weak evidence; `b1_length=:record` gives the
labels closer to Stable32's, which almost never include random-walk FM
there. Below 6 averaged points the ratio says almost nothing, and SigmaTau
reports the label measured at the longest averaging factor that still has 6;
that value comes from simulation, not from the literature. Decimation whitens flicker PM: its lag-1 correlation fades as m
grows, and from m ≈ 16 the test reads it as white PM, whose ADEV has several
times the EDF. A white-PM result is therefore confirmed with R(n), the ratio
of MVAR to AVAR that SP1065 §5.5.2 uses to separate the two PM noises, and
becomes flicker PM when R(n) says so. The result is always one of the five labels above, or for the
Hadamard family one of seven (see below): estimates outside the range are
clamped to its ends, and a τ where the estimate fails takes the type of its
nearest neighbour.

Decimation also moves the FM noises. The estimator of
[riley-2004-lag1-acf](@cite) turns the lag-1 autocorrelation r₁ of the
differenced series into α and rounds to the nearest integer, which is exact
for discrete fractional noise. Phase decimated by m is instead a sampled
continuous process. The second differences of sampled flicker FM and
random-walk FM correlate at r₁ = (9 ln 3 − 16 ln 2)/(8 ln 2) = −0.217 and
r₁ = 1/4, from their generalized autocovariances t² ln|t| and |t|³
[greenhall-2003-edf-stability](@cite), and the estimator reads those as
α = −1.446 and −2.4 where discrete noise gives −1 and −2. Simulated records
reach those values by m ≈ 8. Rounding at −1.5 therefore splits flicker FM
almost evenly between its own label and random-walk FM. For m ≥ 2 SigmaTau
places the limits halfway between the readings 0, −1.446 and −2.4, at −0.72
and −1.92; at m = 1 the record is not decimated and the published rounding
stands. A measured record is itself a sampled process, so flicker FM that
dominates at τ₀ can read near −1.45 there too, close to the rounding limit,
where simulated discrete noise reads −1.

Stable32 does not detrend before the lag-1 estimate. The difference only
matters on records with drift; `identify_noise(...; remove_drift=false)` follows
Stable32, and [Validation](validation.md) quantifies the effect.

The Hadamard family (`hdev`, `mhdev`, `htdev`, `htotdev`, `mhtotdev`) goes one
step further, because it is the family used on drifting clocks. A frequency
drift inflates the B1 ratio, which then reads random-walk FM at long τ
whatever the noise, and the EDF and the HTOTDEV bias factor follow the wrong
label. These estimators therefore identify noise on the record less its
least-squares frequency drift (`remove_drift=true`), as SP1065 §5.5.6
[riley-2008-sp1065](@cite) has deterministic components removed before noise
identification; the deviation itself does not change. The removal also takes
some of the slowest noise with it, so on a record without drift flicker FM
near the longest τ is read as white FM more often. `remove_drift=false`
identifies the noise on the record as given, which is what Stable32 does.

The third difference of the Hadamard family also converges for flicker-walk
FM (`:FWFM`, α = −3) and random-run FM (`:RRFM`, α = −4), and these estimators
identify the two by letting the lag-1 estimator difference a third time
[riley-2004-lag1-acf](@cite). The third differences of sampled random-walk,
flicker-walk and random-run FM correlate at r₁ = −1/3,
(128 ln 2 − 81 ln 3)/(27 ln 3 − 32 ln 2) = −0.035 and 13/33, which the
estimator reads as α = −3.0, −3.93 and −4.57, so for m ≥ 2 the limits sit at
−3.46 and −4.25. HDEV, MHDEV and HTOTDEV then take their EDF at the identified
α [greenhall-2003-edf-stability](@cite) and HTOTDEV its bias factor
[howe-2001-total-estimators](@cite), and MHTOTDEV its own
[bias and EDF](mhtotdev.md) for each. Read as random-walk FM, random-run FM had an HDEV
EDF 26 % too high and an HTOTDEV 6 % too low. The price is paid by
random-walk FM at long τ: its second differences sit at ρ = 0.2, next to the
0.25 at which the estimator differences again, and in simulation it keeps
its label in 91 % of records with 64 averaged points left and 81 % with 32,
where two differences alone keep 98 % and 93 %. A random-walk FM record
labelled flicker-walk FM has an HDEV EDF 2 % lower and an HTOTDEV 3.7 %
higher.

## Equivalent degrees of freedom

An estimated variance from a record with Gaussian noise follows a scaled χ²
distribution. Its number of degrees of freedom, the **EDF**, measures how
much independent information went into the estimate, and it depends on the
estimator, the noise type, m and the record length N. It falls as τ grows,
because a longer τ leaves fewer independent averages; that is why the last
points of a stability plot carry wide intervals. SigmaTau computes it per τ
([`calculate_edf`](@ref)):

- **ADEV, MDEV, HDEV, MHDEV:** the Greenhall–Riley algorithm
  [greenhall-2003-edf-stability](@cite), which evaluates the EDF exactly from
  the noise spectrum rather than from fitted approximations.
- **Total family:** the published EDF coefficients of NIST SP1065
  [riley-2008-sp1065](@cite) and Howe and Peppler
  [howe-2001-total-estimators](@cite). Where none are published the model
  is fitted to the exact EDF of the statistic: for TOTDEV and HTOTDEV on
  white and flicker PM ([Total estimators on PM noise](total_pm.md)) and for
  MHTOTDEV, which is new in SigmaTau ([MHTOTDEV bias and EDF](mhtotdev.md)).
- **Thêo1 and ThêoH:** the empirical per-noise formulas of SP1065 Table 3.
- **PDEV:** the PVAR model of Vernotte, Chen and Rubiola.
- **MTIE and TIE rms** have no χ² model and carry no interval.

`r.edf` holds the values, and `r.neff` the number of terms actually averaged
at each τ.

## Confidence intervals

Given the deviation σ̂ and the EDF ν, the interval at confidence level p is

```math
\sigma_\mathrm{lo} = \hat\sigma \sqrt{\frac{\nu}{\chi^2_{(1+p)/2}(\nu)}}, \qquad
\sigma_\mathrm{hi} = \hat\sigma \sqrt{\frac{\nu}{\chi^2_{(1-p)/2}(\nu)}}
```

([`confidence_intervals`](@ref)). The default level is 68.3 %
(`DEFAULT_CONFIDENCE`), the ±1σ convention of Stable32; pass
`confidence=0.95` for a 95 % interval. The interval is asymmetric about σ̂,
wider above than below, because the χ² distribution is skewed at small ν.
The formula holds for any ν > 0, so a point whose EDF falls below 1 keeps
a χ² interval, a very wide one.
Where an EDF model gives no positive value (SP1065's Thêo1 formula for
random-walk FM beyond τ ≈ 0.63 T), the point has no interval and its bounds
are `NaN`.

## Bias correction

The total estimators and Thêo1 are biased for some noise types: their
expected value differs from the true Allan-family variance by a known factor
that depends on the noise and on τ/T. With `correct_bias=true` (the default)
SigmaTau divides it out using the identified noise type
([`bias_correction`](@ref)):

| Estimator | Correction |
|:---|:---|
| TOTDEV | Greenhall, Howe and Percival: reads low with τ/T on FLFM and RWFM; on PM its exact expectation, high with τ/T |
| MTOTDEV, TTOTDEV | Howe and Peppler: MTOT reads 0.3 % (white PM) to 18 % (random-walk FM) low |
| HTOTDEV | Howe and Peppler on the FM noises; on PM its exact expectation, 25 to 38 % (white) and 4 to 28 % (flicker) high in variance |
| MHTOTDEV | its exact expectation: 10 % (white PM) to 22 % (flicker FM) low in variance |
| Thêo1 | the ThêoBR ratio to ADEV, which needs no noise ID |

The MTOT factors are Table 2 of Howe and Peppler
[howe-2001-total-estimators](@cite): the expected Mod Totvar is 0.995, 0.851,
0.771, 0.717 and 0.679 of MVAR from white PM to random-walk FM, as first
simulated by Howe and Vernotte [howe-1999-modtotvar](@cite), and the test
suite checks the kernel against them. SP1065 Table 11 prints factors of 1.06
to 1.31 under a caption that reads as the same ratio; applied as printed they
lower a statistic that already reads low. At τ = τ₀ MTOT is exactly half of
MVAR for any record, and the correction is that factor of 2.

The published factors hold from τ = 16τ₀. Below it the expectation depends
on the averaging factor, and SigmaTau divides by its exact value at each m
from 2 to 15, computed like the PM factors
([Total estimators on PM noise](total_pm.md)). For MTOT, and for HTOT two
noise types down (white FM to random-run FM):

| Noise (MTOT) | m = 2 | 3 | 4 | 5 | 8 | from m = 16 |
|:---|---:|---:|---:|---:|---:|---:|
| white PM | 0.998 | 0.938 | 0.996 | 0.973 | 0.996 | 0.995 |
| flicker PM | 0.877 | 0.844 | 0.861 | 0.852 | 0.855 | 0.851 |
| white FM | 0.814 | 0.781 | 0.784 | 0.777 | 0.774 | 0.771 |
| flicker FM | 0.735 | 0.722 | 0.722 | 0.719 | 0.718 | 0.717 |
| random-walk FM | 0.687 | 0.681 | 0.681 | 0.680 | 0.679 | 0.679 |

The white PM row dips at odd m, by 6 % at m = 3, and the white FM row is 5 %
above its limit at m = 2. The two flicker rows are the mean of two flicker
models, which are within 2.4 % of it at m = 2 and 1 % from m = 5. MHTOTDEV
has its own table on [its page](mhtotdev.md).

The PM rows of TOTDEV and HTOTDEV have no published source and are derived
in [Total estimators on PM noise](total_pm.md).

`correct_bias=false` returns the raw statistic. Stable32 corrects TOTDEV and
HTOTDEV the same way on the FM noises, applies no factor on PM, and reports
MTOT and Thêo1 raw; the [validation page](validation.md) shows both against
Stable32.

## Drift

A linear frequency drift D adds a term Dτ/√2 to ADEV, which grows as τ¹ and
hides everything else at long τ. The Hadamard family takes a third difference
and does not see it. Otherwise remove it before the analysis:
[`detrend`](@ref) with `method=:quadratic` fits and removes the quadratic
phase term, which is the frequency drift. MHTOTDEV removes a quadratic from
every window before reflecting it, so drift never reaches its reflected
record.

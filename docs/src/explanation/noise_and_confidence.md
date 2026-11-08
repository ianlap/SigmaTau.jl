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
ratios. The result is always one of the five labels above: estimates outside
the range are clamped to WHPM or RWFM, and a τ where the estimate fails takes
the type of its nearest neighbour.

Stable32 does not detrend before the lag-1 estimate. The difference only
matters on records with drift; `identify_noise(...; remove_drift=false)` follows
Stable32, and [Validation](validation.md) quantifies the effect.

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
  [riley-2008-sp1065](@cite) and Howe and Tasset, with the Greenhall–Riley
  value where no coefficients exist (white and flicker PM for TOTDEV and
  HTOTDEV). For MHTOTDEV, which is new in SigmaTau, the coefficients were
  measured by Monte Carlo; see [MHTOTDEV bias and EDF](mhtotdev.md).
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
Where the EDF is below 1 (a single average at the longest τ), SigmaTau falls
back to a symmetric Gaussian interval floored at zero.

## Bias correction

The total estimators and Thêo1 are biased for some noise types: their
expected value differs from the true Allan-family variance by a known factor
that depends on the noise and on τ/T. With `correct_bias=true` (the default)
SigmaTau divides it out using the identified noise type
([`bias_correction`](@ref)):

| Estimator | Correction |
|:---|:---|
| TOTDEV | Howe and Walter: τ/T-dependent, FLFM and RWFM only |
| MTOTDEV, TTOTDEV | SP1065 Table 11 (MTOT overstates MVAR) |
| HTOTDEV | Howe and Tasset, FM noises |
| MHTOTDEV | measured by Monte Carlo, τ/T-dependent |
| Thêo1 | the ThêoBR ratio to ADEV, which needs no noise ID |

`correct_bias=false` returns the raw statistic. Stable32 corrects TOTDEV and
HTOTDEV the same way and reports MTOT and Thêo1 raw; the [validation page](validation.md)
shows both against Stable32.

## Drift

A linear frequency drift D adds a term Dτ/√2 to ADEV, which grows as τ¹ and
hides everything else at long τ. The Hadamard family takes a third difference
and does not see it. Otherwise remove it before the analysis:
[`detrend`](@ref) with `method=:quadratic` fits and removes the quadratic
phase term, which is the frequency drift. MHTOTDEV removes a quadratic from
every window before reflecting it, so drift never reaches its reflected
record.

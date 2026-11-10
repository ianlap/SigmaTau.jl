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
becomes flicker PM when R(n) says so.

R(n) also decides whether the noise is PM at all. With few averaged points
both tests mistake FM for PM often: in simulation the lag-1 test gives white
FM its own label in 91 % of records at 64 averaged points and in 75 % at 32,
and most of the rest come out as PM. That is the costly mistake, because the
EDF of PM noise follows the record length and not the number of averaged
points, so the interval of ADEV or HDEV comes out several times too narrow.
R(n) is 1/2, 27/40 and 33/40 for white, flicker and random-walk FM against
0.26 for flicker PM at m = 16 and less beyond, and it is computed from every
sample where the tests read one in m. From m = 16 a PM result therefore
stands only if R(n) lies nearer its flicker-PM value than 1/2; one that fails
is reported as white FM. On simulated white FM (N = 8192) this raises the
share of 68.3 % ADEV intervals that hold the true value from 56–64 % to
70–73 % at 4 to 32 averaged points, and it leaves PM records as they were.

With fewer than 30 averaged points the rule is stricter, because noise only
reddens with τ and the tests left are weak. A result is PM only if the
result at the longest τ that has 30 points is PM and R(n) at the τ itself
says so: from 6 averaged points it must lie below the midpoint of 1/m and
1/2, the value it takes where white PM and white FM hold equal shares of the
Allan variance, and with fewer it must be white PM outright. Everything else
takes the FM type of the B1 ratio. The B1 ratio alone takes white PM for FM
in one record in five at 16 averaged points and one in three at 8, which put
the FM bias factor of TOTDEV, MTOTDEV and HTOTDEV on PM noise (white PM up
to 33 %, 12 % and 19 % high in variance). In simulation (N = 1024) white PM
is now within 2 % for all three down to 2 averaged points. Two costs remain.
Flicker PM with fewer than 6 averaged points is reported as white FM, so its
interval is wide and the three total variances read 19, 10 and 25 % high at
4 points. And where PM gives way to FM inside the record no label is right
near the crossing: the 68.3 % interval of ADEV, HDEV and TOTDEV then covers
64–96 % at 16 to 4 points, that of MTOTDEV and HTOTDEV 58–86 %, as on main
(54–96 %).

The Hadamard family identifies the noise on the record less its frequency
drift, and removing a drift roughly halves the B1 ratio (0.96, 1.65, 3.37
for white, flicker and random-walk FM at 16 points, against the 1.0, 2.13
and 8.0 its limits expect), which also stops at random-walk FM. So there
the B1 ratio is not used: an FM type with fewer than 30 averaged points is
the one found at the longest τ that has 30. HTOTDEV on flicker FM and redder
noise then reads within 4 % of HDEV at 16 to 3 points, where it read 7–20 %
low in variance; a change from one FM type to another beyond that τ is not
seen, and `noise` states it.
To set the noise type yourself, see [Stating the noise type](@ref).

The result is always one of the five labels above, or for the
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

### Stating the noise type

The tests above read each τ on its own, and they assume power-law noise. On a
measured record a slow periodic term, a step or a stretch of wander can
mislead them where few averaged points remain, and the labels may then run
against the order a sum of power laws allows (whiter noise at a longer τ
than a redder one). Since the label sets the EDF and the bias factor, every
estimator that reports intervals takes `noise` to state it instead (Stable32
lets the user set the noise type too, and allantools' EDF functions take the
exponent α):

```julia
adev(p; noise = :WHFM)                         # one type at every τ
adev(p, [1, 16, 256]; noise = [:WHPM, :WHFM, :WHFM])
hdev(p; noise = powerlaw_fit(p))               # strongest term of a model at each τ
stability(p; noise = :WHFM)
```

A [`PowerLawFit`](@ref) gives, at each τ, the component with the most
variance in the estimator's own kernel, so its labels can only redden as τ
grows. It is as good as the fit: a short record may not resolve the reddest
term. A type redder than an estimator covers counts as its reddest, as in
identification (random-walk FM for the Allan family). The command line has
`--noise auto|fit|TYPE`, the explorer a Noise type menu under Options, and
the Python functions the same `noise` argument with the names as strings.

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

### Matching a reference

The factors SigmaTau derives itself (PM, m < 16, MHTOTDEV) are not in SP1065,
Stable32 or allantools, so `correct_bias` also takes the convention of each:

| `correct_bias` | TOTDEV | MTOTDEV, TTOTDEV | HTOTDEV | Thêo1 | Reproduces |
|:---|:---|:---|:---|:---|:---|
| `true` (default) | all noises | exact below 16τ₀ | all noises, exact below 16τ₀ | ThêoBR, ladder to N/6 (SP1065, IEEE 1139) | |
| `:published` | FM only | published factor at every τ | FM only, published factor at every τ | ThêoBR, ladder to N/30 (Howe 2006) | SP1065, Howe and Peppler |
| `:stable32` | FM only | raw | FM only, published factor at every τ | ThêoBR to N/6 on at most 3000 points | Stable32 |
| `false` | raw | raw | raw | raw | allantools |

MHTOTDEV has no published factor and exists in neither program: it is raw
with `false` and corrected otherwise. The setting changes σ and the interval
that scales with it, not the noise identification or the EDF; where SigmaTau
and Stable32 read a different noise type at a τ, the FM factor still differs
([Coming from Stable32](../guide/from_stable32.md)). The
[validation page](validation.md) shows the default and the raw statistic
against Stable32. The command line names the four values
`--bias exact|published|stable32|none`, and the explorer offers them under
Bias correction.

The two ThêoBR ladders fail in different places. ThêoBR scales Thêo1 at
every τ by one factor, the mean ratio of AVAR to Thêo1 over a ladder of
matched averaging times. Howe's ladder stops at a tenth of the record, so on
a record whose noise changes from PM to FM beyond that point the factor is
the PM one and Thêo1 at long τ reads low: 0.18 to 0.77 of the true variance
in simulation, and 2.4 times low in deviation on the 8192-point validation
record. The ladder of SP1065 and IEEE Std 1139 runs to half the record and
reads 0.56 to 1.12 there, at a price on records of one FM noise: its last
rungs hold one or two AVAR samples each, the factor scatters 2.5 times more
and averages 7 to 9 % low, so bias-removed Thêo1 at short τ carries 15 to
19 % error in σ where ADEV has 4 to 11 %. [`theoh`](@ref) uses ADEV at short
τ and is little affected.

## Drift

A linear frequency drift D adds a term Dτ/√2 to ADEV, which grows as τ¹ and
hides everything else at long τ. The Hadamard family takes a third difference
and does not see it. Otherwise remove it before the analysis:
[`detrend`](@ref) with `method=:quadratic` fits and removes the quadratic
phase term, which is the frequency drift. MHTOTDEV removes a quadratic from
every window before reflecting it, so drift never reaches its reflected
record.

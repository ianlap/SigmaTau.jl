# Conventions

The quantities, symbols and units used by every function of the package and
every page of this documentation. They follow
[SP1065](@cite riley-2008-sp1065) and
[IEEE Std 1139-2022](@cite ieee1139-2022-definitions).

## Phase and frequency

A signal of nominal frequency ``ν_0`` with phase deviation ``φ(t)`` in radians
has the *time deviation* ``x(t)``, in seconds, and the *fractional frequency*
``y(t)``, which is dimensionless
([SP1065 §3.1](@cite riley-2008-sp1065)):

```math
x(t) = \frac{φ(t)}{2πν_0}, \qquad
y(t) = \frac{ν(t) - ν_0}{ν_0} = \frac{dx}{dt} .
```

This documentation calls ``x(t)`` the *phase*. A measurement gives it as the
reading of the clock under test less the reading of the reference. The other
sign changes no deviation and no spectrum, only the sign of a fitted drift.

A record holds samples at a uniform interval ``τ_0``:

| Record | Samples | Unit |
|:---|:---|:---|
| [`PhaseData`](@ref) | ``x_i = x(t_1 + (i-1)τ_0)``, ``i = 1, …, N`` | seconds |
| [`FrequencyData`](@ref) | ``y_i = (x_{i+1} - x_i)/τ_0``, ``i = 1, …, N-1`` | none |

A frequency sample is the average of ``y(t)`` over its interval, with no dead
time between samples, so ``N`` phase samples and ``N - 1`` frequency samples
describe the same record. Every estimator works on phase and integrates a
frequency record first, from ``x_1 = 0``. A missing sample is `NaN`
([Timestamps and gaps](@ref)).

`FrequencyData(p)` and `PhaseData(f)` convert a record to the other kind.
The estimators do not need it.

## Averaging time

A deviation is evaluated at averaging times ``τ = mτ_0``, where the integer
``m`` is the *averaging factor*. Every estimator takes the factors in one of
three forms:

```julia
adev(p)                # octave spacing, m = 1, 2, 4, …
adev(p, Decade)        # a named grid, a TauMode
adev(p, [1, 10, 100])  # the factors themselves
```

A named grid stops at the largest factor the estimator supports on the
record. For ``N`` phase samples that is

| Estimators | Largest ``m`` |
|:---|:---|
| `adev`, `pdev` | ``⌊(N-2)/2⌋`` |
| `totdev` | ``⌊(N-1)/2⌋`` |
| `mdev`, `tdev`, `htotdev` | ``⌊(N-1)/3⌋`` |
| `hdev` | ``⌊(N-2)/3⌋`` |
| `mtotdev`, `ttotdev` | ``⌊N/3⌋`` |
| `mhdev`, `htdev` | ``⌊(N-1)/4⌋`` |
| `mhtotdev` | ``⌊N/4⌋`` |
| `mtie`, `tierms` | ``N-1`` |
| `theo1` | the largest even ``m ≤ N-1`` |

The Allan, Hadamard and parabolic estimators keep at least two terms at
their longest ``τ``; the total estimators keep one subsequence, which the
extension turns into a usable estimate. [`tau_values`](@ref) returns the
factors of a grid, and [`TauMode`](@ref) lists the grids.

Thêo1 is the one estimator whose averaging time is not ``mτ_0``: a result at
factor ``m`` is reported at ``τ = 0.75\,mτ_0``
([SP1065 §5.2.15](@cite riley-2008-sp1065)).

## What an estimator returns

Each estimator returns a [`StabilityResult`](@ref) with one entry per
``τ``:

| Field | Content |
|:---|:---|
| `tau` | averaging time, in seconds |
| `dev` | the deviation, the square root of the variance estimate |
| `noise_type` | the power-law noise type used for the interval and the bias factor |
| `ci` | the confidence interval `(lo, hi)`, at 68.3 % unless `confidence` says otherwise |
| `edf` | the equivalent degrees of freedom behind the interval |
| `neff` | the number of terms averaged |

Deviations of fractional frequency, written ``σ_y(τ)``, are dimensionless.
The time deviations `tdev`, `htdev` and `ttotdev`, written ``σ_x(τ)``, and
`mtie` and `tierms` are in seconds.

[`adev`](@ref) and [`hdev`](@ref) are the overlapping estimators. The
standard recommends the overlapping form of the Allan deviation
([IEEE Std 1139-2022, A.3.2](@cite ieee1139-2022-definitions));
`overlapping=false` gives the non-overlapped forms.

## Noise types

Noise is modeled as a sum of power laws,
``S_y(f) = \sum_α h_α f^α``, and each term has a label:

| Label | Noise | ``α`` |
|:---|:---|---:|
| `:WHPM` | white phase modulation (PM) | 2 |
| `:FLPM` | flicker PM | 1 |
| `:WHFM` | white frequency modulation (FM) | 0 |
| `:FLFM` | flicker FM | −1 |
| `:RWFM` | random-walk FM | −2 |
| `:FWFM` | flicker-walk FM | −3 |
| `:RRFM` | random-run FM | −4 |

The Allan variance converges down to ``α = -2``. The last two labels exist
for the Hadamard family, whose third difference converges down to
``α = -4`` ([SP1065 §3.4](@cite riley-2008-sp1065)).
[Power-law noise](noise.md) gives the slopes and the conversions.

## Spectra

Spectral densities are one-sided, over the Fourier frequency ``f`` in hertz
([IEEE Std 1139-2022, 3.1 and A.2](@cite ieee1139-2022-definitions)):

| Quantity | Symbol | Unit | Relation |
|:---|:---|:---|:---|
| fractional frequency | ``S_y(f)`` | 1/Hz | |
| time deviation | ``S_x(f)`` | s²/Hz | ``S_x = S_y/(2πf)^2`` |
| phase | ``S_φ(f)`` | rad²/Hz | ``S_φ = (2πν_0)^2 S_x`` |
| single-sideband phase noise | ``ℒ(f)`` | dBc/Hz | ``ℒ = S_φ/2``, quoted as ``10\log_{10}ℒ`` |

## Symbols

| Symbol | Meaning |
|:---|:---|
| ``τ_0`` | sample interval, in seconds |
| ``N`` | number of phase samples |
| ``T = (N-1)τ_0`` | duration of the record |
| ``m`` | averaging factor |
| ``τ = mτ_0`` | averaging time |
| ``σ_y(τ)``, ``\mathrm{Mod}\,σ_y(τ)`` | Allan and modified Allan deviation |
| ``σ_x(τ)`` | time deviation |
| ``α`` | exponent of a power-law term of ``S_y(f)`` |
| ``h_α`` | level of that term |
| ``ν`` | equivalent degrees of freedom (EDF) |
| ``B`` | bias factor: the expected value of a variance estimate over that of its reference |
| ``f_h`` | bandwidth of the measurement: the Nyquist frequency ``1/(2τ_0)`` of a sampled record |

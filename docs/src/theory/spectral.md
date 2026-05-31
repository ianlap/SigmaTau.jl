# Theory: Spectral densities

The deviation families ([Allan](allan_family.md), [Total](total_family.md))
describe stability in the *time* domain, as a function of averaging time `τ`.
The same noise is equivalently described in the *frequency* domain by a power
spectral density (PSD). SigmaTau provides the three standard frequency-domain
measures of IEEE Std 1139-2022 [@cite ieee1139-2022-definitions]: the
fractional-frequency PSD `S_y(f)`, the phase PSD `S_x(f)`, and single-sideband
phase noise `ℒ(f)`.

## Definitions

For a clock with phase residual `x(t)` (seconds) and fractional frequency
`y(t) = ẋ(t)`, the one-sided PSDs are related by the time-derivative transfer
function `y = ẋ ⇒ S_y(f) = (2πf)² S_x(f)`:

```math
S_x(f) \quad [\mathrm{s^2/Hz}], \qquad
S_y(f) = (2\pi f)^2\, S_x(f) \quad [1/\mathrm{Hz}].
```

Single-sideband phase noise converts the phase record to radians at the carrier
frequency `ν₀` and reports half the phase PSD on a decibel scale
(IEEE 1139-2022 §3.5):

```math
S_\phi(f) = (2\pi \nu_0)^2\, S_x(f) \quad [\mathrm{rad^2/Hz}], \qquad
\mathscr{L}(f) = \tfrac{1}{2} S_\phi(f) \quad [\mathrm{dBc/Hz}].
```

## Power-law model

Oscillator noise is a sum of power laws. In the fractional-frequency PSD the
exponent is the SP1065 [@cite riley-2008-sp1065] index `α`:

```math
S_y(f) = \sum_{\alpha=-2}^{2} h_\alpha\, f^{\alpha}.
```

| α  | Noise type                | `S_y(f)` | `S_x(f) = S_y/(2πf)²` |
|---:|:--------------------------|:---------|:----------------------|
|  2 | white phase (WPM)         | `f²`     | `f⁰` (flat)           |
|  1 | flicker phase (FPM)       | `f¹`     | `f⁻¹`                 |
|  0 | white frequency (WFM)     | `f⁰`     | `f⁻²`                 |
| −1 | flicker frequency (FFM)   | `f⁻¹`    | `f⁻³`                 |
| −2 | random-walk freq. (RWFM)  | `f⁻²`    | `f⁻⁴`                 |

The `h_α` coefficients are the same ones [`noise_gen`](../reference/stab.md#SigmaTau.noise_gen)
takes as input, which makes synthesis and spectral estimation an exact
round-trip: synthesize a known mixture, estimate the PSD, recover the slope.

## Estimator: Welch's method

`Sy`, `Sx`, and `L` estimate the PSD by Welch's method. The record is split
into overlapping segments of length `nperseg`, each segment is mean-removed and
multiplied by a window (Hann by default), a periodogram is taken with a real
FFT, and the periodograms are averaged. Overlapping segments trade frequency
resolution (set by `nperseg`) for lower variance (set by the segment count).

The normalization is the one-sided "density" convention: power is folded onto
positive frequencies and scaled by the window power and sample rate, so the
estimate is variance-preserving,

```math
\sum_k S_y(f_k)\,\Delta f \;\approx\; \operatorname{var}(y)
```

over the analysis band `[0, f_s/2]` with `f_s = 1/τ₀`. Longer `nperseg` reaches
lower frequencies (`f_min = f_s / nperseg`) at the cost of fewer segments and a
noisier estimate; the defaults (`nperseg = min(N, 256)`, 50 % overlap) follow
the common scipy convention.

!!! note "Discrete S_x ↔ S_y identity"
    `Sx` runs Welch on the phase record and `Sy` on the differenced frequency
    record, so the continuous identity `S_x = S_y/(2πf)²` holds only
    approximately on sampled data — closely at low frequency, where the
    discrete integrator matches the continuous one, and degrading toward the
    Nyquist frequency.

## Usage

```julia
using SigmaTau

p = read_phase("clock.DAT"; time_col=0, value_col=1, tau0=1.0)

sy = Sy(p)                       # fractional-frequency PSD, 1/Hz
sx = Sx(p)                       # phase PSD, s²/Hz
l  = L(p; f_carrier=10e6)        # SSB phase noise at a 10 MHz carrier, dBc/Hz

sy.freq        # one-sided frequency grid (Hz)
sy.psd         # S_y(f)
l.psd          # ℒ(f) in dBc/Hz

# Tune the Welch trade-off: a longer segment reaches lower frequencies.
Sy(p; nperseg=4096, noverlap=2048, window=:hamming)
```

Each call returns a [`SpectralResult`](../reference/types.md#SigmaTau.SpectralResult)
carrying the frequency grid, the spectrum, its units, and the Welch parameters
used. Loading a `Plots` backend renders `S_y` / `S_x` on log-log axes and
`ℒ(f)` on a dB-versus-log-frequency axis:

```julia
using Plots
plot(Sy(p))
plot(L(p; f_carrier=10e6))
```

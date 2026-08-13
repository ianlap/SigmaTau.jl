# Theory: Overview

Time-domain frequency stability analysis describes how a clock's output
varies on different timescales. The classical sample variance diverges
for the dominant noise types in oscillator data (flicker FM, random-walk
FM), so the field uses convergent finite-difference variance estimators
indexed by averaging time `τ`. SigmaTau implements the standard family
of these estimators along with their confidence intervals.

## Phase residuals and fractional frequency

The fractional frequency of an oscillator is its instantaneous frequency error normalized by its nominal frequency. Sampled at an interval $\tau_0$, phase error (in seconds, `PhaseData`) and average fractional frequency (dimensionless, `FrequencyData`) are mathematically equivalent representations of the same data.

For the equations relating phase and frequency, please read NIST SP1065 §4 [riley-2008-sp1065](@cite) and Banerjee & Matsakis [banerjee-2023-timekeeping](@cite) for more details.

## Power-law noise model

Oscillator noise is modeled as a sum of power-law processes indexed by the spectral exponent $\alpha$ of the one-sided PSD ($S_y(f) \propto f^\alpha$). The dominant noise types span from White Phase Modulation ($\alpha = +2$) down to Random Run Frequency Modulation ($\alpha = -4$).

Different estimators handle different subsets of these noise types. For example, the Allan family diverges for $\alpha \le -3$, requiring the Hadamard family (which uses a third difference to converge down to $\alpha = -5$) for records with low-frequency power-law content.

For tables detailing these noise types, variance slopes ($\mu$), and convergence properties, please read NIST SP1065 §4–5 [riley-2008-sp1065](@cite) and IEEE 1139-2022 [ieee1139-2022-definitions](@cite) for more details.

## Estimator family map

| Estimator | Difference order | Phase-averaged? | Boundary-extended? | Drift-insensitive |
|-----------|------------------|-----------------|--------------------|-------------------|
| ADEV      | 2nd              | no              | no                 | no                |
| MDEV      | 2nd              | yes             | no                 | no                |
| HDEV      | 3rd              | no              | no                 | yes               |
| MHDEV     | 3rd              | yes             | no                 | yes               |
| TDEV      | 2nd              | yes             | no                 | no                |
| HTDEV     | 3rd              | yes             | no                 | yes               |
| TOTDEV    | 2nd              | no              | yes                | no                |
| MTOTDEV   | 2nd              | yes             | yes                | no                |
| HTOTDEV   | 3rd              | no              | yes                | yes               |
| MHTOTDEV  | 3rd              | yes             | yes                | yes               |

## MTIE — Maximum Time Interval Error

MTIE complements variance estimators by providing a peak-to-peak phase metric. It characterizes the worst-case maximum phase excursion in any window of length $\tau$ across the record, which is the operative metric for synchronization applications where a single large transient matters (telecom synchronization masks, GNSS holdover budgets).

SigmaTau uses a block-windowed `_mtie_core` kernel for an allocation-free peak/trough scan.

For the MTIE definition, please read NIST SP1065 [riley-2008-sp1065](@cite) and Banerjee & Matsakis [banerjee-2023-timekeeping](@cite) for more details.

!!! note "Implementation status"
    Implemented as [`mtie`](@ref) in `SigmaTau`. It returns a `StabilityResult` with `mtie` units of seconds. There is no published EDF model for this peak statistic, so confidence bounds are not provided.

## Notation used throughout these pages

- `x[i]`: phase residual (seconds), sample interval `τ₀`
- `y[i]`: fractional frequency, `y_i = (x[i+1] − x[i]) / τ₀`
- `m`: averaging factor; `τ = m·τ₀`
- `N`: phase sample count; `M = N − 1`: frequency sample count
- `α`: power-law exponent of `S_y(f)`
- `μ`: corresponding slope of `σ²_y(τ)` on log-log axes
- `σ²_y(τ)`: fractional-frequency variance; `σ_y(τ)`: fractional-frequency
  deviation. Used by ADEV / MDEV / HDEV / MHDEV.
- `σ²_x(τ)`: time variance; `σ_x(τ)`: time deviation, units of seconds.
  Used by TDEV and HTDEV.

## Slope demonstration

On log-log axes the slope of `σ_y(τ)` identifies the dominant noise
type. Pure white PM and pure random-walk FM bracket the table above:

```@example overview
using SigmaTau, Plots, Random
Random.seed!(42)

# Synth: 4096 samples WPM (α=2) and RWFM (α=−2) at τ₀ = 1 s
N = 4096
wpm  = randn(N)
rwfm = cumsum(cumsum(randn(N)))    # double integration of white noise

τs = [2^k for k in 0:9]
r_wpm  = adev(PhaseData(wpm,  1.0), τs; ci=false)
r_rwfm = adev(PhaseData(rwfm, 1.0), τs; ci=false)

plot(r_wpm.tau, r_wpm.dev;
     xscale=:log10, yscale=:log10, marker=:circle,
     xlabel=raw"\(\tau\) (s)", ylabel=raw"\(\sigma_y(\tau)\)",
     label=raw"WPM (\(\alpha = +2\))", legend=:bottomleft)
plot!(r_rwfm.tau, r_rwfm.dev; marker=:diamond, label=raw"RWFM (\(\alpha = -2\))")

# Reference slope guides, anchored at each curve's first point.
plot!(r_wpm.tau, r_wpm.dev[1] .* (r_wpm.tau ./ r_wpm.tau[1]) .^ -1;
      ls=:dash, color=:gray, label=raw"\(\tau^{-1}\) guide")
plot!(r_rwfm.tau, r_rwfm.dev[1] .* (r_rwfm.tau ./ r_rwfm.tau[1]) .^ (1/2);
      ls=:dot, color=:gray, label=raw"\(\tau^{+1/2}\) guide")
```

The WPM record falls along the `τ⁻¹` guide (μ = −1); the RWFM record
rises along `τ^(+1/2)` (μ = +1/2). Read slopes against the guide lines,
not the tick spacing — the τ grid here is octave-spaced while the
ordinate ticks are decades.

## Where to next

- [Allan family](allan_family.md) — ADEV, MDEV, HDEV, TDEV, MHDEV, HTDEV.
- [Total family](total_family.md) — boundary-extended estimators.
- [Confidence](confidence.md) — EDF and χ² intervals.
- [Noise identification](noise_id.md) — lag-1 ACF classification.

## References

- NIST SP1065, *Handbook of Frequency Stability Analysis*, Riley & Howe
  2008. [riley-2008-sp1065](@cite)
- *IEEE Standard Definitions of Physical Quantities for Fundamental
  Frequency and Time Metrology*, IEEE 1139-2022. [ieee1139-2022-definitions](@cite)
- Banerjee & Matsakis, *An Introduction to Modern Timekeeping and Time
  Transfer*, Springer 2023. [banerjee-2023-timekeeping](@cite)

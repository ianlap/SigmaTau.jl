# Theory: Overview

Time-domain frequency stability analysis describes how a clock's output
varies on different timescales. The classical sample variance diverges
for the dominant noise types in oscillator data (flicker FM, random-walk
FM), so the field uses convergent finite-difference variance estimators
indexed by averaging time `τ`. SigmaTau implements the standard family
of these estimators along with their confidence intervals.

## Phase residuals and fractional frequency

Two equivalent representations of the same data:

```math
y_i \;=\; \frac{x_{i+1} - x_i}{\tau_0}, \qquad
x_n \;=\; x_0 + \tau_0 \sum_{i=0}^{n-1} y_i
```

`x(t)` (seconds) is the time-error of the clock relative to a reference;
`y(t)` (dimensionless) is the fractional frequency offset, sampled at
interval `τ₀`. SigmaTau stores phase as `PhaseData` and frequency as
`FrequencyData`; the deviation API accepts either, converting via
prefix-sum where needed.

(Cite [@cite Banerjee2023] for modern presentation; SP1065 §4
[@cite RileyHowe2008] for the conventions.)

## Power-law noise model

Oscillator noise is modeled as a sum of five power-law processes
indexed by `α`, the spectral exponent of the one-sided PSD `S_y(f) ∝ f^α`:

| α  | Noise type | Symbol | ADEV slope (μ_dev) | MDEV slope |
|----|------------|--------|--------------------|------------|
| +2 | White PM   | WPM    | −1                 | −3/2       |
| +1 | Flicker PM | FPM    | −1 (with log term) | −1         |
| 0  | White FM   | WFM    | −1/2               | −1/2       |
| −1 | Flicker FM | FFM    | 0                  | 0          |
| −2 | Random-walk FM | RWFM | +1/2             | +1/2       |

ADEV's degeneracy on WPM/FPM is the historical motivation for MDEV.
HDEV adds drift insensitivity by going to a third difference.

(Cite SP1065 §4–5 [@cite RileyHowe2008]; IEEE 1139-2022
[@cite IEEE1139_2022].)

## Estimator family map

| Estimator | Difference order | Phase-averaged? | Boundary-extended? | Drift-insensitive |
|-----------|------------------|-----------------|--------------------|-------------------|
| ADEV      | 2nd              | no              | no                 | no                |
| MDEV      | 2nd              | yes             | no                 | no                |
| HDEV      | 3rd              | no              | no                 | yes               |
| MHDEV     | 3rd              | yes             | no                 | yes               |
| TDEV      | 2nd              | yes             | no                 | no                |
| LDEV      | 3rd              | yes             | no                 | yes               |
| TOTDEV    | 2nd              | no              | yes                | no                |
| MTOTDEV   | 2nd              | yes             | yes                | no                |
| HTOTDEV   | 3rd              | no              | yes                | yes               |
| MHTOTDEV  | 3rd              | yes             | yes                | yes               |

## Notation used throughout these pages

- `x[i]`: phase residual (seconds), sample interval `τ₀`
- `y[i]`: fractional frequency, `y_i = (x[i+1] − x[i]) / τ₀`
- `m`: averaging factor; `τ = m·τ₀`
- `N`: phase sample count; `M = N − 1`: frequency sample count
- `α`: power-law exponent of `S_y(f)`
- `μ`: corresponding slope of `σ²_y(τ)` on log-log axes
- `σ²_y(τ)`: variance; `σ_y(τ)`: deviation

## Slope demonstration

```@example overview
using SigmaTau, Random
Random.seed!(42)

# Synth: 4096 samples WPM (α=2) and RWFM (α=−2) at τ₀ = 1 s
N = 4096
wpm  = randn(N)
rwfm = cumsum(cumsum(randn(N)))    # double integration of white noise

τs = [1, 4, 16, 64, 256]
σ_wpm  = adev(PhaseData(wpm,  1.0), τs).dev
σ_rwfm = adev(PhaseData(rwfm, 1.0), τs).dev

round.(σ_wpm; sigdigits=3), round.(σ_rwfm; sigdigits=3)
```

The first vector should fall roughly as `τ⁻¹` (μ = −1, WPM); the second
should rise roughly as `τ⁺¹/²` (μ = +1/2, RWFM).

## Where to next

- [Allan family](allan_family.md) — ADEV, MDEV, HDEV, TDEV, MHDEV, LDEV.
- [Total family](total_family.md) — boundary-extended estimators.
- [Confidence](confidence.md) — EDF and χ² intervals.
- [Noise identification](noise_id.md) — lag-1 ACF classification.

## References

- NIST SP1065, *Handbook of Frequency Stability Analysis*, Riley & Howe
  2008. [@cite RileyHowe2008]
- *IEEE Standard Definitions of Physical Quantities for Fundamental
  Frequency and Time Metrology*, IEEE 1139-2022. [@cite IEEE1139_2022]
- Banerjee & Matsakis, *An Introduction to Modern Timekeeping and Time
  Transfer*, Springer 2023. [@cite Banerjee2023]

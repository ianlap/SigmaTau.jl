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

(Cite [@cite banerjee-2023-timekeeping] for modern presentation; SP1065 §4
[@cite riley-2008-sp1065] for the conventions.)

## Power-law noise model

Oscillator noise is modeled as a sum of power-law processes indexed
by `α`, the spectral exponent of the one-sided PSD `S_y(f) ∝ f^α`:

| α  | Noise type | Symbol | ADEV slope (μ_dev) | MDEV slope |
|----|------------|--------|--------------------|------------|
| +2 | White PM         | WPM     | −1                 | −3/2       |
| +1 | Flicker PM       | FPM     | −1 (with log term) | −1         |
|  0 | White FM         | WFM     | −1/2               | −1/2       |
| −1 | Flicker FM       | FFM     | 0                  | 0          |
| −2 | Random-walk FM   | RWFM    | +1/2               | +1/2       |
| −3 | Flicker walk FM  | FWFM    | (Allan diverges)¹  | (diverges)¹ |
| −4 | Random run FM    | RRFM    | (Allan diverges)¹  | (diverges)¹ |

¹ The Allan family's variance integral diverges for `α ≤ −3`. The
Hadamard family (HDEV / MHDEV / HTOTDEV / MHTOTDEV) replaces the
second difference with a third difference and remains finite down to
`α = −4`, so records with very-low-frequency power-law content should
be analyzed with the Hadamard family rather than ADEV / MDEV
[@cite greenhall-1997-third-difference-mvar]. SigmaTau's three-state
clock SDE additionally models random-run FM via the σ₃ Wiener channel
[@cite zucca-2005-clock-model-allan].

ADEV's degeneracy on WPM/FPM is the historical motivation for MDEV.
HDEV adds drift insensitivity by going to a third difference.

(Cite SP1065 §4–5 [@cite riley-2008-sp1065]; IEEE 1139-2022
[@cite ieee1139-2022-definitions].)

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

MTIE complements the variance estimator family with a **peak-to-peak**
phase metric: at each averaging interval `τ`, MTIE reports the largest
phase excursion observed in any window of length `τ` across the entire
record [@cite riley-2008-sp1065]. Where the Allan family characterizes
the *spread* of the phase residual, MTIE characterizes the worst-case
*excursion* — the relevant figure of merit for synchronization
applications where a single large transient is operationally
significant (telecom synchronization masks, GNSS holdover budgets,
financial-trading timestamping). Banerjee & Matsakis frame it as the
canonical time-distribution-network metric
[@cite banerjee-2023-timekeeping].

```math
\mathrm{MTIE}(\tau) \;=\; \max_{1 \le i \le N-m}\;
\Bigl[\max_{0 \le k \le m} x_{i+k} \;-\; \min_{0 \le k \le m} x_{i+k}\Bigr],
\qquad m = \tau / \tau_0.
```

!!! note "Implementation status"
    Implemented as [`mtie`](@ref) in `SigmaTau` (block-windowed
    `_mtie_core` kernel with allocation-free per-window peak/trough
    scan). Returns a `StabilityResult` with `mtie` units of seconds;
    no CI bounds — no published EDF model for the peak statistic.

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

- [Allan family](allan_family.md) — ADEV, MDEV, HDEV, TDEV, MHDEV, HTDEV.
- [Total family](total_family.md) — boundary-extended estimators.
- [Confidence](confidence.md) — EDF and χ² intervals.
- [Noise identification](noise_id.md) — lag-1 ACF classification.

## References

- NIST SP1065, *Handbook of Frequency Stability Analysis*, Riley & Howe
  2008. [@cite riley-2008-sp1065]
- *IEEE Standard Definitions of Physical Quantities for Fundamental
  Frequency and Time Metrology*, IEEE 1139-2022. [@cite ieee1139-2022-definitions]
- Banerjee & Matsakis, *An Introduction to Modern Timekeeping and Time
  Transfer*, Springer 2023. [@cite banerjee-2023-timekeeping]

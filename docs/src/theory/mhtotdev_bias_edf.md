# Theory: MHTOTDEV bias and EDF

MHTOTDEV (modified Hadamard total deviation) is novel to SigmaTau — no other
library computes it, and no published reference gives its bias factor or its
equivalent degrees of freedom (EDF). Every other total-family estimator inherits
a coefficient table from the literature (TOTVAR and HTOT from SP1065 / FCS 2001)
or one reverse-engineered from Stable32 (MTOT). MHTOTDEV has neither, so SigmaTau
**measures** both by Monte Carlo against synthesized known-noise.

This page documents the methodology; it is the reproducible basis for the
shipped `_coeff_mhtot` EDF coefficients and the `bias_correction(:mhtot, …)`
table. The harness is [`tools/mc_mhtotdev.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/tools/mc_mhtotdev.jl).

## The reference: MHDEV

MHTOTDEV is the *total* (boundary-extended) form of the modified-Hadamard
variance MHVAR. The natural unbiased reference is therefore MHDEV (= `mhdev`),
whose EDF SigmaTau already computes from first principles via the
Greenhall–Riley spectral sum (`_calc_edf_core` with difference order `d = 3`,
inner-average `F = 1`). Writing the per-τ variance estimates as

```math
\hat V = \mathrm{MHTOTVAR}(\tau) = \mathtt{mhtotdev}^2, \qquad
\hat W = \mathrm{MHVAR}(\tau)    = \mathtt{mhdev}^2 ,
```

both are computed on the *same* realization at matched averaging factor `m`, so
finite-sample fluctuations are correlated and partly cancel in their ratio.

## Estimators

**EDF.** A variance estimator that is χ²-distributed with `ν` degrees of freedom
satisfies `Var[\hat V] = 2\,E[\hat V]^2/ν`, so from `R` realizations

```math
\widehat{\mathrm{edf}} \;=\; \frac{2\,\bar V^2}{s_V^2},
```

with `\bar V` the sample mean and `s_V^2` the Bessel-corrected sample variance.
A nonparametric bootstrap over the `R` realizations gives the standard error.
The relative precision scales as `SE/edf ≈ √(2/R)`, so `R = 1000` gives ≈4.5 %
and the per-α anchor cell uses `R = 5000`.

**Bias.** The variance-scale bias factor is the ratio of means

```math
B(\alpha) \;=\; \frac{E[\hat V]}{E[\hat W]} \;\approx\; \frac{\bar V}{\bar W},
```

matching the `bias_correction` convention `σ_unbiased = σ_raw / √B` used across
the total family. `B > 1` means MHTOTDEV reads high (correction lowers σ),
`B < 1` means it reads low.

**EDF model.** Per noise type the EDF is fit to the same functional form the
other total estimators use,

```math
\mathrm{edf}(\tau) \;=\; b\,\frac{T}{\tau} - c ,
```

by weighted least squares (weights `1/\mathrm{SE}^2`) over the validity window
`τ ≤ T/10`, where `T = (N-1)\,τ_0`. Coefficients, their standard errors, and
`R²` are reported per α.

## Analytic cross-check

Before any ratio is trusted, the harness verifies that the synthesized noise
reproduces the analytic MHVAR power law `σ²_MH(τ) ∝ τ^{μ(α)}`:

| α  | Noise | μ(α) (variance slope) |
|---:|:------|:----------------------|
|  2 | WHPM  | −3 |
|  1 | FLPM  | −2 |
|  0 | WHFM  | −1 |
| −1 | FLFM  |  0 |
| −2 | RWFM  | +1 |

(`μ(α) = −α − 1` for α ≤ 0; the *modified* estimator rolls WHPM/FLPM off as
`τ^{−3}`/`τ^{−2}` in variance, which is what lets it disambiguate the two PM
noises.) The measured `log MHVAR` vs `log τ` slope must match μ(α) to within a
small tolerance, validating both the `noise_gen` calibration and the MHDEV
reference.

## Monte Carlo design

- **α** ∈ {2, 1, 0, −1, −2}, using the *known* injected α (never `identify_noise`,
  to avoid noise-ID misclassification contaminating the per-α bins).
- **N** ∈ {1025, 2049, 4097, 8193, 16385, 32769} (`2^k + 1`, so `T = 2^k` and
  `T/τ` lands on round values); octave `m`-grid `1 … N÷4`, dropping cells with
  `nsubs = N − 4m + 1 < 16`. Pooling across N gives ≥30 `(T/τ, edf)` points per
  α for the fit.
- **Reproducibility.** Each realization draws from an independent stream
  `Xoshiro(hash((seed, α, N, r)))`, so results are identical regardless of
  thread count. Synthesis takes an `rng` keyword for exactly this purpose.
- **Scale invariance.** `B` and EDF are ratios of like-dimensioned quantities,
  so the arbitrary amplitude of the uncalibrated `_gen_powerlaw_y` draw cancels.

The harness emits a per-cell CSV and a JSON of metadata + fitted coefficients
(git SHA, seed, thread count, N grid, R, the μ(α) check result) under
`tools/artifacts/`, so a shipped coefficient is traceable to the run that
produced it.

!!! note "Provisional vs. authoritative numbers"
    The coefficient tables in `_coeff_mhtot` and `bias_correction(:mhtot, …)`
    are populated from the full sweep (N up to 32769, R = 1000 / 5000), run on a
    workstation. A reduced laptop validation (N ≤ 2049, R = 200–400) reproduces
    the μ(α) slopes and the qualitative bias trend below; the published tables
    use the full-sweep values.

## Result (preliminary, laptop validation)

From the reduced validation run (N ≤ 2049), the measured bias already departs
clearly from the previously-assumed `B = 1`:

| α  | Noise | B(α) | EDF fit `b` | `c` | R² |
|---:|:------|:-----|:------------|:----|:---|
|  2 | WHPM  | ≈ 0.68 | 0.67 | −16.3 | 0.61 |
|  1 | FLPM  | ≈ 0.74 | 0.78 | −6.3  | 0.80 |
|  0 | WHFM  | ≈ 0.81 | 0.94 | −0.24 | 0.91 |
| −1 | FLFM  | ≈ 0.93 | 0.97 |  3.16 | 0.98 |
| −2 | RWFM  | ≈ 1.34 | 0.81 |  3.77 | 0.98 |

MHTOTDEV reads **low** for phase-modulation noise (B < 1) and **high** for
random-walk FM (B > 1) — so treating it as unbiased systematically misstates the
deviation at both ends of the noise range. The full-sweep values, with bootstrap
uncertainties, replace these once the workstation run lands.

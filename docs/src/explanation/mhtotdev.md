# MHTOTDEV bias and EDF

MHTOTDEV (modified Hadamard total deviation) is original to SigmaTau: no other
library computes it, and no published reference gives its bias factor or its
equivalent degrees of freedom (EDF). Every other total-family estimator inherits
a published coefficient table (TOTVAR and MTOT from SP1065 Tables 7–8; HTOT from
the FCS 2001 paper). MHTOTDEV has none, so SigmaTau
**measures** both by Monte Carlo against synthesized known-noise.

This page documents the methodology; it is the reproducible basis for the
shipped `_coeff_mhtot` EDF coefficients and the `bias_correction(:mhtot, …)`
table. The harness is [`tools/mc_mhtotdev.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/tools/mc_mhtotdev.jl).
The estimator under test is the shipped `mhtotdev` kernel, which removes each
window's least-squares quadratic before the reflection, so the coefficients
calibrate exactly the statistic users compute. That per-window detrend makes
MHTOTDEV immune to linear frequency drift on its own and lets it converge for
power laws down to α = −4, like HDEV and MHDEV.

The procedure follows the NIST total-variance literature: Howe, Beard,
Greenhall, Vernotte & Riley, *A Total Estimator of the Hadamard Function Used
for GPS Operations* (32nd PTTI, 2000), for the bias and EDF definitions and the
`τ/τ0 ≥ 16` validity window, and Vernotte & Howe, *Generalization of the Total
Variance Approach to the Different Classes of Structure Functions*, for the
Monte Carlo EDF estimator and the extension methodology.

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
The relative precision scales as `SE/edf ≈ √(2/R)`; the full sweep uses
`R = 3000` per cell, giving ≈2.6 %.

**Bias.** Following the normalized-bias definition of Howe et al. 2000 (eqn 6),
`nbias = E[\hat V]/E[\hat W] - 1`, the variance-scale bias factor is the ratio
of means

```math
B(\alpha) \;=\; \frac{E[\hat V]}{E[\hat W]} \;=\; 1 + \mathrm{nbias} \;\approx\; \frac{\bar V}{\bar W},
```

matching the `bias_correction` convention `σ_unbiased = σ_raw / √B` used across
the total family. `B > 1` means MHTOTDEV reads high (correction lowers σ),
`B < 1` means it reads low.

**EDF model.** The total-variance literature offers two empirical EDF forms:
the Mod-Totvar form (Vernotte & Howe, eqn 2) and the TotHvar form
(Howe et al. 2000, eqn 7),

```math
\mathrm{edf}(\tau) = b\,\frac{T}{\tau} - c
\qquad\text{vs.}\qquad
\mathrm{edf}(\tau) = \frac{T/\tau}{b_0 + b_1\,\tau/T},
```

with `T = (N-1)\,τ_0`. MHTOTDEV is the modified-*Hadamard* total, a hybrid of
the two, so the harness fits **both** by weighted least squares
(weights `1/\mathrm{SE}^2`) and selects the higher-R² form per α. Critically,
the fit is restricted to the validity window `τ/τ_0 ≥ 16`: Howe et al. 2000
note their EDF model "should be used only if … `τ/τ_0 ≥ 16`," below which the
data-extension confers no advantage over the plain estimator and the fit
degrades. Coefficients, standard errors, and `R²` are reported per α for both
forms.

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
- **N** (record length, frequency samples) ∈ {1024, 2048, 4096, 8192, 16384,
  32768}; octave `m`-grid `1 … N÷4`, dropping cells with `nsubs = N − 4m + 1 < 16`.
  Pooling across N gives ≥30 `(T/τ, edf)` points per α for the fit.
- **Both the bias and the EDF use the `τ/τ_0 ≥ 16` window.** The small-`m` cells
  (especially `m = 1`, which is near-degenerate yet carries the tightest standard
  error, hence the largest weight) otherwise skew the bias ratio away from its
  physical value.
- **Reproducibility.** Each realization draws from an independent stream
  `Xoshiro(hash((seed, α, N, r)))`, so results are identical regardless of
  thread count. Synthesis takes an `rng` keyword for exactly this purpose.
- **Scale invariance.** `B` and EDF are ratios of like-dimensioned quantities,
  so the arbitrary amplitude of the uncalibrated `_gen_powerlaw_y` draw cancels.

The harness emits a per-cell CSV and a JSON summary of metadata plus fitted
coefficients (git SHA, seed, thread count, N grid, R, and the μ(α) check
result). The authoritative full-sweep JSON is tracked at
`tools/artifacts/mhtotdev_mc_full.json`; ad-hoc quick-run outputs stay local.

## Result

Full sweep (N = 1024 … 32768, R = 3000, both bias and EDF over `τ/τ_0 ≥ 16`,
per-window quadratic detrend, Kasdin–Walter synthesis, git `9d7a408`). All
EDF fits have R² ≥ 0.9988:

| α  | Noise | EDF `b`, `c` | bias `b0` | `b1` |
|---:|:------|:-------------|:----------|:-----|
|  2 | WHPM  | 1.5429, 3.8504 | 0.9005 | −0.0083 |
|  1 | FLPM  | 1.0409, 2.9661 | 0.8238 | −0.0541 |
|  0 | WHFM  | 0.9355, 2.8108 | 0.7933 | −0.0422 |
| −1 | FLFM  | 0.8780, 2.5258 | 0.7841 | −0.0357 |
| −2 | RWFM  | 0.8580, 2.5361 | 0.7879 | −0.0325 |

MHTOTDEV **reads low** for every noise type: `B ≈ 0.90` for white PM and
`B ≈ 0.79` for the FM noises. Removing a quadratic from each `4m`-sample
window absorbs part of that window's variance before the reflection, and the
loss is nearly the same for every FM noise, so the correction is close to a
constant factor. The bias is modeled τ/T-linearly, `B = b0 + b1·(τ/T)`; the
`b1` term is small and cuts the bias-fit residual by 4 % (WHPM) to 17 %
(FLPM).

Two methodology points matter for reproducing the table: (1) restricting the
bias to the same
`τ/τ_0 ≥ 16` window as the EDF, since the near-degenerate `m = 1` cell otherwise
dominates the weighted ratio and pulls the apparent bias spuriously low; and
(2) the `τ/τ_0 ≥ 16` floor on the EDF fit itself, which holds R² ≥ 0.9988 across
all noise types. The tracked full-sweep artifact records the git SHA, seed,
grid, and per-α fits.

!!! warning "Edge zone at the reach limit"
    The kernel needs `N − 4m + 1 ≥ 1`, so MHTOTDEV's reach limit is `τ = T/4`
    (`T/τ ≥ 4`). The fitted linear form crosses `ν = 1` near
    `T/τ ≈ 3.1–4.1` depending on α (lowest at α = +2, highest at α = −2).
    The χ² interval is still defined below `ν = 1` and is very wide there
    (0.71 σ to 5.9 σ at ν = 0.9 and 68.3 %), as one effective average
    warrants. Treat confidence intervals as valid for `T/τ ≳ 5–6` and
    indicative only in the final octave before the reach limit; the σ values
    themselves are valid all the way to `τ = T/4`.

!!! note "Degrees-of-freedom gain over MHDEV"
    The extension buys the most degrees of freedom for the phase noises. At
    the same τ (N ≥ 4096, 16 ≤ m ≤ N/32) the measured MHTOTDEV EDF is about
    1.36× the Greenhall–Riley MHDEV EDF for WHPM, 1.18× for FLPM, 1.11× for
    WHFM, 1.04× for FLFM and 1.0× for RWFM. For the red FM noises MHTOTDEV's
    advantage is reach (`τ = T/4`) and drift immunity, not precision.

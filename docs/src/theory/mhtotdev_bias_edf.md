# Theory: MHTOTDEV bias and EDF

MHTOTDEV (modified Hadamard total deviation) is novel to SigmaTau — no other
library computes it, and no published reference gives its bias factor or its
equivalent degrees of freedom (EDF). Every other total-family estimator inherits
a published coefficient table (TOTVAR and MTOT from SP1065 Tables 7–8; HTOT from
the FCS 2001 paper). MHTOTDEV has none, so SigmaTau
**measures** both by Monte Carlo against synthesized known-noise.

This page documents the methodology; it is the reproducible basis for the
shipped `_coeff_mhtot` EDF coefficients and the `bias_correction(:mhtot, …)`
table. The harness is [`tools/mc_mhtotdev.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/tools/mc_mhtotdev.jl).
The estimator under test is the shipped `mhtotdev` pipeline, including its
default global least-squares drift removal (`remove_drift=true`; see
[Theory: Total family](total_family.md)), so the coefficients calibrate
exactly the statistic users compute.

The procedure follows the NIST total-variance literature: Howe, Beard,
Greenhall, Vernotte & Riley, *A Total Estimator of the Hadamard Function Used
for GPS Operations* (32nd PTTI, 2000) — for the bias and EDF definitions and the
`τ/τ0 ≥ 16` validity window — and Vernotte & Howe, *Generalization of the Total
Variance Approach to the Different Classes of Structure Functions* — for the
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

with `T = (N-1)\,τ_0`. MHTOTDEV is the modified-*Hadamard* total — a hybrid of
the two — so the harness fits **both** by weighted least squares
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

!!! note "Provisional vs. authoritative numbers"
    The coefficient tables in `_coeff_mhtot` and `bias_correction(:mhtot, …)`
    are populated from the full sweep (N up to 32768, R = 3000), run on a
    workstation. The harness's reduced `quick` mode (N ≤ 2048, R = 200–400)
    reproduces the μ(α) slopes and the qualitative bias trend below; the
    published tables use the full-sweep values.

## Result

Full sweep (N = 1024 … 32768, R = 3000, both bias and EDF over `τ/τ_0 ≥ 16`).
All EDF fits are excellent (R² ≥ 0.998):

| α  | Noise | EDF `b`, `c` | bias `b0` | `b1` |
|---:|:------|:-------------|:----------|:-----|
|  2 | WHPM  | 1.853, 5.482 | 1.064 |  0.017 |
|  1 | FLPM  | 1.219, 3.669 | 0.984 |  0.036 |
|  0 | WHFM  | 1.100, 3.504 | 1.019 | −0.048 |
| −1 | FLFM  | 1.030, 3.387 | 1.213 | −0.321 |
| −2 | RWFM  | 0.813, 2.541 | 1.943 | −3.588 |

MHTOTDEV is **≈ unbiased for white/flicker noise** (`b0 ≈ 1`, `b1 ≈ 0`) and
reads **progressively high for redder FM** — for random-walk FM, `B ≈ 1.9` at
small τ, falling toward 1 as τ → T. The bias is modeled τ/T-linearly,
`B = b0 + b1·(τ/T)`: the `b1` term is negligible for the whiter noises but cuts
the bias-fit residual by ~46 % (FLFM) and ~97 % (RWFM), so the τ/T dependence is
essential at the red end and a constant would over-correct near large τ.

Two methodology points matter for reproducing the table: (1) restricting the
bias to the same
`τ/τ_0 ≥ 16` window as the EDF — the near-degenerate `m = 1` cell otherwise
dominates the weighted ratio and pulls the apparent bias spuriously below 1; and
(2) the `τ/τ_0 ≥ 16` floor on the EDF fit itself, which holds R² ≥ 0.998 across
all noise types. The tracked full-sweep artifact records the git SHA, seed,
grid, and per-α fits.

!!! warning "Edge zone at the reach limit"
    The kernel needs `N − 4m + 1 ≥ 1`, so MHTOTDEV's reach limit is `τ = T/4`
    (`T/τ ≥ 4`). The fitted linear form crosses `ν = 1` near
    `T/τ ≈ 4.1–4.4` depending on α, and a χ² interval computed from `ν < 1`
    is degenerate (the bounds collapse). Treat confidence intervals as valid
    for `T/τ ≳ 5` and indicative only in the final fraction of an octave
    before the reach limit; the σ values themselves are valid all the way to
    `τ = T/4`.

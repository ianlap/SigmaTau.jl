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

Full sweep (N = 1024 … 32768, R = 3000, both bias and EDF over `τ/τ_0 ≥ 16`,
with the wrapper's default `remove_drift=true` global drift removal in the
loop). All EDF fits are excellent (R² ≥ 0.9982):

| α  | Noise | EDF `b`, `c` | bias `b0` | `b1` |
|---:|:------|:-------------|:----------|:-----|
|  2 | WHPM  | 1.6869, 4.3166 | 1.0308 | −0.1197 |
|  1 | FLPM  | 1.1671, 3.2107 | 0.9883 | −0.0686 |
|  0 | WHFM  | 1.0861, 3.3451 | 1.0545 | −0.1706 |
| −1 | FLFM  | 1.0095, 3.3254 | 1.2947 | −0.7193 |
| −2 | RWFM  | 0.7631, 2.4490 | 2.1686 | −5.5884 |

MHTOTDEV is **≈ unbiased for the PM noises** (`b0 ≈ 1`, `b1` small) and reads
**progressively high for redder FM** at small τ — for random-walk FM,
`B ≈ 2.2`. The bias is modeled τ/T-linearly, `B = b0 + b1·(τ/T)`. The `b1`
term now carries two physical effects with the same sign: the
boundary-extension bias falling off toward the record length, and the
suppression of totalized low-frequency power by the global least-squares
drift removal as τ → T (the well-known effect of drift removal on red
noise). Together they cut the bias-fit residual by ~18.5 % (WHPM), ~33.2 %
(FLPM), ~44.7 % (WHFM), ~61.1 % (FLFM), and ~96.0 % (RWFM); even the PM
noises now show a small but nonzero `b1`.

Two methodology points matter for reproducing the table: (1) restricting the
bias to the same
`τ/τ_0 ≥ 16` window as the EDF — the near-degenerate `m = 1` cell otherwise
dominates the weighted ratio and pulls the apparent bias spuriously below 1; and
(2) the `τ/τ_0 ≥ 16` floor on the EDF fit itself, which holds R² ≥ 0.9982 across
all noise types. The tracked full-sweep artifact records the git SHA, seed,
grid, and per-α fits.

!!! warning "Edge zone at the reach limit"
    The kernel needs `N − 4m + 1 ≥ 1`, so MHTOTDEV's reach limit is `τ = T/4`
    (`T/τ ≥ 4`). The fitted linear form crosses `ν = 1` near
    `T/τ ≈ 3.2–4.5` depending on α (lowest at α = +2, highest at α = −2), and a
    χ² interval computed from `ν < 1` is degenerate (the bounds collapse).
    Treat confidence intervals as valid for `T/τ ≳ 5–6` and indicative only
    in the final octave before the reach limit; the σ values themselves are
    valid all the way to `τ = T/4`.

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
default record-level least-squares drift removal (`remove_drift=true`; see
[Theory: Total family](total_family.md)), so the coefficients calibrate
exactly the statistic users compute.

!!! warning "Provisional coefficients"
    The kernel was redefined (per-window **quadratic** least-squares detrend,
    matched to the third difference's annihilation degree), which invalidated
    the previous full-sweep tables. The numbers on this page — and the tables
    currently shipped in `_coeff_mhtot` and `bias_correction(:mhtot, …)` — come
    from the harness's reduced `quick` mode (N ≤ 2048, R = 200/400) and are
    **provisional** pending the replacement full sweep (R = 3000, N up to
    32768, α extended to −3/−4).

The procedure follows the NIST total-variance literature: Howe, Beard,
Greenhall, Vernotte & Riley, *A Total Estimator of the Hadamard Function Used
for GPS Operations* (32nd PTTI, 2000) — for the bias and EDF definitions and the
`τ/τ0 ≥ 16` validity window — and Vernotte & Howe, *Generalization of the Total
Variance Approach to the Different Classes of Structure Functions* — for the
Monte Carlo EDF estimator and the extension methodology.

## The estimator under test

MHTOTDEV follows the Greenhall methodology its MTOTDEV/HTOTDEV siblings use —
overlapping subsegments of length `L = 4m`, a per-window detrend, a 3× time-
reverse (even/mirror) extension, then averaged third differences — with one
deliberate departure: **the per-window detrend is a quadratic least-squares
fit**, not the half-mean slope removal the second-difference totals use.

**The design rule: the detrend degree must equal the kernel's annihilation
degree.** The third difference annihilates phase polynomials of degree ≤ 2,
just as MTOT's second difference annihilates degree ≤ 1 and its linear detrend
matches. With only a slope removal, the stochastic quadratic carried by noise
at α ≤ −3 survives into the time-reverse extension, the seam kink feeds the
boundary-spanning third differences, and the estimator inherits MVAR's
divergence — measured as a log-log variance slope of 2 instead of 3 at α = −4,
with the level growing ∝ T. The matched quadratic detrend restores convergence
over the **full Hadamard range α > −5** (slope 3.00, level N-stable) while
leaving the PM statistics healthy (`B(WHPM) ≈ 0.91`, flat in m and N), verified
by two independent implementations.

The rule was established by a controlled 2×2 study — detrend order (linear vs.
quadratic) × reflection parity (even/mirror vs. odd/point), run on identical
noise realizations. The pathologies that earlier rejected alternative per-window
schemes (`B(WHPM)` of 3.8–13 at m = 16, growing ∝ m, plus χ² EDF collapse)
belong to the **odd (point) reflection** variants and appear there regardless of
detrend order. The quadratic detrend with the unchanged even reflection shows
neither.

**Drift.** Because the quadratic fit removes any deterministic frequency offset
and drift exactly inside every window, **this kernel is drift-immune on its
own**. The wrapper's record-level removal (`remove_drift=true` default) is kept
for noise-identification hygiene — a strongly drifting record otherwise
misclassifies its long-τ points, which would corrupt the EDF and bias lookups —
not because the statistic needs it. The calibration below is nevertheless run
with that default in the loop, so the coefficients match the shipped pipeline
exactly.

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
The relative precision scales as `SE/edf ≈ √(2/R)`; the pending full sweep uses
`R = 3000` per cell, giving ≈2.6 %, while the provisional quick mode below uses
`R = 200/400` (≈10/7 %).

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
(weights `1/\mathrm{SE}^2`) and reports both. The shipped table keeps the
Mod-Totvar form uniformly across α (house convention: one form for the whole
table, so `calculate_edf` has a single branch); the form choice is re-decided
at the full refit. Critically,
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
| −3 | FWFM  | +2 |
| −4 | RRFM  | +3 |

(`μ(α) = −α − 1` across the whole range; the *modified* estimator rolls
WHPM/FLPM off as `τ^{−3}`/`τ^{−2}` in variance, which is what lets it
disambiguate the two PM noises.) The measured `log MHVAR` vs `log τ` slope must
match μ(α) to within a small tolerance, validating both the `noise_gen`
calibration and the MHDEV reference.

The α = −3 and α = −4 rows are the convergence check on the redefined kernel:
MHTOTVAR itself reproduces the same slopes (3.00 at α = −4) with an N-stable
level, where the old half-mean-detrend kernel gave slope 2 with a level ∝ T.

## Monte Carlo design

- **α** ∈ {2, 1, 0, −1, −2, −3, −4}, using the *known* injected α (never
  `identify_noise`, to avoid noise-ID misclassification contaminating the per-α
  bins). The FWFM/RRFM rows became measurable once the quadratic detrend made
  the kernel converge there; they are collected by the harness but are **not
  yet shipped** in `_coeff_mhtot` or `bias_correction` (see the note after the
  result table).
- **N** (record length, frequency samples) ∈ {1024, 2048, 4096, 8192, 16384,
  32768} for the full sweep, {1024, 2048} in quick mode; octave `m`-grid
  `1 … N÷4`. Pooling across N gives ≥30 `(T/τ, edf)` points per α for the fit.
- **Edge cells are measured.** The shipped default grid reaches `m = N÷4`
  (`τ/T ≈ 0.25`), so the calibration must reach it too. MHTOTVAR is defined
  down to a single subsegment and its EDF comes from the ensemble variance
  across realizations, so EDF cells need only `nsubs = N − 4m + 1 ≥ 1`
  (`MIN_NSUBS = 1`). The **bias** ratio additionally needs a well-conditioned
  MHVAR reference, so bias cells keep the `nsubs ≥ 16` floor and require
  `m ≤ (N−1)÷4` (MHVAR's own reach).
- **Both the bias and the EDF fit over the `τ/τ_0 ≥ 16` window.** The small-`m`
  cells (especially `m = 1`, which is near-degenerate yet carries the tightest
  standard error, hence the largest weight) otherwise skew the bias ratio away
  from its physical value.
- **Reproducibility.** Each realization draws from an independent stream
  `Xoshiro(hash((seed, α, N, r)))`, so results are identical regardless of
  thread count. Synthesis takes an `rng` keyword for exactly this purpose.
- **Scale invariance.** `B` and EDF are ratios of like-dimensioned quantities,
  so the arbitrary amplitude of the uncalibrated `_gen_powerlaw_y` draw cancels.

The harness emits a per-cell CSV and a JSON summary of metadata plus fitted
coefficients (git SHA, seed, thread count, N grid, R, and the μ(α) check
result). The tracked artifact path is `tools/artifacts/mhtotdev_mc_full.json`;
ad-hoc quick-run outputs stay local. Reproduce with

```
julia --project=. --threads=auto tools/mc_mhtotdev.jl        # full sweep
julia --project=. --threads=auto tools/mc_mhtotdev.jl quick  # laptop check
```

## Result (PROVISIONAL — quick mode)

Quick mode (`seed = 20260531`, N ∈ {1024, 2048}, R = 200/400, edge cells
included, both bias and EDF fitted over `τ/τ_0 ≥ 16`, with the wrapper's
default `remove_drift=true` in the loop). EDF fits R² ≥ 0.965:

| α  | Noise | EDF `b`, `c` | bias `b0` | `b1` |
|---:|:------|:-------------|:----------|:-----|
|  2 | WHPM  | 1.4449, 1.5717 | 0.9101 | −0.1341 |
|  1 | FLPM  | 0.9951, 1.3987 | 0.8252 | −0.0210 |
|  0 | WHFM  | 0.8929, 1.5941 | 0.7935 | −0.0177 |
| −1 | FLFM  | 0.8298, 1.7928 | 0.7842 | −0.0196 |
| −2 | RWFM  | 0.7722, 1.7331 | 0.7880 | −0.0159 |

These are the values currently shipped in `_coeff_mhtot` and
`bias_correction(:mhtot, …)`. **Replace them from the full-sweep artifact when
it lands.**

**Bias.** With the matched quadratic detrend, MHTOTDEV reads uniformly **LOW**
— `B ≈ 0.78–0.91`, the same sign as its MTOT and TOTH siblings — and `B` is
very nearly a **per-α constant**: the τ/T-linear model `B = b0 + b1·(τ/T)` is
retained for uniformity with `:totvar`, but `b1 ∈ [−0.13, −0.02]` contributes
at most a few percent over the whole grid. The mechanism is simple: the
per-window quadratic fit strips a little genuine noise power at every window
length, roughly independently of τ. There is no longer a "progressively high
for redder FM" trend and no second physical effect riding on `b1` — the
drift-removal-versus-boundary-extension story that the half-mean-detrend
kernel needed is obsolete, because the kernel now removes drift itself and the
record-level removal has nothing left to take out.

One methodology point still matters for reproducing the table: the bias must be
restricted to the same `τ/τ_0 ≥ 16` window as the EDF — the near-degenerate
`m = 1` cell otherwise dominates the weighted ratio and drags the apparent bias
further below 1 than the physical value.

!!! note "α = −3 / −4 not yet shipped"
    The kernel converges for α > −5, and the harness measures FWFM and RRFM,
    but `_coeff_mhtot` and `bias_correction(:mhtot, …)` still stop at α = −2:
    `identify_noise` cannot return `:FWFM`/`:RRFM`, so nothing would look those
    rows up. Extending both the classifier and the tables is a follow-up after
    the full sweep. (In `noise_id_nnls`, where the basis is chosen explicitly
    rather than classified, `:mhtotdev` already carries the full α > −5 basis
    and the α ≤ −3 rows borrow the RWFM EDF for weighting only.)

!!! note "The reach limit is now inside the calibration"
    The kernel needs `N − 4m + 1 ≥ 1`, so MHTOTDEV's reach limit is `τ = T/4`
    (`T/τ ≥ 4`) — exactly where the shipped default `m`-grid ends. Earlier
    calibrations stopped at `T/τ ≥ 8` and had to *extrapolate* both models over
    the last octave; the harness now measures EDF cells down to a single
    subsegment (`MIN_NSUBS = 1`), so `T/τ` down to 4 is fitted data rather than
    extrapolation. Two consequences:

    - The old bias caveat (the model overpredicting `B` by 13–24 % in the
      unmeasured zone) is **resolved by measurement** — the near-constant `B`
      of the quadratic-detrend kernel is fitted across the whole grid,
      including the edge cells.
    - The EDF model now stays `ν > 1` over the entire default grid, so the
      `max(…, 1.0)` floor in `calculate_edf` is a **dormant guard**: it is
      retained so a future coefficient set can never silently trip the
      Gaussian fallback in `confidence_intervals`, but with the shipped
      coefficients it does not bind.

    Confidence intervals are therefore usable across the whole default grid;
    the σ values themselves were always valid to `τ = T/4`. Bias cells still
    stop at `nsubs ≥ 16` (MHVAR, the reference, has no meaningful value below
    that), so the bias model's last octave leans on the fit's near-constancy
    rather than on a directly measured ratio.

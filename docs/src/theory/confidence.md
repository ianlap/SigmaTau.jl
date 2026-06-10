# Theory: Confidence Intervals

Each `StabilityResult.dev` value is a point estimate of the underlying
deviation. The corresponding confidence interval rests on the
chi-squared distribution applied to an *equivalent number of degrees of
freedom* (EDF) that depends on the estimator, the averaging factor `m`,
the record length `N`, the record span `T` (for the total family), and
the identified noise type `α`.

## EDF as a bridge

The variance estimate behaves approximately as a scaled chi-squared
variable:

```math
\mathrm{EDF} \cdot \frac{\hat{\sigma}^2_y(\tau)}{\sigma^2_y(\tau)}
\;\sim\; \chi^2_{\mathrm{EDF}},
\qquad
\mathrm{EDF} = \frac{2\,\bigl(E[\hat\sigma^2_y]\bigr)^2}{\operatorname{Var}[\hat\sigma^2_y]},
```

where `σ̂²_y(τ)` is the variance estimate, `σ²_y(τ)` is the true
variance, and `χ²_EDF` is the chi-squared distribution with `EDF`
degrees of freedom — the defining moment relation is Eq. (1) of
GR03 [greenhall-2003-edf-stability](@cite). For confidence level `CL`
(default `DEFAULT_CONFIDENCE = 0.683`, the 1σ convention shared with
Stable32 and allantools), the lower and upper bounds on `σ²_y(τ)` are:

```math
\sigma^2_{\text{lo}} = \frac{\mathrm{EDF}\cdot\hat\sigma^2}{\chi^2_{(1+\mathrm{CL})/2}},
\quad
\sigma^2_{\text{hi}} = \frac{\mathrm{EDF}\cdot\hat\sigma^2}{\chi^2_{(1-\mathrm{CL})/2}},
```

where `χ²_p` denotes the `p`-quantile of the chi-squared distribution
with `EDF` degrees of freedom: the lower variance bound divides by the
*upper* quantile and vice versa, per SP1065 §5
[riley-2008-sp1065](@cite). Each `StabilityResult.ci` entry carries the
deviations `√σ²_lo` and `√σ²_hi` as its `lo` and `hi` components;
SigmaTau evaluates the quantiles with Distributions.jl rather than
table lookup.

The χ² approximation is empirical and asymptotic in EDF: the interval
narrows as `N/m` grows, and for a fixed record the EDF peaks near white
FM. At small EDF (long τ, short records, or steep red noise) the
interval is wide and asymmetric — the upper bound stretches much
further than the lower bound, reflecting the long right tail of χ² at
low degrees of freedom. When the EDF is non-finite or below 1, SigmaTau
falls back to a normal-approximation envelope

```math
d \;\pm\; K_\alpha \, d \, z / \sqrt{N},
```

where `d` is the reported deviation, `z` is the two-sided normal
quantile at the same confidence level (`z ≈ 1` at `CL = 0.683`), and
`K_α` is the SP1065 simple-interval noise factor —
`0.99, 0.99, 0.87, 0.77, 0.75` for `α = 2, 1, 0, −1, −2`
[riley-2008-sp1065](@cite). The lower limit is floored at zero, since a
deviation is non-negative by construction. The fallback keeps reported
bounds finite at the cost of optimism in the deep red-noise regime.

## Greenhall–Riley 2003 (GR03)

GR03 [greenhall-2003-edf-stability](@cite) gives a single numerical
algorithm — explicitly *not* a set of closed formulas — for the EDF of
the overlapped and non-overlapped ADEV/MDEV/HDEV/MHDEV estimators at
the power-law noise exponents `α ∈ {2, 1, 0, −1, −2, −3, −4}`, subject
to the restriction `α + 2d > 1`. SigmaTau implements the simplified
version of the algorithm (the truncated sum of GR03 Eqs. 10–13) in
`src/edf.jl` (`_calc_edf_core`), always with the overlapped stride
`S = m`. The five noise types SigmaTau's noise identification can
return (`α ∈ {2, …, −2}`) all satisfy the restriction for both
`d = 2` and `d = 3`.

The algorithm models the phase as the τ₀-difference of a pure power-law
process `w(t)` and chains three generalized autocovariance (GACV)
functions: `sw(t, α)`, the GACV of `w`; `sx(t, F, α) = F² Δ_{1/F} Δ_{−1/F} sw`,
the autocovariance of the (averaged) phase; and
`sz(t, F, α, d) = (Δ₁ Δ₋₁)^d sx`, the autocovariance of the `d`-th
difference process whose mean square is the variance being estimated.
The EDF of an estimator with `M` summands is then

```math
\frac{1}{\mathrm{EDF}} = \frac{1}{M\,s_z^2(0)}
\left[ s_z^2(0)
+ 2 \sum_{j=1}^{J-1} \Bigl(1 - \frac{j}{M}\Bigr) s_z^2\!\Bigl(\frac{j}{S}\Bigr)
+ \Bigl(1 - \frac{J}{M}\Bigr) s_z^2\!\Bigl(\frac{J}{S}\Bigr) \right],
```

with the lag sum truncated at `J = min(M, (d+1)·S)`. Symbols: `d` is
the phase-difference order (`d = 2` for the Allan family, `d = 3` for
the Hadamard family); `F` is the filter factor (`F = m` for the
unmodified variances, `F = 1` for the modified ones — this is what
captures the modified-family inner average); `S` is the stride factor
(`S = m`, the overlapped estimator); `M = 1 + ⌊S(N − L)/m⌋` counts the
estimator summands, with filter length `L = m/F + m·d`. SigmaTau
evaluates these sums directly rather than using the asymptotic
coefficient tables of GR03's "full version".

## Total-family EDF

GR03 explicitly excludes the total variances, and no GR03-style
formulas are published for them. SigmaTau uses published coefficient
tables where they exist and documented fallbacks elsewhere:

- **TOTDEV** — `EDF = b(α)·(T/τ) − c(α)` with the SP1065 Table 7
  coefficients for `α ∈ {0, −1, −2}` [riley-2008-sp1065](@cite). For
  WHPM/FLPM (`α = 2, 1`) no totvar-specific value is published, so the
  ADEV-style GR03 EDF (`d = 2`, `F = m`, `S = m`) is used as a
  pragmatic substitute.
- **MTOTDEV** — the same `b(α)·(T/τ) − c(α)` form with the SP1065
  Table 8 coefficients for all five `α ∈ {2, …, −2}`, valid for
  `τ ≥ 16τ₀` [riley-2008-sp1065](@cite). Below that floor (`m < 16`)
  MTOT reduces to MDEV, so the MDEV-style GR03 EDF (`d = 2`, `F = 1`)
  is used instead.
- **HTOTDEV** — `EDF = (T/τ) / (b₀(α) + b₁(α)·τ/T)` with the FCS 2001
  coefficients for `α ∈ {0, −1, −2, −3, −4}`
  [howe-2001-tothvar-steering](@cite); HDEV-style GR03 fallback
  (`d = 3`, `F = m`, `S = m`) at WHPM/FLPM.
- **MHTOTDEV** — novel to SigmaTau; no external reference exists. The
  `b(α)·(T/τ) − c(α)` coefficients were measured by Monte Carlo over
  the `τ/τ₀ ≥ 16` validity window — see
  [Theory: MHTOTDEV bias and EDF](mhtotdev_bias_edf.md).
- **TDEV / HTDEV / TTOTDEV** — rescale MDEV, MHDEV, and MTOTDEV by
  `τ/√3`, which does not change the degrees of freedom; the EDF and the
  χ² bounds are inherited from the wrapped estimator and scaled by the
  same factor.
- **PDEV** — uses the published PVAR EDF model of Vernotte–Chen–Rubiola
  2020 [vernotte-2020-pvar-noninteger](@cite); see
  [Theory: PDEV confidence](pdev_confidence.md).
- **MTIE** — no published EDF model; `ci` and `confidence` are accepted
  for API uniformity but `noise_type`, `edf`, and the CI bounds are
  always returned empty.

## Bias correction summary

| Estimator | B(α) applied? | Notes |
|-----------|---------------|-------|
| ADEV / MDEV / HDEV / MHDEV / TDEV / HTDEV / PDEV | none | Unbiased estimators |
| TOTDEV   | per SP1065 | `B = 1 − a(α)·τ/T` — biased low for FFM/RWFM |
| MTOTDEV / TTOTDEV | per SP1065 | τ-independent; B ≈ 1.27 (variance) under WHFM |
| HTOTDEV  | per FCS 2001 | `B = 1 + a(α)`, `α ∈ {0, …, −4}` — biased low |
| MHTOTDEV | Monte Carlo (SigmaTau) | `B = b₀ + b₁·τ/T`; ≈ 1 for white/flicker noise |

The bias factor `B(α) = E[estimator²] / true variance` is defined on
the variance; the API divides the raw *deviation* by `√B`, which is
equivalent. The correction is applied before the confidence interval is
evaluated, so the reported χ² bounds bracket the corrected deviation.
All four total-family estimators accept `correct_bias` (default
`true`); pass `correct_bias=false` to recover the raw kernel value.

The SP1065 TOTDEV factor `B = 1 − a(α)·(τ/T)` uses `a = 1/(3 ln 2) ≈ 0.481`
for FLFM and `a = 0.75` for RWFM, zero otherwise; `B < 1` and rises
toward 1 as `τ/T → 0` [riley-2008-sp1065](@cite). MTOTDEV uses a
τ-independent variance table `{1.06, 1.17, 1.27, 1.30, 1.31}` for
`α ∈ {2, 1, 0, −1, −2}` [riley-2008-sp1065](@cite)
[riley-2020-r-frequency-stability](@cite) — MTOT is biased *high*, so
the correction lowers σ. HTOTDEV uses the FCS 2001 form `B = 1 + a(α)`
with `a ∈ {−0.005, −0.149, −0.229, −0.283, −0.321}` for
`α ∈ {0, −1, −2, −3, −4}`; all `a < 0`, so `B < 1`, reflecting that
HTOT is biased *low* for the FM noises, while PM noises (`α > 0`) have
no published model and keep `B = 1` [howe-2001-tothvar-steering](@cite).
MHTOTDEV is *not* treated as unbiased: SigmaTau applies a
Monte-Carlo-measured `B = b₀(α) + b₁(α)·(τ/T)` that is ≈ 1 for white
and flicker noise but reads high for redder FM (RWFM `B ≈ 1.9` at small
τ) — see [Theory: MHTOTDEV bias and EDF](mhtotdev_bias_edf.md).

## Implementation contract

`ci=true` is the default for every deviation. The
`StabilityResult.edf` and `ci` vectors are populated only on that path;
pass `ci=false` to skip the χ² evaluation entirely and get them back
empty (`neff`, the per-τ analysis-window count, is populated either
way — it falls out of the kernel loop bounds for free). This is a
deliberate API contract: callers that don't need confidence intervals
pay no noise-identification or χ² cost. (One nuance: the total-family
estimators still populate `noise_type` under `ci=false` whenever
`correct_bias=true`, because the bias factor itself needs the per-τ
noise identification.)

```@example ci
using SigmaTau, Random
Random.seed!(1)
x = randn(2000)
r = adev(PhaseData(x, 1.0), [10, 100]; ci=true)
(r.edf, round.(ci_lower(r); sigdigits=3), round.(ci_upper(r); sigdigits=3))
```

## See also

- [Theory: Allan family](allan_family.md).
- [Theory: Total family](total_family.md).
- [Theory: PDEV confidence](pdev_confidence.md).
- [Theory: MHTOTDEV bias and EDF](mhtotdev_bias_edf.md).
- [API: `SigmaTau`](../reference/stab.md).

## References

- Greenhall & Riley, *Uncertainty of Stability Variances Based on
  Finite Differences*, PTTI 2003 [greenhall-2003-edf-stability](@cite).
- Howe et al., *Total Hadamard Variance*, FCS 2001 [howe-2001-tothvar-steering](@cite).
- SP1065 §5 [riley-2008-sp1065](@cite).
- Riley, *Frequency Stability Analysis Using R*, 2020
  [riley-2020-r-frequency-stability](@cite).
- Vernotte, Chen & Rubiola, PVAR degrees of freedom, 2020
  [vernotte-2020-pvar-noninteger](@cite).

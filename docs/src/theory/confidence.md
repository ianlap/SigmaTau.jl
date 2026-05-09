# Theory: Confidence Intervals

Each `StabilityResult.dev` value is a point estimate of the underlying
deviation. The corresponding confidence interval rests on the
chi-squared distribution applied to an *equivalent number of degrees of
freedom* (EDF) that depends on the estimator, the averaging factor `m`,
the record length `N`, and the noise type `α`.

## EDF as a bridge

The estimator is approximately a sum of `EDF` independent χ²(1)
contributions:

```math
\mathrm{EDF} \cdot \frac{\hat{\sigma}^2_y(\tau)}{\sigma^2_y(\tau)}
\;\sim\; \chi^2_{\mathrm{EDF}}.
```

For confidence level `CL` (default 0.683 ≈ 1σ), the lower and upper
bounds on `σ²_y(τ)` are:

```math
\sigma^2_{\text{lo}} = \frac{\mathrm{EDF}\cdot\hat\sigma^2}{\chi^2_{(1+\mathrm{CL})/2}},
\quad
\sigma^2_{\text{hi}} = \frac{\mathrm{EDF}\cdot\hat\sigma^2}{\chi^2_{(1-\mathrm{CL})/2}}.
```

`StabilityResult.ci_lower` and `ci_upper` are the deviations
`√σ²_{lo,hi}`. (Cite SP1065 §5 [@cite riley-2008-sp1065].)

The χ² approximation is asymptotic in EDF: it tightens as `N/m` grows
and as the noise type approaches white FM. At small EDF (long τ, short
records, or steep red noise) the interval is wide and asymmetric — the
upper bound stretches much further than the lower bound, reflecting
the long right tail of χ² at low degrees of freedom. SigmaTau falls
back to a normal-approximation envelope `d ± Kₙ·d·z/√N` when EDF is
non-finite or below 1, which keeps reported bounds finite at the cost
of optimism in the deep red-noise regime [@cite riley-2008-sp1065].

## Greenhall–Riley 2003 (GR03)

GR03 [@cite greenhall-2003-edf-stability] gives closed-form EDF expressions for
overlapping ADEV/MDEV/HDEV/MHDEV at any α, parameterized by three
dimensionless variances: `sz`, `sx`, `sw`. SigmaTau implements them in
`lib/SigmaTauStability/src/stats/edf.jl`. The formulas span several
pages in the reference; SigmaTau uses the published constants directly
without rederivation.

The decomposition factors the variance of the estimator into a noise
spectrum term (`sw`), a phase-difference shaping term (`sx`), and a
final estimator-difference term (`sz`) of order `d ∈ {2, 3}` for
ADEV-family vs HDEV-family. The averaging factor `F = m` for ADEV/HDEV
and `F = 1` for MDEV/MHDEV captures the modified-family inner average.
The same machinery yields the WPM/FLPM EDF used as a fallback for
TOTDEV and HTOTDEV when no total-specific table value is published.

The total-family estimators do not have published GR03-style EDF
formulas in general. SigmaTau falls back to ADEV/HDEV-style EDF for
WPM/FPM (α = 2, 1) and uses published `a(α)` values
[@cite howe-2001-tothvar-steering] for α ∈ {0, −1, −2}.

## Bias correction summary

| Estimator | B(α) applied? | Notes |
|-----------|---------------|-------|
| ADEV / MDEV / HDEV / MHDEV / TDEV / HTDEV | none | Unbiased estimators |
| TOTDEV   | per SP1065 | bias factor on the variance |
| MTOTDEV  | per SP1065 | ~1.27× under WFM |
| HTOTDEV  | per FCS01  | bias `a(α)` table — α = 0,−1,−2 |
| MHTOTDEV | none      | no published analytic model — HDEV-style fallback (limitation) |

Bias factors are applied to the *variance* before the square root and
before EDF lookup. The SP1065 TOTDEV factor `1 - a(α)·(τ/T)` shrinks
toward 1 as τ/T → 0 [@cite riley-2008-sp1065]. MTOTDEV uses a
τ-independent table `{1.06, 1.17, 1.27, 1.30, 1.31}` for
α ∈ {2, 1, 0, −1, −2} [@cite riley-2020-r-frequency-stability]. HTOTDEV uses the
FCS 2001 form `1 / (1 + a(α))` [@cite howe-2001-tothvar-steering]. MHTOTDEV is treated
as unbiased; Stable32 and AllanLab adopt the same convention.

## Implementation contract

The `StabilityResult.edf` field is a `Vector{Float64}` populated only
when the estimator is called with `calc_ci=true`. When `calc_ci=false`
(the default), it is empty. This is a deliberate API contract: callers
that don't need CI pay no χ² evaluation cost.

```@example ci
using SigmaTau, Random
Random.seed!(1)
x = randn(2000)
r = adev(PhaseData(x, 1.0), [10, 100]; calc_ci=true)
(r.edf, round.(r.ci_lower; sigdigits=3), round.(r.ci_upper; sigdigits=3))
```

## See also

- [Theory: Allan family](allan_family.md).
- [Theory: Total family](total_family.md).
- [API: SigmaTauStability](../reference/stability.md).

## References

- Greenhall & Riley, *Uncertainty of Stability Variances*, PTTI 2003
  [@cite greenhall-2003-edf-stability].
- Howe et al., *Total Hadamard Variance*, FCS 2001 [@cite howe-2001-tothvar-steering].
- SP1065 §5 [@cite riley-2008-sp1065].
- Riley, *Frequency Stability Analysis Using R*, 2020
  [@cite riley-2020-r-frequency-stability].

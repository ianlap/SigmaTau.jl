# Validation

Numerical validation is part of the methodological foundation: a deviation
estimator is only as useful as its agreement with the established
references the time-and-frequency community already trusts. SigmaTau is
cross-validated three ways against:

1. **Stable32** ([Riley & Howe 2008](@cite riley-2008-sp1065)) — desktop
   application, the de facto industry reference for time-and-frequency
   stability analysis. Output precision in the public fixtures is roughly
   five significant figures.
2. **allantools** (A. Wallin) — open-source Python library, second
   independent numerical reference. Output precision is full Float64.
3. **SigmaTau.jl** itself — implementations are pinned to the
   [Riley & Howe 2008](@cite riley-2008-sp1065) and
   [IEEE 1139-2022](@cite ieee1139-2022-definitions) definitions; bias
   corrections and EDF expressions follow
   [Greenhall & Riley 2003](@cite greenhall-2003-edf-stability).

Three-way agreement at the precision floor of the tightest reference is
the bar. Where the three disagree, the disagreement is documented and
attributed to a specific definitional choice (boundary handling, bias
correction policy), not a defect.

## Agreement classes

### Tight agreement (rtol ≤ 1e-5)

The primary kernels — overlapping ADEV, MDEV, TDEV, overlapping HDEV,
MHDEV — agree with both Stable32 and allantools to within the precision
of Stable32's published outputs (~5 significant figures). The
implementations are O(N) prefix-sum forms, exact-by-construction up to
floating-point round-off; against allantools, agreement holds to
rtol ≈ 1e-11 (full Float64 precision).

Representative comparison from `reference/validation/stable32out/`:

| Estimator | τ (s) | Stable32 | allantools | SigmaTau |
|---|---|---|---|---|
| OADEV | 1.0 | 1.00970e+00 | 1.00975e+00 | 1.00975e+00 |
| OADEV | 64  | 1.60850e-02 | 1.60852e-02 | 1.60852e-02 |
| MDEV  | 16  | 1.58180e-02 | 1.58178e-02 | 1.58178e-02 |
| TDEV  | 64  | 9.58680e-02 | 9.58679e-02 | 9.58679e-02 |
| OHDEV | 256 | 4.56720e-03 | 4.56718e-03 | 4.56718e-03 |

(All rows for OADEV, MDEV, TDEV, and OHDEV are tabulated on
[Validation: Stable32](../validation/stable32.md).)

### Documented offsets

The Total family — TOTDEV, HTOTDEV, MTOTDEV — has well-known definitional
choices that produce reproducible offsets relative to Stable32:

- **TOTDEV** at long τ: Stable32 and SigmaTau use different boundary
  reflection conventions for the doubly-extended phase record. Agreement
  is tight at short τ; disagreements at long τ are bounded and tracked
  per-estimator.
- **HTOTDEV**: SigmaTau applies the SP1065 white-FM bias factor B ≈ 1.005
  by default; Stable32 omits it. The remaining ~0.5% offset is exactly
  this bias correction.
- **MTOTDEV**: SigmaTau's biased output is approximately 1.27× larger
  than Stable32's. The factor matches SP1065's tabulated MTOT bias
  B(α=0) ≈ 1.27 (white FM); the underlying unbiased kernel matches
  Stable32 to within ~3%, confirming that the entire offset is
  attributable to the bias correction.

Both biased and unbiased SigmaTau outputs are exposed so users can
reproduce either convention.

### Confidence intervals

SigmaTau's χ² confidence intervals derive equivalent degrees of freedom
from [Greenhall & Riley 2003](@cite greenhall-2003-edf-stability), which
is more conservative than Stable32's legacy approximation. At small EDF
SigmaTau's intervals are slightly wider; the noise-type identification
underneath the EDF formula matches Stable32.

## See also

- [Validation: Methodology](../validation/methodology.md) — three-way
  reference framing and rtol-floor policy.
- [Validation: Stable32](../validation/stable32.md) — full per-estimator
  comparison tables across all τ.

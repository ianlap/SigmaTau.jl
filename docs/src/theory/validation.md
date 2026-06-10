# Validation

Numerical validation is part of the methodological foundation: a deviation
estimator is only as useful as its agreement with the established references
the time-and-frequency community already trusts. SigmaTau is cross-validated
against:

1. **Stable32** ([Riley & Howe 2008](@cite riley-2008-sp1065)) — desktop
   application, the de facto industry reference for time-and-frequency
   stability analysis. Output precision in the public fixtures is roughly
   five significant figures.
2. **allantools** (A. Wallin) — open-source Python library, second
   independent numerical reference. Output precision is full Float64.
3. **Inlined legacy kernels and internal identities** — implementations are
   pinned to the [Riley & Howe 2008](@cite riley-2008-sp1065) and
   [IEEE 1139-2022](@cite ieee1139-2022-definitions) definitions; bias
   corrections and EDF expressions follow
   [Greenhall & Riley 2003](@cite greenhall-2003-edf-stability).

Agreement at the precision floor of the tightest reference is
the bar. Where references disagree, the disagreement is documented and
attributed to a specific definitional choice (bias-correction policy,
noise-identification convention), not a defect.

## Agreement classes

### Tight agreement (fixture precision)

The primary kernels — overlapping ADEV, MDEV, TDEV, overlapping HDEV,
MHDEV — agree with both Stable32 and allantools to within the precision
of Stable32's published outputs (5 significant figures, rtol ≈ 5e-5). The
implementations are O(N) prefix-sum forms, exact-by-construction up to
floating-point round-off; against allantools, agreement holds to
rtol ≈ 1e-11 (full Float64 precision).

Representative comparison from `test/fixtures/validation/stable32out/`:

| Estimator | τ (s) | Stable32 | allantools | SigmaTau |
|---|---|---|---|---|
| OADEV | 1.0 | 1.00970e+00 | 1.00975e+00 | 1.00975e+00 |
| OADEV | 64  | 1.60850e-02 | 1.60852e-02 | 1.60852e-02 |
| MDEV  | 16  | 1.58180e-02 | 1.58178e-02 | 1.58178e-02 |
| TDEV  | 64  | 9.58680e-02 | 9.58679e-02 | 9.58679e-02 |
| OHDEV | 256 | 4.56720e-03 | 4.56718e-03 | 4.56718e-03 |

(All rows for OADEV, MDEV, TDEV, and OHDEV are tabulated on
[Validation: Stable32](../validation/stable32.md).)

### The total family

The Total family — TOTDEV, HTOTDEV, MTOTDEV — agrees with both references
once each program's bias-correction policy is accounted for. The raw
kernels (`correct_bias=false`) match allantools to machine precision
(≤ 4.4e-15 on the shipped fixture); against Stable32 the picture is
per-estimator:

- **TOTDEV**: both programs apply the Howe/Walter `B = 1 − a·τ/T`
  correction for the divergent FM noises. SigmaTau's default matches
  Stable32's σ to ≤ 0.005 % at every τ, including the longest.
- **HTOTDEV**: Stable32 applies the Howe 2005 bias correction by default
  (factors exactly `1/√(1+a)`), and so does SigmaTau — the default
  output matches Stable32 to ≤ 0.003 % at every matching-α τ.
- **MTOTDEV**: Stable32 reports the raw estimator; SigmaTau's raw kernel
  matches it to ≤ 0.004 % (the fixture's 5-significant-figure floor).
  SigmaTau's default additionally divides by `√B(α)` (SP1065 Table 11)
  to unbias toward the parent modified-Allan statistic, which places it
  3 – 14 % below Stable32 depending on the identified noise type. This
  is a deliberate policy difference, not a kernel disagreement; pass
  `correct_bias=false` to reproduce Stable32 (and allantools).

Both corrected and raw outputs are exposed so users can reproduce either
convention.

### Confidence intervals

SigmaTau's χ² confidence intervals derive equivalent degrees of freedom
from [Greenhall & Riley 2003](@cite greenhall-2003-edf-stability) — the
same algorithm Stable32 has used since v1.41 (stated in that paper,
p. 268). Wherever the two programs agree on the identified noise type,
the 68.3 % bounds agree to ≤ 0.05 % for the core estimators. The two
residuals are documented on
[Validation: Stable32](../validation/stable32.md): a noise-identification
detrend convention that flips one τ on the shipped fixture
(`identify_noise(x, m_values; detrend=false)` reproduces Stable32's α at
every τ), and total-family EDF substitutions above α = 0, where the
published TOTVAR/HTOT coefficient tables end.

## See also

- [Validation: Methodology](../validation/methodology.md) — three-way
  reference framing and rtol-floor policy.
- [Validation: Stable32](../validation/stable32.md) — full per-estimator
  comparison tables across all τ.

# Theory: Total Family

Total estimators improve confidence at long τ on short records. They
extend the input series at its boundaries before applying the same
finite-difference operator that defines the corresponding Allan-family
estimator. Four are in standard use: TOTDEV, MTOTDEV, HTOTDEV, MHTOTDEV.

## Why "total"?

Allan-family estimators waste data at long τ: with `m ≈ N/2`, the
overlapping estimator has only a handful of independent windows, so
its EDF (and hence its CI) is poor. Total estimators add reflected
copies of the input at the boundaries so the same number of windows is
available at every τ — buying confidence at long τ at the cost of a
small bias that depends on the noise type.

The reflection is performed on the phase record itself (or, for the
Hadamard total, on the frequency record). The original sequence of
length `N` becomes an extended sequence of length `3N − 4`, after
which the same finite-difference operator that defines ADEV, MDEV,
HDEV, or MHDEV is applied without further modification. This is what
the "total" prefix denotes: same operator, extended record. Confidence
improves because every τ now sees a window count comparable to the
short-τ regime, at the price of a small noise-type-dependent bias
that SP1065 quantifies with the `B(α)` table.

(Cite GHP99 [@cite Greenhall1999] for the original TOTVAR construction.)

## TOTVAR

The total Allan variance, on a phase series extended by reflection
to length `3N − 4`:

```math
\mathrm{TOTVAR}(\tau) \;=\; \frac{1}{2(m\tau_0)^2 (N-2)}
\sum_{i=2}^{N-1} \bigl(x_{i-m}^* - 2 x_i^* + x_{i+m}^*\bigr)^2
```

where `x*` is the reflection-extended series. The construction details
(reflection sign and length) are in `legdocs/equations/total.md` and
GHP99.

```julia
totdev(PhaseData(x, τ₀), τs)
```

## MTOTDEV

A phase-averaged variant of TOTVAR. Howe & Vernotte 1999. For each of
`N − 3m + 1` subsegments of length `3m`, the algorithm half-average
detrends the segment, applies symmetric reflection to length `9m`,
forms cumulative second differences, and averages over the `6m` valid
positions inside the extended segment. The estimator inherits MDEV's
sensitivity separation between WPM and FPM: where ADEV's spectrum is
flat under those two phase-noise regimes, MTOTDEV (like MDEV) resolves
them.

```julia
mtotdev(PhaseData(x, τ₀), τs)
```

(Cite HV99 [@cite Howe1999].)

## HTOTDEV

A third-difference total estimator, drift-insensitive like HDEV.
Originated in Howe 2000 and refined with bias coefficients in the
2001 FCS paper. For `m > 1` the algorithm converts phase to fractional
frequency, splits into `Ny − 3m + 1` segments of length `3m`,
half-average-detrends each, reflects to `9m`, and forms cumulative
Hadamard differences. The `m = 1` branch falls back to ordinary HDEV
as a documented exception — there is no useful reflection at the
shortest τ. Howe 2005 contributed enhancements to the long-τ behaviour
relative to the original 2000 construction.

```julia
htotdev(PhaseData(x, τ₀), τs)
```

(Cite H00 [@cite Howe2000]; the bias `a(α)` table is from FCS01
[@cite Howe2001]; long-τ refinements in [@cite Howe2005].)

## MHTOTDEV

A phase-averaged third-difference total estimator. There is no
dedicated canonical paper for MHTOTDEV; the construction follows the
HV99 modified-total methodology applied to the FCS01 Hadamard total.
For each of `N − 4m + 1` subsegments of phase length `3m + 1`, the
algorithm linearly detrends the segment, reflects symmetrically,
forms cumulative third differences, and applies an `m`-point moving
average. It is the long-τ equivalent of MHDEV (and of LDEV, which is
proportional to MHDEV by `√(3τ²/10)`).

```julia
mhtotdev(PhaseData(x, τ₀), τs)
```

## Bias correction policy

!!! note "Bias correction default"

    SigmaTau applies the SP1065 B(α) bias correction by default for
    MTOTDEV and HTOTDEV. Stable32 reports the *uncorrected* values for
    these estimators. This means SigmaTau's MTOTDEV is approximately
    1.27× higher than Stable32's MTOTDEV under white FM (α = 0). The
    underlying kernel without B(α) matches Stable32 to ~3%.

    For numeric comparison and the per-α bias factors, see
    [Validation: Stable32](../validation/stable32.md).

(Cite SP1065 §5 [@cite RileyHowe2008] for the B(α) tables.)

## Demonstration

```@example total
using SigmaTau, Random
Random.seed!(7)

# Short record, comparing ADEV vs TOTVAR confidence at the largest τ
N = 1024
x = cumsum(randn(N))   # WFM
τs = [1, 4, 16, 64, 256]

a = adev(PhaseData(x, 1.0), τs; calc_ci=true)
t = totdev(PhaseData(x, 1.0), τs; calc_ci=true)

# CI half-width at the largest τ — TOTVAR should be tighter
ci_half(r, i) = (r.ci_upper[i] - r.ci_lower[i]) / 2
last_τ = length(τs)
round.((ci_half(a, last_τ), ci_half(t, last_τ)); sigdigits=3)
```

The TOTVAR half-width at the largest τ should be smaller than ADEV's
on this short record.

## Implementation notes

- The reflection-extension implementation lives in
  `lib/SigmaTauStability/src/core/total.jl`.
- Bias correction is applied in
  `lib/SigmaTauStability/src/api/total.jl` via the `bias_correction`
  helper from `lib/SigmaTauStability/src/stats/edf.jl`.
- The MHTOT EDF model uses an HDEV-style approximation (no published
  analytic form); known limitation tracked as `R-MED-6`.

## See also

- [Theory: Allan family](allan_family.md).
- [Theory: Confidence](confidence.md).
- [Validation: Stable32](../validation/stable32.md).

## References

- TOTVAR: Greenhall, Howe & Percival 1999 [@cite Greenhall1999].
- MTOT: Howe & Vernotte 1999 [@cite Howe1999].
- HTOT: Howe 2000 [@cite Howe2000]; FCS01 [@cite Howe2001];
  Howe 2005 [@cite Howe2005].
- SP1065 §5 [@cite RileyHowe2008].

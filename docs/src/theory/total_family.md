# Theory: Total Family

A *total* estimator is one that operates on a **data-extended** version
of the input record. The same finite-difference operator that defines
its corresponding Allan-family estimator (ADEV, MDEV, HDEV, or MHDEV)
is then applied without modification — only the input changes. The
purpose is statistical: at long τ, an Allan-family estimator runs out
of independent windows and its EDF (and confidence interval) collapse;
extending the record gives every τ a comparable window count.

The four total estimators in standard use — TOTDEV, MTOTDEV, HTOTDEV,
MHTOTDEV — **do not share a single extension method.** Each has its
own scheme, varying in three dimensions: whether the extension is
applied to the whole record or to per-τ subsegments, whether the
segment is detrended first and how, and whether the operation is in
phase or fractional-frequency. Every section below names the specific
scheme for that estimator.

The cost of extension is a small noise-type-dependent bias `B(α)`,
quantified per estimator in SP1065 §5 [riley-2008-sp1065](@cite). SigmaTau
applies it by default; see the bias-correction policy section below.

---

## TOTDEV — total deviation

The Total deviation (TOTDEV) employs a whole-record symmetric reflection at both ends of the original phase series, producing an extended sequence. This reflection is applied exactly once to the full record without prior detrending or subsegmentation. The standard overlapping ADEV second-difference operator is subsequently run over this extended sequence. The primary advantage of this boundary extension is a substantial improvement in long-τ statistical confidence, as the second-difference operator is never starved for data near the record boundaries.

For the formal mathematical definition of the TOTDEV estimator, please read NIST SP1065 §5.2.10 [riley-2008-sp1065](@cite) for more details.

```julia
totdev(PhaseData(x, τ₀), τs)
```

---

## MTOTDEV — modified total deviation

The Modified total deviation (MTOTDEV) utilizes a per-subsegment extension strategy rather than a whole-record reflection. The estimator iterates over overlapping subsegments of the phase record. Within each subsegment, it performs half-average detrending followed by symmetric boundary reflection, subsequently forming cumulative second differences over the extended segment. This methodology effectively inherits MDEV's capacity to separate WPM and FPM noise via phase-averaging, while providing the long-τ confidence benefits of a total estimator.

For the formal mathematical definition of the MTOTDEV estimator and a complete algorithmic breakdown of the per-subsegment extension, please read NIST SP1065 §5.2.11 [riley-2008-sp1065](@cite) for more details.

```julia
mtotdev(PhaseData(x, τ₀), τs)
```

---

## TTOT — time-total deviation

TTOT (Time-Total Deviation) is the time-deviation rescaling of MTOTDEV,
analogous to how TDEV rescales MDEV [riley-2008-sp1065](@cite). It
gives a `σ_x`-units summary (seconds) of long-τ stability for a record
with WPM-dominated noise, leveraging MTOTDEV's per-subsegment extension
to extend the usable τ range beyond TDEV's reach on short records
[banerjee-2023-timekeeping](@cite).

```math
\sigma_{x,\text{TTOT}}(\tau) \;=\; \frac{\tau}{\sqrt{3}}\,
\sigma_{y,\text{MTOT}}(\tau).
```

The `√3` factor matches TDEV's exactly because TTOT inherits MTOTDEV's
modified second-difference operator [riley-2008-sp1065](@cite).

In SigmaTau:

```julia
ttotdev(PhaseData(x, τ₀), τs)
```

[`ttotdev`](@ref) wraps [`mtotdev`](@ref) and rescales the centerline
and CI bounds by `τ / √3`. The `correct_bias`, `ci`, and
`confidence` kwargs pass through unchanged; the `edf` and `noise_type`
columns are reused as-is since a time rescaling does not change the
degrees of freedom.

---

## HTOTDEV — Hadamard total deviation

The Hadamard total deviation (HTOTDEV) diverges from the preceding total estimators by operating in the frequency domain rather than the phase domain. The estimator converts the phase record to fractional frequency before walking over overlapping subsegments. Each subsegment undergoes half-average detrending (which effectively absorbs linear frequency drift) and symmetric reflection, followed by the application of cumulative third differences. This configuration successfully preserves HDEV's intrinsic insensitivity to linear frequency drift while delivering the enhanced long-τ confidence characteristic of total estimators.

Note that at the shortest averaging factor ($m=1$), symmetric reflection fails to produce a statistically useful extension. Consequently, SigmaTau falls back to the ordinary HDEV at $m=1$, strictly adhering to the FCS01 algorithm design. 

For the comprehensive mathematical formulation of the HTOTDEV estimator, its specific drift-handling mechanics, and its documented bias coefficients, please read NIST SP1065 §5.2.12 [riley-2008-sp1065](@cite) for more details.

```julia
htotdev(PhaseData(x, τ₀), τs)
```

---

## MHTOTDEV — modified Hadamard total deviation

The Modified Hadamard total deviation (MHTOTDEV) is an original contribution to SigmaTau. It is the long-τ extension of MHDEV, combining phase-averaged third differences with boundary extension. Because there is no canonical publication defining MHTOTDEV, its construction adapts the modified-total methodology to the Hadamard total framework, integrating the modified third-difference operator.

The estimator uses a per-subsegment **quadratic least-squares** detrend on the phase record. It processes overlapping subsegments by fitting and removing a quadratic in phase (offset + frequency offset + drift), applying symmetric reflection, and computing cumulative third differences with an internal moving average.

The detrend degree matches the kernel's annihilation degree: the third difference annihilates phase polynomials of degree ≤ 2, just as MTOTDEV's second difference annihilates degree ≤ 1 and its linear (half-average) detrend matches. A simple slope removal would leave the stochastic quadratic from α ≤ −3 noise in the window, and the seam kink at the reflection boundary would feed the boundary-spanning third differences, collapsing the kernel's low-frequency response. With the matched quadratic detrend, the estimator converges over the full Hadamard range α > −5.

Because the quadratic fit exactly removes any deterministic frequency offset and drift inside every window, **this kernel is drift-immune on its own**. The [`mhtotdev`](@ref) wrapper still removes the record's least-squares frequency drift (`remove_drift=true` by default) so that the per-τ *noise identification* sees the underlying noise rather than the ramp. The statistic itself does not need it.

The formal definition computed by the kernel is:

```math
\mathrm{MHTOTVAR}(\tau) \;=\; \frac{1}{(m\tau_0)^2 \,(N - 4m + 1)}
\sum_{n} \frac{1}{6\,m^2 \,n_{\text{avg}}}
\sum_{\text{window}} \bigl(\text{third-diff avg}\bigr)^2 ,
```

where the inner sum runs over $n_{\text{avg}}$ valid $m$-point average windows of the cumulated third-difference sequence on the extended segment.

!!! info "Original contribution"
    MHTOTDEV is original to SigmaTau, completing the total estimator family. The authoritative definition lives in the package source (`src/kernels.jl` and `src/deviations.jl`), providing the missing link alongside TOTDEV, MTOTDEV, and HTOTDEV.

```julia
mhtotdev(PhaseData(x, τ₀), τs)
```

**Where it sits.** Across the eight time-domain stability estimators
that SigmaTau implements, MHTOTDEV is the all-three corner of the
cube formed by three independent design axes:

- **Difference order** — second (Allan) vs third (Hadamard).
- **Inner averaging** — none (standard) vs phase-averaged (modified).
- **Boundary handling** — none vs total-style extension.

MHTOTDEV is Hadamard × modified × total, so it carries every
property that distinguishes the family: drift rejection (by the
per-window quadratic detrend, in the kernel itself) and `α ∈ {−4, −3}`
convergence from Hadamard, WPM/FPM disambiguation from the
phase-averaged inner kernel, and tight long-τ confidence from the
per-subsegment extension. It is the right tool when a record has
linear frequency drift, divergent low-frequency noise, ambiguity
between WPM and FPM, *and* a record length too short to let any of
the simpler estimators reach long τ with usable confidence.

!!! note "Noise identification and drift on the other estimators"
    `mhtotdev` feeds its noise identification the drift-removed
    record, so its bias correction and χ² intervals classify the
    underlying noise. Every *other* estimator's CI machinery runs on
    the record as given — on a strongly drifting record the
    identification can classify long-τ points as redder than the true
    noise even for drift-immune deviations like `hdev`/`mhdev`. For
    confidence intervals on a drifting record, remove the
    deterministic drift first (`detrend(fd)` on the frequency record,
    or `detrend(pd; method=:quadratic)` on phase), the same hygiene
    SP1065 recommends generally. See
    [the drifting-clock tutorial](../tutorials/03_drifting_clock_hadamard.md).

---

## Summary: extension scheme by estimator

| Estimator | Where | Pre-extension detrend | Domain | m=1 special case |
|-----------|-------|-----------------------|--------|------------------|
| TOTDEV    | whole record | none           | phase     | no  |
| MTOTDEV   | per subsegment | half-average (offset) | phase     | no  |
| HTOTDEV   | per subsegment | half-average (drift, on `y`) | frequency | yes — falls back to HDEV |
| MHTOTDEV  | per subsegment | quadratic LS (offset + drift)¹ | phase     | no  |

¹ Matched to the third difference's annihilation degree, which makes the
kernel drift-immune on its own. The public wrapper additionally removes the
record's least-squares drift (`remove_drift=true` default) purely as
noise-identification hygiene — see the drift-handling note above.

---

## Bias correction policy

!!! note "Bias correction default"

    SigmaTau applies the published bias corrections by default
    (`correct_bias=true`) for the whole total family. Stable32's policy
    is per-estimator:

    - **TOTDEV**: Stable32 applies the Howe/Walter `B = 1 − a·τ/T`
      correction (FLFM/RWFM only). SigmaTau's default matches it;
      `correct_bias=false` reproduces allantools' raw kernel.
    - **HTOTDEV**: Stable32 applies the Howe 2005 correction (factors
      exactly `1/√(1+a)`). SigmaTau's default matches Stable32 to
      ≤ 0.003 % on the validation fixture; `correct_bias=false`
      reproduces allantools instead.
    - **MTOTDEV**: Stable32 reports the *uncorrected* estimator.
      SigmaTau's raw kernel matches it to ≤ 0.004 %; the default `√B(α)`
      correction (SP1065 Table 11) places SigmaTau's output 3 – 14 %
      below Stable32 depending on noise type. Pass `correct_bias=false`
      to reproduce Stable32 (and allantools).

    In short: `correct_bias=false` reproduces Stable32 for MTOTDEV but
    *not* for TOTDEV or HTOTDEV, where the default already matches
    Stable32 and `correct_bias=false` matches allantools. Side-by-side
    numbers are on [Validation: Stable32](../validation/stable32.md).

(Cite SP1065 §5 [riley-2008-sp1065](@cite) and FCS01 [howe-2001-tothvar-steering](@cite) for
the `a(α)` / `B(α)` tables.)

---

## Demonstration

```@example total
using SigmaTau, Random
Random.seed!(7)

# Short record, comparing ADEV vs TOTVAR confidence at the largest τ
N = 1024
x = cumsum(randn(N))   # WFM
τs = [1, 4, 16, 64, 256]

a = adev(PhaseData(x, 1.0), τs; ci=true)
t = totdev(PhaseData(x, 1.0), τs; ci=true)

# CI half-width at the largest τ — TOTVAR should be tighter
ci_half(r, i) = (r.ci[i].hi - r.ci[i].lo) / 2
last_τ = length(τs)
round.((ci_half(a, last_τ), ci_half(t, last_τ)); sigdigits=3)
```

The TOTVAR half-width at the largest τ should be smaller than ADEV's
on this short record — the data-extension is doing its job.

---

## Implementation notes

- All four total kernels live in `src/kernels.jl`.
- Bias correction is applied in `src/deviations.jl` via the
  `bias_correction` helper from `src/edf.jl`.
- The MHTOTDEV EDF model `edf = b·(T/τ) − c` and the bias
  `B = b₀ + b₁·(τ/T)` are Monte-Carlo-measured (no published form
  exists for MHTOTDEV); the currently shipped coefficients are
  provisional quick-mode values for the quadratic-detrend kernel. See
  [Theory: MHTOTDEV bias and EDF](mhtotdev_bias_edf.md).

---

## See also

- [Theory: Allan family](allan_family.md) — the corresponding
  non-extended estimators.
- [Theory: Confidence](confidence.md) — EDF and χ² intervals.
- [Validation: Stable32](../validation/stable32.md) — numerical
  comparison and the `B(α)` table.

---

## References

- TOTVAR: Greenhall, Howe & Percival 1999 [greenhall-1999-totvar](@cite).
- MTOT: Howe & Vernotte 1999 [howe-1999-modtotvar](@cite).
- HTOT: Howe 2000 [howe-2000-tothvar-ptti](@cite); FCS01 [howe-2001-tothvar-steering](@cite);
  Howe 2005 [howe-2005-tothvar-ieee](@cite).
- SP1065 §5 [riley-2008-sp1065](@cite).

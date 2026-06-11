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

**Extension method.** Whole-record symmetric reflection at both ends
of the phase series `x` of length `N`, producing an extended sequence
`x*` of length `3N − 4`. The reflection is applied **once** to the full
record, with no per-τ subsegmentation and no detrending. The
overlapping ADEV second-difference operator then runs on `x*` exactly
as it would on the original record.

**Definition** (SP1065 §5 Eq. 25 [riley-2008-sp1065](@cite);
GHP99 [greenhall-1999-totvar](@cite)):

```math
\mathrm{TOTVAR}(\tau) \;=\; \frac{1}{2\,(m\tau_0)^2 \,(N-2)}
\sum_{i=2}^{N-1} \bigl(x_{i-m}^{*} - 2\,x_{i}^{*} + x_{i+m}^{*}\bigr)^2 .
```

The sum index runs over the original interior positions; the `x_{i±m}^*`
terms reach into the reflected ends when `i ± m` would otherwise leave
the original record. Long-τ confidence improves substantially over ADEV
because the second difference always has data to operate on, regardless
of how close `i` sits to a boundary.

```julia
totdev(PhaseData(x, τ₀), τs)
```

---

## MTOTDEV — modified total deviation

**Extension method.** *Per-subsegment*, not whole-record. The estimator
walks `N − 3m + 1` overlapping subsegments of length `3m` along the
phase record. For each subsegment the algorithm:

1. **Half-average detrends** the subsegment (subtract the mean of the
   first half from the first half, the mean of the second half from
   the second half — preserves the offset between halves so the inner
   averaging step still sees real phase information).
2. **Symmetrically reflects** the detrended segment at both ends to
   build an extended segment of length `9m` (three copies — reversed,
   original, reversed).
3. Forms cumulative second differences and averages over the `6m`
   valid window positions inside the extension.

The outer sum is over the `N − 3m + 1` subsegments. This is HV99's
construction [howe-1999-modtotvar](@cite); SP1065 §5 reproduces the algorithm.

```math
\mathrm{MTOTVAR}(\tau) \;=\; \frac{1}{2\,(m\tau_0)^2 \,(N - 3m + 1)}
\sum_{n} \frac{1}{6m} \sum_{j} \bigl(a_{j+2} - 2\,a_{j+1} + a_{j}\bigr)^2 ,
```

where `a_j` are the `m`-point averages of the cumsum-reflected extended
segment for subsegment `n`, and the inner sum runs over `6m` valid
positions inside that extended segment.

MTOTDEV inherits MDEV's WPM/FPM separation (slope `−3/2` vs `−1`
under deviation) — phase-averaging the inner second-difference
window splits the degeneracy that ADEV / TOTDEV cannot.

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

**Extension method.** *Per-subsegment in the frequency domain.* This
is the first total estimator that does not operate directly on the
phase record. Two-stage:

1. **Phase → frequency.** Convert `x` to fractional frequency
   `y[i] = (x[i+1] − x[i]) / τ₀`, length `N_y = N − 1`.
2. **Per-subsegment processing.** Walk `N_y − 3m + 1` subsegments of
   length `3m` along `y`. For each:
   - Half-average detrend the segment.
   - Symmetric reflection to length `9m` (`[reverse; segment; reverse]`).
   - Cumulative third differences (Hadamard kernel) summed over `6m`
     positions.

**Documented exception at `m = 1`:** at the shortest averaging factor
there is no useful reflection (the third-difference window is `3m + 1 = 4`
samples; reflection would produce trivially-correlated extensions).
SigmaTau falls back to ordinary HDEV at `m = 1`. This matches the
reference HTOTDEV implementations and is not a numerical defect — it
is the design of the FCS01 algorithm [howe-2001-tothvar-steering](@cite).

The construction originates in Howe 2000 [howe-2000-tothvar-ptti](@cite) with bias
coefficients `a(α)` introduced in the FCS01 paper [howe-2001-tothvar-steering](@cite);
Howe 2005 [howe-2005-tothvar-ieee](@cite) contributed long-τ refinements.

HTOTDEV inherits HDEV's drift insensitivity (the third-difference
kernel annihilates linear-in-`t` terms) so a record with significant
linear frequency drift can be totalized without first detrending the
whole record — drift is absorbed into the per-subsegment half-average
step.

```julia
htotdev(PhaseData(x, τ₀), τs)
```

---

## MHTOTDEV — modified Hadamard total deviation

**Extension method.** *Per-subsegment, half-average detrending* — the
same Greenhall extension MTOTDEV and HTOTDEV use. (MHTOTDEV has no
canonical/published form; SigmaTau adopts the Greenhall methodology for
consistency across the modified-total family.) Walks `N − 4m + 1`
subsegments of phase length `4m` along `x`. For each:

1. **Half-average slope removal** — subtract a linear trend whose slope
   is the difference of the two half-segment means (the Greenhall slope
   estimate shared with MTOTDEV / HTOTDEV), removing the subsegment's
   mean frequency offset before the third-difference kernel.
2. Symmetric reflection.
3. Cumulative third differences plus an `m`-point moving average
   (the modified-style inner averaging, parallel to how MDEV relates
   to ADEV).

**Drift handling.** The kernel above is *not* drift-immune, unlike its
parent MHDEV: the per-window detrend removes only the local frequency
offset, so residual deterministic drift makes the reflected copies
non-smooth at the window boundaries and the differences that span a
boundary pick it back up — a record drifting at `2e-10`/day over a
`1e-11` WHFM floor read 6× high at `T/τ ≈ 10` and 30× high at
`T/τ ≈ 5`. Removing the drift *per window* does not work: the
phase-domain two-parameter detrend inflates the WHPM bias to
`B ≈ 3.8` and collapses the χ² EDF model (R² = 0.32), and the
frequency-domain TOTHVAR scheme gives `B(WHPM)` growing with `m` —
both were measured and rejected. Instead, [`mhtotdev`](@ref) removes
the record's **global least-squares frequency drift** before the
kernel (`remove_drift=true` by default; the per-τ noise identification
also runs on the drift-removed record). The least-squares fit is a
linear projection, so a deterministic drift is removed exactly, and
the bias/EDF coefficients below are measured with this removal in the
loop — the calibration matches what the function computes. Pass
`remove_drift=false` for data that is already detrended.

```math
\mathrm{MHTOTVAR}(\tau) \;=\; \frac{1}{(m\tau_0)^2 \,(N - 4m + 1)}
\sum_{n} \frac{1}{6\,m^2 \,n_{\text{avg}}}
\sum_{\text{window}} \bigl(\text{third-diff avg}\bigr)^2 ,
```

with the per-subsegment inner sum running over `n_{\text{avg}}`
valid `m`-point average windows of the cumulated third-difference
sequence on the extended segment. The exact normalization and window
counts are pinned against a direct reference implementation in
`test/stab/legacy_kernels.jl`.

MHTOTDEV is the long-τ extension of MHDEV: phase-averaged third
differences with boundary extension.

!!! info "Original contribution"
    MHTOTDEV is original to SigmaTau. There is no
    canonical paper for it; the construction follows HV99's
    modified-total methodology [howe-1999-modtotvar](@cite) applied
    to the FCS01 Hadamard total [howe-2001-tothvar-steering](@cite),
    with the modified third-difference operator from Greenhall 1997
    [greenhall-1997-third-difference-mvar](@cite). The authoritative
    definition lives in the package source itself; the kernel sits at
    [`src/kernels.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/src/kernels.jl)
    and the public wrapper at
    [`src/deviations.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/src/deviations.jl).
    Put differently: total forms were published for the Allan
    deviation (TOTDEV), the modified Allan deviation (MTOTDEV), and
    the Hadamard deviation (HTOTDEV). The modified Hadamard deviation
    had none — MHTOTDEV is that missing total form.

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
automatic least-squares drift removal above) and `α ∈ {−4, −3}`
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
| MHTOTDEV  | per subsegment | half-average (offset)¹ | phase     | no  |

¹ Plus a record-level least-squares drift removal in the public wrapper
(`remove_drift=true` default) — see the drift-handling note above.

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
  exists for MHTOTDEV); see
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

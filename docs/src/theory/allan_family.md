# Theory: Allan Family

The Allan family covers the six estimators built from second- and
third-difference operators on phase data: ADEV, MDEV, TDEV, HDEV, MHDEV,
and HTDEV. They share the same boundary handling (no extension) and
differ in difference order, phase-averaging, and final scaling.

## ADEV — overlapping Allan deviation

The Allan deviation (ADEV) is the classical two-sample variance. It is degenerate for White and Flicker Phase Modulation (both yield a variance slope $\mu = -1$), so it cannot independently disambiguate between these two noise types. 

For the mathematical definition and the overlapping estimator form, please read NIST SP1065 §5 [riley-2008-sp1065](@cite) for more details.

```julia
adev(PhaseData(x, τ₀), [1, 2, 4, 8, 16])
```

## MDEV — modified Allan deviation

The Modified Allan deviation (MDEV) adds an inner phase-averaging step to the standard ADEV calculation, which decouples WPM ($\mu = -3/2$) from FPM ($\mu = -1$).

Internally, SigmaTau computes MDEV using Greenhall's third-difference construction via prefix sums [greenhall-1997-third-difference-mvar](@cite). For the canonical definition, please read NIST SP1065 §5.2 [riley-2008-sp1065](@cite) for more details.

```julia
mdev(PhaseData(x, τ₀), τs)
```

## PDEV — parabolic deviation

The parabolic variance PVAR (with `PDEV = √PVAR`) is Vernotte,
Lenczner, Bourgeois & Rubiola's wavelet variance built on a
least-squares fit: instead of the finite-difference kernels of the
Allan family, the estimator weights the phase with a quadratic
(parabolic) window [vernotte-2016-pvar](@cite). The least-squares
construction gives it better long-τ noise-type discrimination than
MVAR at comparable confidence, which is the property that motivated it
[vernotte-2016-pvar](@cite).

```math
\mathrm{PVAR}(\tau) \;=\; \frac{2}{N_e\,m^2\,\tau_0^2}
\sum_{i} \biggl[\sum_{k=0}^{m} c_k\,x_{i+k}\biggr]^2 ,
```

with `c_k` the quadratic weighting kernel (parabolic shape over the
window) constructed so that the operator integrates to zero against
both constants and linear ramps but responds nonzero to quadratic
curvature in `x(t)` [vernotte-2016-pvar](@cite).

!!! note "Implementation status"
    Implemented as [`pdev`](@ref) (parabolic deviation σ = √PVAR);
    `m=1` reproduces ADEV exactly and a linear-trend annihilation test
    pins the kernel. χ² confidence intervals use the
    Vernotte–Chen–Rubiola closed-form EDF
    [vernotte-2020-pvar-noninteger](@cite); see
    [Theory: PDEV confidence](pdev_confidence.md).

## TDEV — time deviation

The Time deviation (TDEV) measures time-error stability in units of seconds (a $\sigma_x$ quantity, not $\sigma_y$). It is a scaled version of MDEV.

TDEV is the standard tool for characterizing the time error of a clock or time-distribution system (e.g., telecommunications synchronization networks, GNSS time transfer, atomic-clock comparison links). For a detailed definition, please read NIST SP1065 §5.2.7 [riley-2008-sp1065](@cite) for more details.

```julia
tdev(PhaseData(x, τ₀), τs)
```

## HDEV — overlapping Hadamard deviation

The Hadamard deviation (HDEV) uses a third-difference kernel with two advantages over ADEV: it automatically filters out linear frequency drift, and its variance integral converges down to $\alpha = -4$. It is the preferred estimator for clocks with known frequency drift (e.g., Cesium, Rubidium) and for characterizing long power-law tails.

For the overlapping Hadamard estimator definition and a demonstration of its drift-insensitivity, please read NIST SP1065 §5.4 [riley-2008-sp1065](@cite) for more details.

```julia
hdev(PhaseData(x, τ₀), τs)
```

## MHDEV — modified Hadamard deviation

The Modified Hadamard deviation (MHDEV) is a phase-averaged third difference. The relationship between MHDEV and HDEV mirrors MDEV and ADEV: the phase-averaging window splits the WPM and FPM degeneracy, while the third-difference kernel preserves drift insensitivity. MHDEV is the right choice when a record contains linear frequency drift and phase noise that ADEV or MDEV cannot disambiguate.

For the exact formulation, please read NIST SP1065 §5.4 [riley-2008-sp1065](@cite) for more details.

```julia
mhdev(PhaseData(x, τ₀), τs)
```

## Theo1 — theoretical variance #1

Theo1 (Theoretical Variance #1) extends the usable averaging-time range of a stability estimator by sampling longer windows than ADEV permits. It extracts stability information from internal cross-products over windows of length `(N − 1)·τ₀`, yielding a factor-of-3 extension in the long-τ reach for a fixed record length. It has a noise-type-dependent bias that is corrected via the Theo1 bias table.

For the Theo1 definition, its bias-correction procedures, and EDF, please read NIST SP1065 §5.2.14 [riley-2008-sp1065](@cite) for more details.

!!! note "Implemented"
    `theo1(data, m_values; ci=true, correct_bias=true, confidence=…)`,
    with the usual `TauMode` / zero-argument / `FrequencyData` dispatch
    chains. `m_values` are even Theo1 averaging factors
    (`tau_values(mode, N, :theo1)` rounds a grid to even and dedupes);
    `r.tau` holds the effective `0.75·m·τ₀`. `correct_bias=true`
    (default) applies the ThêoBR automatic bias removal — the average of 
    the matched-τ `AVAR/Theo1` variance ratios. EDF/CIs use the 
    per-noise-type empirical formulas of SP1065 §5.2.15. The raw kernel 
    is cross-validated against allantools 2024.06 to machine precision.

## ThêoH — composite Theo1 + ADEV estimator

ThêoH (pronounced "theo-H") is a composite estimator combining overlapping ADEV at short averaging factors with bias-corrected Theo1 (ThêoBR) at longer averaging factors. It produces a single continuous curve from `m = 1` out to the Theo1 long-τ limit. The crossover point ensures tight confidence intervals across the entire range that neither ADEV nor Theo1 alone could provide. Banerjee & Matsakis recommend ThêoH as the default long-record stability summary [banerjee-2023-timekeeping](@cite).

For the construction of the composite estimator and the calculation of the crossover point, please read NIST SP1065 §5.2.16 [riley-2008-sp1065](@cite) for more details.

!!! note "Implemented"
    `theoh(data, m_values; ci=true, confidence=…)`, with the usual
    `TauMode` / zero-argument / `FrequencyData` dispatch chains.
    `m_values` are τ-grid factors (target `τ = m·τ₀`, not necessarily
    even): points with `m·τ₀ ≤ k` come from overlapping ADEV at `m`,
    points above from ThêoBR at the even Theo1 factor nearest `4m/3`.
    CIs use the Greenhall–Riley ADEV EDF on the ADEV segment and the 
    SP1065 Table 3 Theo1 EDF on the ThêoBR segment.

## Dynamic Allan deviation — DADEV

The dynamic Allan deviation generalizes σ_y(τ) to a
**time-resolved** stability map: instead of a single scalar curve
versus τ, DADEV produces a 2-D surface σ_y(t, τ) by sliding an
analysis window across the record and computing ADEV inside each
window [mckelvy-2025-telemetrystability](@cite). The original
construction is due to Galleani & Tavella (2009); McKelvy et al.
2025 apply it to telemetry-based stability estimation of high-
precision atomic clocks and demonstrate its ability to localize
non-stationary stability events that a static ADEV averages over.

The estimator structure is

```math
\sigma_y^2(t_c,\,\tau) \;=\; \frac{1}{2(N_w - 2m)\,(m\tau_0)^2}
\sum_{i \in W(t_c)} \bigl(x_{i+2m} - 2 x_{i+m} + x_i\bigr)^2 ,
```

where `W(t_c)` is an analysis window of `N_w` phase samples centered
at coordinate time `t_c`. The window slides across the record to
produce DADEV at every `t_c` of interest
[mckelvy-2025-telemetrystability](@cite).

!!! note "Implemented"
    `dadev(data, m_values; window, step=window÷2)`, with the usual
    `TauMode` / `FrequencyData` dispatch chains (a `TauMode` grid is
    clamped to the *window* length, not the record length). Each window
    is evaluated by the same overlapping-ADEV kernel that `adev` uses,
    so the `window = N` map degenerates to the static estimate exactly.
    Returns a [`DynamicStabilityResult`](@ref) — a `t × τ` matrix with
    `NaN` wherever a window cannot support the averaging factor. No
    confidence intervals are produced: the literature gives no EDF
    model for the time-resolved map.

## Dynamic Hadamard deviation — DHDEV

The dynamic Hadamard deviation is the same sliding-window construction
as DADEV with the overlapping third-difference (Hadamard) kernel inside
each window [mckelvy-2025-telemetrystability](@cite):

```math
\sigma_{H}^2(t_c,\,\tau) \;=\; \frac{1}{6(N_w - 3m)\,(m\tau_0)^2}
\sum_{i \in W(t_c)} \bigl(x_{i+3m} - 3 x_{i+2m} + 3 x_{i+m} - x_i\bigr)^2 .
```

Because the third difference annihilates terms quadratic in `t`, each
window's estimate rejects linear frequency drift *within that window*
(the same property HDEV adds over ADEV; SP1065 §5
[riley-2008-sp1065](@cite)) — a drifting clock's stability map shows
its noise floor rather than the drift ramp, without pre-detrending.

!!! note "Implemented"
    `dhdev(data, m_values; window, step=window÷2)`, parallel to
    `dadev` in every respect (dispatch chains, window-clamped `TauMode`
    grids, `DynamicStabilityResult` return, no CI machinery). The
    window must keep `N_w − 3m ≥ 2`, so the supported τ range at a
    given window is 2/3 of DADEV's.

## HTDEV — Hadamard time deviation

HTDEV is to MHDEV what TDEV is to MDEV: a τ-scaled time-domain
deviation built on the third-difference (Hadamard) kernel instead of
the second-difference (Allan) kernel. Like TDEV, HTDEV has units of
seconds (it is a `\sigma_x` quantity). SigmaTau implements it as
`mhdev` followed by the scaling `τ / √(10/3)`:

```math
\sigma_{x,\mathrm{HT}}(\tau) \;=\; \frac{\tau}{\sqrt{10/3}} \, \sigma_{y,\mathrm{MH}}(\tau).
```

```julia
htdev(PhaseData(x, τ₀), τs)
```

The `\sqrt{10/3}` factor mirrors the `\sqrt{3}` factor in TDEV: each
arises from the integral of the kernel's sampling function squared
against the frequency-to-time-error map. The second-difference (Allan)
kernel produces `1/3`; the third-difference (Hadamard) kernel produces
`3/10`.

**Use case.** HTDEV inherits TDEV's purpose — characterizing the time
error of a clock or time-distribution system — and adds two
properties TDEV cannot offer:

- **Linear frequency-drift insensitivity.** The third-difference
  kernel annihilates terms linear in `t`, so a record from a Cesium,
  Rubidium, or H-maser clock with significant drift can be analyzed
  without first removing the drift. SP1065 §5 makes the equivalent
  point for HDEV vs ADEV: "the Hadamard deviation may be used to
  reject linear frequency drift when a stability analysis is performed"
  [riley-2008-sp1065](@cite). HTDEV carries that benefit into the
  time-domain.
- **Wider noise-type convergence.** The Hadamard family converges
  over `α ∈ {−4, −3}` (frequency walk-walk and random-run FM) where
  the Allan family diverges [greenhall-1997-third-difference-mvar](@cite). HTDEV is the
  time-domain extension of that range.

For records dominated by white-PM noise without drift, TDEV is fine
and slightly more efficient (tighter CI per τ on shorter records).
For records with drift, divergent low-frequency noise, or both, HTDEV
gives a numerically valid time-stability summary where TDEV would be
contaminated.

!!! info "Original contribution"
    HTDEV is original to SigmaTau. The standard
    time-and-frequency references — SP1065
    [riley-2008-sp1065](@cite), IEEE 1139-2022
    [ieee1139-2022-definitions](@cite), NBS-TN-1337
    [sullivan-1990-tn1337](@cite) — do not define it. The
    authoritative definition lives in the package source itself; the
    public wrapper sits at
    [`src/deviations.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/src/deviations.jl)
    and uses the same `τ`-rescale-of-a-Modified-deviation pattern as
    TDEV (SP1065 §5.2.7 [riley-2008-sp1065](@cite)) applied to the
    modified-Hadamard kernel rather than the modified-Allan kernel.
    The earlier name `ldev` was retained temporarily as a deprecated alias
    and removed before the public release.

## Slope vs noise table

For each estimator, the deviation slope `μ_dev` versus the spectral
exponent `α` of `S_y(f)`:

| α  | ADEV | MDEV | HDEV | MHDEV | TDEV | HTDEV |
|----|------|------|------|-------|------|------|
| +2 (WPM)  | −1     | −3/2 | −1     | −3/2 | −1/2  | −1/2 |
| +1 (FPM)  | −1*    | −1   | −1*    | −1   | 0     | 0    |
|  0 (WFM)  | −1/2   | −1/2 | −1/2   | −1/2 | +1/2  | +1/2 |
| −1 (FFM)  |  0     |  0   |  0     |  0   | +1    | +1   |
| −2 (RWFM) | +1/2   | +1/2 | +1/2   | +1/2 | +3/2  | +3/2 |

(*) ADEV/HDEV's response to flicker PM includes a `log(2π·m)` factor
(see SP1065 §5).

## Demonstration

```@example allan
using SigmaTau, Random
Random.seed!(0)

# Pure WPM noise
N = 8192
x = randn(N)
τs = [1, 4, 16, 64, 256]

a  = adev(PhaseData(x, 1.0), τs).dev
m  = mdev(PhaseData(x, 1.0), τs).dev

# Slopes (deviation log-log)
slope(σ) = (log10(σ[end]) - log10(σ[1])) / (log10(τs[end]) - log10(τs[1]))
round.([slope(a), slope(m)]; sigdigits=3)
```

ADEV slope should be near −1; MDEV slope should be near −3/2 — the
characteristic split that makes MDEV able to disambiguate WPM from FPM.

## Implementation notes

- All cores in `src/kernels.jl` take
  `Vector{Float64}` and return raw arrays; the public API in
  `src/deviations.jl` wraps them and returns `StabilityResult`.
- `MDEV/MHDEV` use a prefix-sum form algebraically equivalent to the
  textbook `1/m⁴` form (Greenhall's third-difference construction
  [greenhall-1997-third-difference-mvar](@cite)); the kernels are
  pinned against direct textbook-form reference implementations in
  `test/stab/legacy_kernels.jl`.
- `tdev` and `htdev` are scaling wrappers; they do no extra kernel work.

## See also

- [Theory: Total family](total_family.md) — boundary-extended estimators.
- [Theory: Confidence](confidence.md) — EDF and CI for these estimators.
- [API: `SigmaTau`](../reference/stab.md) — function signatures.

## References

- SP1065 §5 [riley-2008-sp1065](@cite).
- NBS Technical Note 1337 [sullivan-1990-tn1337](@cite).
- Greenhall, *Third-difference approach to MVAR*, IEEE T-IM 1997
  [greenhall-1997-third-difference-mvar](@cite).
- IEEE 1139-2022 [ieee1139-2022-definitions](@cite) for canonical ADEV / MDEV / HDEV
  / TDEV / MHDEV definitions; HTDEV is **not** defined there.

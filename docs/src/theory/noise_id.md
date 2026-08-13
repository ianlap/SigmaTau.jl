# Theory: Noise Identification

Each averaging-factor `m` of an estimator may correspond to a different
dominant noise type. The EDF formula and the χ² CI both depend on the
α exponent of the local power-law process, so SigmaTau classifies the
noise *per τ* before assigning a confidence interval.

## Why classify?

The Greenhall–Riley EDF formulas (see [Confidence](confidence.md))
parameterize on α explicitly. Pretending the noise is white when it is
flicker yields an overconfident CI; pretending it is random-walk when
it is white yields an underconfident one. Classification matters
quantitatively, not just descriptively.

The five power-law families SigmaTau recognises are white phase
modulation (WHPM, α=2), flicker phase modulation (FLPM, α=1), white
frequency modulation (WHFM, α=0), flicker frequency modulation
(FLFM, α=−1), and random-walk frequency modulation (RWFM, α=−2). The
α exponent feeds straight into the EDF and bias-correction tables.

## Lag-1 ACF on differenced phase

The primary method for identifying noise types in SigmaTau relies on the lag-1 autocorrelation $r_1$ of the differenced phase series. The process involves detrending and differencing the input phase, calculating $r_1$, and mapping that value to the dominant spectral exponent $\alpha$. This methodology is robust at boundaries between adjacent power-law families.

For the exact procedure, required thresholds, and theoretical justification, please read Riley & Greenhall 2004 [riley-2004-lag1-acf](@cite) and NIST SP1065 §6 [riley-2008-sp1065](@cite) for more details.

The $\alpha$-to-$r_1$ mapping is implemented in `src/noise.jl`.

```julia
identify_noise(x, [1, 2, 4, 8]; dmin=0, dmax=2)
```

The first argument is the raw phase vector (a `Vector{Float64}`), not a `PhaseData` record; `m_values` is a vector of averaging factors (`τ = m·τ₀`). The result is a `Vector{Symbol}`, one entry per `m`: each entry is one of `:WHPM`, `:FLPM`, `:WHFM`, `:FLFM`, `:RWFM`.

## B1/R(n) fallback

For short records or borderline cases, SigmaTau automatically falls back to the Allan-variance B1 ratio with the R(n) factor for noise disambiguation. 

For a complete explanation of the B1 ratio and R(n) test mechanics, please read NIST SP1065 §6 [riley-2008-sp1065](@cite) for more details.

## NEFF_RELIABLE = 30

The lag-1 method needs a minimum effective sample count for `r₁` to
be a reliable estimator. SigmaTau uses `NEFF_RELIABLE = 30` as the
threshold; below this, classification falls back to B1/R(n) regardless
of `r₁`. The threshold tracks the Riley R 2020 recommendation
[riley-2020-r-frequency-stability](@cite) and the historical SigmaTau policy mandate.

## NNLS power-law decomposition

`identify_noise` answers "which single α dominates at this τ" one τ at a
time. [`noise_id_nnls`](@ref) answers a different question — "how much of
each α is present across the whole record" — with a single fit to the
entire σ̂²(τ) profile:

```math
\min_{a \ge 0} \; \bigl\lVert W^{1/2}\,(K a - \hat{v}) \bigr\rVert^2 ,
\qquad
K_{ki} = (\tau_k/\tau_0)^{\mu(\alpha_i)} , \quad
\hat{v}_k = \hat\sigma^2(\tau_k) .
```

The nonnegativity constraint is what makes the problem well posed: variance
contributions cannot be negative, and without it the near-collinear
power-law columns trade large positive and negative amplitudes. SigmaTau
solves it *exactly* — every one of the `2ⁿ − 1` supports is tried as an
unconstrained least-squares problem restricted to that support and kept only
if feasible, so the returned `a` is the global optimum rather than an
active-set approximation.

```julia
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))
r = noise_id_nnls(p; deviation = :mhdev)
r.amplitudes    # variance of each α-component at τ = τ₀
r.dominant      # dominant noise type at each τ
```

### Response exponents μ(α)

The variance of a pure α power law scales as `σ²_D(τ) ∝ τ^μ`:

| α | type | μ, modified family | μ, unmodified family |
|--:|:-----|-------------------:|---------------------:|
|  2 | WHPM | −3 | −2 |
|  1 | FLPM | −2 | −2 (with a `ln(2π f_h τ)` factor) |
|  0 | WHFM | −1 | −1 |
| −1 | FLFM |  0 |  0 |
| −2 | RWFM |  1 |  1 |
| −3 | FWFM |  2 |  2 |
| −4 | RRFM |  3 |  3 |

The modified family (MDEV/MHDEV) separates the two PM noises, so `μ = −α−1`
holds across the whole range. The unmodified family (ADEV/HDEV) does not:
white PM and flicker PM both give τ⁻², differing only by a slowly varying
logarithm that no power-law column can resolve. A basis for an unmodified
estimator therefore carries **at most one** of α ∈ {1, 2}; asking for both is
an error.

### Which α a given estimator can carry

| estimator | converges for | default basis |
|:----------|:--------------|:--------------|
| `:mdev` | α > −3 | 2, 1, 0, −1, −2 |
| `:mhdev` | α > −5 | 2, 1, 0, −1, −2, −3, −4 |
| `:hdev` | α > −5 | 2, 0, −1, −2, −3, −4 |
| `:mhtotdev` | α > −5 | 2, 1, 0, −1, −2, −3, −4 |

- **Allan family** (`:mdev`) — second differences, convergent for α > −3, so
  the basis stops at RWFM.
- **Hadamard family** (`:mhdev`, `:hdev`) — third differences, convergent for
  α > −5, so FWFM and RRFM are in range.
- **`:mhtotdev`** — Hadamard-family in form *and* in range: its per-window
  detrend is a **quadratic** least-squares fit, matched to the third
  difference's annihilation degree, so no stochastic quadratic survives into
  the time-reverse extension to leave a seam the third difference could pick
  up. The boundary-extended kernel therefore keeps the full third-difference
  low-frequency response and converges for α > −5, like `:mhdev` — verified on
  pure synthesized RRFM (variance slope 3.00, level stable in N). It carries
  the full modified-family basis. (The pre-quadratic kernel, which used a
  half-mean slope removal, had an effectively second-difference low-frequency
  response and excluded α ≤ −3 here; see
  [MHTOTDEV bias and EDF](mhtotdev_bias_edf.md).)

  Two practical notes: the boundary extension still *bends* the response away
  from τ^μ over the last octaves — a bend, not a different power — so prefer
  the empirical-basis route (below) when those octaves matter; and the shipped
  EDF tables stop at α = −2, so the α ≤ −3 rows borrow the RWFM EDF for
  weighting only (this does not enter the fitted amplitudes).

### Caveats

The basis is analytic and asymptotic. Two consequences worth knowing:

- **Small τ.** Every estimator is still converging toward its power law at
  m = 1–2 (the modified family sits 2–3× above the asymptote there, settling
  by m ≈ 8), and the EDF weighting trusts exactly those points most. The
  default grid therefore starts at m = 4; including m = 1 lets a pure WFM
  record buy off its excess with a spurious WPM column.
- **Adjacent columns.** Neighbouring power laws are near-collinear over a
  finite τ range, so individual amplitudes of adjacent α are correlated and
  noisy from realization to realization. `dominant` (and the fitted total)
  are the robust products; the per-α split is not a precision measurement.
  The Monte Carlo study in `tools/nnls_noise_id.jl` quantifies this with an
  *empirical* basis — ensemble-mean σ̂²(τ) of synthesized unit noises — which
  removes the asymptotic-basis error and isolates the sampling error.

## Synthesizing α = −4 (RRFM)

[`noise_gen`](@ref) reaches random-run FM through the `h` path only:

```julia
noise_gen(PhaseData, 8192, 1.0; h = Dict(-4 => 1e-40))
```

`sigma1[-4]` raises an `ArgumentError`: the Allan-variance σ_y(τ₀)
convention has no finite value for α ≤ −3, so there is nothing to target.

The component is built by double discrete integration of white noise rather
than by the f^(α/2) spectral shaper, which zeroes the DC bin of a periodic
record and so caps the lowest Fourier component exactly where an f⁻⁴ law puts
almost all of its power. With `w` white of per-sample variance `q` at τ₀ and
`y = cumsum(cumsum(w))`, the cascade of two discrete integrators has

```math
|H_2(f)|^2 = \frac{1}{16\sin^4(\pi f \tau_0)} ,
```

and white `w` has the one-sided PSD `S_w(f) = 2 q τ₀`, so

```math
S_y(f) = \frac{q\,\tau_0}{8\sin^4(\pi f \tau_0)}
\;\xrightarrow[f \to 0]{}\;
\frac{q}{8\pi^4\tau_0^3}\, f^{-4} .
```

Matching `S_y(f) = h₋₄ f⁻⁴` fixes the driver variance at

```math
q = 8\pi^4 \tau_0^3\, h_{-4} ,
```

which is what the generator uses. Because `Δ²y = √q · w` exactly, the
normalization is verifiable to sampling precision from the synthesized series
itself.

## The `noise_type` field

`adev`, `mdev`, etc. populate `StabilityResult.noise_type` with one
classification per τ. Empty when classification was not requested or
not possible. When a particular m fails (e.g. degenerate variance,
too few points after differencing), the kernel carries forward the
last reliable classification rather than emitting `:unknown` mid-run.

## Implementation notes

- 5σ outlier rejection plus linear detrend run before classification.
- `_noise_id_lag1acf` and `_noise_id_b1rn` are the two internal kernels;
  the public `identify_noise` chooses between them based on `N_eff`.
- Per-m quadratic detrending is on by default to remove frequency
  offset and drift from each decimated subseries; pass `detrend=false`
  to match Stable32's no-detrend convention when comparing fixtures
  point-for-point.
- **Scale-invariance.** The lag-1 ACF degeneracy guard uses a
  relative threshold (`ssx ≤ eps · ‖x‖²`) so classification is
  invariant under positive linear rescaling — phase records in
  seconds (typically ~1e-9..1e-12) classify identically to the same
  shape rescaled into any other units. This guard exists only to
  prevent 0/0 on truly constant input; meaningful tiny-magnitude
  records pass through normally.

## See also

- [Theory: Confidence](confidence.md).
- [API: `SigmaTau`](../reference/stab.md).

## References

- Riley & Greenhall, *Power-law noise identification using the lag-1
  autocorrelation*, PTTI 2004 [riley-2004-lag1-acf](@cite).
- SP1065 §6 [riley-2008-sp1065](@cite).
- Riley, *Frequency Stability Analysis Using R*, 2020
  [riley-2020-r-frequency-stability](@cite).

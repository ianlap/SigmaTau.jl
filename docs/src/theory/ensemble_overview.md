# Theory: Clock State-Space Models

`SigmaTau.Est` models a precision oscillator as a linear stochastic
differential equation (SDE) whose state vector promotes the polynomial
coefficients of the deterministic clock model — phase, frequency, and
(optionally) frequency aging — to states driven by independent Wiener
processes. Discretising the SDE on a uniform sample interval `τ` yields
the propagator `Φ` and the innovation covariance `Q` that drive the
discrete-time [Kalman filter](kalman.md).

This page collects the time-domain models in one place: the
deterministic polynomial skeleton, the two- and three-state SDEs that
generalise it, and the closed-form discretisation `(Φ, Q)` consumed by
the Kalman recursion.

## Polynomial clock model

The polynomial clock model decomposes phase deviation into a low-order
Taylor expansion plus a stochastic residual:

```math
x(t) \;=\; x_0 + y_0 \, t + \tfrac{1}{2}\, D\, t^{2} + \varepsilon(t),
```

with `x₀` the initial offset, `y₀` the initial fractional-frequency
offset, `D` the linear frequency drift, and `ε(t)` the stochastic
noise [Banerjee & Matsakis 2023](@cite banerjee-2023-timekeeping). The
same decomposition appears in IEEE-1139 as the optimal time-prediction
model `x̂(t_p) = x(t₀) + y(t₀) τ_p + ½ D τ²_p`
[IEEE 1139-2022](@cite ieee1139-2022-definitions), and in the HP
Application Note 1289 frequency-stability primer
[Allan, Ashby & Hodge 1997](@cite allan-1997-hp-app1289).

Banerjee and Matsakis warn that the polynomial decomposition is best
read as a local Taylor expansion: real oscillator drifts can be
logarithmic, so the polynomial form is an approximation around the
present epoch [Banerjee & Matsakis 2023](@cite banerjee-2023-timekeeping).
The `½ D t²` phase contribution from frequency drift produces ADEV
and MDEV slopes `σ_y(τ) = D τ / √2` and a TIE of roughly
`1.2 σ_x(τ)` [Allan, Ashby & Hodge 1997](@cite allan-1997-hp-app1289).

The polynomial coefficients become the states of the canonical
three-state SDE below; flicker noises (any non-integer fractional
integration order) cannot be represented as a finite-state SDE driven
by Wiener processes, so the SDE family truncates to the integer-FM
ladder and excludes flicker FM by construction
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan).

## Two-state clock SDE

The two-state SDE models phase $X_1$ and frequency $X_2$ as a pair
of integrated Wiener processes:

```math
\begin{aligned}
dX_1 &= (X_2 + \mu_1)\,dt + \sigma_1\,dW_1, \\
dX_2 &= \mu_2\,dt + \sigma_2\,dW_2.
\end{aligned}
```

`σ₁` drives white-FM phase noise and `σ₂` drives random-walk-FM
frequency noise [Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan).
This is the special case of the general three-state SDE with
`μ₃ = σ₃ = 0`, and is the canonical model when the dominant noises
are white FM and random-walk FM
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan).

The closed-form Allan variance for the stationary two-state clock is

```math
\sigma_y^{2}(\tau) \;=\; \frac{\sigma_1^{2}}{\tau} \;+\; \frac{\sigma_2^{2}\,\tau}{3} \;+\; \frac{\tau^{2}}{2}\,c_3^{2},
```

with the `τ⁻¹` branch from white FM and the `τ¹` branch from
random-walk FM
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan). The
diffusion coefficients are recovered from observed slopes: `σ₁` from
the `τ⁻¹/²` ADEV branch at small `τ` and `σ₂` from the `τ⁺¹/²`
branch at large `τ` [Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan).

In the time-and-frequency community the metrological fractional
frequency `y(t)` is identified with `dX₁/dt` rather than with the
state `X₂` — a distinction that is often misused in the literature
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan). The
phase-noise of a real clock is the sum of direct phase shocks and the
integral of frequency shocks, so the two-state model already requires
both to be modelled separately rather than collapsed into a scalar
process [Stein 2003](@cite stein-2003-timescales). For cesium-class
clocks the two-state model is empirically supported by Allan-deviation
slopes of `τ⁻¹/²` for `τ < 10⁶ s` and `τ⁺¹/²` for `τ > 10⁷ s`
[Stein 2003](@cite stein-2003-timescales). Tryon and Jones formulated
the original cesium-clock state-space estimator with white-noise
increments scaled by `√δ` to keep variances proportional to interval
length [Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).
USNO maser-ensemble Kalman implementations continue to use the same
two-state form
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

In SigmaTau the two-state SDE corresponds to [`TwoStateClock`](@ref).

## Three-state clock SDE

The three-state SDE adds a frequency-aging (drift) state $X_3$
driven by its own Wiener process, with `X₂` now coupled to `X₃`
through an integrating term:

```math
\begin{aligned}
dX_1 &= (X_2 + \mu_1)\,dt + \sigma_1\,dW_1, \\
dX_2 &= (X_3 + \mu_2)\,dt + \sigma_2\,dW_2, \\
dX_3 &= \mu_3\,dt + \sigma_3\,dW_3.
\end{aligned}
```

This is the canonical Wiener-driven model for cesium, hydrogen-maser,
and rubidium clocks and underlies the discrete Kalman filter used in
clock-ensemble timekeeping
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan),
[Stein 2003](@cite stein-2003-timescales),
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).
With `σ₃ > 0` the integrated frequency increments cease to be
stationary, so the resulting Allan variance acquires explicit
dependence on the absolute epoch `t_k`:

```math
\sigma_y^{2}(t_k, \tau) \;=\;
\frac{\sigma_1^{2}}{\tau}
+ \frac{\sigma_2^{2}\,\tau}{3}
+ \frac{\sigma_3^{2}\,\tau^{3}}{20}
+ \sigma_3^{2}\!\left(\frac{\tau^{3}}{3} + \frac{\tau^{2}\, t_k}{2}\right)
+ \frac{\tau^{2}}{2}\!\bigl[c_3 + \mu_3(\tau + t_k)\bigr]^{2}.
```

The Hadamard variance, by contrast, depends on `τ` alone and remains
finite under random-run drift:

```math
H\sigma_y^{2}(\tau) \;=\; \frac{\sigma_1^{2}}{\tau} + \frac{\sigma_2^{2}\,\tau}{6} + \frac{11}{120}\,\sigma_3^{2}\,\tau^{3} + \frac{\mu_3^{2}\,\tau^{4}}{6},
```

which is one operational reason to prefer the
[Hadamard family](allan_family.md) for records in which `σ₃` is
non-negligible
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan).

Tryon and Jones distinguish three frequency-drift hypotheses for the
state $X_3$ — drift-free, deterministic constant drift, and
random-walk drift — and report the first valid statistical test
between them on a 333-day NBS dataset of seven cesium clocks: the
constant-drift model improves `−2 ln L` by 41 over the no-drift model
(highly significant under `χ²(6)`), while adding random-walk drift on
top yields essentially zero further improvement
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).
Frequency drift is statistically significant only for end-of-life or
freshly commissioned cesium clocks; middle-life clocks showed none
over a one-year window
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).

In SigmaTau the three-state SDE corresponds to
[`ThreeStateClock`](@ref).

## State-transition matrix Φ

The state-transition matrix `Φ` is the closed-form discrete-time
propagator that advances the clock state through one sample interval
`τ`. Without process noise the SDE integrates exactly to a Taylor
expansion in `τ`; the three-state propagator is upper-triangular with
polynomial-in-`τ` entries:

```math
\Phi \;=\; \begin{bmatrix}
1 & \tau & \tau^{2}/2 \\
0 & 1     & \tau \\
0 & 0     & 1
\end{bmatrix},
```

applied to the state $[x, y, \omega]^\top$
[Stein 2003](@cite stein-2003-timescales). Frequency aging contributes
the `τ²/2` off-diagonal entry to the phase row because aging
integrates twice into phase
[Stein 2003](@cite stein-2003-timescales). The discrete iterative
solution of the linear SDE is

```math
X(t_{k+1}) \;=\; \Phi_\tau\, X(t_k) + B_\tau M + J_k, \qquad J_k \sim \mathcal{N}(0, Q),
```

with the Gaussian innovation `J_k` carrying all the stochastic
content [Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan).
The discrete-time recursion is exact for the linear Wiener-driven SDE:
there is no discretisation error in the Wiener integration
[Zucca & Tavella 2015](@cite zucca-2015-clock-error-jumps).

For an `m`-clock ensemble of `n`-th-order homogeneous clocks the
global propagator factors as `F = A ⊗ I_m`, where `A` is the per-clock
polynomial transition matrix
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale). Stacking
per-clock states block-diagonally is the standard NBS Kalman
convention
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation),
and is also how Breakiron's USNO maser ensemble assembles its
two-state per-clock blocks
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

In SigmaTau, `Φ` is built by the
[`state_transition`](@ref) hook on each `AbstractClockModel`.

## Process-noise covariance Q

The discrete-time process-noise covariance `Q` is the covariance of
the Gaussian innovation `J_k` driving the clock state forward by one
sample interval `τ`. It is obtained by integrating
`Φ_t · diag(σ_i²) · Φ_tᵀ` over `[0, τ]`. For the canonical three-state
SDE the entries are exact polynomials in `τ` that depend only on the
diffusion coefficients `σ₁, σ₂, σ₃` — independent of the absolute
epoch `t_k`
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan):

```math
Q \;=\;
\begin{pmatrix}
\sigma_1^{2}\tau + \tfrac{1}{3}\sigma_2^{2}\tau^{3} + \tfrac{1}{20}\sigma_3^{2}\tau^{5}
& \tfrac{1}{2}\sigma_2^{2}\tau^{2} + \tfrac{1}{8}\sigma_3^{2}\tau^{4}
& \tfrac{1}{6}\sigma_3^{2}\tau^{3} \\
\tfrac{1}{2}\sigma_2^{2}\tau^{2} + \tfrac{1}{8}\sigma_3^{2}\tau^{4}
& \sigma_2^{2}\tau + \tfrac{1}{3}\sigma_3^{2}\tau^{3}
& \tfrac{1}{2}\sigma_3^{2}\tau^{2} \\
\tfrac{1}{6}\sigma_3^{2}\tau^{3}
& \tfrac{1}{2}\sigma_3^{2}\tau^{2}
& \sigma_3^{2}\tau
\end{pmatrix}.
```

The `t_k` dependence visible in the three-state Allan variance
therefore lives in the joint `X₂, X₃` propagation, not in the
innovation covariance itself
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan). The
two-state restriction (`σ₃ = 0`) collapses to the standard Breakiron
form

```math
Q(\Delta t) \;=\;
\begin{pmatrix}
s_\eta^{2}\,\Delta t + \tfrac{1}{3}s_\nu^{2}\,(\Delta t)^{3} & \tfrac{1}{2}s_\nu^{2}\,(\Delta t)^{2} \\
\tfrac{1}{2}s_\nu^{2}\,(\Delta t)^{2} & s_\nu^{2}\,\Delta t
\end{pmatrix},
```

with random-walk-FM density `s_η²` and random-run-FM density `s_ν²`
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

Tryon and Jones use the small-`δ` diagonal approximation
`M(t) = δ(t) Q` with diagonal entries
`(σ_ε², σ_η², σ_α²)` — the off-diagonal cross terms from the full
integration are dropped, and the two are not equivalent at finite `δ`
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).
Breakiron calibrates the spectral densities `s_ε², s_η², s_ν²` by
least-squares fitting the Hadamard variance, which is insensitive to
deterministic drift
[Breakiron 2001](@cite breakiron-2001-kalman-timescales). The
practical recipe is unchanged: identify the dominant noise from the
ADEV slope, read off `σ_i` at the dominant `τ`, then assemble `Q`
[Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan).

When `Φ`, `Q`, and `R` are block-diagonal in a multi-clock ensemble,
`P` is also block-diagonal, so the Kalman recursion factors into
per-clock 2×2 inversions and cost is linear in clock count
[Breakiron 2001](@cite breakiron-2001-kalman-timescales). When `Q` is
unknown a priori it can be calibrated from observed Allan- or
Hadamard-deviation slopes.

In SigmaTau, `Q` is built by the [`process_noise`](@ref) hook on
each `AbstractClockModel`.

## See also

- [Theory: Kalman Filter](kalman.md) — discrete-time estimator that
  consumes `Φ` and `Q`.
- [Theory: Clock Ensembles](clock_ensembles.md) — joint stacked-state
  Kalman model built on top of the per-clock SDE.
- [Theory: Clock Steering with PID Controllers](steering.md) —
  closes the loop between the estimator and the oscillator's
  frequency control input.
- [API: `SigmaTau.Est`](../reference/est.md) — exported types and
  functions.

## References

- [IEEE 1139-2022](@cite ieee1139-2022-definitions) — IEEE Standard
  Definitions of Physical Quantities for Fundamental Frequency and
  Time Metrology.
- [Banerjee & Matsakis 2023](@cite banerjee-2023-timekeeping) —
  *An Introduction to Modern Timekeeping and Time Transfer*.
- [Allan, Ashby & Hodge 1997](@cite allan-1997-hp-app1289) —
  HP Application Note 1289.
- [Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation) —
  state-space MLE for NBS cesium clocks.
- [Stein 2003](@cite stein-2003-timescales) — time-scale-algorithm
  unification.
- [Breakiron 2001](@cite breakiron-2001-kalman-timescales) — USNO
  maser-ensemble Kalman filter.
- [Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan) — the
  canonical three-state clock SDE and the integrated `Q` matrix.
- [Yan et al. 2023](@cite yan-2023-structured-kf-timescale) —
  homogeneous-ensemble factorisation `F = A ⊗ I_m`.

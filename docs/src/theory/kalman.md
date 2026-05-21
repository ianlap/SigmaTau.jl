# Theory: Kalman Filter

The Kalman filter is the recursive linear-Gaussian state estimator
used throughout atomic-clock timekeeping to fuse a clock's state-space
dynamics — phase, frequency, and (optionally) drift — with noisy
phase-difference measurements. For the discrete-time SDE clock model
introduced in [Clock state-space models](ensemble_overview.md),
discretised via the propagator `Φ` and innovation covariance `Q`, the
Kalman recursion is the optimal linear estimator under Gaussian noise
and the maximum-likelihood estimator under independent Gaussian
innovations
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation),
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

This page collects the standard predict / update recursion, the
innovation and Kalman gain that drive it, the unconditional covariance
propagator that `SigmaTau.Est` adds alongside `predict!` for
holdover-band projection, and a pedagogical reading list. The
implementation lives at [`KalmanFilter`](@ref) in `src/est/filters.jl`.

## Standard Kalman filter

Each Kalman step has two phases. The prediction step propagates the
previous posterior estimate `x̂` and covariance `P` through the
discrete transition `Φ` and adds `Q`. The update step corrects the
prediction with a measurement-residual term weighted by the Kalman
gain `K`
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation):

```math
\begin{aligned}
\hat x_k^- &= \Phi\, \hat x_{k-1}^+, \\
P_k^-      &= \Phi\, P_{k-1}^+ \Phi^\top + Q, \\
C_k        &= H\, P_k^- H^\top + R, \\
K_k        &= P_k^- H^\top C_k^{-1}, \\
\hat x_k^+ &= \hat x_k^- + K_k (z_k - H\, \hat x_k^-), \\
P_k^+      &= (I - K_k H)\, P_k^-.
\end{aligned}
```

Innovations are Gaussian white noise iff the model and parameters are
correct, which makes them simultaneously the likelihood ingredient for
parameter estimation and the structural diagnostic for model adequacy
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).

A Kalman filter that includes phase as a state is **not detectable**
when only phase differences are observed, so phase-covariance entries
grow without bound
[Breakiron 2001](@cite breakiron-2001-kalman-timescales). The
atomic-clock-ensemble system inherits this undetectability; the
conventional Kalman covariance therefore grows on the non-observable
common-mode component while observable linear combinations of states
remain tight.

In `SigmaTau` the standard recursion is implemented by
[`KalmanFilter`](@ref) with [`predict!`](@ref) and
[`update!`](@ref). The filter is mutable but uses out-of-place
`StaticArrays` math for AD-cleanness.

## Innovation

The innovation is the measurement residual computed at each update —
the difference between the actual measurement and its model-predicted
value `H Φ x̂`:

```math
\nu_k \;=\; z_k - H\, \hat x_k^-.
```

It serves jointly as the update term and as a model-fit diagnostic
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation),
[Stein 2003](@cite stein-2003-timescales). When the model and its
parameters are correct, the innovation sequence is Gaussian white
noise with covariance `C = H P Hᵀ + R`, which makes it the natural
ingredient for both maximum-likelihood parameter estimation and
goodness-of-fit testing
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).

In a multi-clock ensemble the innovation drives the random-shock
estimate for each pair of clocks:

```math
\hat s_{ij}(t_k) \;=\; K_{ij}(t_{k+1})\,\bigl[\,z_{ij}(t_{k+1}) - H\Phi\,\hat x_{ij}(t_k)\,\bigr]
```

[Stein 2003](@cite stein-2003-timescales).

## Kalman gain

The Kalman gain is the matrix that weights the innovation in the
update step, balancing trust between model prediction and new
measurement:

```math
K \;=\; P H^\top C^{-1}, \qquad C \;=\; H P H^\top + R.
```

For a clock-difference 3-state system the gain is a 3-element column
vector `K_{ij} = [k^x, k^y, k^ω]^\top` that distributes the innovation
across phase, frequency, and frequency-aging shock estimates
[Stein 2003](@cite stein-2003-timescales). Inverse-variance weighting
of clocks within a Kalman ensemble mean systematically underestimates
each clock's variance because each clock contributes to the ensemble it
is compared to; the bias is `σ_u² = σ² / (1 − w)` and motivates
referencing each clock to a mean of the rest
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

## Unconditional covariance propagation

`SigmaTau.Est` exposes a [`prop!`](@ref) method alongside `predict!`
that advances `x ← Φ(dt) x` and `P ← Φ(dt) P Φ(dt)' + Q(dt)`
unconditionally — without the `est.k > 0` gate that `predict!`
inherits from the legacy MATLAB filter, and without bumping the step
counter. The intended use is producing a 1-σ covariance band around a
deterministic forward projection (the shaded holdover-bound pattern in
[`tutorials/05_holdover_comparison.md`](../tutorials/05_holdover_comparison.md))
without disturbing the live filter's update sequencing.

The dt-aware [`state_transition`](@ref) / [`process_noise`](@ref)
overloads on `AbstractClockModel` are what make this practical: `prop!`
calls them at the caller-supplied horizon `h·τ₀` instead of being
locked to `model.tau`, so a single converged filter can re-seed
side-channel projections at arbitrary horizons. `Φ` has the group
property and `Q` is additive under it
(`Q(dt₁ + dt₂) = Φ(dt₂) Q(dt₁) Φ(dt₂)' + Q(dt₂)`), so two `prop!`s of
`dt₁` then `dt₂` produce the same state and covariance as one `prop!`
of `dt₁ + dt₂` — a property locked at `rtol = 1e-14` in
[`test/est/runtests.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/test/est/runtests.jl).

## See also

- [Theory: Clock State-Space Models](ensemble_overview.md) — supplies
  `Φ`, `Q`, and the underlying SDE the Kalman filter inverts.
- [Theory: Clock Steering with PID Controllers](steering.md) — closes
  the loop with a controller acting on the Kalman estimates.
- [API: `SigmaTau.Est`](../reference/est.md) —
  [`KalmanFilter`](@ref), [`predict!`](@ref),
  [`update!`](@ref), [`prop!`](@ref).

## References

- [Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation) —
  state-space MLE for NBS cesium clocks.
- [Breakiron 2001](@cite breakiron-2001-kalman-timescales) — USNO
  maser-ensemble two-state Kalman filter.
- [Stein 2003](@cite stein-2003-timescales) — innovation-driven
  random-shock estimates in clock ensembles.
- [Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan) —
  three-state SDE underlying the Kalman model.
- [Chaudhari 2022](@cite chaudhari-2022-uva-kalman-chapter) —
  pedagogical chapter on the Kalman filter and its variants
  (EKF / UKF / particle filter); useful as a from-first-principles
  refresher of the recursion this module implements.

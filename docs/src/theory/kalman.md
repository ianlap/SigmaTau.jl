# Theory: Kalman Filter and Variants

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
innovation and Kalman gain that drive it, three implementation
variants (U-D factorisation, ALS noise tuning, online adaptive tuning),
the structured-Kalman variant that handles the undetectability of a
clock-only ensemble, and Wu's LTI-equivalence framework for
characterising the filter's steady-state performance.

## Standard Kalman filter

Each Kalman step has two phases. The prediction step propagates the
previous posterior estimate `x̂` and covariance `P` through the
discrete transition `Φ` and adds `Q`. The update step corrects the
prediction with a measurement-residual term weighted by the Kalman
gain `K`
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation),
[Kubczak et al. 2019](@cite kubczak-2019-fast-sync-rubidium):

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

The recursion is the dynamic-system generalisation of the
Gauss–Plackett recursive least squares, repurposed by Tryon and Jones
so that model parameters anywhere in the state-space form can be
estimated by maximum likelihood
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).
Innovations are Gaussian white noise iff the model and parameters are
correct, which makes them simultaneously the likelihood ingredient for
parameter estimation and the structural diagnostic for model adequacy
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).

The discrete algebraic Riccati equation gives the steady-state error
covariance:

```math
P_p \;=\; A P_p A^\top + G Q_w G^\top - (A P_p C^\top + G S_{wv})(C P_p C^\top + R_v)^{-1}(C P_p A^\top + S_{wv}^\top G^\top),
```

with optional cross-covariance `S_{wv}` between process and
measurement noise
[Åkesson et al. 2008](@cite akesson-2008-generalized-als).

A Kalman filter that includes phase as a state is **not detectable**
when only phase differences are observed, so phase-covariance entries
grow without bound
[Breakiron 2001](@cite breakiron-2001-kalman-timescales). The
atomic-clock-ensemble system inherits this undetectability; the
conventional Kalman covariance therefore grows without bound and is
not strictly optimal
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale). Reducing
to a two-state (frequency, drift) filter avoids the unbounded phase
covariance and is reported to outperform the three-state Tryon-style
filter on simulated cesium data
[Breakiron 2001](@cite breakiron-2001-kalman-timescales). When `Φ`,
`Q`, and `R` are block-diagonal in a multi-clock ensemble, `P` is also
block-diagonal, so the recursion factors into per-clock 2×2 inversions
and cost is linear in clock count
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

Practical Kalman filtering reduces to a covariance-estimation problem
because `Q` and `R` are typically unknown, and ad-hoc manual tuning is
unreliable
[Åkesson et al. 2008](@cite akesson-2008-generalized-als). In the
rubidium-disciplining context the prediction-step covariance is often
parameterised as `P_pred = A P Aᵀ + e Q` with a single scalar tuning
factor `e` that interpolates between trust-prediction (low `e`) and
trust-measurement (high `e`)
[Kubczak et al. 2019](@cite kubczak-2019-fast-sync-rubidium).

In SigmaTau the standard recursion is implemented by
[`StandardKalmanFilter`](@ref) with [`predict!`](@ref) and
[`update!`](@ref).

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

[Stein 2003](@cite stein-2003-timescales). The Tryon–Jones MLE
strategy minimises `−2 ln L` over the residuals from a single Kalman
pass:

```math
L \;=\; \sum_t \bigl[\ln\lvert C(t)\rvert + I^\top(t)\, C^{-1}(t)\, I(t)\bigr],
```

an outer-optimiser-driven procedure that estimates parameters anywhere
in the state-space form
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation). A
Cholesky factorisation of an augmented matrix containing `C`, `H P`,
and `I` avoids explicit inversion of the innovation covariance
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).

The innovation autocovariance carries identifiable information about
`Q`, `R`, and any cross-correlation `S_{wv}`; the Autocovariance
Least-Squares (ALS) machinery below solves for these as a linear
least-squares problem
[Åkesson et al. 2008](@cite akesson-2008-generalized-als). For exact
innovation generation in jump-augmented simulation, Cholesky-factor
`Q = A Aᵀ` and draw `J_k = A Z` with `Z ∼ N(0, I)`
[Zucca & Tavella 2015](@cite zucca-2015-clock-error-jumps).

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
[Stein 2003](@cite stein-2003-timescales). For a hydrogen maser with
diffusion coefficients `(σ₀², σ₁², σ₂²) = (1.3 × 10⁻²², 3.7 × 10⁻²⁶,
7.0 × 10⁻³⁴)` the steady-state Kalman gains are
`(K_{s1}, K_{s2}) = (0.1294, 1.6771 × 10⁻⁵ s⁻¹)`
[Wu et al. 2023](@cite wu-2023-kf-performance-lti).

The steady-state gain is determined by `Q` and `R` via the discrete
algebraic Riccati equation; observations enter only as inputs to the
LTI structure
[Wu et al. 2023](@cite wu-2023-kf-performance-lti),
[Åkesson et al. 2008](@cite akesson-2008-generalized-als).
Inverse-variance weighting of clocks within a Kalman ensemble mean
systematically underestimates each clock's variance because each
clock contributes to the ensemble it is compared to; the bias is
`σ_u² = σ² / (1 − w)`
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

## U-D factorized Kalman filter

The U-D factorised Kalman filter (the "Bierman–Thornton" form)
propagates and updates the covariance factors `U` (unit
upper-triangular) and `D` (diagonal) under the decomposition
`P = U D Uᵀ`, instead of working with `P` directly
[Ramos et al. 2022](@cite ramos-2022-ud-kalman-filter).

The standard Kalman covariance propagation
`P⁻ = F_{k−1} P_{k−1}⁺ F_{k−1}ᵀ + G_{k−1} Q_{k−1} G_{k−1}ᵀ` is rewritten
in factored form as

```math
P^{-} \;=\; \bar U\,\bar D\,\bar U^{\top}
       \;=\; [\,F U^{+} \;\; G\,]
              \begin{bmatrix} D^{+} & 0 \\ 0 & Q \end{bmatrix}
              [\,F U^{+} \;\; G\,]^{\top}
       \;=\; W \breve D W^{\top}.
```

The Weighted Modified Gram–Schmidt (WMGS) recursion extracts `Ū` and
`D̄`:

```math
v_n = w_n, \qquad
v_k = w_k - \sum_{j=k+1}^{n} u(k,j)\, v_j, \qquad
u(k,j) = \frac{w_k\,\breve D\, v_j^{\top}}{v_j\,\breve D\, v_j^{\top}}
```

[Ramos et al. 2022](@cite ramos-2022-ud-kalman-filter). The
measurement update in Bierman/Carlson rank-one form is

```math
P^{+} = \bar U\!\bigl[\bar D - \bar D \bar w\, a_i\, \bar w^{\top}\bar D\bigr]\bar U^{\top},\quad
\bar w = \bar U^{\top} H_i^{\top},\quad
a_i = (H_i \bar U \bar D \bar U^{\top} H_i^{\top} + r_i)^{-1},
```

with the bracketed term symmetric positive-definite so that a `U D Uᵀ`
decomposition exists
[Ramos et al. 2022](@cite ramos-2022-ud-kalman-filter). The modified
Agee–Turner / Carlson recursion implements the rank-one update with
stable pivoting; the unmodified Agee–Turner is unstable when the
rank-one update enters with a negative sign, so the Carlson form is
preferred [Ramos et al. 2022](@cite ramos-2022-ud-kalman-filter).

For correlated measurement noise `R_c`, the filter pre-decorrelates
via `R_c = U_R D_R U_Rᵀ` and `z = U_R⁻¹ y`; sequential scalar
measurement updates assume uncorrelated measurements
[Ramos et al. 2022](@cite ramos-2022-ud-kalman-filter).

The U-D form avoids forming or storing `P` directly, propagating only
the `U` and `D` factors. All square-root operations live inside `D`,
so the propagation and update steps perform no square roots, and
positive semi-definiteness can be enforced by sign-checking the
entries of `D`
[Ramos et al. 2022](@cite ramos-2022-ud-kalman-filter).

!!! info "What U-D factorisation actually buys you"
    Compared with Potter's square-root filter `P = S Sᵀ`, the U-D
    factorisation does **not** necessarily improve numerical
    precision — efficiency and symmetry preservation are the actual
    wins. The Ramos summary explicitly contradicts the casual claim
    that UD is more numerically precise; treatments that inherit a
    "UD = better numerics" framing should be revised
    [Ramos et al. 2022](@cite ramos-2022-ud-kalman-filter).

!!! note "Planned implementation"
    [`UDFactorizedFilter`](@ref) is currently an empty struct stub in
    SigmaTauEnsemble. The mathematical recipe above is documented;
    no Bierman/Thornton predict / update loop is wired up yet.

## Generalized ALS noise tuning

Generalised Autocovariance Least-Squares (ALS) noise tuning is an
**offline** procedure that estimates a Kalman filter's process-noise
covariance `Q_w`, measurement-noise covariance `R_v`, and (optionally)
cross term `S_{wv}` from the autocovariance of innovations recorded
under an existing suboptimal filter
[Åkesson et al. 2008](@cite akesson-2008-generalized-als).

The Åkesson 2008 generalisation extends Odelson's uncorrelated-noise
ALS to systems with correlated process and measurement noise, casts
the estimation problem as a symmetric semidefinite least-squares
problem with an L-curve-tuned regularisation toward an initial guess,
and solves it with a Mehrotra predictor–corrector interior-point
method
[Åkesson et al. 2008](@cite akesson-2008-generalized-als):

```math
\min_X \;\tfrac{1}{2}\,\|A_{LS}\,\mathrm{vec}(X) - \mathrm{vec}(\hat{\mathcal{R}}_e(L))\|^{2} + \tfrac{1}{2}\,\lambda\,\|\mathrm{vec}(X - X_0)\|^{2}
\quad \text{s.t.} \; X \succeq 0.
```

Discretising a continuous-time clock SDE with white noise generally
produces correlated process and measurement noise in discrete time,
motivating the cross term `S_{wv}` that this generalisation handles
[Åkesson et al. 2008](@cite akesson-2008-generalized-als).
Regularisation toward an initial guess `X₀` via the parameter `λ`
trades bias for variance in the estimated covariance; `λ` is selected
via the L-curve method
[Åkesson et al. 2008](@cite akesson-2008-generalized-als). When the
noise-shaping matrix `Φ_w` is unknown, only the product
`Φ_w Q_w Φ_wᵀ` is identifiable; a trace-surrogate rank minimisation
plus SVD recovers a minimal disturbance representation
[Åkesson et al. 2008](@cite akesson-2008-generalized-als). The
interior-point solver uses a Mehrotra predictor–corrector with
centring parameter `σ = (μ̂_aff / μ̂)³` and step lengths preserving
positive semi-definiteness via generalised eigenvalue problems
[Åkesson et al. 2008](@cite akesson-2008-generalized-als).

The vectorised steady-state Lyapunov form linear in the noise vector
`q` and measurement variance `R` is

```math
P_s \;=\; (I - A \otimes A)^{-1}\bigl[(G \otimes G)\,M\,q + (A K \otimes A K)\,R\bigr],
```

solved jointly for the noise parameters from the autocorrelation
matrix of the innovation sequence as

```math
[q\;\; R]^{\top} \;=\; (A_{LS}^{\top} A_{LS})^{-1} A_{LS}^{\top}\, R_1(N)_s
```

[Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium).

On the Van der Vusse CSTR benchmark over 200 Monte-Carlo realisations,
ALS-tuned filtering achieves 17.7 % performance improvement vs 25.0 %
for the omniscient ideal Kalman — about 71 % of the ideal gap closed
[Åkesson et al. 2008](@cite akesson-2008-generalized-als). For a
PRS10 rubidium GPSDO, ALS converges in 50 iterations from a
deliberately bad initial guess `(q₁, q₂, q₃, R) = (1, 1, 1, 5)` to
physically meaningful values `q₁ = 4.70 × 10⁻¹⁶`,
`q₂ = 1.23 × 10⁻¹⁸`, `q₃ = 1.68 × 10⁻²⁰`, `R = 1.86 × 10⁻¹⁴`
[Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium).

## Adaptive Kalman filter

An adaptive Kalman filter retunes its noise covariances `Q` and `R`
**online** from the observed innovation sequence rather than relying
on a fixed a priori specification, addressing the fact that practical
Kalman filtering reduces to a covariance-estimation problem because
`Q` and `R` are typically unknown
[Åkesson et al. 2008](@cite akesson-2008-generalized-als),
[Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium).

In the rubidium-disciplining context Liu et al. couple the standard
predict / update Kalman recursion with an iterative ALS loop that
updates `(q₁, q₂, q₃, R)` from the running innovation autocorrelation,
then refreshes the `Q` matrix without requiring multiple Hadamard
differencing of the clock-difference data
[Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium):

```math
P_s \;=\; (I - A \otimes A)^{-1}\bigl[(G \otimes G)\,M\,q + (A K \otimes A K)\,R\bigr],
\qquad
[q\;\; R]^{\top} \;=\; (A_{LS}^{\top} A_{LS})^{-1} A_{LS}^{\top}\, R_1(N)_s.
```

Initial state-vector entries are constructed from three consecutive
clock-difference samples by symmetric finite differencing of phase to
recover frequency and drift estimates
[Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium). The resulting
GPSDO loop disciplines an SRS PRS10 rubidium against a GPS 1PPS
reference and reaches long-term stability comparable to a cesium
clock: after full discipline the PRS10 reaches mean frequency
deviation `1.85 × 10⁻¹⁵` and daily drift `4.39 × 10⁻¹⁵`
[Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium). Disciplined
ADEV reaches `4.12 × 10⁻¹³` at 86 400 s and `3.06 × 10⁻¹³` at
100 000 s; medium-term ADEV (1 000–20 000 s) degrades versus
free-running due to GPS-receiver jitter
[Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium). Locked
clock-difference statistics: mean 0.068 ns, standard deviation
2.568 ns, peak-to-peak 11.358 ns
[Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium).

## Structured Kalman filter for timescales

The structured Kalman filter is an alternative implementation of the
conventional Kalman filter (CKF) for atomic-clock ensembles that
decomposes the system using observable Kalman canonical decomposition
parameterised by a transformation matrix `Γ`, runs the Kalman update
on the **observable** subsystem only, and propagates the **unobservable**
subsystem open-loop
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale).

It addresses the long-known undetectability problem in atomic-clock-
ensemble Kalman filters: because phase differences are the only
observable, the conventional Kalman covariance grows without bound and
incurs numerical instability
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale). The
unbounded covariance growth lives only in the unobservable subsystem,
which is propagated open-loop, so the structured implementation
avoids the numerical instability seen in CKF realisations
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale).

The atomic time prediction error is the weighted average across the
`m` clocks:

```math
TA[k] \;:=\; D\,(x[k] - \hat x[k]), \qquad
D \;=\; \tfrac{1}{m}\, C\,(I_n \otimes \mathbf{1}_m^{\top}).
```

The discrete process-noise covariance for the homogeneous `n`-th-order
`m`-clock ensemble used in the structured filter is

```math
W \;=\; \int_{0}^{\tau} A_t\, \mathrm{diag}(q_1^{2}, \ldots, q_n^{2})\, A_t^{\top}\, dt \,\otimes\, I_m
```

[Yan et al. 2023](@cite yan-2023-structured-kf-timescale). The
transformation `Γ` is selected by minimising the convex quadratic cost

```math
J(\Gamma) \;=\; \sum_{k=0}^{T} \bigl[\delta_1\,(\mathbb{E}[TA[k]])^{2} + \delta_2\,V[TA[k]]\bigr],
```

with global minimum at `Γ = 0` under symmetric initial conditions
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale).

Setting `Γ = 0` with matched initial covariance
`P̂₀ = (I_n ⊗ V) P₀ (I_n ⊗ V)^\top` recovers CKF asymptotically
(with improved numerical stability); the structured KF differs from
CKF only by a `Γ L̂_k(H ε + w)` term in the unobservable prediction
error, which vanishes when `Γ = 0`
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale).
Optimising `Γ` over the convex objective combining the squared mean
and variance of the atomic-time prediction error can further narrow
the confidence interval, and the structured KF shows better
robustness in Allan deviation across choices of initial covariance
`P̂₀` than Greenhall's covariance-correction CKF
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale).

The symmetry conditions guaranteeing `Γ = 0` optimality are realistic
when all clocks share a UTC[0] receiver in the same lab, but break
when the receivers are independent — making nonzero `Γ` potentially
preferable in heterogeneous deployments
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale).

## Kalman filter performance bounds

Wu's LTI-equivalence framework treats the steady-state two-state
atomic-clock Kalman filter as three equivalent discrete-time linear
time-invariant systems with transfer functions `H_x(z)`, `H_y(z)`,
and `H_{y_int}(z)` derived from the converged Kalman gains
`(K_{s1}, K_{s2})`
[Wu et al. 2023](@cite wu-2023-kf-performance-lti):

```math
H_x(z) \;=\; \frac{K_{s1}(1 - z^{-1}) + K_{s2}\, T\, z^{-1}}{(1 - K_{s1})(1 - z^{-1})^{2} + K_{s1}(1 - z^{-1}) + K_{s2}\, T\, z^{-1}}.
```

```math
H_{y_{\text{int}}}(z) \;=\; \frac{K_{s2}\, T\, z^{-1}}{(1 - K_{s1})(1 - z^{-1})^{2} + K_{s1}(1 - z^{-1}) + K_{s2}\, T\, z^{-1}}.
```

The integrated-frequency output `ŷ_int,k` stays in time units and is
low-pass filtered, useful for downstream stability analysis
[Wu et al. 2023](@cite wu-2023-kf-performance-lti).

By exploiting LTI superposition, each component of the actual
observation — polynomial trend, periodic fluctuation, stochastic
noise — can be analysed independently, even when that component is
not part of the KF model
[Wu et al. 2023](@cite wu-2023-kf-performance-lti). The transfer-
function structure depends on the KF model (via `Φ`, `H`) and the
coefficients depend on `Q`, `R` via the steady-state gains;
observations enter only as inputs
[Wu et al. 2023](@cite wu-2023-kf-performance-lti).

The Allan variance for a clock with white PM, white FM, and
random-walk FM is

```math
\sigma_y^{2}(\tau) \;=\; \frac{3\sigma_0^{2}}{\tau^{2}} + \frac{\sigma_1^{2}}{\tau} + \frac{\sigma_2^{2}\,\tau}{3},
```

used to recover the diffusion coefficients by slope fits
[Wu et al. 2023](@cite wu-2023-kf-performance-lti). Spectral-density
relationships connect the time-domain `σ_i²` to one-sided `h_α`
coefficients: `σ₀² = h₂ f_h / (4π²)`, `σ₁² = h₀ / 2`,
`σ₂² = h₋₂ · 2π²`
[Wu et al. 2023](@cite wu-2023-kf-performance-lti). The two-state
Kalman `Q` matrix used by the analysis is

```math
Q \;=\;
\begin{pmatrix}
\sigma_1^{2}\,T + \tfrac{1}{3}\sigma_2^{2}\,T^{3} & \tfrac{1}{2}\sigma_2^{2}\,T^{2} \\
\tfrac{1}{2}\sigma_2^{2}\,T^{2} & \sigma_2^{2}\,T
\end{pmatrix}, \qquad R \;=\; \sigma_0^{2}.
```

When a periodic fluctuation `A₀ sin(2π f_p t + φ₀)` is present in
observations, it adds a term
`(A₀ · 2π f_p)² sin⁴(π f_p τ) / (π τ)²` to the Allan variance
[Wu et al. 2023](@cite wu-2023-kf-performance-lti).

This lets the analyst predict periodic-fluctuation magnitude/phase
response and noise-driven estimate variance in closed form from the
steady-state gains alone, without rerunning the filter — a practical
diagnostic for tuning `Q` and `R` against an observed Allan-deviation
ladder
[Wu et al. 2023](@cite wu-2023-kf-performance-lti).

## See also

- [Theory: Clock State-Space Models](ensemble_overview.md) — supplies
  `Φ`, `Q`, and the underlying SDE the Kalman filter inverts.
- [Theory: Clock Steering with PID Controllers](steering.md) — closes
  the loop with a controller acting on the Kalman estimates.
- [Theory: Time-Scale Algorithms and Oscillator Networks](ensembles_and_oscillator_networks.md) —
  multi-clock ensemble timescales built on top of the Kalman recursion.
- [API: `SigmaTau.Est`](../reference/est.md) —
  [`StandardKalmanFilter`](@ref), [`UDFactorizedFilter`](@ref),
  [`predict!`](@ref), [`update!`](@ref).

## References

- [Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation) —
  state-space MLE for NBS cesium clocks.
- [Breakiron 2001](@cite breakiron-2001-kalman-timescales) — USNO
  maser-ensemble two-state Kalman filter.
- [Zucca & Tavella 2005](@cite zucca-2005-clock-model-allan) —
  three-state SDE underlying the Kalman model.
- [Zucca & Tavella 2015](@cite zucca-2015-clock-error-jumps) —
  innovation generation for jump-augmented simulation.
- [Åkesson et al. 2008](@cite akesson-2008-generalized-als) —
  generalised ALS noise tuning.
- [Kubczak et al. 2019](@cite kubczak-2019-fast-sync-rubidium) —
  rubidium fast-sync KF with single-scalar `Q`-tuning.
- [Ramos et al. 2022](@cite ramos-2022-ud-kalman-filter) — U-D
  factorised Kalman filter survey.
- [Wu et al. 2023](@cite wu-2023-kf-performance-lti) — LTI-equivalent
  steady-state KF performance bounds.
- [Yan et al. 2023](@cite yan-2023-structured-kf-timescale) —
  structured Kalman filter for clock ensembles.
- [Liu et al. 2024](@cite liu-2024-adaptive-kf-rubidium) — adaptive
  ALS-driven KF for rubidium GPSDO.
- [Stein 2003](@cite stein-2003-timescales) — innovation-driven
  random-shock estimates in clock ensembles.
- [Chaudhari 2022](@cite chaudhari-2022-uva-kalman-chapter) —
  pedagogical chapter on the Kalman filter and its variants
  (EKF / UKF / particle filter); useful as a from-first-principles
  refresher of the recursion this module implements.

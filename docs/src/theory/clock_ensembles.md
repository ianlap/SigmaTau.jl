# Theory: Clock Ensembles

A clock ensemble combines multiple physical oscillators into a single
realised timescale that is more uniform than any of its members. Two
analytic frameworks are relevant to `SigmaTau`:

1. The **stacked-state Kalman ensemble** of
   [Stein 2003](@cite stein-2003-timescales) and
   [Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation) —
   the joint state vector concatenates the per-clock states, `Φ` and
   `Q` are block-diagonal, and `H` selects the `N − 1` phase
   differences against a reference clock. This is what
   [`ClockEnsemble`](@ref) implements.
2. The **three-cornered-hat** estimator (and its `M`-clock
   generalisation, the `N`-cornered hat) for separating individual-
   clock variances from pairwise stability measurements when no fourth,
   more-stable reference is available
   [Riley & Howe 2008](@cite riley-2008-sp1065). This is a
   post-processing layer on top of the per-pair stability deviations
   from [`SigmaTau.Stab`](../reference/stab.md); a worked example lives
   at [`tutorials/06_three_cornered_hat.md`](../tutorials/06_three_cornered_hat.md).

## Basic time-scale equation

A time-scale algorithm estimates the time error of each clock in an
ensemble; the corrected ensemble time is more uniform than the time of
any individual clock [Stein 2003](@cite stein-2003-timescales). Because
only pairwise difference measurements are available — only
`z_{ij}(t_k) = x_i(t_k) − x_j(t_k) + v(t_k)` is observable —
individual clock corrections are formally unobservable, and the
algorithm picks one solution from an infinite ambiguous family by
imposing closure constraints
[Stein 2003](@cite stein-2003-timescales). The standard closure is a
weighted-sum-zero on the estimated phase, frequency, and frequency-
aging shocks across the ensemble:

```math
\sum_{i=1}^{N} a_i(t_k)\,\hat\varepsilon_i(t_k) \;=\; 0, \qquad
\sum_{i=1}^{N} b_i(t_k)\,\hat\eta_i(t_k) \;=\; 0, \qquad
\sum_{i=1}^{N} c_i(t_k)\,\hat\alpha_i(t_k) \;=\; 0,
```

with weights $a_i, b_i, c_i$ chosen to reflect each clock's noise
level [Stein 2003](@cite stein-2003-timescales). The same observability
obstruction motivated Tryon and Jones's original NBS Kalman ensemble,
which constrains drifts to sum to zero across the seven-clock ensemble
because clock readings are differential
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).
USNO maser-ensemble Kalman timescales follow the same pattern with
inverse-variance weights and an upper bound on any single clock's
contribution
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

In `SigmaTau`, the joint state-space realisation of this framework is
[`ClockEnsemble`](@ref): a homogeneous tuple of `AbstractClockModel`
members with block-diagonal `Φ`/`Q`, an `H` matrix that selects
phase differences against a reference clock, and a 1-σ reference-shared
measurement-noise covariance `R`. The ensemble is itself an
`AbstractClockModel`, so the existing
[`predict!`](@ref) / [`update!`](@ref) / [`prop!`](@ref) loop on
[`StandardKalmanFilter`](@ref) consumes it unchanged.

## Inverse-noise weights

The weight vectors `a_i, b_i, c_i` in the closure constraints above are
not unique; choosing them is what picks a single solution from the
otherwise under-determined ensemble problem
[Stein 2003](@cite stein-2003-timescales). Stein §VI–VII proposes the
inverse-noise specialisation

```math
a_i \;\propto\; 1/q_{1,i}, \qquad
b_i \;\propto\; 1/q_{2,i}, \qquad
c_i \;\propto\; 1/q_{3,i},
```

normalised so each vector sums to one. Splitting the weight choice —
phase weights from random-walk-phase-noise levels and frequency
weights from random-walk-frequency-noise levels — produces a
time-scale 3 % better than the best clock at short `τ` and 8 % better
at long `τ` [Stein 2003](@cite stein-2003-timescales).

`SigmaTau.Est` derives these weights automatically from each clock's
diffusion coefficients (`q1`, `q2`, optional `q3`), surfaces them on
[`ClockEnsemble`](@ref) as [`EnsembleWeights`](@ref), and lets the
caller pass an explicit `weights::EnsembleWeights{N}` to override the
auto-derivation when desired. The joint stacked Kalman filter itself
does not consume the weights — they are an interpretation layer for
recovering individual-clock shock estimates from the joint state
(Stein §VI–VII eqs. 6.2 / 6.3 / 7.2 / 7.3), and shipping them as a
field of `ClockEnsemble` keeps the API stable for the future
per-pair shock-recombination pass.

## Observable vs non-observable modes

Even with closure-constraint weights chosen, the reference clock's
absolute phase is not observable in a Kalman ensemble that sees only
phase differences. Its covariance diagonal grows unbounded; every
observable linear combination of states (`x_i − x_ref` for `i ≠ ref`)
stays tight [Stein 2003](@cite stein-2003-timescales),
[Yan et al. 2023](@cite yan-2023-structured-kf-timescale). On the
`ClockEnsemble` output this manifests as a monotonically growing
`P[1, 1]` entry alongside a Riccati-bounded differential-mode
covariance — see the
"ClockEnsemble observable vs non-observable modes" testset in
[`test/est/runtests.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/test/est/runtests.jl)
for an explicit numerical demonstration.

## Three-cornered hat

The three-cornered hat is the classical technique for separating
individual-clock variances from a triple of pairwise frequency-
stability measurements when no fourth, more-stable reference is
available [Riley & Howe 2008](@cite riley-2008-sp1065). It is the
`M = 3` instance of the more general `M`-cornered closure that
recovers each clock's variance from the full set of `\binom{M}{2}`
pairwise variances under the assumption that the clocks are
statistically independent.

Frequency stability cannot be assessed with a lone oscillator using
classical phase techniques; verification requires pairwise comparisons
or three-cornered-hat comparisons among `≥ 3` clocks
[Riley & Howe 2008](@cite riley-2008-sp1065). The NIST SP1065 analysis
workflow places the three-cornered hat as a late stage: data precision
check, gap/outlier/jump preprocessing, drift analysis, variance
analysis, spectral analysis, outlier recognition, plotting, variance
selection, then the three-cornered hat for separating the oscillator
under test from the reference
[Riley & Howe 2008](@cite riley-2008-sp1065).

Tutorial [`tutorials/06_three_cornered_hat.md`](../tutorials/06_three_cornered_hat.md)
synthesises three independent free-running clocks via
`_gen_powerlaw_phase`, builds the three pairwise difference records,
runs [`adev`](@ref) on each, and solves the classical TCH linear system
to recover each clock's individual σ_y(τ). End-to-end the recovered σ
tracks ground truth to within ~5 % at small τ and illustrates the
"TCH break-points" (negative-variance recoveries clamped to zero) at
long τ. The tutorial also includes a prose callout for the two real-
world failure modes — correlated noises and one clock dominating —
which break the independence assumption underlying the closure.

## See also

- [Theory: Clock State-Space Models](ensemble_overview.md) — per-clock
  SDE underlying the stacked-state ensemble.
- [Theory: Kalman Filter and Variants](kalman.md) — supplies the
  recursion the ensemble model is consumed by.
- [Theory: Clock Steering with PID Controllers](steering.md) — closes
  the loop on a single clock or on a steered ensemble realisation.
- [`tutorials/06_three_cornered_hat.md`](../tutorials/06_three_cornered_hat.md) —
  worked TCH on synthetic data.
- [`tutorials/07_clock_ensemble.md`](../tutorials/07_clock_ensemble.md) —
  two-clock paper time scale on a Stein 2003 Figure 4 fixture.
- [API: `SigmaTau.Est`](../reference/est.md) — [`ClockEnsemble`](@ref),
  [`EnsembleWeights`](@ref).

## References

- [Stein 2003](@cite stein-2003-timescales) — basic time-scale equation
  and the closure-ambiguity unification.
- [Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation) —
  NBS Kalman-ensemble drift-sum-zero constraint.
- [Breakiron 2001](@cite breakiron-2001-kalman-timescales) — USNO
  maser-ensemble inverse-variance weighting and bias correction.
- [Riley & Howe 2008](@cite riley-2008-sp1065) — NIST SP1065
  three-cornered-hat workflow.
- [Yan et al. 2023](@cite yan-2023-structured-kf-timescale) —
  homogeneous-ensemble factorisation `F = A ⊗ I_m` and observable
  Kalman decomposition.

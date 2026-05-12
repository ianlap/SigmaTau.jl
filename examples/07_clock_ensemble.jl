# # Clock ensembles: building a paper-clock time scale from N clocks
#
# A *time scale* is a paper clock — a virtual reference assembled
# from the measurements of an ensemble of real clocks, whose σ_y(τ)
# is below the σ_y(τ) of any single member. This is how UTC, GPS
# time, and lab time scales are built.
#
# This tutorial follows Stein 2003, *"Time Scales Demystified"*
# (Proc. 2003 IFCS), §V–VII. Stein's stacked-state Kalman formulation
# concatenates the per-clock states into one long state vector, makes
# Φ and Q block-diagonal, and selects `N-1` phase differences against
# a reference clock as the measurement. The reference clock's
# absolute phase is non-observable (Stein §II), but every *observable*
# linear combination of states remains tight.
#
# This tutorial:
#
# 1. Builds two heterogeneous `ThreeStateClock`s (a Cs-like and an
#    Rb-like clock with different noise floors).
# 2. Combines them via `clockA + clockB` into a `ClockEnsemble` —
#    the operator constructs the joint state-space model.
# 3. Synthesises a noisy phase-difference record.
# 4. Streams the differences through `StandardKalmanFilter`
#    unmodified — the existing `predict!`/`update!` loop consumes the
#    ensemble model directly because `ClockEnsemble <: AbstractClockModel`.
# 5. Forms a weighted paper-clock time scale from the filter's joint
#    state estimate and compares its σ_y(τ) to each free-running
#    clock — qualitatively reproducing Stein Figure 4.
#
# See [Three-cornered hat](06_three_cornered_hat.md) for the
# complementary noise-separation problem on the *same* pairwise
# difference data.

using SigmaTau
using SigmaTau.Stab: _gen_powerlaw_phase
using LinearAlgebra
using Random
using Statistics: mean
using FFTW                          # AbstractFFTs backend for noise synth

# ## 1. Two clocks of different grade
#
# - **Clock A** — Cs-beam-like. Larger WFM (`q1`) — the dominant
#   short-τ noise of a commercial caesium standard.
# - **Clock B** — Rb-like. Smaller WFM but larger RWFM (`q2`) — the
#   Rb stalwart trade-off (better short-term, worse long-term).
#
# Both clocks tick at `τ = 1 s` and report phase only (`q0` is the
# WPM measurement-noise floor of the comparator).

Random.seed!(20260512)

τ  = 1.0
clockA = ThreeStateClock(tau=τ, q0=1e-22, q1=1e-22, q2=1e-32, q3=0.0)
clockB = ThreeStateClock(tau=τ, q0=1e-22, q1=2e-23, q2=1e-30, q3=0.0)

# ## 2. Build the ensemble
#
# The ergonomic ensemble constructor is `Base.:+`. The result is a
# `ClockEnsemble{2, ThreeStateClock}` carrying both clocks, Stein
# auto-weights, and a reference-clock index (default 1).

ensemble = clockA + clockB

println("ensemble type:    ", typeof(ensemble))
println("joint nstates:    ", nstates(ensemble))
println("reference clock:  ", ensemble.ref)
println("Stein weights a:  ", ensemble.weights.a)
println("Stein weights b:  ", ensemble.weights.b)

# Take a look at the joint Φ — block-diagonal, with each 3×3 block
# being the polynomial integrator of one member clock.

Φ = state_transition(ensemble, τ)
println("\nJoint Φ (6×6 block-diagonal):")
display(Φ)

# And the measurement matrix H — a `(N-1) × (N·3) = 1 × 6` row
# vector picking the phase difference `x_B − x_A`.

H = measurement_matrix(ensemble)
println("\nMeasurement matrix H (selects x_B − x_A):")
display(H)

# ## 3. Synthesise a phase-difference record
#
# Independent random seeds for each clock guarantee statistically
# independent phase records — Stein's precondition for the joint
# ensemble update. We use the package's `_gen_powerlaw_phase`
# generator with α = 0 (white FM) so the dominant noise is the WFM
# channel, matching the Cs/Rb-like noise floors we set on `clockA`
# and `clockB`.

N = 4096

# WFM-dominated records, scaled to each clock's q1 floor. The
# `_gen_powerlaw_phase` helper returns unit-scale phase; we multiply
# by `sqrt(q1)` so the WFM σ_y at τ=1 matches the clock's diffusion.
xA = sqrt(clockA.q1) .* _gen_powerlaw_phase(0, N; tau0=τ)
xB = sqrt(clockB.q1) .* _gen_powerlaw_phase(0, N; tau0=τ)

# Measurement: noisy phase difference (`z = xB − xA + WPM noise`),
# where the WPM is the quadrature sum of both clocks' floor.
z = (xB .- xA) .+ sqrt(clockA.q0 + clockB.q0) .* randn(N)

# ## 4. Run the joint Kalman filter
#
# `StandardKalmanFilter` consumes the ensemble unmodified.
# Initial state seeds the phase blocks from the first sample of each
# clock; initial covariance is generous on the absolute-phase
# diagonal entries (we know those are non-observable) and tight on
# the velocity/aging entries.

ns = nstates(ensemble)        # = 6
x0 = zeros(ns)
x0[1] = xA[1]                 # phase block of clock A
x0[4] = xB[1]                 # phase block of clock B
P0 = Matrix(1e-14 * I(ns))
P0[1, 1] = 1.0                # absolute phase A: free to wander
P0[4, 4] = 1.0                # absolute phase B: free to wander

kf = StandardKalmanFilter(x0, P0)

xA_est = zeros(N)
xB_est = zeros(N)
yA_est = zeros(N)
yB_est = zeros(N)
for k in 1:N
    predict!(kf, ensemble, τ)
    update!(kf, ensemble, z[k])

    xA_est[k] = kf.x[1]
    yA_est[k] = kf.x[2]
    xB_est[k] = kf.x[4]
    yB_est[k] = kf.x[5]
end

# Sanity: the filter's predicted difference matches the observed
# difference to within the WPM noise floor.

diff_est = xB_est .- xA_est
println("\nResidual std (z − Hx̂): ",
        sqrt(mean((z .- diff_est).^2)),
        "    expected ≈ ", sqrt(clockA.q0 + clockB.q0))

# ## 5. Form the paper-clock time scale
#
# The classical Stein construction assembles a paper clock as the
# Stein-weighted sum of individual phase estimates:
#
# ```math
# x_\mathrm{TS}(t_k) = \sum_i a_i \, \hat x_i(t_k)
# ```
#
# Because the reference-clock phase is non-observable, the *absolute*
# value of `x_TS` is arbitrary, but its *increments* are well-defined
# and that's what σ_y(τ) measures.

x_TS = ensemble.weights.a[1] .* xA_est .+ ensemble.weights.a[2] .* xB_est

# ## 6. Compare ADEV of the time scale to each clock
#
# `adev` computes the Allan deviation curve at the requested integer
# `m` multiples of `τ₀`. We expect:
# - At every τ the ensemble σ_y sits at or below the better of the
#   two individual clocks at that τ — Stein Figure 4.

m_grid  = unique(round.(Int, exp10.(range(0, log10(N ÷ 4); length = 12))))
adev_A  = adev(PhaseData(xA,   τ), m_grid; calc_ci=false)
adev_B  = adev(PhaseData(xB,   τ), m_grid; calc_ci=false)
adev_TS = adev(PhaseData(x_TS, τ), m_grid; calc_ci=false)

# A few sample points for the prose record (full curves in the
# plot below):

for (i, τi) in enumerate(adev_A.tau)
    println("τ = ", τi,
            "    σ_y(A) = ",  adev_A.dev[i],
            "    σ_y(B) = ",  adev_B.dev[i],
            "    σ_y(TS) = ", adev_TS.dev[i])
end

# ## 7. Plot (Stein Figure 4 reproduction)
#
# When run interactively, the recipe in `SigmaTauRecipesBaseExt`
# renders each `StabilityResult` as a τ-vs-σ_y log-log curve.

# ```julia
# using Plots
# plot(adev_A,  label="Clock A (Cs-like)")
# plot!(adev_B,  label="Clock B (Rb-like)")
# plot!(adev_TS, label="Time scale (paper clock)", linewidth=2)
# xlabel!("τ (s)")
# ylabel!("σ_y(τ)")
# title!("Two-clock time scale (Stein 2003 Fig. 4)")
# ```
#
# The ensemble curve sits below the better of the two free-running
# clocks across the τ band — the paper clock is more stable than
# either of its components, exactly as Stein's algorithm promises.

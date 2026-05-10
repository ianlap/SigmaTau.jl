# estimators/filters.jl — Estimator Definitions and Standard Loop

"""
    AbstractEstimator

Supertype for state estimators consumed by the standardized
`predict!` / `update!` loop. Concrete subtypes hold the running state
mean and covariance (or an equivalent factorization) and must overload
`predict!` and `update!` against an [`AbstractClockModel`](@ref).

Shipped subtypes: [`StandardKalmanFilter`](@ref),
[`UDFactorizedFilter`](@ref) (stub), and [`KuramotoOscillator`](@ref)
(stub).
"""
abstract type AbstractEstimator end

# ── Legacy-compat helpers ────────────────────────────────────────────────────
#
# The legacy MATLAB/Julia filter applies safe_sqrt(Pm[i,i])^2 to every
# diagonal element of P after each update.  This effectively clamps any
# diagonal value whose abs is < 1e-10 to exactly 0.0 and takes abs() of
# the rest.  This is numerically significant when noise diffusion
# coefficients are very small (q ~ 1e-23).
#
# We expose this as an opt-in flag so researchers get an AD-clean path by
# default while legacy parity remains testable.

"""
    safe_sqrt_sq(x) → Float64

Legacy MATLAB artifact: `safe_sqrt(x)^2`.
If `|x| < 1e-10` → `0.0`, else `|x|`.
"""
@inline safe_sqrt_sq(x::Float64) = abs(x) < 1e-10 ? 0.0 : abs(x)

"""
    clamp_covariance_diag(P::SMatrix{2,2}) → SMatrix{2,2}

Apply `safe_sqrt_sq` to diagonal entries of a 2×2 static covariance matrix.
Returns a new (non-mutating) SMatrix.
"""
function clamp_covariance_diag(P::SMatrix{2,2,Float64})
    @SMatrix [safe_sqrt_sq(P[1,1]) P[1,2];
              P[2,1]               safe_sqrt_sq(P[2,2])]
end

"""
    clamp_covariance_diag(P::SMatrix{3,3}) → SMatrix{3,3}

Apply `safe_sqrt_sq` to diagonal entries of a 3×3 static covariance matrix.
Returns a new (non-mutating) SMatrix.
"""
function clamp_covariance_diag(P::SMatrix{3,3,Float64})
    @SMatrix [safe_sqrt_sq(P[1,1]) P[1,2]               P[1,3];
              P[2,1]               safe_sqrt_sq(P[2,2])  P[2,3];
              P[3,1]               P[3,2]                safe_sqrt_sq(P[3,3])]
end

# ── StandardKalmanFilter ─────────────────────────────────────────────────────

"""
    StandardKalmanFilter{V,M}

Mutable discrete-time Kalman filter state.

- `x::V`  — state vector  (SVector)
- `P::M`  — covariance    (SMatrix / Symmetric wrapper)
- `k::Int` — step counter
- `legacy_compat::Bool` — when `true`, apply `safe_sqrt_sq` diagonal
  clamping after every update to reproduce MATLAB-era behavior.
  Default `false` for AD-clean operation.
"""
mutable struct StandardKalmanFilter{V<:AbstractVector{Float64}, M<:AbstractMatrix{Float64}} <: AbstractEstimator
    x::V
    P::M
    k::Int
    legacy_compat::Bool
end

"""
    StandardKalmanFilter(x0, P0; legacy_compat=false)

Construct a `StandardKalmanFilter` from initial state `x0` and covariance
`P0`, converting to `SVector`/`SMatrix` internally for zero-allocation
dispatch.
"""
function StandardKalmanFilter(x0::AbstractVector{Float64}, P0::AbstractMatrix{Float64};
                              legacy_compat::Bool = false)
    n = length(x0)
    x = SVector{n, Float64}(x0...)
    P = SMatrix{n, n, Float64}(P0...)
    return StandardKalmanFilter(x, P, 0, legacy_compat)
end

"""
    UDFactorizedFilter <: AbstractEstimator

Bierman/Thornton U-D factorized Kalman filter intended for
low-observability scenarios such as lunar distance, where direct
covariance propagation can lose positive-definiteness. Propagates the
upper-triangular U and diagonal D factors via Weighted Modified
Gram–Schmidt and rank-one Bierman/Carlson updates.

!!! note "Stub implementation"
    This type is exported but no fields, no `predict!`, and no
    `update!` methods are defined for it yet. The intended
    implementation follows Ramos 2022 (WMGS time update + modified
    Agee–Turner rank-one measurement update with `U_R^{-1}`
    pre-decorrelation for correlated measurements).
"""
struct UDFactorizedFilter <: AbstractEstimator end # For low-observability lunar distance

"""
    KuramotoOscillator <: AbstractEstimator

Phase-coupled oscillator network framed as a distributed estimator,
targeted at pLEO SWaP-constrained nearest-neighbor clock
synchronization rather than centralized Kalman estimation.

!!! note "Stub implementation"
    This type is exported but no fields, no `predict!`, and no
    `update!` methods are defined for it yet. The intended
    implementation realizes the Kuramoto-style phase coupling on a
    nearest-neighbor topology suited to small spacecraft.
"""
struct KuramotoOscillator <: AbstractEstimator end # pLEO SWaP constrained nearest-neighbor

# ── PID steering controller ──────────────────────────────────────────────────

"""
    PIDController(; g_p=0.1, g_i=0.01, g_d=0.05)

Discrete PID controller for clock steering, matching the legacy
`filter.jl` controller. Holds running phase-error sum and the last
emitted steer for fold-in via `predict!(…; steering=…)`.
"""
Base.@kwdef mutable struct PIDController
    g_p::Float64 = 0.1
    g_i::Float64 = 0.01
    g_d::Float64 = 0.05
    sumx::Float64 = 0.0
    last_steer::Float64 = 0.0
end

"""
    step!(pid::PIDController, x::AbstractVector{<:Real}) → Float64

Compute and store the next steer value from the current Kalman state
estimate. Sign convention: drives phase (and frequency, when present) toward
zero. Mirrors legacy `step!(::PIDController, ::Vector{Float64})`.
"""
function step!(pid::PIDController, x::AbstractVector{<:Real})
    pid.sumx += x[1]
    steer = -pid.g_p * x[1] - pid.g_i * pid.sumx
    length(x) >= 2 && (steer -= pid.g_d * x[2])
    pid.last_steer = steer
    return steer
end

"""
    steer_to_correction(steer::Float64, ns::Int, dt::Float64) → SVector{ns,Float64}

Build the steering correction vector that `predict!(…; steering=…)`
expects. Phase component is `steer·dt`, frequency component is `steer`,
higher-order states are zero.
"""
function steer_to_correction(steer::Float64, ns::Int, dt::Float64)
    if ns == 1
        return SVector{1,Float64}(steer * dt)
    elseif ns == 2
        return SVector{2,Float64}(steer * dt, steer)
    else
        return SVector{ns,Float64}(steer * dt, steer, ntuple(_ -> 0.0, ns - 2)...)
    end
end

# ── The Standardized Update Loop ─────────────────────────────────────────────

"""
    predict!(est::StandardKalmanFilter, model::AbstractClockModel, dt::Real;
             steering::Union{Nothing,AbstractVector{Float64}} = nothing)

Propagate the estimator state forward in time by `dt`. Φ and Q are
re-derived from `model` for the supplied `dt` via the dt-aware
[`state_transition`](@ref) / [`process_noise`](@ref) overloads, so
`predict!(est, model, dt)` covers arbitrary horizons — not only the
discretization step `model.tau`. Passing `dt = model.tau` reproduces
the historical behaviour bit-exactly.

Optionally adds a `steering` correction vector to the predicted state
mean — matches the legacy `filter_step!` semantics where a PID
controller's last steer is folded into the state after the Φ
propagation (phase = +`u·dt`, frequency = +`u`).

On the first step (k == 0) the prediction is skipped — the initial
state is used directly, matching the legacy `filter_step!` convention
where prediction only fires when `s.k > 1` (after the first
increment). For an unconditional propagation that ignores this gate
(e.g. covariance-band projection from a fresh side-channel filter)
use [`prop!`](@ref).
"""
function predict!(est::StandardKalmanFilter, model::AbstractClockModel, dt::Real;
                  steering::Union{Nothing,AbstractVector{Float64}} = nothing)
    Phi = state_transition(model, dt)
    Q   = process_noise(model, dt)

    if est.k > 0
        x_pred = Phi * est.x
        if steering !== nothing
            n = length(x_pred)
            s = SVector{n, Float64}(ntuple(i -> i <= length(steering) ? steering[i] : 0.0, n))
            x_pred = x_pred + s
        end
        est.x = x_pred
        est.P = Phi * est.P * Phi' + Q
    end

    return est
end

"""
    prop!(est::StandardKalmanFilter, model::AbstractClockModel, dt::Real;
          steering::Union{Nothing,AbstractVector{Float64}} = nothing)

Unconditional covariance propagation. Advances `est.x ← Φ(dt) x` and
`est.P ← Φ(dt) P Φ(dt)' + Q(dt)` regardless of `est.k`. Unlike
[`predict!`](@ref) — which gates on `est.k > 0` to match the legacy
MATLAB `if k > 1` convention — `prop!` always propagates, which is
what you want for producing a 1σ covariance band around a
deterministic forward projection (e.g. shaded ±1σ holdover bounds).

Does not increment `est.k` (only [`update!`](@ref) does), so a sequence
`prop!(est, model, h·τ)` is safe to call from a freshly constructed
filter or in a side-channel "what-if" copy without disturbing the live
filter's update sequencing.

The optional `steering` correction is added to the predicted state mean
exactly as in `predict!`: phase = `+u·dt`, frequency = `+u`, higher
states zero.
"""
function prop!(est::StandardKalmanFilter, model::AbstractClockModel, dt::Real;
               steering::Union{Nothing,AbstractVector{Float64}} = nothing)
    Phi = state_transition(model, dt)
    Q   = process_noise(model, dt)

    x_pred = Phi * est.x
    if steering !== nothing
        n = length(x_pred)
        s = SVector{n, Float64}(ntuple(i -> i <= length(steering) ? steering[i] : 0.0, n))
        x_pred = x_pred + s
    end
    est.x = x_pred

    Pm   = SMatrix(est.P)
    Pnew = Phi * Pm * Phi' + Q
    est.P = Symmetric((Pnew + Pnew') ./ 2.0)

    return est
end

"""
    update!(est::StandardKalmanFilter, model::AbstractClockModel, z)

Scalar or vector measurement update.

Computes innovation, Kalman gain, and posterior covariance using
out-of-place StaticArrays math (AD-friendly — no in-place mutation).

When `est.legacy_compat == true`, the MATLAB-era `safe_sqrt_sq`
diagonal clamping is applied to P after the standard update.
"""
function update!(est::StandardKalmanFilter, model::AbstractClockModel, z::Union{Real, AbstractVector})
    est.k += 1

    H = measurement_matrix(model)
    R = measurement_noise(model)

    z_vec = z isa Real ? SVector{1, Float64}(z) : SVector{length(z), Float64}(z...)

    # Innovation
    ν = z_vec - H * est.x

    # Innovation covariance & Kalman gain
    Pm = SMatrix(est.P)   # ensure plain SMatrix for arithmetic
    S  = H * Pm * H' + R
    K  = Pm * H' / S

    # State update
    est.x = est.x + K * ν

    # Covariance update  (standard form, then symmetrize)
    n     = length(est.x)
    I_mat = SMatrix{n, n, Float64}(I)
    Pnew  = (I_mat - K * H) * Pm
    Pnew  = (Pnew + Pnew') ./ 2.0

    # Optional legacy diagonal clamping
    if est.legacy_compat
        Pnew = clamp_covariance_diag(SMatrix(Pnew))
    end

    est.P = Symmetric(Pnew)

    return est
end



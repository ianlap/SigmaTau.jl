# estimators/filters.jl — Estimator Definitions and Standard Loop

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

struct UDFactorizedFilter <: AbstractEstimator end # For low-observability lunar distance
struct KuramotoOscillator <: AbstractEstimator end # pLEO SWaP constrained nearest-neighbor

# ── The Standardized Update Loop ─────────────────────────────────────────────

"""
    predict!(est::StandardKalmanFilter, model::AbstractClockModel, dt::Float64)

Propagate the estimator state forward in time by `dt`.

On the first step (k == 0) the prediction is skipped — the initial
state is used directly, matching the legacy `filter_step!` convention
where prediction only fires when `s.k > 1` (after the first
increment).
"""
function predict!(est::StandardKalmanFilter, model::AbstractClockModel, dt::Float64)
    Phi = state_transition(model)
    Q   = process_noise(model)

    if est.k > 0
        est.x = Phi * est.x
        est.P = Phi * est.P * Phi' + Q
    end

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



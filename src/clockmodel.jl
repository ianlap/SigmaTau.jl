# Clock state-space model: phase, frequency and (optionally) drift driven by
# white noise, the model the Kalman filters and timescales run on, and its
# link to the power-law coefficients h_α.

"""
    ClockModel(; tau0, R=0, q1=0, q2=0, q3=0, drift=q3 > 0)

Discrete clock model with state `[phase, frequency]`, or
`[phase, frequency, drift]` when `drift = true`. Each state is driven by
white noise with the diffusion coefficients of Zucca & Tavella (2005):

| field | noise | S_y(f) | unit |
|:------|:------|:-------|:-----|
| `R`  | white PM, phase variance per sample at `tau0` | h₂ = 8π² τ₀ R | s² |
| `q1` | white FM, phase diffusion | h₀ = 2 q1 | s |
| `q2` | random-walk FM, frequency diffusion | h₋₂ = q2 / 2π² | 1/s |
| `q3` | random-run FM, drift diffusion | h₋₄ = q3 / 8π⁴ | 1/s³ |

`R` is the noise of a reading, taken as white; `q1`–`q3` are process noise,
see [`process_noise`](@ref). Build one from a noise fit with
[`clock_model`](@ref), and pass it to [`powerlaw_sigma`](@ref) or
[`powerlaw_psd`](@ref) to see the σ(τ) or spectrum it stands for.

$(TYPEDFIELDS)
"""
struct ClockModel
    "Sample interval in seconds, to which `R` refers."
    tau0::Float64
    "White PM phase variance per sample (s²)."
    R::Float64
    "White FM diffusion coefficient (s)."
    q1::Float64
    "Random-walk FM diffusion coefficient (1/s)."
    q2::Float64
    "Random-run FM diffusion coefficient (1/s³)."
    q3::Float64
    "Whether the state carries a frequency drift."
    drift::Bool
    function ClockModel(tau0, R, q1, q2, q3, drift)
        tau0 > 0 || throw(ArgumentError("ClockModel: tau0 must be positive, got $tau0"))
        min(R, q1, q2, q3) >= 0 || throw(ArgumentError("ClockModel: noise levels must be ≥ 0"))
        drift || q3 == 0 || throw(ArgumentError("ClockModel: q3 > 0 needs drift = true"))
        return new(tau0, R, q1, q2, q3, drift)
    end
end

ClockModel(; tau0::Real, R::Real=0.0, q1::Real=0.0, q2::Real=0.0, q3::Real=0.0,
           drift::Bool=q3 > 0) = ClockModel(Float64(tau0), Float64(R), Float64(q1),
                                            Float64(q2), Float64(q3), drift)

function Base.show(io::IO, ::MIME"text/plain", m::ClockModel)
    print(io, m.drift ? "3" : "2", "-state ClockModel, τ₀ = ", m.tau0, " s")
    @printf(io, "\n  R  %.4g s²\n  q1 %.4g s\n  q2 %.4g 1/s", m.R, m.q1, m.q2)
    m.drift && @printf(io, "\n  q3 %.4g 1/s³", m.q3)
end

"""
    clock_model(fit::PowerLawFit; drift=…) → ClockModel
    clock_model(h::AbstractDict; tau0, drift=…) → ClockModel

The [`ClockModel`](@ref) of a power-law model: white FM (α = 0),
random-walk FM (−2) and random-run FM (−4) drive the states, and the PM is
the reading noise `R`. `drift` defaults to whether h₋₄ is positive. `h`
takes α or noise names as keys, as in [`powerlaw_sigma`](@ref).

White PM (α = 2) is `R = h₂/8π²τ₀`. Flicker PM (α = 1) is not white and has
no state; the filter averages it like white PM over its own response time
and follows the slower part as if it were the clock. It is folded into `R`
as the measurement variance of
[van-dierendonck-1984-allan-kalman](@cite), the PM at the frequencies the
tracker does not follow: `R` is the white PM with the model's PM Allan
variance at the lag where PM hands over to FM. Leaving it out would make
the filter trust its readings too much.

Flicker FM and flicker walk (α = −1, −3) have no finite-state form and are
left out. Compare `powerlaw_sigma(model, τ)` with `powerlaw_sigma(fit, τ)`
to see what that costs over the τ you care about.

# Examples

```julia
using SigmaTau
p = noise_gen(PhaseData, 65536, 1.0; h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30))
m = clock_model(powerlaw_fit(p; alphas = [2, 0, -2]))
process_noise(m, 10.0)        # Q over a 10 s step
```
"""
function clock_model(fit::PowerLawFit; drift::Union{Bool,Nothing}=nothing)
    isfinite(fit.tau0) || throw(ArgumentError(
        "clock_model: the fit has no sample interval; use clock_model(Dict(zip(fit.alphas, fit.h)); tau0)"))
    return clock_model(_powerlaw_h(fit); tau0=fit.tau0, drift=drift)
end

function clock_model(h::AbstractDict; tau0::Real, drift::Union{Bool,Nothing}=nothing)
    hd = _powerlaw_h(h)
    hα(α) = get(hd, α, 0.0)
    q3 = 8π^4 * hα(-4)
    return ClockModel(; tau0, R=_reading_noise(hd, Float64(tau0)), q1=hα(0) / 2, q2=2π^2 * hα(-2),
                      q3, drift=drift === nothing ? q3 > 0 : drift)
end

# Variance `dev` at lag m of the components `αs` of a power-law model.
function _model_var(dev::Symbol, hd::Dict{Int,Float64}, αs::NTuple{N,Int}, m::Int, τ₀::Float64) where {N}
    σ² = 0.0
    for α in αs
        h = get(hd, α, 0.0)
        h > 0 && (σ² += h * _powerlaw_var(dev, α, m, τ₀))
    end
    return σ²
end

# Last lag, up to `_PM_HANDOVER_MAX`, at which the PM is still the larger part
# of the model's variance, and 1 when it never is: the Allan variance, or the
# Hadamard where a term steeper than random-walk FM keeps the Allan from
# converging. PM over FM falls steadily with the lag, so bisection finds it.
const _PM_HANDOVER_MAX = 2^16
function _pm_handover(hd::Dict{Int,Float64}, τ₀::Float64)
    dev = get(hd, -3, 0.0) > 0 || get(hd, -4, 0.0) > 0 ? :hdev : :adev
    leads(m) = _model_var(dev, hd, (2, 1), m, τ₀) > _model_var(dev, hd, (0, -1, -2, -3, -4), m, τ₀)
    leads(1) || return 1
    lo, hi = 1, _PM_HANDOVER_MAX
    leads(hi) && return hi
    while hi - lo > 1
        mid = (lo + hi) ÷ 2
        leads(mid) ? (lo = mid) : (hi = mid)
    end
    return lo
end

# R for white plus flicker PM: the white PM with the same Allan variance,
# 3R/τ², at the handover lag. White PM alone gives h₂/8π²τ₀ at every lag.
function _reading_noise(hd::Dict{Int,Float64}, τ₀::Float64)
    get(hd, 1, 0.0) > 0 || return get(hd, 2, 0.0) / (8π^2 * τ₀)
    m = _pm_handover(hd, τ₀)
    return _model_var(:adev, hd, (2, 1), m, τ₀) * (m * τ₀)^2 / 3
end

_powerlaw_tau0(m::ClockModel) = m.tau0

function _powerlaw_h(m::ClockModel)
    h = Dict(2 => 8π^2 * m.tau0 * m.R, 0 => 2m.q1, -2 => m.q2 / (2π^2))
    m.drift && (h[-4] = m.q3 / (8π^4))
    return h
end

"""
    state_transition(m::ClockModel, δ) → Matrix

The matrix Φ(δ) that carries the state over a step of δ seconds:
`[1 δ; 0 1]`, or `[1 δ δ²/2; 0 1 δ; 0 0 1]` with drift.
"""
function state_transition(m::ClockModel, δ::Real)
    # Filled in place: matrix literals are not trim-safe.
    Φ = Matrix{Float64}(I, m.drift ? 3 : 2, m.drift ? 3 : 2)
    Φ[1, 2] = δ
    m.drift && (Φ[1, 3] = δ^2 / 2; Φ[2, 3] = δ)
    return Φ
end

"""
    process_noise(m::ClockModel, δ) → Matrix

Covariance Q(δ) of the noise the clock picks up over a step of δ seconds,
the white-noise drives integrated over the step (Greenhall 2003, §2):

    Q₁₁ = q1 δ + q2 δ³/3 + q3 δ⁵/20    Q₁₂ = q2 δ²/2 + q3 δ⁴/8    Q₁₃ = q3 δ³/6
    Q₂₂ = q2 δ + q3 δ³/3               Q₂₃ = q3 δ²/2              Q₃₃ = q3 δ

Without drift, Q is the upper-left 2×2 block. It composes over steps:
`Q(a + b) = Φ(b) Q(a) Φ(b)ᵀ + Q(b)`.
"""
function process_noise(m::ClockModel, δ::Real)
    (; q1, q2, q3) = m
    q11 = q1 * δ + q2 * δ^3 / 3 + q3 * δ^5 / 20
    q12 = q2 * δ^2 / 2 + q3 * δ^4 / 8
    q22 = q2 * δ + q3 * δ^3 / 3
    Q = zeros(m.drift ? 3 : 2, m.drift ? 3 : 2)
    Q[1, 1], Q[1, 2], Q[2, 1], Q[2, 2] = q11, q12, q12, q22
    if m.drift
        Q[1, 3] = Q[3, 1] = q3 * δ^3 / 6
        Q[2, 3] = Q[3, 2] = q3 * δ^2 / 2
        Q[3, 3] = q3 * δ
    end
    return Q
end

# Clock state-space model: phase, frequency and (optionally) drift driven by
# white noise, with optional Markov states for flicker FM (in the frequency)
# and flicker PM (in the readings): the model the Kalman filters and
# timescales run on, and its link to the power-law coefficients h_α.

"""
    ClockModel(; tau0, R=0, q1=0, q2=0, q3=0, drift=q3 > 0, flicker_fm=0, flicker_pm=0, …)

Discrete clock model with state `[phase, frequency]`, or
`[phase, frequency, drift]` when `drift = true`. Each state is driven by
white noise with the diffusion coefficients of Zucca & Tavella (2005):

| field | noise | S_y(f) | unit |
|:------|:------|:-------|:-----|
| `R`  | white PM, phase variance per sample at `tau0` | h₂ = 8π² τ₀ R | s² |
| `q1` | white FM, phase diffusion | h₀ = 2 q1 | s |
| `q2` | random-walk FM, frequency diffusion | h₋₂ = q2 / 2π² | 1/s |
| `q3` | random-run FM, drift diffusion | h₋₄ = q3 / 8π⁴ | 1/s³ |
| `flicker_fm` | flicker FM, floor σ_y = √(2 ln 2 · h₋₁) | h₋₁ | – |
| `flicker_pm` | flicker PM of the readings | h₁ | s² |

`R` is the noise of a reading, taken as white, and `flicker_pm` the flicker
in it; `q1`–`q3` and `flicker_fm` are process noise, see
[`process_noise`](@ref). Build one from a noise fit with
[`clock_model`](@ref), and pass it to [`powerlaw_sigma`](@ref) or
[`powerlaw_psd`](@ref) to see the σ(τ) or spectrum it stands for.

With `flicker_fm > 0` the state gains Markov frequency processes that add
to the frequency, the model of flicker FM of
[davis-greenhall-stacey-2005-flicker](@cite) (its eqs 2.1–2.10, 3.3–3.4):
independent single-pole processes of equal variance whose time constants
rise geometrically. The paper spaces them by 8 or 16; here they are a
decade apart, and their sum is flicker FM within about 3 % in σ(τ) across
the τ they span. With `flicker_pm > 0` it gains Markov phase processes in
the same way, which add to each reading and not to the phase: measurement
noise with memory, which the filter estimates and keeps out of the phase.
That is the same construction one integration up and is ours, with no
published clock filter behind it.

The states sit only where they are needed: `fm_band` and `pm_band` are the
first and last τ (s) at which that flicker holds at least a quarter of the
model's Allan variance, looked for from τ₀ to `tau_max` (default 10⁶ τ₀),
and the time constants run a decade apart from the first to at least the
last. That is typically three to five states each. A flicker that never
reaches a quarter gets no states. What the states leave out at τ₀ is added
as white noise, to `q1` for flicker FM and to `R` for flicker PM, so the
model's Allan variance at τ₀ is that of the sampled flicker. Pass a band to
place the states yourself, `(0, 0)` for none. See the clock model page for
the levels and the accuracy.

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
    "Flicker FM level h₋₁."
    flicker_fm::Float64
    "First and last τ (s) the Markov frequency states span; `(0, 0)` for none."
    fm_band::NTuple{2,Float64}
    "Flicker PM level h₁ (s²)."
    flicker_pm::Float64
    "First and last τ (s) the Markov phase states span; `(0, 0)` for none."
    pm_band::NTuple{2,Float64}
    function ClockModel(tau0, R, q1, q2, q3, drift, flicker_fm=0.0, fm_band=(0.0, 0.0), flicker_pm=0.0, pm_band=(0.0, 0.0))
        tau0 > 0 || throw(ArgumentError("ClockModel: tau0 must be positive, got $tau0"))
        min(R, q1, q2, q3, flicker_fm, flicker_pm) >= 0 || throw(ArgumentError("ClockModel: noise levels must be ≥ 0"))
        drift || q3 == 0 || throw(ArgumentError("ClockModel: q3 > 0 needs drift = true"))
        for (h, b) in ((flicker_fm, fm_band), (flicker_pm, pm_band))
            b == (0, 0) || (h > 0 && tau0 <= b[1] <= b[2] < Inf) || throw(ArgumentError(
                "ClockModel: a band runs from at least tau0 upward and needs its flicker level, got $b"))
        end
        return new(tau0, R, q1, q2, q3, drift, flicker_fm, fm_band, flicker_pm, pm_band)
    end
end

function ClockModel(; tau0::Real, R::Real=0.0, q1::Real=0.0, q2::Real=0.0, q3::Real=0.0, drift::Bool=q3 > 0,
                    flicker_fm::Real=0.0, flicker_pm::Real=0.0, tau_max::Real=1e6 * tau0,
                    fm_band::Union{Nothing,NTuple{2,Real}}=nothing, pm_band::Union{Nothing,NTuple{2,Real}}=nothing)
    fm, pm = _flicker_bands(Float64(tau0), Float64(R), Float64(q1), Float64(q2), Float64(flicker_fm), Float64(flicker_pm), Float64(tau_max))
    return ClockModel(Float64(tau0), Float64(R), Float64(q1), Float64(q2), Float64(q3), drift,
                      Float64(flicker_fm), fm_band === nothing ? fm : Float64.(fm_band),
                      Float64(flicker_pm), pm_band === nothing ? pm : Float64.(pm_band))
end

# Where each flicker is needed: the first and last τ, on a grid of eight per
# decade from τ₀ to tau_max, at which it holds at least a quarter of the
# model's Allan variance (3 × flicker ≥ the rest). SP1065's long-τ forms:
# 3R/τ², h₁(1.038 + 3 ln(πτ/τ₀))/(4π²τ²), q1/τ, 2 ln 2 h₋₁, q2 τ/3.
function _flicker_bands(tau0::Float64, R::Float64, q1::Float64, q2::Float64, hfm::Float64, hpm::Float64, tau_max::Float64)
    fm_lo = fm_hi = pm_lo = pm_hi = 0.0
    0 < tau0 <= tau_max < Inf || return (0.0, 0.0), (0.0, 0.0)
    for j in 0:floor(Int, 8 * log10(tau_max / tau0) + 1e-9)
        τ = tau0 * 10.0^(j / 8)
        white = 3R / τ^2 + q1 / τ + q2 * τ / 3
        fm = 2log(2) * hfm
        pm = hpm * (1.038 + 3log(π * τ / tau0)) / (4π^2 * τ^2)
        if fm > 0 && 3fm >= white + pm
            fm_lo == 0 && (fm_lo = τ)
            fm_hi = τ
        end
        if pm > 0 && 3pm >= white + fm
            pm_lo == 0 && (pm_lo = τ)
            pm_hi = τ
        end
    end
    return (fm_lo, fm_hi), (pm_lo, pm_hi)
end

# Time constants of the Markov states on a band, a decade apart.
_band_taus(b::NTuple{2,Float64}) =
    b[1] > 0 ? Float64[b[1] * 10.0^j for j in 0:ceil(Int, log10(b[2] / b[1]) - 1e-9)] : Float64[]
_fm_taus(m::ClockModel) = _band_taus(m.fm_band)
_pm_taus(m::ClockModel) = _band_taus(m.pm_band)

# Variance of each Markov state: a sum of processes whose rates are a factor
# ρ apart has spectrum σ²/(f ln ρ) inside its band. For the frequency states
# that is S_y = h₋₁/f, so σ² = h₋₁ ln 10; for the phase states
# S_x = h₁/(4π² f), so σ² = h₁ ln 10/4π².
_fm_var(m::ClockModel) = m.flicker_fm * log(10)
_pm_var(m::ClockModel) = m.flicker_pm * log(10) / (4π^2)

# g(x) = 2x − 3 + 4e⁻ˣ − e⁻²ˣ: the phase variance an integrated Markov process
# of rate 1/T gains over δ = xT is σ²T² g(x), and its Allan variance at
# τ = xT is σ² g(x)/x². The series avoids the cancellation at small x.
function _markov_g(x::Float64)
    x > 0.05 && return 2x - 3 + 4exp(-x) - exp(-2x)
    return x^3 * (2 / 3 + x * (-1 / 2 + x * (7 / 30 + x * (-1 / 12 + x * (31 / 1260 + x * (-1 / 160 + x * 127 / 90720))))))
end

# White FM that stands for the flicker FM the Markov states leave out at τ₀:
# it brings the model's Allan variance there to that of sampled flicker FM,
# 2 h₋₁.
function _fm_q1(m::ClockModel)
    σ², x = _fm_var(m), 0.0
    for T in _fm_taus(m)
        x += σ² * _markov_g(m.tau0 / T) / (m.tau0 / T)^2
    end
    return m.tau0 * (2m.flicker_fm - x)
end

# White PM that does the same for the flicker PM: a stationary Markov phase
# process has Allan variance σ²(1 − a)(3 − a)/τ², a = exp(−τ/T).
function _pm_R(m::ClockModel)
    m.flicker_pm > 0 || return 0.0
    σ², x = _pm_var(m), 0.0
    for T in _pm_taus(m)
        a = exp(-m.tau0 / T)
        x += σ² * (1 - a) * (3 - a)
    end
    return max(m.flicker_pm * _powerlaw_var(:adev, 1, 1, m.tau0) * m.tau0^2 - x, 0.0) / 3
end

# Variance of the white noise on a reading.
_reading_var(m::ClockModel) = m.R + _pm_R(m)

# A reading in state terms: the phase plus the Markov phase states.
function _reading_row(m::ClockModel)
    h = zeros(_nstates(m))
    h[1] = 1.0
    h[end-length(_pm_taus(m))+1:end] .= 1.0
    return h
end

# The model with other white levels.
_with(m::ClockModel; R=m.R, q1=m.q1, q2=m.q2) =
    ClockModel(m.tau0, R, q1, q2, m.q3, m.drift, m.flicker_fm, m.fm_band, m.flicker_pm, m.pm_band)

_npoly(m::ClockModel) = m.drift ? 3 : 2
_nstates(m::ClockModel) = _npoly(m) + length(_fm_taus(m)) + length(_pm_taus(m))

function Base.show(io::IO, ::MIME"text/plain", m::ClockModel)
    print(io, _nstates(m), "-state ClockModel, τ₀ = ", m.tau0, " s")
    @printf(io, "\n  R  %.4g s²\n  q1 %.4g s\n  q2 %.4g 1/s", m.R, m.q1, m.q2)
    m.drift && @printf(io, "\n  q3 %.4g 1/s³", m.q3)
    for (name, h, T) in (("FM", m.flicker_fm, _fm_taus(m)), ("PM", m.flicker_pm, _pm_taus(m)))
        h > 0 || continue
        @printf(io, "\n  flicker %s %.4g, ", name, h)
        isempty(T) ? print(io, "no states") : @printf(io, "%d Markov states, %.4g to %.4g s", length(T), first(T), last(T))
    end
end

"""
    clock_model(fit::PowerLawFit; drift=…, flicker=false, …) → ClockModel
    clock_model(h::AbstractDict; tau0, drift=…, flicker=false, …) → ClockModel

The [`ClockModel`](@ref) of a power-law model: white FM (α = 0),
random-walk FM (−2) and random-run FM (−4) drive the states, and the PM is
the reading noise `R`. `drift` defaults to whether h₋₄ is positive. `h`
takes α or noise names as keys, as in [`powerlaw_sigma`](@ref).

White PM (α = 2) is `R = h₂/8π²τ₀`. Flicker PM (α = 1) is not white and by
default has no state; the filter averages it like white PM over its own
response time and follows the slower part as if it were the clock. It is
folded into `R` as the measurement variance of
[van-dierendonck-1984-allan-kalman](@cite), the PM at the frequencies the
tracker does not follow: `R` is the white PM with the model's PM Allan
variance at the lag where PM hands over to FM. Leaving it out would make
the filter trust its readings too much.

Flicker FM and flicker walk (α = −1, −3) have no finite-state form and are
left out. Compare `powerlaw_sigma(model, τ)` with `powerlaw_sigma(fit, τ)`
to see what that costs over the τ you care about.

With `flicker = true` flicker FM and flicker PM are kept instead, as the
Markov states of [`ClockModel`](@ref), and `R` is the white PM alone.
`fm_band`, `pm_band` and `tau_max` pass through; from a fit of a deviation
`tau_max` defaults to its longest τ.

# Examples

```julia
using SigmaTau
p = noise_gen(PhaseData, 65536, 1.0; h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30))
m = clock_model(powerlaw_fit(p; alphas = [2, 0, -2]))
process_noise(m, 10.0)        # Q over a 10 s step
```
"""
function clock_model(fit::PowerLawFit; tau_max::Real=fit.source in POWERLAW_DEVIATIONS ? maximum(fit.x) : 1e6 * fit.tau0, kw...)
    isfinite(fit.tau0) || throw(ArgumentError(
        "clock_model: the fit has no sample interval; use clock_model(Dict(zip(fit.alphas, fit.h)); tau0)"))
    return clock_model(_powerlaw_h(fit); tau0=fit.tau0, tau_max, kw...)
end

function clock_model(h::AbstractDict; tau0::Real, drift::Union{Bool,Nothing}=nothing, flicker::Bool=false, kw...)
    hd = _powerlaw_h(h)
    hα(α) = get(hd, α, 0.0)
    q3 = 8π^4 * hα(-4)
    return ClockModel(; tau0, R=flicker ? hα(2) / (8π^2 * tau0) : _reading_noise(hd, Float64(tau0)),
                      q1=hα(0) / 2, q2=2π^2 * hα(-2), q3, drift=drift === nothing ? q3 > 0 : drift,
                      flicker_fm=flicker ? hα(-1) : 0.0, flicker_pm=flicker ? hα(1) : 0.0, kw...)
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
    m.flicker_fm > 0 && (h[-1] = m.flicker_fm)
    m.flicker_pm > 0 && (h[1] = m.flicker_pm)
    return h
end

"""
    state_transition(m::ClockModel, δ) → Matrix

The matrix Φ(δ) that carries the state over a step of δ seconds:
`[1 δ; 0 1]`, or `[1 δ δ²/2; 0 1 δ; 0 0 1]` with drift. Each Markov state
of time constant T adds a row and column after these, the frequency states
first: it decays by `a = exp(−δ/T)`, and a frequency state moves the phase
by `T (1 − a)` times itself
([davis-greenhall-stacey-2005-flicker](@cite), eqs 2.7 and 3.3).
"""
function state_transition(m::ClockModel, δ::Real)
    # Filled in place: matrix literals are not trim-safe.
    np, n = _npoly(m), _nstates(m)
    Φ = Matrix{Float64}(I, n, n)
    Φ[1, 2] = δ
    m.drift && (Φ[1, 3] = δ^2 / 2; Φ[2, 3] = δ)
    for (j, T) in enumerate(_fm_taus(m))
        Φ[1, np+j] = -T * expm1(-δ / T)
        Φ[np+j, np+j] = exp(-δ / T)
    end
    nf = np + length(_fm_taus(m))
    for (j, T) in enumerate(_pm_taus(m))
        Φ[nf+j, nf+j] = exp(-δ / T)
    end
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

A Markov frequency state m, with time constant T, variance σ² and
`a = exp(−δ/T)`, adds its own drive integrated over the step
([davis-greenhall-stacey-2005-flicker](@cite), eqs 2.8–2.9 and 3.4), and a
Markov phase state p only its own:

    Q₁₁ += σ² T² (2δ/T − 3 + 4a − a²)    Q₁ₘ = σ² T (1 − a)²    Qₘₘ = σ² (1 − a²)
    Qₚₚ = σ² (1 − a²)

The flicker FM that the states leave out at τ₀ adds white FM to `q1`.
"""
function process_noise(m::ClockModel, δ::Real)
    (; q2, q3) = m
    q1 = m.q1 + _fm_q1(m)
    np, n = _npoly(m), _nstates(m)
    q11 = q1 * δ + q2 * δ^3 / 3 + q3 * δ^5 / 20
    q12 = q2 * δ^2 / 2 + q3 * δ^4 / 8
    q22 = q2 * δ + q3 * δ^3 / 3
    Q = zeros(n, n)
    Q[1, 1], Q[1, 2], Q[2, 1], Q[2, 2] = q11, q12, q12, q22
    if m.drift
        Q[1, 3] = Q[3, 1] = q3 * δ^3 / 6
        Q[2, 3] = Q[3, 2] = q3 * δ^2 / 2
        Q[3, 3] = q3 * δ
    end
    σ² = _fm_var(m)
    for (j, T) in enumerate(_fm_taus(m))
        Q[1, 1] += σ² * T^2 * _markov_g(δ / T)
        Q[1, np+j] = Q[np+j, 1] = σ² * T * expm1(-δ / T)^2
        Q[np+j, np+j] = -σ² * expm1(-2δ / T)
    end
    nf = np + length(_fm_taus(m))
    for (j, T) in enumerate(_pm_taus(m))
        Q[nf+j, nf+j] = -_pm_var(m) * expm1(-2δ / T)
    end
    return Q
end

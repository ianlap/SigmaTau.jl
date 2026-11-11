# Clock state-space model: phase, frequency and (optionally) drift driven by
# white noise, with optional Markov states for flicker FM (in the frequency)
# and flicker PM (in the readings): the model the Kalman filters and time
# scales run on, and its link to the power-law coefficients h_α.

"""
    ClockModel(; tau0, R=0, q1=0, q2=0, q3=0, drift=q3 > 0, flicker_fm=0,
               flicker_pm=0, tau_max=1e6 * tau0, fm_band=nothing, pm_band=nothing)

State-space model of a clock read every `tau0` seconds: phase, frequency
and, with `drift = true`, frequency drift driven by white noises, with
optional Markov states for flicker noise.

The phase ``x``, frequency ``y`` and drift ``d`` follow the three-state
model of [Zucca and Tavella (2015), Eq. (1)](@cite zucca-2015-clock-jumps)
without its deterministic terms,

```math
dx = y\\,dt + \\sqrt{q_1}\\,dW_1, \\qquad dy = d\\,dt + \\sqrt{q_2}\\,dW_2,
\\qquad dd = \\sqrt{q_3}\\,dW_3 ,
```

with independent Wiener processes ``W_i``; ``q_1``, ``q_2`` and ``q_3``
are the squares of that paper's diffusion coefficients. A reading is the
phase plus white noise of variance ``R``. The levels stand for the
coefficients of ``S_y(f) = \\sum_α h_α f^α`` as follows, with the white PM
taken up to the Nyquist frequency ``1/(2τ_0)``:

| Field        | Noise                   | Level                      | Unit |
|:-------------|:------------------------|:---------------------------|:-----|
| `R`          | white PM of a reading   | ``h_2 = 8π^2 τ_0 R``       | s²   |
| `q1`         | white FM                | ``h_0 = 2 q_1``            | s    |
| `q2`         | random-walk FM          | ``h_{-2} = q_2/(2π^2)``    | 1/s  |
| `q3`         | random-run FM           | ``h_{-4} = q_3/(8π^4)``    | 1/s³ |
| `flicker_fm` | flicker FM              | ``h_{-1}``                 | –    |
| `flicker_pm` | flicker PM of a reading | ``h_1``                    | s²   |

`R` and `flicker_pm` describe the readings; `q1` to `q3` and `flicker_fm`
are the clock's own noise, the process noise of [`process_noise`](@ref).
[`clock_model`](@ref) builds a model from a noise fit, and
[`powerlaw_sigma`](@ref) and [`powerlaw_psd`](@ref) take a model and return
the deviation or the spectrum of the power law it stands for.

With `flicker_fm > 0` the state gains Markov frequency processes that add
to the frequency: the flicker FM model of
[Davis et al. (2005), Eqs. (2.1)–(2.10)](@cite davis-greenhall-stacey-2005-flicker),
independent first-order processes of equal variance whose time constants
rise geometrically. The paper spaces them by a factor of 8 or 16; here
they are a decade apart, each of variance ``h_{-1} \\ln 10``. With
`flicker_pm > 0` the state gains Markov phase processes of variance
``h_1 \\ln 10/(4π^2)`` that add to each reading and not to the phase:
measurement noise with memory, which a filter estimates and keeps out of
the phase. The paper treats flicker FM only; using its construction on the
readings is this package's choice.

The states are placed where the flicker matters. By default `fm_band` and
`pm_band` are the first and last τ, between τ₀ and `tau_max`, at which
that flicker holds at least a quarter of the model's Allan variance. The
time constants start at the lower end of the band and rise by decades
until one reaches the upper end, which makes
``1 + ⌈\\log_{10}(τ_2/τ_1)⌉`` states for a band from ``τ_1`` to ``τ_2``.
A flicker that never holds a quarter gets no states. What the states leave
out at τ₀ is added as white noise, to `q1` for flicker FM and to `R` for
flicker PM, so that the model's Allan variance at τ₀ is that of the
sampled flicker. [Flicker states](@ref) gives the accuracy of the
approximation.

# Keywords

- `tau0`: sample interval in seconds.
- `R`: white PM variance of a reading, in s².
- `q1`, `q2`, `q3`: white FM, random-walk FM and random-run FM levels.
- `drift`: whether the state carries a frequency drift; `q3 > 0` needs it.
- `flicker_fm`, `flicker_pm`: flicker levels ``h_{-1}`` and ``h_1``.
- `tau_max`: longest τ, in seconds, at which the default bands are looked
  for.
- `fm_band`, `pm_band`: first and last τ, in seconds, that the Markov
  states span, in place of the default; `(0, 0)` for no states.

# Fields

$(TYPEDFIELDS)

# Examples

A hydrogen maser with white FM of ``10^{-12}`` at 1 s and a flicker floor
of ``1.5 × 10^{-13}``:

```jldoctest
julia> using SigmaTau

julia> ClockModel(; tau0 = 1.0, q1 = 5e-25, q2 = 1e-33, flicker_fm = 1.6e-26)
8-state ClockModel, τ₀ = 1.0 s
  R  0 s²
  q1 5e-25 s
  q2 1e-33 1/s
  flicker FM 1.6e-26, 6 Markov states, 10 to 1e+06 s
```
"""
struct ClockModel
    "Sample interval in seconds, to which `R` refers."
    tau0::Float64
    "White PM variance of a reading (s²)."
    R::Float64
    "White FM level ``q_1`` (s)."
    q1::Float64
    "Random-walk FM level ``q_2`` (1/s)."
    q2::Float64
    "Random-run FM level ``q_3`` (1/s³)."
    q3::Float64
    "Whether the state carries a frequency drift."
    drift::Bool
    "Flicker FM level ``h_{-1}``."
    flicker_fm::Float64
    "First and last τ (s) that the Markov frequency states span; `(0, 0)` for none."
    fm_band::NTuple{2,Float64}
    "Flicker PM level ``h_1`` (s²)."
    flicker_pm::Float64
    "First and last τ (s) that the Markov phase states span; `(0, 0)` for none."
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

# No reading noise and no process noise: a filter's gain is 0/0 once its start is used up.
_noiseless(m::ClockModel) = _reading_var(m) == 0 && iszero(process_noise(m, m.tau0))

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

# The model with every noise level times s: the same filter gains and
# predictions, and s times the covariance.
_scaled(m::ClockModel, s::Real) =
    ClockModel(m.tau0, s * m.R, s * m.q1, s * m.q2, s * m.q3, m.drift, s * m.flicker_fm, m.fm_band, s * m.flicker_pm, m.pm_band)

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
    clock_model(fit; drift, flicker=true, tau_max, fm_band, pm_band) → ClockModel
    clock_model(h::AbstractDict; tau0, drift, flicker=true, tau_max, fm_band,
                pm_band) → ClockModel

The [`ClockModel`](@ref) of a power-law noise model.

White FM, random-walk FM and random-run FM (α = 0, −2, −4) become `q1`,
`q2` and `q3`, and white PM (α = 2) becomes `R`, by the relations listed
under [`ClockModel`](@ref). The flicker noises have no finite state, and
`flicker` decides what becomes of them.

By default, flicker FM and flicker PM (α = −1, 1) are kept as the Markov
states of [`ClockModel`](@ref), across the τ where each holds a quarter of
the model's Allan variance, and `R` is the white PM alone. A fit without
flicker gives the same model with either setting.

With `flicker = false`, flicker FM is left out and flicker PM is folded
into `R`. The model then has two states, or three with drift, as the
`:stein` time scale of [`timescale`](@ref) needs, and it has less noise
than the fit: in simulation the filter of a clock on its flicker floor
claimed several times less than the error it made
([What the states buy](@ref)). The fold rests on what a filter does with
PM: it averages the PM that is faster than its own response and follows
the rest as if it were the clock. `R` stands for the part it averages, as
the measurement variance of
[Van Dierendonck et al. (1984), Eqs. (68)–(69)](@cite van-dierendonck-1984-allan-kalman)
is the phase noise that a tracking loop does not follow:

```math
R = \\frac{(m τ_0)^2\\, σ_{\\mathrm{PM}}^2(m τ_0)}{3} ,
```

with ``σ_{\\mathrm{PM}}^2`` the Allan variance of the white and flicker PM
together and ``m`` the last lag, up to 65 536, at which it exceeds that of
the FM terms (the Hadamard variances are compared when the model has
α = −3 or −4). For white PM alone this is ``h_2/(8π^2 τ_0)`` at every lag.
The paper gives the principle; the formula is this package's.

Flicker-walk FM (α = −3) is left out with either setting.
`powerlaw_sigma(model, τ)` beside `powerlaw_sigma(fit, τ)` shows what a
model leaves out at the τ of interest.

# Arguments

- `fit`: a [`PowerLawFit`](@ref) that has a sample interval, or a
  [`NoiseFit`](@ref).
- `h`: a `Dict` of ``h_α`` keyed by ``α`` or by noise type, as in
  [`powerlaw_sigma`](@ref).

# Keywords

- `tau0`: sample interval in seconds; a fit carries its own.
- `drift`: whether the model has a drift state; by default, whether
  ``h_{-4}`` is positive.
- `flicker`: keep flicker FM and flicker PM as Markov states; `false`
  leaves flicker FM out and folds flicker PM into `R`.
- `tau_max`, `fm_band`, `pm_band`: passed to [`ClockModel`](@ref). From a
  [`PowerLawFit`](@ref) of a deviation curve, `tau_max` defaults to its
  longest τ.

# Examples

```jldoctest
julia> using SigmaTau

julia> h = Dict(2 => 8π^2 * 1e-22, 0 => 2e-24, -2 => 1e-30 / 2π^2);

julia> m = clock_model(h; tau0 = 1.0);

julia> (m.R, m.q1, m.q2) .≈ (1e-22, 1e-24, 1e-30)
(true, true, true)
```

From a record, with the flicker that the fit finds kept as states:

```julia
p = noise_gen(PhaseData, 65536, 1.0;
              h = Dict(0 => 1e-24, -1 => 1.6e-26, -2 => 5e-34))
m = clock_model(powerlaw_fit(p; alphas = [2, 1, 0, -1, -2]))
```
"""
function clock_model(fit::PowerLawFit; tau_max::Real=fit.source in POWERLAW_DEVIATIONS ? maximum(fit.x) : 1e6 * fit.tau0, kw...)
    isfinite(fit.tau0) || throw(ArgumentError(
        "clock_model: the fit has no sample interval; use clock_model(Dict(zip(fit.alphas, fit.h)); tau0)"))
    return clock_model(_powerlaw_h(fit); tau0=fit.tau0, tau_max, kw...)
end

function clock_model(h::AbstractDict; tau0::Real, drift::Union{Bool,Nothing}=nothing, flicker::Bool=true, kw...)
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
    state_transition(m::ClockModel, δ) → Matrix{Float64}

Matrix ``Φ(δ)`` that carries the state of `m` over a step of `δ` seconds.

For phase, frequency and drift,

```math
Φ(δ) = \\begin{pmatrix} 1 & δ & δ^2/2 \\\\ 0 & 1 & δ \\\\ 0 & 0 & 1 \\end{pmatrix} ,
```

or its upper-left 2 × 2 block without drift. Each Markov state of time
constant ``T`` adds a row and a column after these, the frequency states
first: it decays by ``a = e^{-δ/T}``, and a frequency state moves the
phase by ``T(1 - a)`` times itself
([Davis et al. (2005), Eqs. (2.7) and (3.3)](@cite davis-greenhall-stacey-2005-flicker)).

# Examples

```jldoctest
julia> using SigmaTau

julia> state_transition(ClockModel(; tau0 = 1.0, q1 = 1e-24, q3 = 1e-40), 10.0)
3×3 Matrix{Float64}:
 1.0  10.0  50.0
 0.0   1.0  10.0
 0.0   0.0   1.0
```
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
    process_noise(m::ClockModel, δ) → Matrix{Float64}

Covariance ``Q(δ)`` of the noise that the state of `m` gathers over a step
of `δ` seconds.

The white noises of the model, integrated over the step, give
([Greenhall (2003), Eq. (5)](@cite greenhall-2003-jones-tryon);
[Zucca and Tavella (2015), Eq. (3)](@cite zucca-2015-clock-jumps))

```math
Q(δ) = \\begin{pmatrix}
q_1 δ + q_2 δ^3/3 + q_3 δ^5/20 & q_2 δ^2/2 + q_3 δ^4/8 & q_3 δ^3/6 \\\\
q_2 δ^2/2 + q_3 δ^4/8 & q_2 δ + q_3 δ^3/3 & q_3 δ^2/2 \\\\
q_3 δ^3/6 & q_3 δ^2/2 & q_3 δ
\\end{pmatrix} ,
```

or its upper-left 2 × 2 block without drift. Steps compose:
``Q(a + b) = Φ(b)\\,Q(a)\\,Φ(b)^{\\mathsf T} + Q(b)``, with ``Φ`` from
[`state_transition`](@ref).

A Markov frequency state ``j`` of time constant ``T`` and variance
``σ^2``, with ``a = e^{-δ/T}``, adds
([Davis et al. (2005), Eqs. (2.8)–(2.9) and (3.4)](@cite davis-greenhall-stacey-2005-flicker))

```math
Q_{11} \\mathrel{+}= σ^2 T^2 \\left(2δ/T - 3 + 4a - a^2\\right), \\qquad
Q_{1j} = σ^2 T (1 - a)^2, \\qquad Q_{jj} = σ^2 (1 - a^2) ,
```

and a Markov phase state only its own ``Q_{jj} = σ^2 (1 - a^2)``. The
white FM that stands for the flicker FM the states leave out at ``τ_0``
is added to ``q_1``.

# Examples

```jldoctest
julia> using SigmaTau

julia> m = ClockModel(; tau0 = 1.0, q1 = 1e-24, q2 = 3e-30);

julia> process_noise(m, 10.0) ≈ [1e-23+1e-27 1.5e-28; 1.5e-28 3e-29]
true

julia> Φ, Q = state_transition(m, 5.0), process_noise(m, 5.0);

julia> process_noise(m, 10.0) ≈ Φ * Q * Φ' + Q
true
```
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

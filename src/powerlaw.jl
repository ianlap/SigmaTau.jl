# Power-law noise model S_y(f) = Σ h_α f^α: fit the h_α to a measured σ(τ) or
# spectrum, and convert a model to σ(τ), S_y(f), S_x(f) or ℒ(f).

const POWERLAW_DEVIATIONS = (:adev, :mdev, :hdev, :mhdev)

"""
    PowerLawFit

Result of [`powerlaw_fit`](@ref): the coefficients ``h_α`` of the power-law
model ``S_y(f) = \\sum_α h_α f^α``, with ``S_y`` in 1/Hz, and the fit they
give. [`powerlaw_sigma`](@ref) and [`powerlaw_psd`](@ref) convert it to a
deviation or to another spectrum.

# Fields

$(TYPEDFIELDS)
"""
struct PowerLawFit
    "Exponents ``α`` of the components, highest first."
    alphas::Vector{Int}
    "Coefficients ``h_α ≥ 0``, one per ``α``."
    h::Vector{Float64}
    "What was fitted: a deviation (`:adev`, `:mdev`, `:hdev`, `:mhdev`) or a spectrum (`:Sy`, `:Sx`, `:Sphi`, `:L`)."
    source::Symbol
    "Sample interval of the record, in seconds; `NaN` for a spectrum given as points."
    tau0::Float64
    "Abscissa: ``τ`` in seconds for a deviation, ``f`` in Hz for a spectrum."
    x::Vector{Float64}
    "Measured values: ``σ(τ)`` for a deviation, ``S_y(f)`` in 1/Hz for a spectrum."
    y::Vector{Float64}
    "The model at each `x`, in the units of `y`."
    fit::Vector{Float64}
end

function Base.show(io::IO, ::MIME"text/plain", f::PowerLawFit)
    print(io, "PowerLawFit to ", f.source, " (", length(f.x), " points)")
    for (α, h) in zip(f.alphas, f.h)
        @printf(io, "\n  h%-3s %-5s %.4g", string(α), _nnls_symbol(α), h)
    end
end

# n-point Gauss–Legendre rule on [-1, 1], by Newton iteration on Pₙ.
function _gauss_legendre(n::Int)
    x = zeros(n); w = zeros(n)
    for i in 1:n
        z = cospi((i - 0.25) / (n + 0.5))
        dp = 0.0
        for _ in 1:100
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p1, p2 = ((2j - 1) * z * p1 - (j - 1) * p2) / j, p1
            end
            dp = n * (z * p1 - p2) / (z^2 - 1)
            dz = p1 / dp
            z -= dz
            abs(dz) < 1e-15 && break
        end
        x[i] = z
        w[i] = 2 / ((1 - z^2) * dp^2)
    end
    return x, w
end
const _GL_X, _GL_W = _gauss_legendre(16)

function _gl(f, a::Float64, b::Float64)
    c, r = (a + b) / 2, (b - a) / 2
    s = 0.0
    @inbounds for k in eachindex(_GL_X)
        s += _GL_W[k] * f(c + r * _GL_X[k])
    end
    return s * r
end

# Lobes of sin^p(mθ) integrated one by one; past this many, sin^p(mθ) is
# replaced by its mean, which moves the integral by about 1/_PL_LOBES².
const _PL_LOBES = 2048

# (c, p, q) of an estimator's response c·sin^p(πfτ)/(πfτ)^q to S_y(f); with
# the phase in steps of τ₀ the modified pair has m²sin²(πfτ₀) for two of
# the powers of πfτ.
_powerlaw_response(dev::Symbol) =
    dev === :adev ? (2.0, 4, 2) : dev === :hdev ? (8 / 3, 6, 2) :
    dev === :mdev ? (2.0, 6, 4) : (8 / 3, 8, 4)

# Variance σ²(τ = m·τ₀) that the unit component S_y(f) = (sin(πfτ₀)/(πτ₀))^α
# on [0, 1/(2τ₀)] gives in each estimator. That is the spectrum of a record
# sampled at τ₀ (and of `noise_gen`), which is f^α well below the Nyquist
# frequency. With θ = πfτ₀ and s = sin θ, the transfer functions are
#   AVAR 2 sin⁴(mθ)/(m² s²)    HVAR  (8/3) sin⁶(mθ)/(m² s²)
#   MVAR 2 sin⁶(mθ)/(m⁴ s⁴)    MHVAR (8/3) sin⁸(mθ)/(m⁴ s⁴)
function _powerlaw_var(dev::Symbol, α::Int, m::Int, τ₀::Float64)
    c, p, q = _powerlaw_response(dev)
    k = c / (Float64(m)^q * (π * τ₀)^(α + 1))   # (π τ₀)⁻¹ from df = dθ/(πτ₀)
    env(θ) = (s = sin(θ); s^(α - q))
    full(θ) = sin(m * θ)^p * env(θ)
    lobe = π / m
    nl = min(_PL_LOBES, floor(Int, (π / 2) / lobe))
    total = 0.0
    for j in 0:nl-1
        total += _gl(full, j * lobe, (j + 1) * lobe)
    end
    θc = nl * lobe
    if nl < _PL_LOBES
        θc < π / 2 && (total += _gl(full, θc, π / 2))
    else
        # Mean of sin^p over a period: 3/8, 5/16, 35/128 for p = 4, 6, 8.
        mean_p = p == 4 ? 3 / 8 : p == 6 ? 5 / 16 : 35 / 128
        edges = exp.(range(log(θc), log(π / 2); length=65))
        for j in 1:64
            total += mean_p * _gl(env, edges[j], edges[j+1])
        end
    end
    return k * total
end

function _powerlaw_alphas(alphas, deviation::Symbol)
    all(isinteger, alphas) || throw(ArgumentError("powerlaw: α must be whole numbers, got $(collect(alphas))"))
    αs = reverse!(sort!(unique(Int.(collect(alphas)))))
    isempty(αs) && throw(ArgumentError("powerlaw: `alphas` is empty"))
    αmin = deviation === :adev || deviation === :mdev ? -2 : -4
    for α in αs
        αmin <= α <= 2 || throw(ArgumentError(
            "powerlaw: α = $α is outside $αmin…2, where $deviation converges"))
    end
    return αs
end

_powerlaw_h(h::AbstractDict) = Dict(_noise_alpha(k) => Float64(v) for (k, v) in h)
_powerlaw_h(f::PowerLawFit) = Dict(zip(f.alphas, f.h))
_powerlaw_tau0(h) = NaN
_powerlaw_tau0(f::PowerLawFit) = f.tau0

"""
    powerlaw_fit(data, m_values=Int[]; deviation=:mhdev, alphas=2:-1:-2) → PowerLawFit
    powerlaw_fit(s::SpectralResult; alphas=2:-1:-2, f_carrier=nothing) → PowerLawFit
    powerlaw_fit(f, S; quantity=:Sy, alphas=2:-1:-2, f_carrier=nothing) → PowerLawFit

Coefficients ``h_α`` of the power-law model ``S_y(f) = \\sum_α h_α f^α``,
fitted by nonnegative least squares to the deviation of a record or to a
spectrum.

# From a record

The deviation of `data`, a `PhaseData` or `FrequencyData` record that may
have gaps, is computed without bias correction on the octave grid, or at
`m_values`. Each component enters with the variance it has in that
estimator on a record sampled at ``τ_0`` (see [`powerlaw_sigma`](@ref)).
That response is exact from ``m = 1``, so every point is fitted, and it
separates white PM from flicker PM even under ADEV, where both fall close
to ``τ^{-2}``.

A point is weighted by ``ν/(2σ^4)``, the inverse variance of a
``χ^2``-distributed estimate, with ``σ^2`` from the model and the EDF ``ν``
from the point's dominant component and its number of terms. An averaging
time at which gaps cost terms therefore counts for less.

# From a spectrum

The spectrum is a result of [`psd_y`](@ref), [`psd_x`](@ref) or
[`ssb_phase_noise`](@ref), or points `f`, `S` such as the phase-noise table
of a data sheet. A Welch estimate is fitted with the spectrum of a sampled
record, ``h_α (\\sin(πfτ_0)/(πτ_0))^α``, and points with ``h_α f^α``. Every
point has the same fractional weight, and a point at ``f = 0`` is left out.

# Keywords

- `deviation=:mhdev`: the estimator fitted, `:adev`, `:mdev`, `:hdev` or
  `:mhdev`.
- `alphas=2:-1:-2`: the exponents fitted. The Allan pair takes ``α`` from 2
  to −2. The Hadamard pair takes it to −4 and is insensitive to linear
  frequency drift.
- `quantity=:Sy`: what `S` holds: `:Sy` (1/Hz), `:Sx` (s²/Hz), `:Sphi`
  (rad²/Hz) or `:L` (dBc/Hz).
- `f_carrier`: the carrier frequency in Hz, needed for ``S_φ(f)`` and
  ``ℒ(f)``.

!!! note "Reading the result"
    A coefficient cannot be negative, so a component the record does not
    have comes out small and positive, at the expense of its neighbors. In
    simulation (white PM, white FM and random-walk FM, 65 536 samples) the
    components present were recovered about 2 % low with a scatter of 2 to
    3 %, except ``h_{-2}``, which rests on the few longest averaging times
    and scattered by 20 %. Compare `fit` with `y` before trusting a small
    coefficient.

[`noise_id_nnls`](@ref) fits the asymptotic power laws of the variance to
name the noise types present; this fit gives the calibrated ``h_α``.

# Examples

```julia
p = noise_gen(PhaseData, 65536, 1.0; h=Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30))
f = powerlaw_fit(p)                 # from MHDEV
powerlaw_sigma(f, [1.0, 100.0])     # the ADEV the model gives at 1 s and 100 s
```
"""
function powerlaw_fit(data::AbstractTimingData, m_values::AbstractVector{<:Integer}=Int[];
                      deviation::Symbol=:mhdev, alphas=2:-1:-2)
    deviation in POWERLAW_DEVIATIONS || throw(ArgumentError(
        "powerlaw_fit: deviation must be one of $(join(POWERLAW_DEVIATIONS, ", ")), got :$deviation"))
    αs = _powerlaw_alphas(alphas, deviation)
    τ₀ = data.tau0
    N  = _nphase(data)
    ms = isempty(m_values) ? _default_m_values(N, deviation) : collect(Int, m_values)
    # The estimators skip gaps, phase or frequency alike.
    res = deviation === :adev ? adev(data, ms; ci=false) :
          deviation === :mdev ? mdev(data, ms; ci=false) :
          deviation === :hdev ? hdev(data, ms; ci=false) : mhdev(data, ms; ci=false)
    return _powerlaw_fit_curve(res, deviation, αs, τ₀, ms)
end

"""
    powerlaw_fit(r::StabilityResult; alphas=…) → PowerLawFit

The same fit to an existing ADEV, MDEV, HDEV or MHDEV result `r`, such as one
clock's Allan deviation from [`groslambert`](@ref).
"""
function powerlaw_fit(r::StabilityResult; alphas=2:-1:-2)
    deviation = r.deviation_type
    deviation in POWERLAW_DEVIATIONS || throw(ArgumentError(
        "powerlaw_fit: a result of $deviation cannot be fitted; use :adev, :mdev, :hdev or :mhdev"))
    isfinite(r.tau0) && r.tau0 > 0 || throw(ArgumentError("powerlaw_fit: the result has no sample interval"))
    αs = _powerlaw_alphas(alphas, deviation)
    ms = round.(Int, r.tau ./ r.tau0)
    return _powerlaw_fit_curve(r, deviation, αs, r.tau0, ms)
end

function _powerlaw_fit_curve(res::StabilityResult, deviation::Symbol, αs::Vector{Int}, τ₀::Real, ms::Vector{Int})
    keep = findall(k -> isfinite(res.dev[k]) && res.dev[k] > 0, eachindex(res.dev))
    length(keep) >= length(αs) || throw(ArgumentError(
        "powerlaw_fit: $(length(keep)) usable τ points for $(length(αs)) components; " *
        "pass fewer `alphas` or a longer record"))
    mk, τ, v = ms[keep], res.tau[keep], res.dev[keep] .^ 2
    K = [_powerlaw_var(deviation, α, m, τ₀) for m in mk, α in αs]
    # The record length each τ's term count stands for: N without gaps, less
    # where gaps cost terms.
    Ns = [_gap_equivalent_n(deviation, mk[k], res.neff[keep[k]]) for k in eachindex(keep)]
    # χ² weights νₖ/(2σ²ₖ) taken at the model, not the measured σ̂²ₖ, which
    # would favor the points that happened to come out low.
    halfν = _nnls_edf_weights(deviation, K, v, αs, mk, τ, Ns, τ₀) .* v .^ 2
    h, _ = _nnls_weighted(K, v, halfν ./ v .^ 2)
    for _ in 1:3
        model = max.(K * h, 1e-3 .* v)
        h, _ = _nnls_weighted(K, v, halfν ./ model .^ 2)
    end
    return PowerLawFit(αs, h, deviation, τ₀, τ, sqrt.(v), sqrt.(K * h))
end

function powerlaw_fit(s::SpectralResult; alphas=2:-1:-2, f_carrier=nothing)
    τ₀ = 1 / (s.nperseg * (s.freq[2] - s.freq[1]))
    return _powerlaw_fit_psd(s.freq, s.psd, s.spectral_type, alphas, f_carrier, τ₀)
end

powerlaw_fit(f::AbstractVector{<:Real}, S::AbstractVector{<:Real};
             quantity::Symbol=:Sy, alphas=2:-1:-2, f_carrier=nothing) =
    _powerlaw_fit_psd(Float64.(f), Float64.(S), quantity, alphas, f_carrier, NaN)

# S_y from S_x, S_φ or ℒ. A sampled phase record is differenced, not
# differentiated, to frequency: |2 sin(πfτ₀)/τ₀|² instead of (2πf)².
function _powerlaw_to_sy(quantity::Symbol, f::Vector{Float64}, S::Vector{Float64},
                         f_carrier, τ₀::Float64)
    quantity === :Sy && return S
    d2 = isnan(τ₀) ? (2π .* f) .^ 2 : (2 .* sin.(π .* f .* τ₀) ./ τ₀) .^ 2
    quantity === :Sx && return d2 .* S
    quantity === :Sphi || quantity === :L ||
        throw(ArgumentError("powerlaw: quantity must be :Sy, :Sx, :Sphi or :L, got :$quantity"))
    f_carrier === nothing && throw(ArgumentError("powerlaw: S_φ(f) and ℒ(f) need `f_carrier`"))
    Sφ = quantity === :L ? 2 .* exp10.(S ./ 10) : S
    return d2 .* Sφ ./ (2π * f_carrier)^2
end

function _powerlaw_fit_psd(f::Vector{Float64}, S::Vector{Float64}, quantity::Symbol,
                           alphas, f_carrier, τ₀::Float64)
    length(f) == length(S) || throw(DimensionMismatch(
        "powerlaw_fit: $(length(f)) frequencies but $(length(S)) values"))
    αs = _powerlaw_alphas(alphas, :mhdev)
    Sy = _powerlaw_to_sy(quantity, f, S, f_carrier, τ₀)
    keep = findall(k -> f[k] > 0 && isfinite(Sy[k]) && Sy[k] > 0, eachindex(f))
    length(keep) >= length(αs) || throw(ArgumentError(
        "powerlaw_fit: $(length(keep)) usable points for $(length(αs)) components"))
    fk, v = f[keep], Sy[keep]
    g = isnan(τ₀) ? fk : sin.(π .* fk .* τ₀) ./ (π * τ₀)
    K = [g[k]^α for k in eachindex(fk), α in αs]
    # Equal fractional weight, against the measurement and then the model.
    h, _ = _nnls_weighted(K, v, 1.0 ./ v .^ 2)
    model = max.(K * h, 1e-3 .* v)
    h, _ = _nnls_weighted(K, v, 1.0 ./ model .^ 2)
    return PowerLawFit(αs, h, quantity, τ₀, fk, v, K * h)
end

"""
    powerlaw_sigma(h, tau; deviation=:adev, tau0) → Vector{Float64}

Deviation ``σ(τ)`` that the power-law model ``S_y(f) = \\sum_α h_α f^α``
gives in an estimator, for a record sampled at `tau0`.

The model is taken as the spectrum of a sampled record,
``h_α (\\sin(πfτ_0)/(πτ_0))^α`` up to the Nyquist frequency
``f_h = 1/(2τ_0)``. That is the spectrum [`noise_gen`](@ref) produces, and
it is ``h_α f^α`` well below ``f_h``. It is integrated through the transfer
function of the estimator at ``τ = mτ_0``:

```math
σ^2(τ) = \\int_0^{f_h} S_y(f)\\,\\lvert H_m(f) \\rvert^2\\,df .
```

| `deviation` | ``\\lvert H_m(f) \\rvert^2``, with ``θ = πfτ_0`` |
|:------------|:-------------------------------------------------|
| `:adev`     | ``2\\sin^4(mθ)/(m^2\\sin^2θ)``                     |
| `:mdev`     | ``2\\sin^6(mθ)/(m^4\\sin^4θ)``                     |
| `:hdev`     | ``(8/3)\\sin^6(mθ)/(m^2\\sin^2θ)``                 |
| `:mhdev`    | ``(8/3)\\sin^8(mθ)/(m^4\\sin^4θ)``                 |

For ADEV at ``m ≫ 1`` this is the usual conversion
([SP1065 §7.1, Eq. (67)](@cite riley-2008-sp1065);
[IEEE Std 1139-2022, Eq. (B.4)](@cite ieee1139-2022-definitions)):
``h_0/(2τ)`` for white FM, ``2\\ln 2\\,h_{-1}`` for flicker FM and
``(2π^2/3)\\,h_{-2}\\,τ`` for random-walk FM. The white FM term and the
white PM term, ``3 f_h h_2/(4π^2 τ^2)``, hold at every ``m``. At ``τ_0``
the sampled record is above the conversion by 44 % in variance on flicker
FM and 50 % on random-walk FM, and within 1 % from ``16τ_0``. On flicker
PM it is 19 % above at ``τ_0`` and still 4 % above at ``64τ_0``.

# Arguments

- `h`: a [`PowerLawFit`](@ref), a [`NoiseFit`](@ref), a
  [`ClockModel`](@ref), or a `Dict` of ``h_α`` keyed by ``α`` or by noise
  type, as in [`noise_gen`](@ref).
- `tau`: averaging times in seconds, each a whole multiple of `tau0`.

# Keywords

- `deviation=:adev`: `:adev`, `:mdev`, `:hdev` or `:mhdev`.
- `tau0`: the sample interval in seconds. A fit or a clock model carries
  its own.

# Examples

White FM with random-walk FM, at 1 s, 100 s and 10⁴ s:

```jldoctest
julia> using SigmaTau

julia> h = Dict(:WHFM => 2e-24, :RWFM => 1e-30);

julia> print(round.(powerlaw_sigma(h, [1.0, 100.0, 1e4]; tau0=1.0); sigdigits=4))
[1.0e-12, 1.032e-13, 2.567e-13]
```
"""
function powerlaw_sigma(h, tau::AbstractVector{<:Real}; deviation::Symbol=:adev,
                        tau0::Real=_powerlaw_tau0(h))
    deviation in POWERLAW_DEVIATIONS || throw(ArgumentError(
        "powerlaw_sigma: deviation must be one of $(join(POWERLAW_DEVIATIONS, ", ")), got :$deviation"))
    τ₀ = Float64(tau0)
    isfinite(τ₀) && τ₀ > 0 || throw(ArgumentError("powerlaw_sigma: pass `tau0`, the sample interval in seconds"))
    hd = _powerlaw_h(h)
    _powerlaw_alphas(keys(hd), deviation)
    return map(tau) do τ
        m = round(Int, τ / τ₀)
        m >= 1 && isapprox(m * τ₀, τ; rtol=1e-9) || throw(ArgumentError(
            "powerlaw_sigma: τ = $τ is not a whole multiple of tau0 = $τ₀"))
        sqrt(sum(hα * _powerlaw_var(deviation, α, m, τ₀) for (α, hα) in hd))
    end
end

"""
    powerlaw_psd(h, f; quantity=:Sy, f_carrier=nothing) → Vector{Float64}

The power-law model `h` as a spectrum at the Fourier frequencies `f`, in
Hz. `h` is a [`PowerLawFit`](@ref), a [`NoiseFit`](@ref), a
[`ClockModel`](@ref), or a `Dict` of ``h_α`` keyed by ``α`` or by noise
type. `quantity` selects the spectrum
([IEEE Std 1139-2022, Eqs. (A.13), (A.15) and (A.16)](@cite ieee1139-2022-definitions)):

- `:Sy`: ``S_y(f) = \\sum_α h_α f^α``, in 1/Hz.
- `:Sx`: ``S_x(f) = S_y(f)/(2πf)^2``, in s²/Hz.
- `:Sphi`: ``S_φ(f) = (2πν_0)^2 S_x(f)``, in rad²/Hz, at the carrier
  frequency ``ν_0`` given as `f_carrier`, in Hz.
- `:L`: ``ℒ(f) = \\tfrac{1}{2} S_φ(f)``, in dBc/Hz.

# Examples

White FM with ``σ_y = 10^{-12}`` at 1 s, on a 10 MHz carrier:

```jldoctest
julia> using SigmaTau

julia> L = powerlaw_psd(Dict(:WHFM => 2e-24), [1.0, 10.0, 100.0];
                        quantity=:L, f_carrier=10e6);

julia> print(round.(L; digits=2))
[-100.0, -120.0, -140.0]
```
"""
function powerlaw_psd(h, f::AbstractVector{<:Real}; quantity::Symbol=:Sy, f_carrier=nothing)
    hd = _powerlaw_h(h)
    Sy = [sum(hα * Float64(fk)^α for (α, hα) in hd) for fk in f]
    quantity === :Sy && return Sy
    Sx = Sy ./ (2π .* f) .^ 2
    quantity === :Sx && return Sx
    quantity === :Sphi || quantity === :L ||
        throw(ArgumentError("powerlaw_psd: quantity must be :Sy, :Sx, :Sphi or :L, got :$quantity"))
    f_carrier === nothing && throw(ArgumentError("powerlaw_psd: S_φ(f) and ℒ(f) need `f_carrier`"))
    Sφ = (2π * f_carrier)^2 .* Sx
    return quantity === :Sphi ? Sφ : 10 .* log10.(0.5 .* Sφ)
end

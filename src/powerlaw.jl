# Power-law noise model S_y(f) = Σ h_α f^α: fit the h_α to a measured σ(τ) or
# spectrum, and convert a model between σ(τ), S_y(f), S_x(f) and ℒ(f).

const POWERLAW_DEVIATIONS = (:adev, :mdev, :hdev, :mhdev)

"""
    PowerLawFit

Result of [`powerlaw_fit`](@ref): the coefficients of
`S_y(f) = Σ h_α f^α` (S_y in 1/Hz) and the fit they give. Pass it to
[`powerlaw_sigma`](@ref) or [`powerlaw_psd`](@ref) to convert it to another
domain.

$(TYPEDFIELDS)
"""
struct PowerLawFit
    "Exponents α of the components, highest first."
    alphas::Vector{Int}
    "Coefficients h_α ≥ 0, one per α."
    h::Vector{Float64}
    "What was fitted: a deviation (`:adev`, `:mdev`, `:hdev`, `:mhdev`) or a spectrum (`:Sy`, `:Sx`, `:L`)."
    source::Symbol
    "Sample interval of the record in seconds; `NaN` for a spectrum given as points."
    tau0::Float64
    "Abscissa: τ in seconds for a deviation, f in Hz for a spectrum."
    x::Vector{Float64}
    "Measured values: σ(τ) for a deviation, S_y(f) in 1/Hz for a spectrum."
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

# Variance σ²(τ = m·τ₀) that the unit component S_y(f) = (sin(πfτ₀)/(πτ₀))^α
# on [0, 1/(2τ₀)] gives in each estimator. That is the spectrum of a record
# sampled at τ₀ (and of `noise_gen`), which is f^α well below the Nyquist
# frequency. With θ = πfτ₀ and s = sin θ, the transfer functions are
#   AVAR 2 sin⁴(mθ)/(m² s²)    HVAR  (8/3) sin⁶(mθ)/(m² s²)
#   MVAR 2 sin⁶(mθ)/(m⁴ s⁴)    MHVAR (8/3) sin⁸(mθ)/(m⁴ s⁴)
function _powerlaw_var(dev::Symbol, α::Int, m::Int, τ₀::Float64)
    c, p, q = dev === :adev ? (2.0, 4, 2) : dev === :hdev ? (8 / 3, 6, 2) :
              dev === :mdev ? (2.0, 6, 4) : (8 / 3, 8, 4)
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
    αs = sort!(unique(Int.(collect(alphas))); rev=true)
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

"""
    powerlaw_fit(data [, m_values]; deviation=:mhdev, alphas=2:-1:-2) → PowerLawFit
    powerlaw_fit(s::SpectralResult; alphas=2:-1:-2, f_carrier=nothing) → PowerLawFit
    powerlaw_fit(f, S; quantity=:Sy, alphas=2:-1:-2, f_carrier=nothing) → PowerLawFit

Fit the power-law coefficients h_α of `S_y(f) = Σ h_α f^α` by nonnegative
least squares, either to a record's σ(τ) or to a spectrum.

**From a record** (`PhaseData` or `FrequencyData`): the deviation (`:adev`,
`:mdev`, `:hdev` or `:mhdev`, without bias correction) is computed on the
octave grid, or at `m_values`, and each α contributes its exact σ²(τ) for a
record sampled at τ₀ (see [`powerlaw_sigma`](@ref)). The response is exact
down to m = 1, so every point is fitted, and white and flicker PM stay
apart even under ADEV. The weights are the χ² ones, νₖ/(2σ⁴ₖ), with the EDF
from each point's dominant component and σ²ₖ from the model. The Allan pair takes α from −2 to 2;
the Hadamard pair, which also tolerates a linear frequency drift, from −4
to 2.

**From a spectrum:** a [`SpectralResult`](@ref) of [`psd_y`](@ref),
[`psd_x`](@ref) or [`ssb_phase_noise`](@ref) (which needs `f_carrier`), or
points `f`, `S` such as a datasheet's phase noise table, where `quantity`
is `:Sy` (1/Hz), `:Sx` (s²/Hz) or `:L` (dBc/Hz). A Welch estimate is fitted
with the sampled spectrum `h_α (sin(πfτ₀)/(πτ₀))^α`, points with `h_α f^α`.
Every point has the same fractional weight. The DC bin is left out.

The coefficients cannot go negative, so a component the record does not have
comes out small but positive, at a few percent of its neighbours' expense.
The long-τ octaves carry few samples: h₋₂ from a σ(τ) fit scatters by
about 20 % on 65 536 samples. Compare `fit` with `y` before trusting a
small h_α. [`noise_id_nnls`](@ref) fits pure power laws to name the
dominant noise type; this fit gives the calibrated h_α.

# Examples

```julia
using SigmaTau
p = noise_gen(PhaseData, 65536, 1.0; h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30))
f = powerlaw_fit(p)                 # from MHDEV
powerlaw_sigma(f, [1.0, 100.0])     # ADEV the model predicts at 1 s and 100 s
```
"""
function powerlaw_fit(data::PhaseData, m_values::Vector{Int}=Int[];
                      deviation::Symbol=:mhdev, alphas=2:-1:-2)
    deviation in POWERLAW_DEVIATIONS || throw(ArgumentError(
        "powerlaw_fit: deviation must be one of $(join(POWERLAW_DEVIATIONS, ", ")), got :$deviation"))
    αs = _powerlaw_alphas(alphas, deviation)
    x  = _f64(data.x)
    τ₀ = data.tau0
    ms = isempty(m_values) ? _default_m_values(length(x), deviation) : m_values
    res = deviation === :adev ? adev(data, ms; ci=false) :
          deviation === :mdev ? mdev(data, ms; ci=false) :
          deviation === :hdev ? hdev(data, ms; ci=false) : mhdev(data, ms; ci=false)
    keep = findall(k -> isfinite(res.dev[k]) && res.dev[k] > 0, eachindex(res.dev))
    length(keep) >= length(αs) || throw(ArgumentError(
        "powerlaw_fit: $(length(keep)) usable τ points for $(length(αs)) components; " *
        "pass fewer `alphas` or a longer record"))
    mk, τ, v = ms[keep], res.tau[keep], res.dev[keep] .^ 2
    K = [_powerlaw_var(deviation, α, m, τ₀) for m in mk, α in αs]
    # χ² weights νₖ/(2σ²ₖ) taken at the model, not the measured σ̂²ₖ, which
    # would favour the points that happened to come out low.
    halfν = _nnls_edf_weights(deviation, K, v, αs, mk, τ, length(x), τ₀) .* v .^ 2
    h, _ = _nnls_weighted(K, v, halfν ./ v .^ 2)
    for _ in 1:3
        model = max.(K * h, 1e-3 .* v)
        h, _ = _nnls_weighted(K, v, halfν ./ model .^ 2)
    end
    return PowerLawFit(αs, h, deviation, τ₀, τ, sqrt.(v), sqrt.(K * h))
end

powerlaw_fit(data::FrequencyData, m_values::Vector{Int}=Int[]; kwargs...) =
    powerlaw_fit(_freq_to_phase(data), m_values; kwargs...)

function powerlaw_fit(s::SpectralResult; alphas=2:-1:-2, f_carrier=nothing)
    τ₀ = 1 / (s.nperseg * (s.freq[2] - s.freq[1]))
    return _powerlaw_fit_psd(s.freq, s.psd, s.spectral_type, alphas, f_carrier, τ₀)
end

powerlaw_fit(f::AbstractVector{<:Real}, S::AbstractVector{<:Real};
             quantity::Symbol=:Sy, alphas=2:-1:-2, f_carrier=nothing) =
    _powerlaw_fit_psd(Float64.(f), Float64.(S), quantity, alphas, f_carrier, NaN)

# S_y from S_x or ℒ. A sampled phase record is differenced, not
# differentiated, to frequency: |2 sin(πfτ₀)/τ₀|² instead of (2πf)².
function _powerlaw_to_sy(quantity::Symbol, f::Vector{Float64}, S::Vector{Float64},
                         f_carrier, τ₀::Float64)
    quantity === :Sy && return S
    d2 = isnan(τ₀) ? (2π .* f) .^ 2 : (2 .* sin.(π .* f .* τ₀) ./ τ₀) .^ 2
    quantity === :Sx && return d2 .* S
    quantity === :L || throw(ArgumentError("powerlaw: quantity must be :Sy, :Sx or :L, got :$quantity"))
    f_carrier === nothing && throw(ArgumentError("powerlaw: ℒ(f) needs `f_carrier`"))
    return d2 .* 2 .* exp10.(S ./ 10) ./ (2π * f_carrier)^2
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

σ(τ) that the power-law model `h` gives in `deviation` (`:adev`, `:mdev`,
`:hdev` or `:mhdev`) for a record sampled at `tau0`. `h` is a
[`PowerLawFit`](@ref) (whose `tau0` is the default) or a `Dict` of h_α keyed
by α or noise type, as in [`noise_gen`](@ref). Every τ must be a whole
multiple of `tau0`.

The model is the spectrum of a sampled record, `h_α (sin(πfτ₀)/(πτ₀))^α` up
to the Nyquist frequency f_h = 1/(2τ₀), integrated through the estimator's
transfer function. Well above m = 1 this gives SP1065's closed forms, such
as `h₀/(2τ)` for white FM in ADEV, `2 ln 2 · h₋₁` for flicker FM and
`(2π²/3) h₋₂ τ` for random-walk FM. White and flicker PM depend on the
bandwidth, here f_h; `3 f_h h₂/(4π²τ²)` holds for white PM at every m.

# Examples

```julia
using SigmaTau
powerlaw_sigma(Dict(:WHFM => 2e-24), [1.0, 10.0, 100.0]; tau0=1.0)   # ≈ 1e-12 τ^-1/2
```
"""
function powerlaw_sigma(h, tau::AbstractVector{<:Real}; deviation::Symbol=:adev,
                        tau0::Real=(h isa PowerLawFit ? h.tau0 : NaN))
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

The power-law model `h` as a spectrum at the frequencies `f` (Hz):
`:Sy`, `S_y(f) = Σ h_α f^α` in 1/Hz; `:Sx`, `S_x(f) = S_y(f)/(2πf)²` in
s²/Hz; or `:L`, `ℒ(f) = ½ (2π f_carrier)² S_x(f)` in dBc/Hz. `h` is a
[`PowerLawFit`](@ref) or a `Dict` of h_α keyed by α or noise type.
"""
function powerlaw_psd(h, f::AbstractVector{<:Real}; quantity::Symbol=:Sy, f_carrier=nothing)
    hd = _powerlaw_h(h)
    Sy = [sum(hα * Float64(fk)^α for (α, hα) in hd) for fk in f]
    quantity === :Sy && return Sy
    Sx = Sy ./ (2π .* f) .^ 2
    quantity === :Sx && return Sx
    quantity === :L || throw(ArgumentError("powerlaw_psd: quantity must be :Sy, :Sx or :L, got :$quantity"))
    f_carrier === nothing && throw(ArgumentError("powerlaw_psd: ℒ(f) needs `f_carrier`"))
    return 10 .* log10.(0.5 .* (2π * f_carrier)^2 .* Sx)
end

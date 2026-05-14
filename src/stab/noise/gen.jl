# noise/gen.jl — Public power-law clock-noise generator.
#
# Wraps the internal `_gen_powerlaw_y` spectral shaper with a magnitude-aware
# API. Users specify the noise mixture by α exponent and either a target
# Allan deviation σ_y(τ=τ₀) per noise type or the corresponding fractional-
# frequency PSD coefficient h_α. Contributions for different α are drawn
# independently and summed.

# σ_y(τ=τ₀) from h_α via SP1065 Table 3 with the Nyquist bandwidth
# f_h = 1/(2τ₀) for the bandwidth-dependent WPM (α=2) and FPM (α=1) cases.
function _h_to_sigma1(alpha::Int, h::Float64, tau0::Float64)
    if alpha == 2
        return sqrt(3 * h / (8 * π^2 * tau0^3))
    elseif alpha == 1
        # 2π·f_h·τ = π at τ=τ₀, f_h=1/(2τ₀); SP1065 Table 3 numerator
        # collapses to 1.038 + 3·ln(π).
        coef = 1.038 + 3 * log(π)
        return sqrt(coef * h / (4 * π^2 * tau0^2))
    elseif alpha == 0
        return sqrt(h / (2 * tau0))
    elseif alpha == -1
        return sqrt(2 * log(2) * h)
    elseif alpha == -2
        return sqrt(2 * π^2 * h * tau0 / 3)
    else
        throw(ArgumentError(
            "noise_gen: unsupported α = $alpha; expected ∈ {-2, -1, 0, 1, 2}"))
    end
end

# Empirical σ_y(τ=τ₀) of a fractional-frequency vector via overlapping ADEV
# at m=1. Used to rescale a raw spectral-shaper realization to a known level.
function _measure_sigma1(y::Vector{Float64}, tau0::Float64)
    x = cumsum(y) .* tau0
    return _adev_core(x, [1], tau0)[1]
end

function _noise_gen_y(N::Int, tau0::Float64,
                     sigma1::AbstractDict, h::AbstractDict)
    !isempty(sigma1) && !isempty(h) &&
        throw(ArgumentError("noise_gen: pass either `sigma1` or `h`, not both"))
    isempty(sigma1) && isempty(h) &&
        throw(ArgumentError("noise_gen: noise mixture is empty; pass `sigma1` or `h`"))
    N >= 4 || throw(ArgumentError("noise_gen: N must be ≥ 4 (got $N)"))
    tau0 > 0 || throw(ArgumentError("noise_gen: tau0 must be > 0 (got $tau0)"))

    targets = Dict{Int, Float64}()
    if !isempty(sigma1)
        for (α, s) in sigma1
            α isa Integer ||
                throw(ArgumentError("noise_gen: α keys must be Int, got $(typeof(α))"))
            -2 <= α <= 2 ||
                throw(ArgumentError("noise_gen: α must be ∈ {-2,…,2}, got $α"))
            s >= 0 ||
                throw(ArgumentError("noise_gen: σ values must be ≥ 0, got $s for α=$α"))
            targets[Int(α)] = Float64(s)
        end
    else
        for (α, hα) in h
            α isa Integer ||
                throw(ArgumentError("noise_gen: α keys must be Int, got $(typeof(α))"))
            -2 <= α <= 2 ||
                throw(ArgumentError("noise_gen: α must be ∈ {-2,…,2}, got $α"))
            hα >= 0 ||
                throw(ArgumentError("noise_gen: h values must be ≥ 0, got $hα for α=$α"))
            targets[Int(α)] = _h_to_sigma1(Int(α), Float64(hα), tau0)
        end
    end

    y_total = zeros(Float64, N)
    for (α, σ_target) in targets
        σ_target == 0 && continue
        y_raw = _gen_powerlaw_y(α, N)
        σ_raw = _measure_sigma1(y_raw, tau0)
        σ_raw > 0 || continue
        y_total .+= y_raw .* (σ_target / σ_raw)
    end
    return y_total
end

"""
    noise_gen(::Type{PhaseData},     N, tau0; sigma1=Dict(), h=Dict()) → PhaseData
    noise_gen(::Type{FrequencyData}, N, tau0; sigma1=Dict(), h=Dict()) → FrequencyData

Synthesize a length-`N` clock record whose fractional-frequency power spectrum
is a sum of user-specified power laws `S_y(f) = h_α f^α`. Returns either a
`PhaseData` (integrated to phase with sample interval `tau0`) or a
`FrequencyData` (raw `y` sequence) depending on the first argument.

# α index

| α  | Noise type                | Phase PSD slope |
|---:|:--------------------------|:----------------|
|  2 | WPM — white phase         | flat            |
|  1 | FPM — flicker phase       | 1/f             |
|  0 | WFM — white frequency     | 1/f²            |
| -1 | FFM — flicker frequency   | 1/f³            |
| -2 | RWFM — random-walk freq.  | 1/f⁴            |

# Specifying the mixture

Pass *either* `sigma1` *or* `h` (not both):

- `sigma1[α] = σ` — target Allan deviation σ_y(τ=τ₀) of the α-component.
  Convenient when working from a clock spec sheet that quotes
  σ_y(1 s) directly.
- `h[α] = h_α` — fractional-frequency PSD coefficient in
  `S_y(f) = h_α f^α`. Unambiguous everywhere.

For WPM (α=2) and FPM (α=1) the SP1065 σ_y(τ) ↔ h_α relation depends on a
measurement bandwidth `f_h`. This function uses the Nyquist convention
`f_h = 1/(2τ₀)` — the highest frequency representable on a record sampled
at `τ₀`.

Each per-α component is rescaled to hit its requested amplitude exactly for
the drawn realization, so the empirical σ_y(τ=τ₀) of an isolated component
equals the requested value to numerical precision; the natural slope of the
power-law carries through to larger τ. Components from different α are
statistically independent.

Seed the global RNG (`Random.seed!`) before calling for reproducible output.

# Examples

```julia
using SigmaTau, Random
Random.seed!(0)

# Composite: WFM at σ_y(1 s) = 1e-12 plus RWFM at σ_y(1 s) = 1e-14
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))

# Same record but as fractional-frequency
y = noise_gen(FrequencyData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))

# Specify by h_α (here pure WFM with h_0 = 2e-24 → σ_y(1 s) ≈ 1e-12)
noise_gen(PhaseData, 8192, 1.0; h = Dict(0 => 2e-24))
```
"""
function noise_gen(::Type{PhaseData}, N::Int, tau0::Real;
                   sigma1::AbstractDict = Dict{Int,Float64}(),
                   h::AbstractDict      = Dict{Int,Float64}())
    τ₀ = Float64(tau0)
    y  = _noise_gen_y(N, τ₀, sigma1, h)
    return PhaseData(cumsum(y) .* τ₀, τ₀)
end

function noise_gen(::Type{FrequencyData}, N::Int, tau0::Real;
                   sigma1::AbstractDict = Dict{Int,Float64}(),
                   h::AbstractDict      = Dict{Int,Float64}())
    τ₀ = Float64(tau0)
    y  = _noise_gen_y(N, τ₀, sigma1, h)
    return FrequencyData(y, τ₀)
end

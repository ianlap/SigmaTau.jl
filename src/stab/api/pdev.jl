# api/pdev.jl — User wrapper for parabolic deviation

"""
$(SIGNATURES)

Parabolic deviation σ_PDEV(τ) for a phase record, per Vernotte–Lenczner–
Bourgeois–Rubiola (IEEE T-UFFC 63(4), 2016) and Vernotte–Chen–Rubiola 2020. PDEV
is built from a least-squares parabolic fit and is the recommended estimator when
evaluating the uncertainty of an ω-averaged frequency measurement. At `τ = τ₀`
(i.e. `m = 1`) PDEV reduces to overlapping ADEV.

When `ci=true`, the result populates per-τ noise type, EDF, and χ²-based
confidence bounds using the PVAR EDF model of Vernotte–Chen–Rubiola 2020
(arXiv:2005.13631), `ν ≈ 35 / (A(α)·(m/M) − 12·(m/M)²)` with `M = N − 2m`. At
`m = 1` the EDF is ADEV's (the PVAR(τ₀) ≡ AVAR(τ₀) identity); the large-τ region
(`m ≳ N/4`) uses the paper's semi-log interpolation down to `ν = 1` at `m = N/2`.
The approximation is accurate to ≈ ±10 % except at `m = 2`, where it reads high.
PVAR is an unbiased wavelet variance, so no bias correction is applied.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData(collect(1.0:10.0), 1.0);

julia> r = pdev(p, [1, 2]; ci=false);

julia> round.(r.dev; sigdigits=4)
2-element Vector{Float64}:
 0.0
 0.0
```
"""
function pdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _pdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0

    if !ci
        return StabilityResult(:pdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:pdev, raw_devs, noises, m_values, taus, length(x), (length(x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(x), confidence)

    return StabilityResult(:pdev, taus, raw_devs, noises, lower, upper, edfs)
end

pdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = pdev(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`.
pdev(data::PhaseData,     taus::TauMode; kwargs...) = pdev(data, tau_values(taus, length(data.x), :pdev); kwargs...)
pdev(data::FrequencyData, taus::TauMode; kwargs...) = pdev(data, tau_values(taus, length(data.y), :pdev); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to PDEV's algorithmic
# m-max (`(N − 1) ÷ 2`, see `_default_m_values`).
pdev(data::PhaseData;     kwargs...) = pdev(data, _default_m_values(length(data.x), :pdev); kwargs...)
pdev(data::FrequencyData; kwargs...) = pdev(data, _default_m_values(length(data.y), :pdev); kwargs...)

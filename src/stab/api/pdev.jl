# api/pdev.jl — User wrapper for parabolic deviation

"""
$(SIGNATURES)

Parabolic deviation σ_PDEV(τ) for a phase record, per Vernotte–Lenczner–
Bourgeois–Rubiola (IEEE T-UFFC 63(4), 2016) and Vernotte 2020. PDEV is
built from a least-squares parabolic fit and is the recommended
estimator when evaluating the uncertainty of an ω-averaged frequency
measurement. At `τ = τ₀` (i.e. `m = 1`) PDEV reduces to overlapping
ADEV.

No standard EDF / χ² confidence model is published for PDEV, so the
returned `noise_type`, `ci_lower`, `ci_upper`, and `edf` vectors are
empty. The `calc_ci` and `confidence` kwargs are accepted for API
uniformity.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData(collect(1.0:10.0), 1.0);

julia> r = pdev(p, [1, 2]);

julia> round.(r.dev; sigdigits=4)
2-element Vector{Float64}:
 0.0
 0.0
```
"""
function pdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    raw_devs = _pdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    return StabilityResult(:pdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
end

pdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = pdev(_freq_to_phase(data), m_values; kwargs...)

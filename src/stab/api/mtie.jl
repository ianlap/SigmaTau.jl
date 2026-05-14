# api/mtie.jl — User wrapper for Maximum Time Interval Error

"""
$(SIGNATURES)

Maximum Time Interval Error (MTIE) for a phase record, per ITU-T G.810.
For each averaging factor `m`, returns the largest peak-to-peak phase
excursion observed in any sliding window of `m+1` samples (spanning
`τ = m·τ₀`).

MTIE is a σ_x quantity (units of seconds), reported as a single
deterministic envelope rather than a statistical σ — there is no
published EDF / χ² confidence model for it, so `noise_type`,
`ci_lower`, `ci_upper`, and `edf` are returned empty even when
`calc_ci=true`. The kwarg is accepted for API uniformity with the
other deviations.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData([0.0, 1.0, 0.5, 2.0, 1.5], 1.0);

julia> r = mtie(p, [1, 3]);

julia> r.dev
2-element Vector{Float64}:
 1.5
 2.0
```
"""
function mtie(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    raw_devs = _mtie_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    return StabilityResult(:mtie, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
end

mtie(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mtie(_freq_to_phase(data), m_values; kwargs...)

# Zero-arg convenience: octave-spaced m_values up to MTIE's algorithmic
# m-max (`N − 1`, see `_default_m_values`).
mtie(data::PhaseData;     kwargs...) = mtie(data, _default_m_values(length(data.x), :mtie); kwargs...)
mtie(data::FrequencyData; kwargs...) = mtie(data, _default_m_values(length(data.y), :mtie); kwargs...)

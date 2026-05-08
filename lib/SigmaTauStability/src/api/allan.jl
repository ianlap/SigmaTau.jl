# api/allan.jl — User wrappers for stability calculations

"""
$(SIGNATURES)

Overlapping Allan deviation σ_y(τ) for a phase record, per IEEE 1139-2022 §C.2.
EDF for the χ²-based CI uses the closed-form approximation of [Greenhall2003](@cite).

`m_values` selects the analysis-interval factors (τ = m·τ₀). When
`calc_ci=true`, the result populates per-τ noise type, χ²-based confidence
bounds, and equivalent degrees of freedom.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData(collect(1.0:10.0), 1.0);

julia> r = adev(p, [1, 2]; calc_ci=false);

julia> round.(r.dev; sigdigits=4)
2-element Vector{Float64}:
 0.0
 0.0
```
"""
function adev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _adev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0

    if !calc_ci
        return StabilityResult(:adev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:adev, raw_devs, noises, m_values, taus, length(data.x), (length(data.x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:adev, taus, raw_devs, noises, lower, upper, edfs)
end

"""
$(SIGNATURES)

Modified Allan deviation Mod σ_y(τ) for a phase record, per IEEE 1139-2022 §C.3.

Uses a phase-averaged second difference; better than `adev` at separating
white-PM from flicker-PM noise.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData(collect(1.0:10.0), 1.0);

julia> r = mdev(p, [1, 2]; calc_ci=false);

julia> round.(r.dev; sigdigits=4)
2-element Vector{Float64}:
 0.0
 0.0
```
"""
function mdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _mdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0

    if !calc_ci
        return StabilityResult(:mdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:mdev, raw_devs, noises, m_values, taus, length(data.x), (length(data.x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:mdev, taus, raw_devs, noises, lower, upper, edfs)
end

"""
    tdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Time Deviation (TDEV) for the given PhaseData.

TDEV is defined as `τ · MDEV / √3`. Confidence-interval bounds inherit MDEV's
χ²/Gaussian limits scaled by the same `τ/√3` factor.
"""
function tdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    res = mdev(data, m_values; calc_ci=calc_ci, confidence=confidence)
    factor = res.tau ./ sqrt(3.0)

    if !calc_ci
        return StabilityResult(:tdev, res.tau, res.dev .* factor, Symbol[], Float64[], Float64[], Float64[])
    end

    return StabilityResult(:tdev, res.tau, res.dev .* factor, res.noise_type,
                           res.ci_lower .* factor, res.ci_upper .* factor, res.edf)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
adev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = adev(_freq_to_phase(data), m_values; kwargs...)
mdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mdev(_freq_to_phase(data), m_values; kwargs...)
tdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = tdev(_freq_to_phase(data), m_values; kwargs...)

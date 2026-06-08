# api/hadamard.jl — User wrappers for Hadamard stability calculations

"""
    hdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Overlapping Hadamard Deviation for the given PhaseData.
"""
function hdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _hdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0

    if !ci
        return StabilityResult(:hdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=3)
    edfs = calculate_edf(:hdev, raw_devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(x), confidence)

    return StabilityResult(:hdev, taus, raw_devs, noises, lower, upper, edfs)
end

"""
    mhdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Modified Hadamard Deviation for the given PhaseData.
"""
function mhdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _mhdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0

    if !ci
        return StabilityResult(:mhdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=3)
    edfs = calculate_edf(:mhdev, raw_devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(x), confidence)

    return StabilityResult(:mhdev, taus, raw_devs, noises, lower, upper, edfs)
end

"""
    htdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Hadamard time deviation (HTDEV) for the given `PhaseData`
— a third-difference time deviation, related to MHDEV by
`σ_x,HT(τ) = τ / √(10/3) · σ_y,MH(τ)`. HTDEV has units of seconds
(it is a σ_x quantity); HTDEV is to MHDEV what TDEV is to MDEV. The
construction is original to this package; SP1065, IEEE 1139-2022,
and NBS-TN-1337 do not define it.
"""
function htdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    res = mhdev(data, m_values; ci=ci, confidence=confidence)
    factor = res.tau ./ sqrt(10.0 / 3.0)

    if !ci
        return StabilityResult(:htdev, res.tau, res.dev .* factor, Symbol[], Float64[], Float64[], Float64[])
    end

    return StabilityResult(:htdev, res.tau, res.dev .* factor, res.noise_type,
                           res.ci_lower .* factor, res.ci_upper .* factor, res.edf)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
hdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = hdev(_freq_to_phase(data),  m_values; kwargs...)
mhdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mhdev(_freq_to_phase(data), m_values; kwargs...)
htdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = htdev(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`.
hdev(data::PhaseData,      taus::TauMode; kwargs...) = hdev(data,  tau_values(taus, length(data.x), :hdev);  kwargs...)
hdev(data::FrequencyData,  taus::TauMode; kwargs...) = hdev(data,  tau_values(taus, length(data.y), :hdev);  kwargs...)
mhdev(data::PhaseData,     taus::TauMode; kwargs...) = mhdev(data, tau_values(taus, length(data.x), :mhdev); kwargs...)
mhdev(data::FrequencyData, taus::TauMode; kwargs...) = mhdev(data, tau_values(taus, length(data.y), :mhdev); kwargs...)
htdev(data::PhaseData,     taus::TauMode; kwargs...) = htdev(data, tau_values(taus, length(data.x), :htdev); kwargs...)
htdev(data::FrequencyData, taus::TauMode; kwargs...) = htdev(data, tau_values(taus, length(data.y), :htdev); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to each kernel's
# algorithmic m-max (see `_default_m_values`). All kwargs pass through.
hdev(data::PhaseData;      kwargs...) = hdev(data,  _default_m_values(length(data.x), :hdev);  kwargs...)
hdev(data::FrequencyData;  kwargs...) = hdev(data,  _default_m_values(length(data.y), :hdev);  kwargs...)
mhdev(data::PhaseData;     kwargs...) = mhdev(data, _default_m_values(length(data.x), :mhdev); kwargs...)
mhdev(data::FrequencyData; kwargs...) = mhdev(data, _default_m_values(length(data.y), :mhdev); kwargs...)
htdev(data::PhaseData;     kwargs...) = htdev(data, _default_m_values(length(data.x), :htdev); kwargs...)
htdev(data::FrequencyData; kwargs...) = htdev(data, _default_m_values(length(data.y), :htdev); kwargs...)

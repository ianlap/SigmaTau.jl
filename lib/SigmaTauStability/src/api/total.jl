# api/total.jl — User wrappers for Total stability calculations

"""
    totdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Total Deviation for the given PhaseData.
"""
function totdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _totdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    if !calc_ci
        return StabilityResult(:totdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    B = bias_correction(noises, :totvar, taus, T)
    biased_devs = raw_devs ./ B

    edfs = calculate_edf(:totdev, biased_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(biased_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:totdev, taus, biased_devs, noises, lower, upper, edfs)
end

"""
    mtotdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Modified Total Deviation for the given PhaseData.
"""
function mtotdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _mtotdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    if !calc_ci
        return StabilityResult(:mtotdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    B = bias_correction(noises, :mtot, taus, T)
    biased_devs = raw_devs ./ B

    edfs = calculate_edf(:mtotdev, biased_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(biased_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:mtotdev, taus, biased_devs, noises, lower, upper, edfs)
end

"""
    htotdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Hadamard Total Deviation for the given PhaseData.
"""
function htotdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _htotdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    if !calc_ci
        return StabilityResult(:htotdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    B = bias_correction(noises, :htot, taus, T)
    biased_devs = raw_devs ./ B

    edfs = calculate_edf(:htotdev, biased_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(biased_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:htotdev, taus, biased_devs, noises, lower, upper, edfs)
end

"""
    mhtotdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Modified Hadamard Total Deviation for the given PhaseData.
"""
function mhtotdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _mhtotdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    if !calc_ci
        return StabilityResult(:mhtotdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    # MHTOTDEV doesn't have a known bias correction model in FCS 2001/SP1065.

    edfs = calculate_edf(:mhtotdev, raw_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:mhtotdev, taus, raw_devs, noises, lower, upper, edfs)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
totdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)   = totdev(_freq_to_phase(data),   m_values; kwargs...)
mtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = mtotdev(_freq_to_phase(data),  m_values; kwargs...)
htotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = htotdev(_freq_to_phase(data),  m_values; kwargs...)
mhtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mhtotdev(_freq_to_phase(data), m_values; kwargs...)

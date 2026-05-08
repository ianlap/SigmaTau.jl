# api/total.jl — User wrappers for Total stability calculations

"""
    totdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:howe, calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Total Deviation for the given PhaseData. See `_totdev_core` for
the meaning of `detrend`.
"""
function totdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:howe, calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _totdev_core(data.x, m_values, data.tau0; detrend=detrend)
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
    mtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Modified Total Deviation for the given PhaseData. See `_mtotdev_core` for
the meaning of `detrend`.
"""
function mtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _mtotdev_core(data.x, m_values, data.tau0; detrend=detrend)
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
    htotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Hadamard Total Deviation for the given PhaseData. See `_htotdev_core` for
the meaning of `detrend`.
"""
function htotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _htotdev_core(data.x, m_values, data.tau0; detrend=detrend)
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
    mhtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, confidence::Float64=0.95)

Modified Hadamard Total Deviation. See `_mhtotdev_core` for
the meaning of `detrend`.

No bias correction is applied. FCS 2001 and NIST SP1065 publish no
bias-correction model for MHTOTDEV; the estimator is treated as
unbiased (B = 1) by policy, matching Stable32 and AllanLab.
`bias_correction(:mhtot, …)` returns ones for the same reason. EDF
uses the empirical SP1065 fit coefficients (`_coeff_mhtot`).
"""
function mhtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _mhtotdev_core(data.x, m_values, data.tau0; detrend=detrend)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    if !calc_ci
        return StabilityResult(:mhtotdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:mhtotdev, raw_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:mhtotdev, taus, raw_devs, noises, lower, upper, edfs)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
totdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)   = totdev(_freq_to_phase(data),   m_values; kwargs...)
mtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = mtotdev(_freq_to_phase(data),  m_values; kwargs...)
htotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = htotdev(_freq_to_phase(data),  m_values; kwargs...)
mhtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mhtotdev(_freq_to_phase(data), m_values; kwargs...)

# api/hadamard.jl — User wrappers for Hadamard stability calculations

"""
    hdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Overlapping Hadamard Deviation for the given PhaseData.
"""
function hdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _hdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0
    
    if !calc_ci
        return StabilityResult(:hdev, taus, raw_devs, Symbol[], Float64[], Float64[])
    end
    
    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:hdev, raw_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)
    
    return StabilityResult(:hdev, taus, raw_devs, noises, lower, upper)
end

"""
    mhdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Modified Hadamard Deviation for the given PhaseData.
"""
function mhdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _mhdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0
    
    if !calc_ci
        return StabilityResult(:mhdev, taus, raw_devs, Symbol[], Float64[], Float64[])
    end
    
    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:mhdev, raw_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)
    
    return StabilityResult(:mhdev, taus, raw_devs, noises, lower, upper)
end

"""
    ldev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Hadamard Time Deviation (LDEV) for the given PhaseData.
"""
function ldev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    res = mhdev(data, m_values; calc_ci=calc_ci, confidence=confidence)
    factor = res.tau ./ sqrt(10.0 / 3.0)
    
    if !calc_ci
        return StabilityResult(:ldev, res.tau, res.dev .* factor, Symbol[], Float64[], Float64[])
    end
    
    return StabilityResult(:ldev, res.tau, res.dev .* factor, res.noise_type, res.ci_lower .* factor, res.ci_upper .* factor)
end

# api/allan.jl — User wrappers for stability calculations

"""
    adev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Allan Deviation for the given PhaseData.
"""
function adev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _adev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    
    if !calc_ci
        return StabilityResult(:adev, taus, raw_devs, Symbol[], Float64[], Float64[])
    end
    
    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:adev, raw_devs, noises, m_values, taus, length(data.x), (length(data.x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)
    
    return StabilityResult(:adev, taus, raw_devs, noises, lower, upper)
end

"""
    mdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)

Computes the Modified Allan Deviation for the given PhaseData.
"""
function mdev(data::PhaseData, m_values::Vector{Int}; calc_ci::Bool=true, confidence::Float64=0.95)
    raw_devs = _mdev_core(data.x, m_values, data.tau0)
    taus = m_values .* data.tau0
    
    if !calc_ci
        return StabilityResult(:mdev, taus, raw_devs, Symbol[], Float64[], Float64[])
    end
    
    noises = identify_noise(data.x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:mdev, raw_devs, noises, m_values, taus, length(data.x), (length(data.x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)
    
    return StabilityResult(:mdev, taus, raw_devs, noises, lower, upper)
end

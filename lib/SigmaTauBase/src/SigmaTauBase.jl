module SigmaTauBase

export AbstractTimingData, PhaseData, FrequencyData, StabilityResult

abstract type AbstractTimingData end

struct PhaseData{T<:AbstractFloat} <: AbstractTimingData
    x::Vector{T}      # Phase residuals in seconds
    tau0::Float64     # Base sample interval
end

struct FrequencyData{T<:AbstractFloat} <: AbstractTimingData
    y::Vector{T}      # Fractional frequency
    tau0::Float64
end

# The unified return struct for all stability calculations
struct StabilityResult
    deviation_type::Symbol  # e.g., :adev, :mdev
    tau::Vector{Float64}
    dev::Vector{Float64}
    noise_type::Vector{Symbol}
    ci_lower::Vector{Float64}
    ci_upper::Vector{Float64}
end

end # module

"""
    FrequencyData{T<:AbstractFloat}

Fractional-frequency samples `y(t)` at uniform interval `tau0` (default `1.0` s).

$(TYPEDFIELDS)
"""
struct FrequencyData{T<:AbstractFloat} <: AbstractTimingData
    "Fractional-frequency samples (dimensionless)."
    y::Vector{T}
    "Base sample interval τ₀ in seconds."
    tau0::Float64
    function FrequencyData(y::Vector{T}, tau0::Real = 1.0) where {T<:AbstractFloat}
        tau0 > 0 || throw(ArgumentError("FrequencyData: tau0 must be positive, got $tau0"))
        length(y) >= 2 || throw(ArgumentError(
            "FrequencyData: need at least 2 frequency samples, got $(length(y))"))
        return new{T}(y, Float64(tau0))
    end
end

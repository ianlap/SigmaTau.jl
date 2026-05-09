"""
    FrequencyData{T<:AbstractFloat}

Fractional-frequency samples `y(t)` at uniform interval `tau0`.

$(TYPEDFIELDS)
"""
struct FrequencyData{T<:AbstractFloat} <: AbstractTimingData
    "Fractional-frequency samples (dimensionless)."
    y::Vector{T}
    "Base sample interval τ₀ in seconds."
    tau0::Float64
end

"""
    PhaseData{T<:AbstractFloat}

Phase residuals `x(t)` sampled at uniform interval `tau0`.

$(TYPEDFIELDS)
"""
struct PhaseData{T<:AbstractFloat} <: AbstractTimingData
    "Phase samples in seconds."
    x::Vector{T}
    "Base sample interval τ₀ in seconds."
    tau0::Float64
end

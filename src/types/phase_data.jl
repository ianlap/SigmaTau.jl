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
    function PhaseData(x::Vector{T}, tau0::Real) where {T<:AbstractFloat}
        tau0 > 0 || throw(ArgumentError("PhaseData: tau0 must be positive, got $tau0"))
        length(x) >= 2 || throw(ArgumentError(
            "PhaseData: need at least 2 phase samples, got $(length(x))"))
        return new{T}(x, Float64(tau0))
    end
end

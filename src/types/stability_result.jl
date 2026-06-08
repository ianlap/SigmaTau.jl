"""
    StabilityResult

Unified return type for every stability calculation.

The `noise_type`, `ci_lower`, `ci_upper`, and `edf` vectors are empty when
the calculation was invoked with `ci=false`.

$(TYPEDFIELDS)
"""
struct StabilityResult
    "Which deviation produced this result (e.g. `:adev`, `:mdev`)."
    deviation_type::Symbol
    "Analysis intervals τ in seconds."
    tau::Vector{Float64}
    "Stability deviation σ_y(τ) per interval."
    dev::Vector{Float64}
    "Noise-type symbol identified at each τ (empty unless `ci=true`)."
    noise_type::Vector{Symbol}
    "Lower CI bound (empty unless `ci=true`)."
    ci_lower::Vector{Float64}
    "Upper CI bound (empty unless `ci=true`)."
    ci_upper::Vector{Float64}
    "Equivalent degrees of freedom (empty unless `ci=true`)."
    edf::Vector{Float64}
end

function Base.show(io::IO, r::StabilityResult)
    n = length(r.tau)
    ci = isempty(r.edf) ? "no CI" : "with CI"
    if n == 0
        print(io, "StabilityResult(:", r.deviation_type, ", 0 pts, ", ci, ")")
    else
        print(io, "StabilityResult(:", r.deviation_type, ", ", n, " pts, ",
                  "τ∈[", r.tau[1], ", ", r.tau[end], "] s, ", ci, ")")
    end
end

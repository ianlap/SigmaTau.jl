"""
    StabilityResult

Unified return type for every stability calculation.

The `noise_type`, `ci_lower`, `ci_upper`, and `edf` vectors are empty when
the calculation was invoked with `calc_ci=false`.

$(TYPEDFIELDS)
"""
struct StabilityResult
    "Which deviation produced this result (e.g. `:adev`, `:mdev`)."
    deviation_type::Symbol
    "Analysis intervals τ in seconds."
    tau::Vector{Float64}
    "Stability deviation σ_y(τ) per interval."
    dev::Vector{Float64}
    "Noise-type symbol identified at each τ (empty unless `calc_ci=true`)."
    noise_type::Vector{Symbol}
    "Lower CI bound (empty unless `calc_ci=true`)."
    ci_lower::Vector{Float64}
    "Upper CI bound (empty unless `calc_ci=true`)."
    ci_upper::Vector{Float64}
    "Equivalent degrees of freedom (empty unless `calc_ci=true`)."
    edf::Vector{Float64}
end

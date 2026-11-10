# `stability`: a suite of deviations in one call.
"""
    DEFAULT_DEVIATIONS

The estimators [`stability`](@ref) runs when `devs` is not given: `:adev`,
`:mdev`, `:hdev` and `:mhdev`. All four are fractional-frequency
deviations, so the default suite shares one ordinate. The time deviations,
in seconds, are run on request.
"""
const DEFAULT_DEVIATIONS = (:adev, :mdev, :hdev, :mhdev)

function _stability(data, N::Int, kind::Symbol, tau0::Float64;
                    devs = DEFAULT_DEVIATIONS, m = Octave,
                    ci::Bool = true, confidence::Real = DEFAULT_CONFIDENCE,
                    correct_bias::Union{Bool,Symbol} = true, remove_drift::Bool = true,
                    noise::_NoiseChoice = :auto)
    results = StabilityResult[]
    for sym in devs
        name = Symbol(sym)
        ms = m isa TauMode ? _resolve_m(name, m, N) : collect(Int, m)
        push!(results, deviation(name, data, ms; ci, confidence, correct_bias, remove_drift, noise))
    end
    tau_mode = m isa TauMode ? Symbol(string(m)) : :explicit
    return StabilitySuite(results, kind, tau0, N, ci ? Float64(confidence) : NaN, tau_mode)
end

"""
    stability(data; devs=DEFAULT_DEVIATIONS, m=Octave, ci=true,
              confidence=DEFAULT_CONFIDENCE, correct_bias=true, remove_drift=true,
              noise=:auto) → StabilitySuite

Several deviations of one record in a single call. Each name in `devs` is
run through [`deviation`](@ref), and the results come back as a
[`StabilitySuite`](@ref), indexed by name.

# Arguments

- `data`: a `PhaseData` or `FrequencyData` record.

# Keywords

- `devs`: the estimators to run, as symbols or strings from
  [`estimators`](@ref).
- `m`: a [`TauMode`](@ref) grid, cut for each estimator at the largest
  averaging factor it supports, or an integer vector that every estimator
  takes as it stands. `theo1` reads such a vector as even Thêo1 factors (an
  odd entry gives `NaN`) and `theoh` as factors of ``τ = mτ_0``, so prefer a
  `TauMode` when the suite mixes them with other estimators.
- `ci`, `confidence`, `noise`: as in [`adev`](@ref), for every estimator.
- `correct_bias`: reaches the total estimators and `theo1`.
- `remove_drift`: reaches the Hadamard family.

# Examples

```jldoctest
julia> using SigmaTau

julia> pd = PhaseData(collect(1.0:64.0), 1.0);

julia> suite = stability(pd; devs=(:adev, :hdev), m=Octave, ci=false);

julia> keys(suite)
2-element Vector{Symbol}:
 :adev
 :hdev

julia> suite[:adev].deviation_type
:adev
```
"""
stability(data::PhaseData; kwargs...) =
    _stability(data, length(data.x), :phase, data.tau0; kwargs...)

stability(data::FrequencyData; kwargs...) =
    _stability(data, length(data.y) + 1, :frequency, data.tau0; kwargs...)   # N + 1 phase points

# `stability`: a suite of deviations in one call.
"""
Default deviation set computed by [`stability`](@ref) when `devs` is not given:
overlapping Allan, modified Allan, overlapping Hadamard, and modified Hadamard.
All four are fractional-frequency (σ_y) quantities, so the default suite
overlays on one ordinate; time deviations (`tdev`, `htdev`) are σ_x quantities
in seconds and must be requested explicitly.
"""
const DEFAULT_DEVIATIONS = (:adev, :mdev, :hdev, :mhdev)

function _stability(data, N::Int, kind::Symbol, tau0::Float64;
                    devs = DEFAULT_DEVIATIONS, m = Octave, taus = nothing,
                    ci::Bool = true, confidence::Real = DEFAULT_CONFIDENCE,
                    correct_bias::Bool = true, remove_drift::Bool = true)
    if taus !== nothing
        Base.depwarn("stability: `taus` is deprecated, use `m`", :stability)
        m = taus
    end
    results = StabilityResult[]
    for sym in devs
        name = Symbol(sym)
        ms = m isa TauMode ? _resolve_m(name, m, N) : collect(Int, m)
        push!(results, deviation(name, data, ms; ci, confidence, correct_bias, remove_drift))
    end
    tau_mode = m isa TauMode ? Symbol(string(m)) : :explicit
    return StabilitySuite(results, kind, tau0, N, ci ? Float64(confidence) : NaN, tau_mode)
end

"""
    stability(data::PhaseData; devs=DEFAULT_DEVIATIONS, m=Octave,
              ci=true, confidence=DEFAULT_CONFIDENCE, correct_bias=true,
              remove_drift=true) → StabilitySuite

Compute a suite of stability deviations from one record in a single call, the
batch analog of `adev`, `mdev`, …. Each symbol in `devs` is run through
[`deviation`](@ref) with a per-deviation averaging-factor grid derived from
`m`, so every deviation gets a grid valid for its own algorithmic m-max.
Returns a [`StabilitySuite`](@ref). Also accepts `FrequencyData`.

`devs` is any iterable of estimator names from [`estimators`](@ref)
(symbols or strings). `m` is a [`TauMode`](@ref) (e.g. `Octave`, `Decade`)
or an integer vector of averaging factors applied to every deviation, as in
the estimator functions.
`correct_bias` reaches only the estimators that take it (the total family and
`theo1`), and `remove_drift` only the Hadamard family. Note that an explicit `m`
vector is passed to every deviation verbatim: `:theo1` interprets it as even
Thêo1 factors (odd entries give NaN rows) and `:theoh` as τ-grid factors, so
prefer a `TauMode` when mixing the Theo family with others.

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

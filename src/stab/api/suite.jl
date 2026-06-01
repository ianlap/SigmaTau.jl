# api/suite.jl — compute-all entry point.

"""
Default deviation set computed by [`stability`](@ref) when `devs` is not given —
a Stable32-like core: overlapping Allan, modified Allan, overlapping Hadamard,
and time deviation.
"""
const DEFAULT_DEVIATIONS = (:adev, :mdev, :hdev, :tdev)

# Per-symbol routing with keyword filtering. The Allan/Hadamard families forward
# only `calc_ci`/`confidence`; the total family also forwards any extra kwargs
# (`correct_bias`), letting each total function apply it; `pdev` forwards
# `calc_ci`/`confidence`; `mtie` forwards only `calc_ci`. The fallback rejects
# unknown symbols.
_dispatch_dev(::Val{:adev},  data, m; calc_ci, confidence, kwargs...) = adev(data, m;  calc_ci, confidence)
_dispatch_dev(::Val{:mdev},  data, m; calc_ci, confidence, kwargs...) = mdev(data, m;  calc_ci, confidence)
_dispatch_dev(::Val{:tdev},  data, m; calc_ci, confidence, kwargs...) = tdev(data, m;  calc_ci, confidence)
_dispatch_dev(::Val{:hdev},  data, m; calc_ci, confidence, kwargs...) = hdev(data, m;  calc_ci, confidence)
_dispatch_dev(::Val{:mhdev}, data, m; calc_ci, confidence, kwargs...) = mhdev(data, m; calc_ci, confidence)
_dispatch_dev(::Val{:htdev}, data, m; calc_ci, confidence, kwargs...) = htdev(data, m; calc_ci, confidence)

_dispatch_dev(::Val{:totdev},   data, m; calc_ci, confidence, kwargs...) = totdev(data, m;   calc_ci, confidence, kwargs...)
_dispatch_dev(::Val{:mtotdev},  data, m; calc_ci, confidence, kwargs...) = mtotdev(data, m;  calc_ci, confidence, kwargs...)
_dispatch_dev(::Val{:ttotdev},  data, m; calc_ci, confidence, kwargs...) = ttotdev(data, m;  calc_ci, confidence, kwargs...)
_dispatch_dev(::Val{:htotdev},  data, m; calc_ci, confidence, kwargs...) = htotdev(data, m;  calc_ci, confidence, kwargs...)
_dispatch_dev(::Val{:mhtotdev}, data, m; calc_ci, confidence, kwargs...) = mhtotdev(data, m; calc_ci, confidence, kwargs...)

_dispatch_dev(::Val{:mtie}, data, m; calc_ci, confidence, kwargs...) = mtie(data, m; calc_ci)
_dispatch_dev(::Val{:pdev}, data, m; calc_ci, confidence, kwargs...) = pdev(data, m; calc_ci, confidence)

_dispatch_dev(::Val{S}, data, m; kwargs...) where {S} =
    throw(ArgumentError("stability: unknown deviation :$S"))

function _stability(data, N::Int, kind::Symbol, tau0::Float64;
                    devs = DEFAULT_DEVIATIONS, taus = Octave,
                    calc_ci::Bool = true, confidence::Float64 = DEFAULT_CONFIDENCE,
                    kwargs...)
    results = StabilityResult[]
    for sym in devs
        m = taus isa TauMode ? tau_values(taus, N, sym) : collect(Int, taus)
        push!(results, _dispatch_dev(Val(sym), data, m; calc_ci, confidence, kwargs...))
    end
    tau_mode = taus isa TauMode ? Symbol(string(taus)) : :explicit
    return StabilitySuite(results, kind, tau0, N, calc_ci ? confidence : NaN, tau_mode)
end

"""
    stability(data::PhaseData; devs=DEFAULT_DEVIATIONS, taus=Octave,
              calc_ci=true, confidence=DEFAULT_CONFIDENCE, kwargs...) → StabilitySuite

Compute a suite of stability deviations from one record in a single call — the
batch analog of `adev`, `mdev`, …. Each symbol in `devs` is dispatched to its
public deviation function with a per-deviation averaging-factor grid derived from
`taus`, so every deviation gets a grid valid for its own algorithmic m-max.
Returns a [`StabilitySuite`](@ref). Also accepts `FrequencyData`.

`devs` is any iterable of deviation symbols drawn from `(:adev, :mdev, :tdev,
:hdev, :mhdev, :htdev, :totdev, :mtotdev, :ttotdev, :htotdev, :mhtotdev, :mtie,
:pdev)`. `taus` is a [`TauMode`](@ref) (e.g. `Octave`, `Decade`) or an explicit
`Vector{Int}` of averaging factors applied to every deviation. The extra keyword
argument `correct_bias` is forwarded only to the deviations that accept it (the
total family); `pdev` receives `calc_ci`/`confidence`; `mtie` receives only
`calc_ci`.

# Examples

```jldoctest
julia> using SigmaTau

julia> pd = PhaseData(collect(1.0:64.0), 1.0);

julia> suite = stability(pd; devs=(:adev, :hdev), taus=Octave, calc_ci=false);

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
    _stability(data, length(data.y), :frequency, data.tau0; kwargs...)

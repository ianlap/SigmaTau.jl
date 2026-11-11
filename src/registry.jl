# The estimator registry and the by-name entry points shared by `stability`,
# the explorer's compiled core and the command line: `estimators`, `deviation`,
# `tau_mode` and `summarize`.

const _ESTIMATORS = (
    (name=:adev,     label="ADEV",     family=:allan, axis=:sigma_y, correct_bias=false, description="Overlapping Allan deviation, the standard"),
    (name=:mdev,     label="MDEV",     family=:allan, axis=:sigma_y, correct_bias=false, description="Modified Allan deviation; separates white and flicker PM"),
    (name=:hdev,     label="HDEV",     family=:allan, axis=:sigma_y, correct_bias=false, description="Overlapping Hadamard deviation; ignores linear frequency drift"),
    (name=:mhdev,    label="MHDEV",    family=:allan, axis=:sigma_y, correct_bias=false, description="Modified Hadamard deviation"),
    (name=:pdev,     label="PDEV",     family=:allan, axis=:sigma_y, correct_bias=false, description="Parabolic deviation, for averaging (Ω) counters"),
    (name=:totdev,   label="TOTDEV",   family=:total, axis=:sigma_y, correct_bias=true,  description="Total Allan deviation; tighter intervals at long τ"),
    (name=:mtotdev,  label="MTOTDEV",  family=:total, axis=:sigma_y, correct_bias=true,  description="Total modified Allan deviation"),
    (name=:htotdev,  label="HTOTDEV",  family=:total, axis=:sigma_y, correct_bias=true,  description="Total Hadamard deviation"),
    (name=:mhtotdev, label="MHTOTDEV", family=:total, axis=:sigma_y, correct_bias=true,  description="Total modified Hadamard deviation"),
    (name=:theo1,    label="THEO1",    family=:theo,  axis=:sigma_y, correct_bias=true,  description="Thêo1, τ up to 3/4 of the record"),
    (name=:theoh,    label="THEOH",    family=:theo,  axis=:sigma_y, correct_bias=false, description="ThêoH: ADEV at short τ joined to bias-corrected Thêo1"),
    (name=:tdev,     label="TDEV",     family=:time,  axis=:sigma_x, correct_bias=false, description="Time deviation, σ_x in seconds"),
    (name=:htdev,    label="HTDEV",    family=:time,  axis=:sigma_x, correct_bias=false, description="Hadamard time deviation"),
    (name=:ttotdev,  label="TTOTDEV",  family=:time,  axis=:sigma_x, correct_bias=true,  description="Total time deviation"),
    (name=:mtie,     label="MTIE",     family=:time,  axis=:sigma_x, correct_bias=false, description="Maximum time interval error"),
    (name=:tierms,   label="TIE rms",  family=:time,  axis=:sigma_x, correct_bias=false, description="RMS time interval error"),
)

# ── Call forms every estimator shares ────────────────────────────────────
# deviations.jl defines each estimator on phase data and a `Vector{Int}` of
# averaging factors.
const _ESTIMATOR_NAMES = map(e -> e.name, _ESTIMATORS)

for f in (_ESTIMATOR_NAMES..., :dadev, :dhdev)
    # Frequency data is integrated to phase. The estimators that skip gaps
    # take a gapped frequency record as it is, because a missing frequency
    # sample loses the phase step across it (gaps.jl).
    if f === :adev
        # ADEV alone corrects for dead time between the samples (deadtime.jl).
        @eval adev(data::FrequencyData, m_values::Vector{Int}; dead_time::Real=0, kwargs...) =
            dead_time != 0 ? _dead_time_adev(data, m_values, Float64(dead_time); kwargs...) :
            _has_gaps(data) ? _gap_result(:adev, _gap_record(data), m_values; kwargs...) :
                              adev(_freq_to_phase(data), m_values; kwargs...)
    elseif f in (:mdev, :tdev, :hdev, :mhdev, :htdev, :mtie, :tierms)
        @eval $f(data::FrequencyData, m_values::Vector{Int}; kwargs...) = _has_gaps(data) ?
            _gap_result($(QuoteNode(f)), _gap_record(data), m_values; kwargs...) :
            $f(_freq_to_phase(data), m_values; kwargs...)
    else
        @eval $f(data::FrequencyData, m_values::Vector{Int}; kwargs...) =
            $f(_freq_to_phase(data), m_values; kwargs...)
    end
    # Any integer vector (`adev(pd, 1:10)`), forwarded as a `Vector{Int}` copy.
    @eval $f(data::AbstractTimingData, m_values::AbstractVector{<:Integer}; kwargs...) =
        $f(data, collect(Int, m_values); kwargs...)
end

# A `TauMode` grid, or none for octaves up to the estimator's largest m (the
# :theo1 grid is rounded to even averaging factors in `tau_values`).
for f in _ESTIMATOR_NAMES
    @eval $f(data::AbstractTimingData, taus::TauMode=Octave; kwargs...) =
        $f(data, tau_values(taus, _nphase(data), _grid_kernel($(QuoteNode(f)), kwargs)); kwargs...)
end

# Two non-overlapped ADEV (HDEV) terms span (d + 1)·m + 1 samples, the reach
# of MDEV (MHDEV), so the grid stops where theirs does.
_grid_kernel(f::Symbol, kwargs) =
    get(kwargs, :overlapping, true) ? f : f === :adev ? :mdev : f === :hdev ? :mhdev : f

"""
    estimators() → Vector{NamedTuple}

The estimators that [`deviation`](@ref) and [`stability`](@ref) run, one
entry each, with the fields

- `name::Symbol`: the function name, such as `:adev`.
- `label::String`: the display name, such as `"ADEV"`.
- `family::Symbol`: `:allan`, `:total`, `:theo` or `:time`.
- `axis::Symbol`: `:sigma_y` (fractional frequency) or `:sigma_x` (seconds).
- `correct_bias::Bool`: whether the `correct_bias` keyword applies.
- `description::String`: one line for menus and help.

The dynamic deviations [`dadev`](@ref) and [`dhdev`](@ref) return maps over
time and τ, not σ(τ) curves, and are not listed.

# Examples

```jldoctest
julia> using SigmaTau

julia> e = estimators();

julia> e[1].name, e[1].label, e[1].family, e[1].axis
(:adev, "ADEV", :allan, :sigma_y)

julia> [x.name for x in e if x.axis === :sigma_x]
5-element Vector{Symbol}:
 :tdev
 :htdev
 :ttotdev
 :mtie
 :tierms
```
"""
estimators() = collect(_ESTIMATORS)

"""
    deviation(name::Symbol, data, m=Octave; ci=true, confidence=DEFAULT_CONFIDENCE,
              correct_bias=true, remove_drift=true, noise=:auto,
              dead_time=0) → StabilityResult

Run the estimator called `name`, any name from [`estimators`](@ref), on a
`PhaseData` or `FrequencyData` record. `deviation(:adev, data, m)` is
`adev(data, m)`; use it when the estimator is chosen at run time, from a
menu or a command-line option.

`m` is a [`TauMode`](@ref) grid or an integer vector of averaging factors.
The keywords are those of [`adev`](@ref). `correct_bias` reaches only the
estimators that take it (the total estimators and `theo1`) and
`remove_drift` only the Hadamard family, and a `dead_time` other than zero
is an error for any estimator but `:adev`. `mtie` and `tierms` accept `ci`
and `confidence` and report no interval. The non-overlapped forms are
reached through [`adev`](@ref) and [`hdev`](@ref) directly.

The function is written as plain branches, so it compiles under
`juliac --trim`.

# Examples

```jldoctest
julia> using SigmaTau

julia> pd = PhaseData(collect(1.0:64.0) .^ 2, 1.0);

julia> deviation(:hdev, pd, 1:3; ci=false).dev == hdev(pd, 1:3; ci=false).dev
true
```
"""
deviation(name::Symbol, data::AbstractTimingData, m::TauMode = Octave; kwargs...) =
    deviation(name, data, _resolve_m(name, m, _nphase(data)); kwargs...)

deviation(name::Symbol, data::AbstractTimingData, m::AbstractVector{<:Integer}; kwargs...) =
    _deviation(name, data, collect(Int, m); kwargs...)

function _resolve_m(name::Symbol, m::TauMode, N::Int)
    name in (:dadev, :dhdev) && throw(_dynamic_error(name))
    name in _ESTIMATOR_NAMES || throw(_unknown_error(name))
    return tau_values(m, N, name)
end

_dynamic_error(name::Symbol) = ArgumentError(
    "$name is a dynamic deviation; call $name(data; window, step) instead")
_unknown_error(name::Symbol) = ArgumentError(
    "unknown deviation :$name; see estimators() for the names")

function _deviation(name::Symbol, data::AbstractTimingData, ms::Vector{Int};
                    ci::Bool = true, confidence::Real = DEFAULT_CONFIDENCE,
                    correct_bias::Union{Bool,Symbol} = true, remove_drift::Bool = true,
                    noise::_NoiseChoice = :auto, dead_time::Real = 0)
    c = Float64(confidence)
    b = correct_bias
    name === :adev     && return adev(data, ms; ci=ci, confidence=c, noise=noise, dead_time=Float64(dead_time))
    dead_time == 0 || throw(ArgumentError(
        "dead_time applies to adev alone: the bias functions are those of the two-sample variance"))
    name === :mdev     && return mdev(data, ms; ci=ci, confidence=c, noise=noise)
    name === :tdev     && return tdev(data, ms; ci=ci, confidence=c, noise=noise)
    name === :hdev     && return hdev(data, ms; ci=ci, confidence=c, remove_drift=remove_drift, noise=noise)
    name === :mhdev    && return mhdev(data, ms; ci=ci, confidence=c, remove_drift=remove_drift, noise=noise)
    name === :htdev    && return htdev(data, ms; ci=ci, confidence=c, remove_drift=remove_drift, noise=noise)
    name === :totdev   && return totdev(data, ms; ci=ci, confidence=c, correct_bias=b, noise=noise)
    name === :mtotdev  && return mtotdev(data, ms; ci=ci, confidence=c, correct_bias=b, noise=noise)
    name === :ttotdev  && return ttotdev(data, ms; ci=ci, confidence=c, correct_bias=b, noise=noise)
    name === :htotdev  && return htotdev(data, ms; ci=ci, confidence=c, correct_bias=b, remove_drift=remove_drift, noise=noise)
    name === :mhtotdev && return mhtotdev(data, ms; ci=ci, confidence=c, correct_bias=b, remove_drift=remove_drift, noise=noise)
    name === :theo1    && return theo1(data, ms; ci=ci, confidence=c, correct_bias=b, noise=noise)
    name === :theoh    && return theoh(data, ms; ci=ci, confidence=c, noise=noise)
    name === :pdev     && return pdev(data, ms; ci=ci, confidence=c, noise=noise)
    name === :mtie     && return mtie(data, ms; ci=ci, confidence=c)
    name === :tierms   && return tierms(data, ms; ci=ci, confidence=c)
    name in (:dadev, :dhdev) && throw(_dynamic_error(name))
    throw(_unknown_error(name))
end

"""
    tau_mode(s::AbstractString) → TauMode

The [`TauMode`](@ref) named by `s`, as typed on a command line or chosen in
a menu. Case and separators are ignored: `"half-octave"`, `"half_octave"`,
`"halfoctave"` and `"HalfOctave"` all give `HalfOctave`, and `"all"` or
`"alltaus"` gives `AllTaus`. Any other name throws an `ArgumentError`.

# Examples

```jldoctest
julia> using SigmaTau

julia> tau_mode("half-decade")
HalfDecade::TauMode = 5
```
"""
function tau_mode(s::AbstractString)
    k = replace(lowercase(strip(s)), "-" => "", "_" => "", " " => "")
    k == "octave"        && return Octave
    k == "halfoctave"    && return HalfOctave
    k == "quarteroctave" && return QuarterOctave
    k == "decade"        && return Decade
    k == "halfdecade"    && return HalfDecade
    (k == "all" || k == "alltaus") && return AllTaus
    throw(ArgumentError("unknown τ grid \"$s\"; use octave, half-octave, " *
                        "quarter-octave, decade, half-decade or all"))
end

"""
    summarize(data::PhaseData) → NamedTuple
    summarize(data::FrequencyData) → NamedTuple

A quick look at a record before analysis, as the explorer's Data page
shows it:

- `n`: number of samples.
- `tau0`: sample interval, in seconds.
- `span`: `n·tau0`, in seconds.
- `gaps`: number of `NaN` samples.
- `outliers`: the indices [`find_outliers`](@ref) returns at `nsigma = 5`.
- `drift_per_day`: least-squares slope of the fractional frequency against
  time, per day; `NaN` with fewer than three finite frequency samples.

# Examples

```jldoctest
julia> using SigmaTau

julia> s = summarize(PhaseData(0.5e-9 .* (0:99) .^ 2, 1.0));

julia> (s.n, s.gaps, round(s.drift_per_day; sigdigits=3))
(100, 0, 8.64e-5)
```
"""
function summarize(data::PhaseData)
    x = data.x
    n = length(x)
    y = Float64[(Float64(x[i + 1]) - Float64(x[i])) / data.tau0 for i in 1:n-1]
    return _summary(n, data.tau0, count(isnan, x), find_outliers(data; nsigma=5), y)
end

function summarize(data::FrequencyData)
    n = length(data.y)
    return _summary(n, data.tau0, count(isnan, data.y), find_outliers(data; nsigma=5),
                    Vector{Float64}(data.y))
end

function _summary(n::Int, tau0::Float64, gaps::Int, outliers::Vector{Int}, y::Vector{Float64})
    st = 0.0; sy = 0.0; sxy = 0.0; sxx = 0.0; k = 0
    for i in eachindex(y)
        isfinite(y[i]) || continue
        t = (i - 1) * tau0
        st += t; sy += y[i]; sxy += t * y[i]; sxx += t * t; k += 1
    end
    drift = k > 2 ? (k * sxy - st * sy) / (k * sxx - st * st) * 86400.0 : NaN
    return (n = n, tau0 = tau0, span = n * tau0, gaps = gaps,
            outliers = outliers, drift_per_day = drift)
end

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

"""
    estimators() → Vector{NamedTuple}

The estimators [`deviation`](@ref) and [`stability`](@ref) know, one entry
per estimator with fields

- `name::Symbol`: the function name, e.g. `:adev`;
- `label::String`: the display name, e.g. `"ADEV"`;
- `family::Symbol`: `:allan`, `:total`, `:theo` or `:time`;
- `axis::Symbol`: `:sigma_y` (fractional frequency) or `:sigma_x` (seconds);
- `correct_bias::Bool`: whether the `correct_bias` keyword applies;
- `description::String`: one line for menus and help.

The dynamic deviations `dadev` / `dhdev` return maps, not σ(τ) curves, and
are not listed.

# Examples

```jldoctest
julia> using SigmaTau

julia> e = estimators(); e[1]
(name = :adev, label = "ADEV", family = :allan, axis = :sigma_y, correct_bias = false, description = "Overlapping Allan deviation, the standard")

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
              correct_bias=true) → StabilityResult

Run the estimator called `name` (any `name` from [`estimators`](@ref)) on a
`PhaseData` or `FrequencyData` record. `m` is a [`TauMode`](@ref), resolved
against the estimator's own m-max, or any integer vector of averaging factors.
`correct_bias` is passed only to the estimators that take it (the total family
and `theo1`); `mtie` and `tierms` accept `ci` and `confidence` but report no
intervals.

`deviation(:adev, pd, m)` is `adev(pd, m)`; use it when the estimator is
chosen at run time (a menu, a command-line flag, a loop). It is written as
plain branches, so it compiles under `juliac --trim`.

# Examples

```jldoctest
julia> using SigmaTau

julia> pd = PhaseData(collect(1.0:64.0) .^ 2, 1.0);

julia> deviation(:hdev, pd, 1:3; ci=false).dev == hdev(pd, 1:3; ci=false).dev
true
```
"""
deviation(name::Symbol, data::AbstractTimingData, m::TauMode = Octave; kwargs...) =
    deviation(name, data, _resolve_m(name, m, length(data)); kwargs...)

deviation(name::Symbol, data::AbstractTimingData, m::AbstractVector{<:Integer}; kwargs...) =
    _deviation(name, data, collect(Int, m); kwargs...)

function _resolve_m(name::Symbol, m::TauMode, N::Int)
    name in (:dadev, :dhdev) && throw(_dynamic_error(name))
    any(e -> e.name === name, _ESTIMATORS) || throw(_unknown_error(name))
    return tau_values(m, N, name)
end

_dynamic_error(name::Symbol) = ArgumentError(
    "$name is a dynamic deviation; call $name(data; window, step) instead")
_unknown_error(name::Symbol) = ArgumentError(
    "unknown deviation :$name; see estimators() for the names")

function _deviation(name::Symbol, data::AbstractTimingData, ms::Vector{Int};
                    ci::Bool = true, confidence::Real = DEFAULT_CONFIDENCE,
                    correct_bias::Bool = true)
    c = Float64(confidence)
    b = correct_bias
    name === :adev     && return adev(data, ms; ci=ci, confidence=c)
    name === :mdev     && return mdev(data, ms; ci=ci, confidence=c)
    name === :tdev     && return tdev(data, ms; ci=ci, confidence=c)
    name === :hdev     && return hdev(data, ms; ci=ci, confidence=c)
    name === :mhdev    && return mhdev(data, ms; ci=ci, confidence=c)
    name === :htdev    && return htdev(data, ms; ci=ci, confidence=c)
    name === :totdev   && return totdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    name === :mtotdev  && return mtotdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    name === :ttotdev  && return ttotdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    name === :htotdev  && return htotdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    name === :mhtotdev && return mhtotdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    name === :theo1    && return theo1(data, ms; ci=ci, confidence=c, correct_bias=b)
    name === :theoh    && return theoh(data, ms; ci=ci, confidence=c)
    name === :pdev     && return pdev(data, ms; ci=ci, confidence=c)
    name === :mtie     && return mtie(data, ms; ci=ci, confidence=c)
    name === :tierms   && return tierms(data, ms; ci=ci, confidence=c)
    name in (:dadev, :dhdev) && throw(_dynamic_error(name))
    throw(_unknown_error(name))
end

"""
    tau_mode(s::AbstractString) → TauMode

Parse a τ-grid name, as typed on a command line or chosen in a menu, into a
[`TauMode`](@ref). Case and the separator are ignored, so `"half-octave"`,
`"half_octave"`, `"halfoctave"` and `"HalfOctave"` all give `HalfOctave`;
`"all"` and `"alltaus"` give `AllTaus`.

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

A quick check of a record before analysis, as the explorer's Data page shows
it:

- `n`: number of samples;
- `tau0`: sample interval in seconds;
- `span`: `n · tau0`, the record length in seconds;
- `gaps`: number of `NaN` samples;
- `outliers`: indices from [`find_outliers`](@ref) at `nsigma = 5`;
- `drift_per_day`: least-squares slope of the fractional frequency against
  time, per day (`NaN` with fewer than three finite frequency samples).

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

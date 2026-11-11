# Data records (PhaseData, FrequencyData) and result types.
"""
    AbstractTimingData

Supertype of the data records [`PhaseData`](@ref) and
[`FrequencyData`](@ref).
"""
abstract type AbstractTimingData end


"""
    PhaseData(x, tau0=1.0; source="user", start=NaN)

Phase samples ``x(t)`` in seconds at the uniform sample interval `tau0`: the
time deviation of a clock from its reference
([IEEE Std 1139-2022, 3.1](@cite ieee1139-2022-definitions)). Every
estimator takes this record, so the sample interval travels with the data.

`x` can be any real vector, such as a view, a range or integers, and is
copied to a `Vector` of floats. `NaN` marks a missing sample (see
[`fillgaps`](@ref)). A record holds at least two samples, and
`length(data)` returns their number.

[`read_phase`](@ref) builds a record from a file, and `PhaseData(f)`
converts a [`FrequencyData`](@ref) record to phase.

# Examples

```jldoctest
julia> using SigmaTau

julia> PhaseData([0.0, 1.5e-9, 2.9e-9, 4.6e-9], 10.0)
PhaseData{Float64}(4 samples, τ₀=10 s)
```

# Fields

$(TYPEDFIELDS)
"""
struct PhaseData{T<:AbstractFloat} <: AbstractTimingData
    "Phase samples, in seconds."
    x::Vector{T}
    "Sample interval τ₀, in seconds."
    tau0::Float64
    """Origin of the record: `"user"` when constructed directly, the file
    path when read by [`read_phase`](@ref). [`detrend`](@ref),
    [`fillgaps`](@ref) and [`remove_outliers`](@ref) keep it, and
    [`save`](@ref) writes it to the file header."""
    source::String
    """UTC time of the first sample, in seconds since 2000-01-01T00:00Z, or
    `NaN` when unknown. The readers set it from ISO 8601 or MJD timestamps
    and from a `# start:` header line; [`start_time`](@ref) returns it as a
    `DateTime`."""
    start::Float64
    function PhaseData(x::Vector{T}, tau0::Real = 1.0;
                       source::AbstractString = "user", start::Real = NaN) where {T<:AbstractFloat}
        tau0 > 0 || throw(ArgumentError("PhaseData: tau0 must be positive, got $tau0"))
        length(x) >= 2 || throw(ArgumentError(
            "PhaseData: need at least 2 phase samples, got $(length(x))"))
        return new{T}(x, Float64(tau0), String(source), Float64(start))
    end
end

# Any real vector (a view, a range, integers): copied into a `Vector` of its
# float type, so `PhaseData(view(x, 1:100), 1.0)` and `PhaseData(1:10)` work.
PhaseData(x::AbstractVector{<:Real}, tau0::Real = 1.0; source::AbstractString = "user", start::Real = NaN) =
    PhaseData(Vector{float(eltype(x))}(x), tau0; source=source, start=start)


"""
    FrequencyData(y, tau0=1.0; source="user", start=NaN)

Fractional-frequency samples ``y(t) = (ν(t) - ν_0)/ν_0``, dimensionless, at
the uniform sample interval `tau0`
([IEEE Std 1139-2022, 3.1](@cite ieee1139-2022-definitions)). Each sample
is taken as the average over its interval, with no dead time. Every
estimator accepts the record and integrates it to phase before the
analysis.

Like [`PhaseData`](@ref), the record takes any real vector, marks a missing
sample with `NaN`, holds at least two samples and supports `length`.

[`read_frequency`](@ref) builds a record from a file, and `FrequencyData(p)`
converts a [`PhaseData`](@ref) record to frequency. The two conversions are
described below.

# Examples

```jldoctest
julia> using SigmaTau

julia> FrequencyData([1.2e-11, -0.4e-11, 0.7e-11], 10.0)
FrequencyData{Float64}(3 samples, τ₀=10 s)
```

# Fields

$(TYPEDFIELDS)
"""
struct FrequencyData{T<:AbstractFloat} <: AbstractTimingData
    "Fractional-frequency samples, dimensionless."
    y::Vector{T}
    "Sample interval τ₀, in seconds."
    tau0::Float64
    """Origin of the record: `"user"` when constructed directly, the file
    path when read by [`read_frequency`](@ref). [`detrend`](@ref),
    [`fillgaps`](@ref) and [`remove_outliers`](@ref) keep it, and
    [`save`](@ref) writes it to the file header."""
    source::String
    """UTC time of the first sample, in seconds since 2000-01-01T00:00Z, or
    `NaN` when unknown. The readers set it from ISO 8601 or MJD timestamps
    and from a `# start:` header line; [`start_time`](@ref) returns it as a
    `DateTime`."""
    start::Float64
    function FrequencyData(y::Vector{T}, tau0::Real = 1.0;
                           source::AbstractString = "user", start::Real = NaN) where {T<:AbstractFloat}
        tau0 > 0 || throw(ArgumentError("FrequencyData: tau0 must be positive, got $tau0"))
        length(y) >= 2 || throw(ArgumentError(
            "FrequencyData: need at least 2 frequency samples, got $(length(y))"))
        return new{T}(y, Float64(tau0), String(source), Float64(start))
    end
end

FrequencyData(y::AbstractVector{<:Real}, tau0::Real = 1.0; source::AbstractString = "user", start::Real = NaN) =
    FrequencyData(Vector{float(eltype(y))}(y), tau0; source=source, start=start)

"""
    length(data::PhaseData) → Int
    length(data::FrequencyData) → Int

Number of samples in the record.
"""
Base.length(p::PhaseData) = length(p.x)
Base.length(f::FrequencyData) = length(f.y)

# Phase samples a record stands for: N frequency samples integrate to N + 1.
_nphase(p::PhaseData) = length(p.x)
_nphase(f::FrequencyData) = length(f.y) + 1

"""
    start_time(data) → Union{DateTime, Nothing}

UTC time of the first sample of a [`PhaseData`](@ref) or
[`FrequencyData`](@ref) record, to the millisecond, or `nothing` when the
record does not carry it. The readers take it from ISO 8601 or MJD
timestamps and from the `# start:` line that [`save`](@ref) writes. Sample
`i` is at `start_time(data) + (i − 1)·τ₀`.
"""
start_time(d::AbstractTimingData) =
    isnan(d.start) ? nothing : Dates.DateTime(2000) + Dates.Millisecond(round(Int64, d.start * 1000))

# ISO 8601 UTC text of seconds since 2000-01-01, to the nanosecond. By hand:
# Dates.format and padding with a Char are not resolved under juliac --trim.
function _pad(v::Int64, n::Int)
    t = string(v)
    while length(t) < n
        t = "0" * t
    end
    return t
end

function _iso_text(sec::Float64)
    whole = floor(Int64, sec)
    ns = round(Int64, (sec - whole) * 1e9)
    ns == 1_000_000_000 && (whole += 1; ns = 0)
    day, s = fldmod(whole, 86400)
    date = Dates.Date(2000, 1, 1) + Dates.Day(day)
    y = Int64(Dates.year(date)); mo = Int64(Dates.month(date)); d = Int64(Dates.day(date))
    t = _pad(y, 4) * "-" * _pad(mo, 2) * "-" * _pad(d, 2) * "T" * _pad(s ÷ 3600, 2) * ":" *
        _pad(s % 3600 ÷ 60, 2) * ":" * _pad(s % 60, 2)
    ns == 0 && return t * "Z"
    return t * "." * rstrip(_pad(ns, 9), '0') * "Z"
end

# One-line displays of the data records: sample count, τ₀, and the origin
# when it is not the "user" default. Printing every sample would swamp the
# REPL at realistic record lengths.
function Base.show(io::IO, p::PhaseData)
    print(io, "PhaseData{", eltype(p.x), "}(", length(p.x), " samples, τ₀=",
              _show4g(p.tau0), " s")
    p.source == "user" || print(io, ", source=", repr(p.source))
    isnan(p.start) || print(io, ", start=", _iso_text(p.start))
    print(io, ")")
end

function Base.show(io::IO, f::FrequencyData)
    print(io, "FrequencyData{", eltype(f.y), "}(", length(f.y), " samples, τ₀=",
              _show4g(f.tau0), " s")
    f.source == "user" || print(io, ", source=", repr(f.source))
    isnan(f.start) || print(io, ", start=", _iso_text(f.start))
    print(io, ")")
end


"""
    StabilityResult

Result of every stability estimator: one deviation per averaging time, with
the noise type, the equivalent degrees of freedom and the confidence
interval when the estimator was called with `ci=true`.

`ci` and `edf` are empty for a call with `ci=false`, and always for `mtie`
and `tierms`. `noise_type` is empty too, unless a bias correction needed
the noise type. `neff` is always filled. The interval at point `i` runs
from `r.ci[i].lo` to `r.ci[i].hi`; [`ci_lower`](@ref) and
[`ci_upper`](@ref) return the limits as vectors. Displaying a result prints
it as a table.

# Fields

$(TYPEDFIELDS)
"""
struct StabilityResult
    "Estimator that produced the result, such as `:adev`."
    deviation_type::Symbol
    "Averaging times τ, in seconds."
    tau::Vector{Float64}
    "Deviation at each τ: dimensionless for the σ_y estimators, in seconds for `tdev`, `htdev`, `ttotdev`, `mtie` and `tierms`."
    dev::Vector{Float64}
    "Noise type at each τ, identified or stated: `:WHPM`, `:FLPM`, `:WHFM`, `:FLFM`, `:RWFM`, `:FWFM` or `:RRFM`."
    noise_type::Vector{Symbol}
    "Confidence interval `(lo, hi)` at each τ."
    ci::Vector{@NamedTuple{lo::Float64, hi::Float64}}
    "Equivalent degrees of freedom at each τ."
    edf::Vector{Float64}
    "Number of terms averaged at each τ; 0 where the deviation is `NaN`."
    neff::Vector{Int}
    "Sample interval τ₀ of the record, in seconds; `NaN` if unknown."
    tau0::Float64
    "Number of phase samples in the record; 0 if unknown."
    n::Int
    "Confidence level of `ci`; `NaN` when there is no interval."
    confidence::Float64
end

# Without the record metadata (files written by v0.4, results made by hand).
StabilityResult(deviation_type::Symbol, tau, dev, noise_type, ci, edf, neff) =
    StabilityResult(deviation_type, tau, dev, noise_type, ci, edf, neff, NaN, 0, NaN)

# Element type of `StabilityResult.ci`.
const _CIBound = @NamedTuple{lo::Float64, hi::Float64}

# Lower and upper limits, as `confidence_intervals` returns them, in the form
# stored in `ci`.
_zip_ci(lower::Vector{Float64}, upper::Vector{Float64}) =
    _CIBound[(lo = lo, hi = hi) for (lo, hi) in zip(lower, upper)]

# Each limit times a factor per τ, for the time deviations.
_scale_ci(ci::Vector{_CIBound}, factor::Vector{Float64}) =
    _CIBound[(lo = c.lo * f, hi = c.hi * f) for (c, f) in zip(ci, factor)]

"""
    ci_lower(r::StabilityResult) → Vector{Float64}

Lower confidence limits of `r`, one per averaging time; empty when the
result carries no interval.
"""
ci_lower(r::StabilityResult) = getfield.(r.ci, :lo)

"""
    ci_upper(r::StabilityResult) → Vector{Float64}

Upper confidence limits of `r`, one per averaging time; empty when the
result carries no interval.
"""
ci_upper(r::StabilityResult) = getfield.(r.ci, :hi)

# The show methods print 4 significant digits (%.4g). The stored vectors keep
# full precision.
_show4g(x::Real) = @sprintf("%.4g", x)

function Base.show(io::IO, r::StabilityResult)
    n = length(r.tau)
    ci = isempty(r.edf) ? "no CI" : "with CI"
    if n == 0
        print(io, "StabilityResult(:", r.deviation_type, ", 0 pts, ", ci, ")")
    else
        print(io, "StabilityResult(:", r.deviation_type, ", ", n, " pts, ",
                  "τ∈[", _show4g(r.tau[1]), ", ", _show4g(r.tau[end]), "] s, ", ci, ")")
    end
end

# REPL display: the one-line summary, then a table with one row per τ. The
# interval columns appear when the result carries them. A result longer than
# 12 rows shows 11 and an ellipsis line.
function Base.show(io::IO, ::MIME"text/plain", r::StabilityResult)
    show(io, r)
    n = length(r.tau)
    n == 0 && return
    has_ci = !isempty(r.edf)
    maxrows = 12
    nshow = n > maxrows ? maxrows - 1 : n
    header = has_ci ? ["tau", "dev", "neff", "ci_lo", "ci_hi", "edf", "noise"] :
                      ["tau", "dev", "neff"]
    rows = map(1:nshow) do i
        has_ci ?
            [_show4g(r.tau[i]), _show4g(r.dev[i]), string(r.neff[i]),
             _show4g(r.ci[i].lo), _show4g(r.ci[i].hi), _show4g(r.edf[i]),
             isempty(r.noise_type) ? "" : string(r.noise_type[i])] :
            [_show4g(r.tau[i]), _show4g(r.dev[i]), string(r.neff[i])]
    end
    widths = map(length, header)
    for row in rows, j in eachindex(widths)
        widths[j] = max(widths[j], length(row[j]))
    end
    # Numeric columns right-align; the trailing noise column prints as is.
    pad(s, j) = (has_ci && j == length(header)) ? s : lpad(s, widths[j])
    print(io, "\n  ", join((pad(header[j], j) for j in eachindex(header)), "  "))
    for row in rows
        print(io, "\n  ", join((pad(row[j], j) for j in eachindex(header)), "  "))
    end
    n > maxrows && print(io, "\n  … ", n, " rows total")
    return
end


"""
    StabilitySuite

The results of one [`stability`](@ref) call: a [`StabilityResult`](@ref)
per estimator, in the order requested, with the record and grid settings
that produced them.

`suite[:adev]` returns a result by estimator name and `suite[1]` by
position. The suite iterates in order and supports `keys`, `haskey`,
`length` and `pairs`. [`save`](@ref) writes it to one file.

# Fields

$(TYPEDFIELDS)
"""
struct StabilitySuite
    "One result per estimator, in the order requested."
    results::Vector{StabilityResult}
    "Kind of the input record: `:phase` or `:frequency`."
    data_kind::Symbol
    "Sample interval τ₀, in seconds."
    tau0::Float64
    "Number of phase samples in the record."
    n::Int
    "Confidence level of the intervals; `NaN` for `ci=false`."
    confidence::Float64
    "Grid the suite was run on: the `TauMode` name, such as `:Octave`, or `:explicit`."
    tau_mode::Symbol
end

Base.length(s::StabilitySuite) = length(s.results)
Base.iterate(s::StabilitySuite, state...) = iterate(s.results, state...)
Base.eltype(::Type{StabilitySuite}) = StabilityResult
Base.keys(s::StabilitySuite) = Symbol[r.deviation_type for r in s.results]
Base.haskey(s::StabilitySuite, k::Symbol) = any(r -> r.deviation_type === k, s.results)
Base.pairs(s::StabilitySuite) = (r.deviation_type => r for r in s.results)

Base.getindex(s::StabilitySuite, i::Integer) = s.results[i]

function Base.getindex(s::StabilitySuite, k::Symbol)
    idx = findfirst(r -> r.deviation_type === k, s.results)
    idx === nothing && throw(KeyError(k))
    return s.results[idx]
end

function Base.show(io::IO, s::StabilitySuite)
    devs = join(("$d" for d in keys(s)), ", ")
    ci = isnan(s.confidence) ? "no CI" : "CI@$(s.confidence)"
    print(io, "StabilitySuite(", length(s), " devs [", devs, "], ",
              s.data_kind, ", N=", s.n, ", τ₀=", _show4g(s.tau0), " s, ", ci, ")")
end

function Base.show(io::IO, ::MIME"text/plain", s::StabilitySuite)
    ci = isnan(s.confidence) ? "no CI" : "CI @ $(s.confidence)"
    print(io, "StabilitySuite: ", length(s), " deviation(s), ", s.data_kind,
              " data, N=", s.n, ", τ₀=", _show4g(s.tau0), " s, m=", s.tau_mode, ", ", ci)
    for r in s.results
        n = length(r.tau)
        rng = n == 0 ? "" : "  τ∈[$(_show4g(r.tau[1])), $(_show4g(r.tau[end]))] s"
        print(io, "\n  ", rpad(string(r.deviation_type), 9), lpad(n, 3), " pts", rng)
    end
end


"""
    SpectralResult

Result of the spectral-density estimates [`psd_y`](@ref), [`psd_x`](@ref)
and [`ssb_phase_noise`](@ref): the one-sided frequency grid, the estimated
density and the Welch settings that produced it.

# Fields

$(TYPEDFIELDS)
"""
struct SpectralResult
    "Quantity estimated: `:Sy` for ``S_y(f)``, `:Sx` for ``S_x(f)`` or `:L` for ``ℒ(f)``."
    spectral_type::Symbol
    "Fourier frequencies ``f``, in Hz."
    freq::Vector{Float64}
    "Spectral density at each ``f``, in the unit named by `units`."
    psd::Vector{Float64}
    "Unit of `psd`: `:per_Hz` (``S_y``, 1/Hz), `:s2_per_Hz` (``S_x``, s²/Hz) or `:dBc_per_Hz` (``ℒ``)."
    units::Symbol
    "Welch segment length, in samples."
    nperseg::Int
    "Overlap of successive segments, in samples."
    noverlap::Int
    "Taper applied to each segment: `:hann`, `:hamming` or `:rectangular`."
    taper::Symbol
end

function Base.show(io::IO, r::SpectralResult)
    n = length(r.freq)
    rng = n == 0 ? "" : ", f∈[$(r.freq[1]), $(r.freq[end])] Hz"
    print(io, "SpectralResult(:", r.spectral_type, ", ", n, " bins", rng,
              ", ", r.units, ", nperseg=", r.nperseg, ")")
end


"""
    DynamicStabilityResult

Result of the dynamic deviations [`dadev`](@ref) and [`dhdev`](@ref): a map
``σ_y(t, τ)`` of the deviation computed inside a window of `window` phase
samples that slides along the record.

Row `i` of `dev` is the window centered at `t[i]`, and column `j` the
averaging time `tau[j]`. An entry is `NaN` where the window is too short
for the averaging factor.

The result carries no noise type, EDF or interval.
[Galleani (2011)](@cite galleani-2011-davar-confidence) derives confidence
surfaces for the dynamic Allan variance; they are not implemented here.

# Fields

$(TYPEDFIELDS)
"""
struct DynamicStabilityResult
    "Estimator that produced the result: `:dadev` or `:dhdev`."
    deviation_type::Symbol
    "Center time of each window, in seconds from the first sample; one per row of `dev`."
    t::Vector{Float64}
    "Averaging times τ, in seconds; one per column of `dev`."
    tau::Vector{Float64}
    "Deviation map, `length(t)` by `length(tau)`."
    dev::Matrix{Float64}
    "Window length, in samples."
    window::Int
    "Sample interval τ₀, in seconds."
    tau0::Float64
end

function Base.show(io::IO, r::DynamicStabilityResult)
    nt, nτ = size(r.dev)
    print(io, "DynamicStabilityResult(:", r.deviation_type, ", ", nt, "×", nτ, " map")
    nt > 0 && print(io, ", t∈[", _show4g(r.t[1]), ", ", _show4g(r.t[end]), "] s")
    nτ > 0 && print(io, ", τ∈[", _show4g(r.tau[1]), ", ", _show4g(r.tau[end]), "] s")
    print(io, ", window=", r.window, ")")
end

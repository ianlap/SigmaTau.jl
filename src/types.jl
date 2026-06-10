# types.jl — shared data records and result types.
#
# Merged from types/{abstract, phase_data, frequency_data, stability_result,
# stability_suite, spectral_result}.jl. Order matters: PhaseData/FrequencyData
# subtype AbstractTimingData, and the show methods follow their structs.


# ──────────────────────────────────────────────────────────────────────
# ── types/abstract.jl ───────────────────────────────────────────────────

"""
    AbstractTimingData

Supertype of timing data records (`PhaseData`, `FrequencyData`).
"""
abstract type AbstractTimingData end


# ──────────────────────────────────────────────────────────────────────
# ── types/phase_data.jl ─────────────────────────────────────────────────

"""
    PhaseData{T<:AbstractFloat}

Phase residuals `x(t)` sampled at uniform interval `tau0` (default `1.0` s).

$(TYPEDFIELDS)
"""
struct PhaseData{T<:AbstractFloat} <: AbstractTimingData
    "Phase samples in seconds."
    x::Vector{T}
    "Base sample interval τ₀ in seconds."
    tau0::Float64
    function PhaseData(x::Vector{T}, tau0::Real = 1.0) where {T<:AbstractFloat}
        tau0 > 0 || throw(ArgumentError("PhaseData: tau0 must be positive, got $tau0"))
        length(x) >= 2 || throw(ArgumentError(
            "PhaseData: need at least 2 phase samples, got $(length(x))"))
        return new{T}(x, Float64(tau0))
    end
end


# ──────────────────────────────────────────────────────────────────────
# ── types/frequency_data.jl ─────────────────────────────────────────────

"""
    FrequencyData{T<:AbstractFloat}

Fractional-frequency samples `y(t)` at uniform interval `tau0` (default `1.0` s).

$(TYPEDFIELDS)
"""
struct FrequencyData{T<:AbstractFloat} <: AbstractTimingData
    "Fractional-frequency samples (dimensionless)."
    y::Vector{T}
    "Base sample interval τ₀ in seconds."
    tau0::Float64
    function FrequencyData(y::Vector{T}, tau0::Real = 1.0) where {T<:AbstractFloat}
        tau0 > 0 || throw(ArgumentError("FrequencyData: tau0 must be positive, got $tau0"))
        length(y) >= 2 || throw(ArgumentError(
            "FrequencyData: need at least 2 frequency samples, got $(length(y))"))
        return new{T}(y, Float64(tau0))
    end
end


# ──────────────────────────────────────────────────────────────────────
# ── types/stability_result.jl ───────────────────────────────────────────

"""
    StabilityResult

Unified return type for every stability calculation.

The `noise_type`, `ci`, and `edf` vectors are empty when the calculation was
invoked with `ci=false` (and always for `mtie`); `neff` is always populated.
Per-τ confidence bounds read as `r.ci[i].lo` / `r.ci[i].hi`; the
[`ci_lower`](@ref) / [`ci_upper`](@ref) accessors return them as plain
vectors.

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
    "Confidence interval `(lo, hi)` per τ (empty unless `ci=true`)."
    ci::Vector{@NamedTuple{lo::Float64, hi::Float64}}
    "Equivalent degrees of freedom (empty unless `ci=true`)."
    edf::Vector{Float64}
    "Number of analysis windows averaged at each τ (0 where the kernel returns NaN)."
    neff::Vector{Int}
end

# Element type of `StabilityResult.ci`.
const _CIBound = @NamedTuple{lo::Float64, hi::Float64}

# Zip separate lower/upper bound vectors (as returned by
# `confidence_intervals`) into the tuple form stored in `ci`.
_zip_ci(lower::Vector{Float64}, upper::Vector{Float64}) =
    _CIBound[(lo = lo, hi = hi) for (lo, hi) in zip(lower, upper)]

# Multiply each CI bound by a per-τ factor (tdev/htdev/ttotdev rescaling).
_scale_ci(ci::Vector{_CIBound}, factor::Vector{Float64}) =
    _CIBound[(lo = c.lo * f, hi = c.hi * f) for (c, f) in zip(ci, factor)]

"""
    ci_lower(r::StabilityResult) → Vector{Float64}

Lower confidence bounds of `r` as a plain vector (`getfield.(r.ci, :lo)`).
Empty when the result carries no CIs.
"""
ci_lower(r::StabilityResult) = getfield.(r.ci, :lo)

"""
    ci_upper(r::StabilityResult) → Vector{Float64}

Upper confidence bounds of `r` as a plain vector (`getfield.(r.ci, :hi)`).
Empty when the result carries no CIs.
"""
ci_upper(r::StabilityResult) = getfield.(r.ci, :hi)

# Display rounding for the show methods below: 4 significant digits via %.4g.
# The stored vectors keep full precision — this is presentation only.
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

# Full REPL display: the one-line summary followed by an aligned per-τ table.
# CI columns appear only when the result carries them. Results longer than
# 12 rows truncate to 11 plus an ellipsis line.
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
             string(r.noise_type[i])] :
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


# ──────────────────────────────────────────────────────────────────────
# ── types/stability_suite.jl ────────────────────────────────────────────

"""
    StabilitySuite

Ordered collection of [`StabilityResult`](@ref)s computed from one timing record
in a single [`stability`](@ref) call, plus the session metadata needed to
reproduce and serialize the analysis.

Index by deviation symbol to recover a result (`suite[:adev]`) or by position
(`suite[1]`). The suite is iterable in request order and supports `keys`,
`haskey`, `length`, and `pairs`, so a CLI or GUI can render it as an ordered
table.

$(TYPEDFIELDS)
"""
struct StabilitySuite
    "Per-deviation results, in the order they were requested."
    results::Vector{StabilityResult}
    "Input record kind: `:phase` or `:frequency`."
    data_kind::Symbol
    "Base sample interval τ₀ in seconds."
    tau0::Float64
    "Number of input samples N."
    n::Int
    "Confidence level used for the CI bounds (`NaN` when `ci=false`)."
    confidence::Float64
    "Tau-grid selector, recorded as a Symbol (e.g. `:Octave`, or `:explicit`)."
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
              " data, N=", s.n, ", τ₀=", _show4g(s.tau0), " s, taus=", s.tau_mode, ", ", ci)
    for r in s.results
        n = length(r.tau)
        rng = n == 0 ? "" : "  τ∈[$(_show4g(r.tau[1])), $(_show4g(r.tau[end]))] s"
        print(io, "\n  ", rpad(string(r.deviation_type), 9), lpad(n, 3), " pts", rng)
    end
end


# ──────────────────────────────────────────────────────────────────────
# ── types/spectral_result.jl ────────────────────────────────────────────

"""
    SpectralResult

Unified return type for every spectral-density estimate (`Sy`, `Sx`, `L`).

Mirrors [`StabilityResult`](@ref): a flat, non-parametric record carrying the
one-sided frequency grid, the estimated spectrum, and the Welch parameters that
produced it, so a saved or plotted result is self-describing.

$(TYPEDFIELDS)
"""
struct SpectralResult
    "Which estimate produced this result: `:Sy`, `:Sx`, or `:L`."
    spectral_type::Symbol
    "One-sided frequency grid `f` in Hz."
    freq::Vector{Float64}
    "Spectral density at each `f`; interpret per `units`."
    psd::Vector{Float64}
    "Units of `psd`: `:per_Hz` (S_y, 1/Hz), `:s2_per_Hz` (S_x, s²/Hz), or `:dBc_per_Hz` (ℒ(f))."
    units::Symbol
    "Welch segment length (samples)."
    nperseg::Int
    "Welch segment overlap (samples)."
    noverlap::Int
    "Window applied to each segment: `:hann`, `:hamming`, or `:rectangular`."
    window::Symbol
end

function Base.show(io::IO, r::SpectralResult)
    n = length(r.freq)
    rng = n == 0 ? "" : ", f∈[$(r.freq[1]), $(r.freq[end])] Hz"
    print(io, "SpectralResult(:", r.spectral_type, ", ", n, " bins", rng,
              ", ", r.units, ", nperseg=", r.nperseg, ")")
end

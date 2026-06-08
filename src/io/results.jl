# io/results.jl — Save/load StabilityResult and StabilitySuite as
# tab-separated text (stdlib only).

const _IO_VERSION = "2"

# Shared row writer: the `# calc_ci=` flag, the column header, and one tab-
# delimited row per τ. Used by both `save_result` (v1) and `save_suite` (v2).
# The on-disk token stays `# calc_ci=` (not `# ci=`) on purpose: it is the v1/v2
# file-format key, independent of the public `ci` kwarg, so result files written
# by earlier versions keep loading. Do not rename it without bumping _IO_VERSION.
function _write_result_table(io::IO, r::StabilityResult)
    has_ci = !isempty(r.ci_lower)
    println(io, "# calc_ci=", has_ci)
    println(io, "tau\tdev\tnoise_type\tci_lower\tci_upper\tedf")
    for k in 1:length(r.tau)
        noise_s = has_ci ? string(r.noise_type[k]) : ""
        ci_lo_s = has_ci ? string(r.ci_lower[k])   : "NaN"
        ci_hi_s = has_ci ? string(r.ci_upper[k])   : "NaN"
        edf_s   = has_ci ? string(r.edf[k])        : "NaN"
        println(io, string(r.tau[k]), '\t', string(r.dev[k]), '\t',
                    noise_s, '\t', ci_lo_s, '\t', ci_hi_s, '\t', edf_s)
    end
end

# Shared row parser: turn tab-delimited data rows back into a StabilityResult.
function _parse_result_rows(deviation_type::Symbol, has_ci::Bool,
                            data_lines::AbstractVector{<:AbstractString})
    N        = length(data_lines)
    tau      = Vector{Float64}(undef, N)
    dev      = Vector{Float64}(undef, N)
    noise_s  = Vector{String}(undef, N)
    ci_lower = Vector{Float64}(undef, N)
    ci_upper = Vector{Float64}(undef, N)
    edf_vec  = Vector{Float64}(undef, N)

    for (k, line) in enumerate(data_lines)
        parts = split(line, '\t')
        length(parts) == 6 || error("malformed row $k (expected 6 columns, got $(length(parts)))")
        tau[k]      = parse(Float64, parts[1])
        dev[k]      = parse(Float64, parts[2])
        noise_s[k]  = parts[3]
        ci_lower[k] = parse(Float64, parts[4])
        ci_upper[k] = parse(Float64, parts[5])
        edf_vec[k]  = parse(Float64, parts[6])
    end

    has_ci || return StabilityResult(deviation_type, tau, dev,
                                     Symbol[], Float64[], Float64[], Float64[])
    return StabilityResult(deviation_type, tau, dev,
                           Symbol.(noise_s), ci_lower, ci_upper, edf_vec)
end

"""
    save_result(path::AbstractString, r::StabilityResult) → path

Write a single `StabilityResult` to a tab-delimited text file at `path`.

The file is self-describing: comment lines encode the format version, deviation
type, and whether confidence intervals were computed (`calc_ci`), followed by a
column-name header and one data row per τ value. All six data columns are always
present; `noise_type`, `ci_lower`, `ci_upper`, and `edf` are written as `""` /
`NaN` when `calc_ci=false`.

Returns `path` so the call can be chained. For a whole [`stability`](@ref)
session (multiple deviations + metadata) use [`save_suite`](@ref).

# Examples

```julia
pd = PhaseData(randn(256), 1.0)
r  = adev(pd, [1, 2, 4, 8, 16])
save_result("my_adev.tsv", r)
r2 = load_result("my_adev.tsv")
r2.dev ≈ r.dev   # true
```

See also [`load_result`](@ref), [`save_suite`](@ref).
"""
function save_result(path::AbstractString, r::StabilityResult)
    open(path, "w") do io
        println(io, "# SigmaTau StabilityResult v1")
        println(io, "# deviation_type=", r.deviation_type)
        _write_result_table(io, r)
    end
    return path
end

"""
    load_result(path::AbstractString) → StabilityResult

Read a single `StabilityResult` previously written by [`save_result`](@ref).

The `# calc_ci=` header line determines whether confidence-interval fields are
reconstructed. When `calc_ci=false`, `noise_type`, `ci_lower`, `ci_upper`, and
`edf` are returned as empty vectors, matching the standard deviation API
contract. Files written by [`save_suite`](@ref) are rejected — use
[`load_suite`](@ref) for those.

See also [`save_result`](@ref).
"""
function load_result(path::AbstractString)
    lines = readlines(path)
    (!isempty(lines) && occursin("StabilitySuite", lines[1])) &&
        error("load_result: \"$path\" is a StabilitySuite file — use `load_suite` instead.")

    deviation_type = :unknown
    has_ci         = true   # default: assume CI present if header is absent (v1 compat)

    for line in lines
        if startswith(line, "# deviation_type=")
            deviation_type = Symbol(line[length("# deviation_type=") + 1 : end])
        elseif startswith(line, "# calc_ci=")
            has_ci = strip(line[length("# calc_ci=") + 1 : end]) == "true"
        end
    end

    data_lines = filter(l -> !startswith(l, '#') && !startswith(l, "tau"), lines)
    isempty(data_lines) && error("load_result: no data rows found in \"$path\"")

    return _parse_result_rows(deviation_type, has_ci, data_lines)
end

"""
    save_suite(path::AbstractString, suite::StabilitySuite; source_file="") → path

Write a whole analysis session (multiple deviations + session metadata) to a
self-describing tab-delimited text file at `path` (format v2). Round-trips via
[`load_suite`](@ref).

A key=value metadata header records the package version, an ISO-8601 timestamp,
the (optional) `source_file`, the input kind/τ₀/N, the confidence level, the tau
mode, and the deviation set; each result then follows as a `# --- result <sym>
---` block reusing the same six-column table as [`save_result`](@ref). Returns
`path`.
"""
function save_suite(path::AbstractString, suite::StabilitySuite; source_file::AbstractString="")
    open(path, "w") do io
        println(io, "# SigmaTau StabilitySuite v", _IO_VERSION)
        println(io, "# package_version=", pkgversion(@__MODULE__))
        println(io, "# timestamp=", Dates.now())
        println(io, "# source_file=", source_file)
        println(io, "# data_kind=", suite.data_kind)
        println(io, "# tau0=", suite.tau0)
        println(io, "# n=", suite.n)
        println(io, "# confidence=", suite.confidence)
        println(io, "# tau_mode=", suite.tau_mode)
        println(io, "# deviations=", join((string(r.deviation_type) for r in suite.results), ","))
        for r in suite.results
            println(io, "# --- result ", r.deviation_type, " ---")
            _write_result_table(io, r)
        end
    end
    return path
end

"""
    load_suite(path::AbstractString) → StabilitySuite

Read a `StabilitySuite` previously written by [`save_suite`](@ref) (format v2),
reconstructing every result block and the session metadata (`data_kind`, `tau0`,
`n`, `confidence`, `tau_mode`). Files written by [`save_result`](@ref) are
rejected — use [`load_result`](@ref) for those.
"""
function load_suite(path::AbstractString)
    lines = readlines(path)
    (!isempty(lines) && occursin("StabilitySuite", lines[1])) ||
        error("load_suite: \"$path\" is not a StabilitySuite file — use `load_result` for single results.")

    meta       = Dict{String,String}()
    delim_idxs = Int[]
    for (i, ln) in enumerate(lines)
        if startswith(ln, "# --- result ")
            push!(delim_idxs, i)
        elseif isempty(delim_idxs)
            m = match(r"^# (\w+)=(.*)$", ln)
            m === nothing || (meta[m.captures[1]] = String(strip(m.captures[2])))
        end
    end
    isempty(delim_idxs) && error("load_suite: no result blocks found in \"$path\"")

    results = StabilityResult[]
    for (b, start) in enumerate(delim_idxs)
        stop     = b < length(delim_idxs) ? delim_idxs[b + 1] - 1 : length(lines)
        block    = @view lines[start:stop]
        dev_type = Symbol(match(r"^# --- result (\S+) ---", block[1]).captures[1])
        ci_lines = filter(l -> startswith(l, "# calc_ci="), block)
        has_ci   = !isempty(ci_lines) && strip(split(ci_lines[1], "=")[2]) == "true"
        data     = filter(l -> !startswith(l, '#') && !startswith(l, "tau") && !isempty(strip(l)), block)
        push!(results, _parse_result_rows(dev_type, has_ci, data))
    end

    return StabilitySuite(results,
                          Symbol(get(meta, "data_kind", "phase")),
                          parse(Float64, get(meta, "tau0", "NaN")),
                          parse(Int,     get(meta, "n", "0")),
                          parse(Float64, get(meta, "confidence", "NaN")),
                          Symbol(get(meta, "tau_mode", "explicit")))
end

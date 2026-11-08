# Save and load StabilityResult and StabilitySuite as tab-separated text.

const _IO_VERSION = "2"

# Shared row writer: the `# calc_ci=` flag, the column header, and one tab-
# delimited row per τ. Used by both result files (v1) and suite files (v2).
# The on-disk token stays `# calc_ci=` (not `# ci=`) on purpose: it is the v1/v2
# file-format key, independent of the public `ci` kwarg, so result files written
# by earlier versions keep loading. Do not rename it without bumping _IO_VERSION.
function _write_result_table(io::IO, r::StabilityResult)
    has_ci = !isempty(r.ci)
    println(io, "# calc_ci=", has_ci)
    println(io, "tau\tdev\tnoise_type\tci_lower\tci_upper\tedf\tneff")
    for k in 1:length(r.tau)
        noise_s = has_ci ? string(r.noise_type[k]) : ""
        ci_lo_s = has_ci ? string(r.ci[k].lo)      : "NaN"
        ci_hi_s = has_ci ? string(r.ci[k].hi)      : "NaN"
        edf_s   = has_ci ? string(r.edf[k])        : "NaN"
        println(io, string(r.tau[k]), '\t', string(r.dev[k]), '\t',
                    noise_s, '\t', ci_lo_s, '\t', ci_hi_s, '\t', edf_s, '\t',
                    string(r.neff[k]))
    end
end

# Shared row parser: turn tab-delimited data rows back into a StabilityResult.
# Accepts 7-column rows and the 6-column rows of v0.4 files (neff loads as 0).
function _parse_result_rows(deviation_type::Symbol, has_ci::Bool,
                            data_lines::AbstractVector{<:AbstractString},
                            tau0::Float64 = NaN, n::Int = 0, confidence::Float64 = NaN)
    N        = length(data_lines)
    tau      = Vector{Float64}(undef, N)
    dev      = Vector{Float64}(undef, N)
    noise_s  = Vector{String}(undef, N)
    ci_lo    = Vector{Float64}(undef, N)
    ci_hi    = Vector{Float64}(undef, N)
    edf_vec  = Vector{Float64}(undef, N)
    neff     = zeros(Int, N)

    for (k, line) in enumerate(data_lines)
        parts = split(line, '\t')
        length(parts) in (6, 7) || error("malformed row $k (expected 6 or 7 columns, got $(length(parts)))")
        tau[k]     = parse(Float64, parts[1])
        dev[k]     = parse(Float64, parts[2])
        noise_s[k] = parts[3]
        ci_lo[k]   = parse(Float64, parts[4])
        ci_hi[k]   = parse(Float64, parts[5])
        edf_vec[k] = parse(Float64, parts[6])
        length(parts) == 7 && (neff[k] = parse(Int, parts[7]))
    end

    has_ci || return StabilityResult(deviation_type, tau, dev,
                                     Symbol[], _CIBound[], Float64[], neff, tau0, n, NaN)
    return StabilityResult(deviation_type, tau, dev,
                           Symbol.(noise_s), _zip_ci(ci_lo, ci_hi), edf_vec, neff,
                           tau0, n, confidence)
end

function _save_result(path::AbstractString, r::StabilityResult)
    open(path, "w") do io
        println(io, "# SigmaTau StabilityResult v1")
        println(io, "# deviation_type=", r.deviation_type)
        println(io, "# tau0=", r.tau0)
        println(io, "# n=", r.n)
        println(io, "# confidence=", r.confidence)
        _write_result_table(io, r)
    end
    return path
end

function _load_result(path::AbstractString)
    lines = readlines(path)
    (!isempty(lines) && occursin("StabilitySuite", lines[1])) &&
        error("load: \"$path\" is a StabilitySuite file")

    deviation_type = :unknown
    has_ci         = true   # default: assume CI present if header is absent (v1 compat)
    tau0, n, confidence = NaN, 0, NaN   # absent in v0.4 files

    for line in lines
        if startswith(line, "# deviation_type=")
            deviation_type = Symbol(line[length("# deviation_type=") + 1 : end])
        elseif startswith(line, "# calc_ci=")
            has_ci = strip(line[length("# calc_ci=") + 1 : end]) == "true"
        elseif startswith(line, "# tau0=")
            tau0 = parse(Float64, strip(line[length("# tau0=") + 1 : end]))
        elseif startswith(line, "# n=")
            n = parse(Int, strip(line[length("# n=") + 1 : end]))
        elseif startswith(line, "# confidence=")
            confidence = parse(Float64, strip(line[length("# confidence=") + 1 : end]))
        end
    end

    data_lines = filter(l -> !startswith(l, '#') && !startswith(l, "tau"), lines)
    isempty(data_lines) && error("load: no data rows found in \"$path\"")

    return _parse_result_rows(deviation_type, has_ci, data_lines, tau0, n, confidence)
end

function _save_suite(path::AbstractString, suite::StabilitySuite; source_file::AbstractString="")
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

function _load_suite(path::AbstractString)
    lines = readlines(path)
    (!isempty(lines) && occursin("StabilitySuite", lines[1])) ||
        error("load: \"$path\" is not a StabilitySuite file")

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
    isempty(delim_idxs) && error("load: no result blocks found in \"$path\"")

    tau0 = parse(Float64, get(meta, "tau0", "NaN"))
    n    = parse(Int,     get(meta, "n", "0"))
    conf = parse(Float64, get(meta, "confidence", "NaN"))
    results = StabilityResult[]
    for (b, start) in enumerate(delim_idxs)
        stop     = b < length(delim_idxs) ? delim_idxs[b + 1] - 1 : length(lines)
        block    = @view lines[start:stop]
        dev_type = Symbol(match(r"^# --- result (\S+) ---", block[1]).captures[1])
        ci_lines = filter(l -> startswith(l, "# calc_ci="), block)
        has_ci   = !isempty(ci_lines) && strip(split(ci_lines[1], "=")[2]) == "true"
        data     = filter(l -> !startswith(l, '#') && !startswith(l, "tau") && !isempty(strip(l)), block)
        push!(results, _parse_result_rows(dev_type, has_ci, data, tau0, n, conf))
    end

    return StabilitySuite(results,
                          Symbol(get(meta, "data_kind", "phase")),
                          tau0, n, conf,
                          Symbol(get(meta, "tau_mode", "explicit")))
end

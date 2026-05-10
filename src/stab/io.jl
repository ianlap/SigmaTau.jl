# stab/io.jl — Save/load StabilityResult as tab-separated text (stdlib only)

const _IO_VERSION = "1"

"""
    save_result(path::AbstractString, r::StabilityResult) → path

Write a `StabilityResult` to a tab-delimited text file at `path`.

The file is self-describing: three comment lines encode the format version,
deviation type, and whether confidence intervals were computed (`calc_ci`),
followed by a column-name header and one data row per τ value. All six data
columns are always present; `noise_type`, `ci_lower`, `ci_upper`, and `edf`
are written as `""` / `NaN` when `calc_ci=false`.

Returns `path` so the call can be chained.

# Examples

```julia
pd = PhaseData(randn(256), 1.0)
r  = adev(pd, [1, 2, 4, 8, 16])
save_result("my_adev.tsv", r)
r2 = load_result("my_adev.tsv")
r2.dev ≈ r.dev   # true
```

See also [`load_result`](@ref).
"""
function save_result(path::AbstractString, r::StabilityResult)
    N      = length(r.tau)
    has_ci = !isempty(r.ci_lower)

    open(path, "w") do io
        println(io, "# SigmaTau StabilityResult v", _IO_VERSION)
        println(io, "# deviation_type=", r.deviation_type)
        println(io, "# calc_ci=", has_ci)
        println(io, "tau\tdev\tnoise_type\tci_lower\tci_upper\tedf")
        for k in 1:N
            noise_s  = has_ci ? string(r.noise_type[k]) : ""
            ci_lo_s  = has_ci ? string(r.ci_lower[k])   : "NaN"
            ci_hi_s  = has_ci ? string(r.ci_upper[k])   : "NaN"
            edf_s    = has_ci ? string(r.edf[k])         : "NaN"
            println(io, string(r.tau[k]), '\t', string(r.dev[k]), '\t',
                        noise_s, '\t', ci_lo_s, '\t', ci_hi_s, '\t', edf_s)
        end
    end
    return path
end

"""
    load_result(path::AbstractString) → StabilityResult

Read a `StabilityResult` previously written by [`save_result`](@ref).

The `# calc_ci=` header line determines whether confidence interval fields are
reconstructed. When `calc_ci=false`, `noise_type`, `ci_lower`, `ci_upper`, and
`edf` are returned as empty vectors, matching the standard deviation API contract.

See also [`save_result`](@ref).
"""
function load_result(path::AbstractString)
    lines = readlines(path)

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

    N        = length(data_lines)
    tau      = Vector{Float64}(undef, N)
    dev      = Vector{Float64}(undef, N)
    noise_s  = Vector{String}(undef, N)
    ci_lower = Vector{Float64}(undef, N)
    ci_upper = Vector{Float64}(undef, N)
    edf_vec  = Vector{Float64}(undef, N)

    for (k, line) in enumerate(data_lines)
        parts = split(line, '\t')
        length(parts) == 6 || error("load_result: malformed row $k in \"$path\" (expected 6 columns, got $(length(parts)))")
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

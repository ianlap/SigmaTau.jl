# io/save.jl — `save`: the canonical writer for timing data, results, and
# suites. Timing records go out as commented two-column text that
# `read_phase` / `read_frequency` round-trip with no keywords; results and
# suites delegate to the existing `save_result` / `save_suite` writers.

"""
    save(path::AbstractString, data::PhaseData)      → path
    save(path::AbstractString, data::FrequencyData)  → path
    save(path::AbstractString, r::StabilityResult)   → path
    save(path::AbstractString, s::StabilitySuite; source_file="") → path

Write a SigmaTau object to a plain text file at `path`. Returns `path` so the
call can be chained.

For `PhaseData` / `FrequencyData` the file is two delimited columns — sample
time `(i-1)·τ₀` in seconds and the sample value — preceded by a comment
header recording the provenance:

```
# SigmaTau phase data
# source: user
# tau0: 1.0
0.0     -1.7881e-10
1.0     2.1406e-10
…
```

`source` is `"user"` for records constructed directly from a vector and the
originating file path for records that came from [`read_phase`](@ref) /
[`read_frequency`](@ref). The delimiter follows the extension convention of
the readers: `,` for `.csv`, tab otherwise. [`read_phase`](@ref) /
[`read_frequency`](@ref) skip the comment header automatically and pick up
`# tau0:` when no `tau0` keyword is given, so `read_phase(save(path, pd))`
reproduces `pd` exactly — and any external tool can consume the file as a
plain two-column table. Typical use: generate noise with [`noise_gen`](@ref)
and hand it to software outside SigmaTau.

For `StabilityResult` / `StabilitySuite` the call delegates to the
tab-delimited writers ([`save_result`](@ref) / [`save_suite`](@ref), which
remain available under their original names); read those files back with
[`load_result`](@ref) / [`load_suite`](@ref).

!!! note
    FileIO.jl also exports a `save` function. If both packages are loaded,
    qualify the call as `SigmaTau.save(path, data)`.

# Examples

```julia
pd = noise_gen(PhaseData, 4096, 1.0; sigma1=Dict(0 => 1e-12))
save("white_pm.txt", pd)
pd2 = read_phase("white_pm.txt")     # tau0 and values round-trip

r = adev(pd)
save("white_pm_adev.tsv", r)         # same file save_result writes
```

See also [`read_phase`](@ref), [`read_frequency`](@ref),
[`load_result`](@ref), [`load_suite`](@ref).
"""
save(path::AbstractString, data::PhaseData) =
    _save_data(path, data.x, data.tau0, data.source, "phase")
save(path::AbstractString, data::FrequencyData) =
    _save_data(path, data.y, data.tau0, data.source, "frequency")
save(path::AbstractString, r::StabilityResult) = save_result(path, r)
save(path::AbstractString, s::StabilitySuite; source_file::AbstractString="") =
    save_suite(path, s; source_file=source_file)

function _save_data(path::AbstractString, v::AbstractVector{<:Real},
                    tau0::Float64, source::String, kind::String)
    d = lowercase(splitext(path)[2]) == ".csv" ? ',' : '\t'
    open(path, "w") do io
        println(io, "# SigmaTau ", kind, " data")
        println(io, "# source: ", source)
        println(io, "# tau0: ", tau0)
        for (i, val) in enumerate(v)
            println(io, (i - 1) * tau0, d, val)
        end
    end
    return path
end

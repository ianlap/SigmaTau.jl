# `save` and `load`: one writer and one reader for records, results and
# suites. A record is written as commented two-column text, which
# `read_phase` and `read_frequency` read back without keywords; results and
# suites are written as tab-separated tables (io/results.jl).

"""
    save(path, data::PhaseData) → path
    save(path, data::FrequencyData) → path
    save(path, r::StabilityResult) → path
    save(path, s::StabilitySuite; source_file="") → path

Write a record, a result or a suite to a text file and return `path`.
[`load`](@ref) reads the file back.

# Records

A `PhaseData` or `FrequencyData` record is written as two columns, the
sample time ``(i-1)τ_0`` in seconds and the sample, after a comment header:

```
# SigmaTau phase data
# source: user
# tau0: 1.0
0.0     -1.7881e-10
1.0     2.1406e-10
…
```

The columns are separated by a comma when `path` ends in `.csv`, and by a
tab otherwise. `source` is `"user"` for a record made from a vector, and
the path of the file for a record that was read. A record that knows the
UTC time of its first sample writes it as well, as
`# start: 2026-11-08T21:41:03Z`.

[`read_phase`](@ref) and [`read_frequency`](@ref) skip the header, take
``τ_0`` and the start time from it, and return the record as it was saved.
Any other program can read the file as a two-column table.

# Results

A `StabilityResult` is written as a tab-separated table with one row for
each ``τ`` and the columns `tau`, `dev`, `noise_type`, `ci_lower`,
`ci_upper`, `edf` and `neff`. Comment lines before it give the format
version, the estimator, ``τ_0``, the number of samples, the confidence
level and whether the result has intervals. Without intervals, the
interval and EDF columns hold `NaN`.

A `StabilitySuite` is written as one such table for each result, after a
header with the package version, the time of writing, `source_file`, the
kind of data, ``τ_0``, the number of samples, the confidence level, the
``τ`` grid and the estimators.

!!! note "FileIO.jl"
    FileIO.jl exports `save` and `load` as well. With both packages loaded,
    write `SigmaTau.save(path, data)` and `SigmaTau.load(path)`.

# Examples

```jldoctest
julia> using SigmaTau

julia> x = PhaseData([0.0, 1.5e-9, 1.0e-9, 2.5e-9], 10.0);

julia> y = read_phase(save(joinpath(mktempdir(), "x.txt"), x));

julia> y.x == x.x && y.tau0 == x.tau0
true

julia> r = adev(x, [1]; ci=false);

julia> load(save(joinpath(mktempdir(), "r.tsv"), r)).dev == r.dev
true
```
"""
save(path::AbstractString, data::PhaseData) =
    _save_data(path, data.x, data.tau0, data.source, data.start, "phase")
save(path::AbstractString, data::FrequencyData) =
    _save_data(path, data.y, data.tau0, data.source, data.start, "frequency")
save(path::AbstractString, r::StabilityResult) = _save_result(path, r)
save(path::AbstractString, s::StabilitySuite; source_file::AbstractString="") =
    _save_suite(path, s; source_file=source_file)

"""
    load(path) → StabilityResult, StabilitySuite, PhaseData or FrequencyData

Read a file written by [`save`](@ref). Its first line says what the file
holds: a result, a suite, or a phase or frequency record. Other data files
are read with [`read_phase`](@ref), [`read_frequency`](@ref) or
[`read_data`](@ref).

A result file written by version 0.4 carries no ``τ_0``, number of samples,
confidence level or `neff`. It loads with `tau0 = NaN`, `n = 0`,
`confidence = NaN` and `neff` zero.
"""
function load(path::AbstractString)
    first = open(readline, path)
    occursin("StabilitySuite", first) && return _load_suite(path)
    startswith(first, "# SigmaTau phase data") && return read_phase(path)
    startswith(first, "# SigmaTau frequency data") && return read_frequency(path)
    startswith(first, "# SigmaTau") || throw(ArgumentError(
        "load: $path was not written by save; read a data file with read_phase, read_frequency or read_data"))
    return _load_result(path)
end

function _save_data(path::AbstractString, v::AbstractVector{<:Real},
                    tau0::Float64, source::String, start::Float64, kind::String)
    d = lowercase(splitext(path)[2]) == ".csv" ? ',' : '\t'
    open(path, "w") do io
        println(io, "# SigmaTau ", kind, " data")
        println(io, "# source: ", source)
        println(io, "# tau0: ", tau0)
        isnan(start) || println(io, "# start: ", _iso_text(start))
        for (i, val) in enumerate(v)
            println(io, (i - 1) * tau0, d, val)
        end
    end
    return path
end

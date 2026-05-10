# io/read.jl — Read raw phase / frequency data files into PhaseData /
# FrequencyData with optional preprocessing (scaling, detrend, gap fill).

"""
    _read_columns(path; header=0, delim=nothing) → Matrix{Float64}

stdlib `DelimitedFiles.readdlm` wrapper. Auto-detects the delimiter from
the file extension unless `delim` is given:

- `.csv` → `,`
- `.tsv` → `\\t`
- everything else → whitespace (`readdlm` default)

The result must be coercible to `Matrix{Float64}`. For CSV files with
string headers, set `header` to the number of header rows to skip.
"""
function _read_columns(path::AbstractString; header::Int=0, delim::Union{Char,Nothing}=nothing)
    isfile(path) || throw(ArgumentError("read_columns: file not found: $path"))
    ext = lowercase(splitext(path)[2])
    d   = if delim !== nothing
        delim
    elseif ext == ".csv"
        ','
    elseif ext == ".tsv"
        '\t'
    else
        nothing
    end
    raw = if d === nothing
        readdlm(path; skipstart=header)
    else
        readdlm(path, d; skipstart=header)
    end
    M = try
        Matrix{Float64}(raw)
    catch err
        throw(ArgumentError("read_columns: failed to coerce $path to Float64 matrix " *
                            "(non-numeric content?). Original error: $err"))
    end
    return M
end

"""
    _ingest(path; tau0, time_col, value_col, header, delim, scaling, detrend, fillgaps)
        → (values::Vector{Float64}, tau0::Float64)

Internal helper shared by `read_phase` and `read_frequency`. Returns the
post-preprocessing value vector and the resolved `tau0`.
"""
function _ingest(path::AbstractString;
                 tau0::Union{Real,Nothing},
                 time_col::Int,
                 value_col::Int,
                 header::Int,
                 delim::Union{Char,Nothing},
                 scaling::Real,
                 detrend::Symbol,
                 fillgaps::Bool)
    M = _read_columns(path; header=header, delim=delim)
    ncols = size(M, 2)
    value_col >= 1 || throw(ArgumentError("read: value_col must be ≥ 1"))
    value_col <= ncols ||
        throw(ArgumentError("read: value_col=$value_col exceeds available columns ($ncols)"))

    have_time = time_col >= 1
    if have_time
        time_col <= ncols ||
            throw(ArgumentError("read: time_col=$time_col exceeds available columns ($ncols)"))
        t = @view M[:, time_col]
    end
    v = Vector{Float64}(M[:, value_col])

    if scaling != 1
        v .*= scaling
    end

    # Resolve tau0
    resolved_tau0 = if tau0 !== nothing
        Float64(tau0)
    elseif have_time && length(t) > 1
        Float64(median(diff(t)))
    else
        throw(ArgumentError("read: tau0 must be supplied when no time column is present"))
    end

    # Gap filling (requires a time column)
    if fillgaps
        have_time ||
            throw(ArgumentError("read: fillgaps=true requires a time column (got time_col=0)"))
        _, v_eq = _make_equispaced(t, v)
        v, _    = _howe_fillgaps_core(v_eq)
    end

    if detrend !== :none
        _detrend_core!(v, detrend)
    end

    return v, resolved_tau0
end

"""
    read_phase(path; tau0=nothing, time_col=1, value_col=2, header=0,
                     delim=nothing, scaling=1.0,
                     detrend=:none, fillgaps=false) → PhaseData

Read a phase-data file and return a [`PhaseData`](@ref) ready for stability
analysis.

# Keyword arguments

- `tau0`     — sample interval in seconds. Auto-inferred from the time column
               (median of `diff(t)`) when omitted; required when `time_col=0`.
- `time_col` — 1-based column index for timestamps (`0` = no time column).
- `value_col` — 1-based column index for phase samples.
- `header`   — number of header lines to skip.
- `delim`    — field delimiter (`Char`). Auto-detected from extension when
               `nothing` (`,` for `.csv`, tab for `.tsv`, whitespace otherwise).
- `scaling`  — multiply phase samples by this factor (e.g. `1e-12` for ps→s).
- `detrend`  — one of `:none`, `:mean`, `:endpoint`, `:linear`. See [`detrend`](@ref).
- `fillgaps` — when `true`, equispace the record on the minimum time spacing
               and impute `NaN`s with [`fillgaps`](@ref) (Howe's algorithm).
               Requires a time column.

# Examples

```julia
pd = read_phase("phase.csv")                          # 2-col CSV, auto tau0
pd = read_phase("p.tsv"; detrend=:linear)             # remove linear trend
pd = read_phase("p.txt"; time_col=0, tau0=1.0)        # no time column
pd = read_phase("g.csv"; fillgaps=true, scaling=1e-9) # ns→s + Howe fill
```
"""
function read_phase(path::AbstractString;
                    tau0::Union{Real,Nothing}=nothing,
                    time_col::Int=1,
                    value_col::Int=2,
                    header::Int=0,
                    delim::Union{Char,Nothing}=nothing,
                    scaling::Real=1.0,
                    detrend::Symbol=:none,
                    fillgaps::Bool=false)
    v, t0 = _ingest(path;
                    tau0=tau0, time_col=time_col, value_col=value_col,
                    header=header, delim=delim,
                    scaling=scaling, detrend=detrend, fillgaps=fillgaps)
    return PhaseData(v, t0)
end

"""
    read_frequency(path; tau0=nothing, time_col=1, value_col=2, header=0,
                         delim=nothing, scaling=1.0,
                         detrend=:none, fillgaps=false) → FrequencyData

Read a fractional-frequency file. Same keyword arguments as
[`read_phase`](@ref); only the return type differs.
"""
function read_frequency(path::AbstractString;
                        tau0::Union{Real,Nothing}=nothing,
                        time_col::Int=1,
                        value_col::Int=2,
                        header::Int=0,
                        delim::Union{Char,Nothing}=nothing,
                        scaling::Real=1.0,
                        detrend::Symbol=:none,
                        fillgaps::Bool=false)
    v, t0 = _ingest(path;
                    tau0=tau0, time_col=time_col, value_col=value_col,
                    header=header, delim=delim,
                    scaling=scaling, detrend=detrend, fillgaps=fillgaps)
    return FrequencyData(v, t0)
end

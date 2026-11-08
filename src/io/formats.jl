# Instrument file formats: detection, and readers for the files that are not
# plain numeric columns (TimeLab .tim, Stable32 .phd/.frd, frequency-counter
# logs, logs with named columns). Plain columns go through read_phase /
# read_frequency.

"""
Input formats [`read_data`](@ref) understands, as [`detect_format`](@ref)
names them.
"""
const INPUT_FORMATS = (:columns, :stable32, :timelab, :counter, :table)

# What a format reader returns: the samples, τ₀, whether they are phase, and
# which samples fillgaps put back (empty when none were).
struct _Record
    v::Vector{Float64}
    tau0::Float64
    phase::Bool
    filled::BitVector
end
_Record(v::Vector{Float64}, tau0::Float64, phase::Bool) = _Record(v, tau0, phase, falses(0))
_Record(vf::Tuple{Vector{Float64},BitVector}, tau0::Float64, phase::Bool) = _Record(vf[1], tau0, phase, vf[2])

_split_fields(line::AbstractString, d::Union{Char,Nothing}) =
    d === nothing ? split(line; keepempty=false) : map(strip, split(line, d))

function _head_lines(path::AbstractString, n::Int)
    out = String[]
    io = open(path, "r")
    try
        for line in eachline(io)
            push!(out, rstrip(line, '\r'))
            length(out) >= n && break
        end
    finally
        close(io)
    end
    return out
end

_isnumfield(f::AbstractString) = tryparse(Float64, f) !== nothing

# Keysight / Agilent SCPI real format: "+1.00000000012345E+07".
const _SCPI_REAL = r"^[+-]\d\.\d+E[+-]\d+$"
# SCPI's "not a number" reading (9.91E+37), written for a missed gate.
_scpi_nan(v::Float64) = v == 9.91e37 ? NaN : v

const _TIMELAB_FIELD = r"^(STR|DBL|S32|BLN|S64|U32) 0x[0-9A-Fa-f]+ \"([^\"]*)\"\s?(.*)$"

function _delim_for(path::AbstractString, line::AbstractString, delim::Char='\0')
    delim == '\0' || return delim
    ext = lowercase(splitext(path)[2])
    ext == ".csv" && return ','
    ext == ".tsv" && return '\t'
    occursin(',', line) && return ','
    occursin('\t', line) && return '\t'
    return nothing
end

"""
    detect_format(path; delim=nothing) → Symbol

Name the layout of a data file from its extension and first lines, one of
[`INPUT_FORMATS`](@ref):

- `:timelab`: a TimeLab `.tim` measurement log (also what the Microchip
  53100A, 5120A/5125A and TimePod write through TimeLab): typed header
  fields such as `DBL 0x… "Sample Interval" 0.1`, then a `TIC n` block of
  phase samples in seconds.
- `:stable32`: a Stable32 `.phd` (phase) or `.frd` (frequency) file.
- `:counter`: frequency-counter readings in SCPI real format
  (`+1.00000000012345E+07`, as the Keysight 53230A and 53220A log them),
  one per line, optionally after a timestamp.
- `:table`: numeric columns under a line of column names (a CSV export,
  most GPSDO and logging-script output).
- `:columns`: plain numeric columns, read by [`read_phase`](@ref) and
  [`read_frequency`](@ref).

`delim` is the field delimiter, by default comma for `.csv`, tab for
`.tsv`, else whitespace.
"""
detect_format(path::AbstractString; delim::Union{AbstractChar,Nothing}=nothing) =
    _detect_format(path, delim === nothing ? '\0' : Char(delim))

function _detect_format(path::AbstractString, delim::Char)
    isfile(path) || throw(ArgumentError("detect_format: file not found: " * String(path)))
    lines = _head_lines(path, 200)
    any(l -> startswith(l, "TIC "), lines) &&
        any(l -> match(_TIMELAB_FIELD, l) !== nothing, lines) && return :timelab
    lowercase(splitext(path)[2]) in (".phd", ".frd") && return :stable32
    # Any other Stable32 export (a text header closed by "# Header End") is
    # plain columns to read_phase.
    any(l -> strip(l) == "# Header End", lines) && return :columns
    data = String[]
    for l in lines
        s = strip(l)
        (isempty(s) || startswith(s, '#')) && continue
        push!(data, String(s))
        length(data) >= 8 && break
    end
    isempty(data) && return :columns
    d = _delim_for(path, data[1], delim)
    first_fields = _split_fields(data[1], d)
    if length(data) >= 2 && !all(_isnumfield, first_fields)
        next = _split_fields(data[2], d)
        length(next) == length(first_fields) && any(_isnumfield, next) && return :table
        return :columns
    end
    all(data) do l
        f = _split_fields(l, d)
        1 <= length(f) <= 2 && occursin(_SCPI_REAL, f[end]) && all(_isnumfield, f)
    end && return :counter
    return :columns
end

# τ₀: the explicit value, else the header's, else the median time step.
function _resolve_tau0(tau0::Float64, header::Float64, t::Vector{Float64}, what::String)
    isnan(tau0) || return tau0
    isnan(header) || return header
    length(t) >= 2 && return _median_step(t)
    throw(ArgumentError("read_data: " * what * " has no sample interval; pass tau0"))
end

# With fillgaps, samples go on an even grid of the time column (or, without
# one, stay where they are) and NaN samples are filled. Returns the samples
# and the mask of filled ones.
function _equispace(t::Vector{Float64}, v::Vector{Float64}, fillgaps::Bool)
    fillgaps || return v, falses(0)
    v_eq = length(t) >= 2 ? _make_equispaced(t, v)[2] : v
    filled, mask = _howe_fillgaps_core(v_eq)
    return filled, BitVector(mask)
end

function _read_timelab(path::AbstractString, tau0::Float64)
    fields = Dict{String,String}()
    v = Float64[]
    block = ""
    io = open(path, "r")
    try
        for raw in eachline(io)
            line = strip(raw)
            (isempty(line) || startswith(line, ';')) && continue
            if isempty(block)
                m = match(_TIMELAB_FIELD, line)
                if m !== nothing
                    k = m.captures[2]; val = m.captures[3]
                    (k === nothing || val === nothing) || (fields[String(k)] = String(strip(val)))
                    continue
                end
                block = String(first(split(line)))
                block == "TIC" || throw(ArgumentError(
                    "read_data: TimeLab block \"" * block * "\" is not supported (only TIC phase blocks)"))
                continue
            end
            x = tryparse(Float64, line)
            x === nothing && throw(ArgumentError("read_data: non-numeric TimeLab sample \"" * String(line) * "\" in " * String(path)))
            push!(v, x)
        end
    finally
        close(io)
    end
    isempty(block) && throw(ArgumentError("read_data: no TIC data block in " * String(path)))
    scale = tryparse(Float64, get(fields, "Scale Factor", "1"))
    scale === nothing || scale == 1 || (v .*= scale)
    si = tryparse(Float64, get(fields, "Sample Interval", ""))
    return _Record(v, _resolve_tau0(tau0, si === nothing ? NaN : si, Float64[], "this TimeLab file"), true)
end

function _read_stable32(path::AbstractString, tau0::Float64, hint::Int, fillgaps::Bool)
    nskip, htau = _scan_comment_header(path)
    # Record type: the header's "Type:" line, else the extension.
    phase = lowercase(splitext(path)[2]) != ".frd"
    for l in _head_lines(path, nskip)
        m = match(r"^Type:\s*(phase|freq)"i, strip(l))
        m === nothing || (phase = lowercase(something(m.captures[1], "phase")) == "phase")
    end
    hint == 0 || (phase = hint > 0)
    M = _read_columns(path; header=nskip)
    size(M, 1) >= 2 || throw(ArgumentError("read_data: fewer than 2 samples in " * String(path)))
    v = Vector{Float64}(M[:, end])
    t = size(M, 2) >= 2 ? Vector{Float64}(M[:, 1]) : Float64[]
    # A Stable32 time tag is an MJD (days); a large value with sub-day steps.
    length(t) >= 2 && t[1] > 30000 && _median_step(t) < 1 && (t .*= 86400)
    t0 = _resolve_tau0(tau0, htau === nothing ? NaN : htau, t, "this Stable32 file")
    return _Record(_equispace(t, v, fillgaps), t0, phase)
end

# Nominal frequency of a counter log: the mean reading to 4 significant
# figures (10 MHz, 5 MHz, 10.23 MHz, 1 Hz ...), unless given.
function _nominal(f::Vector{Float64}, nominal::Float64)
    isnan(nominal) || return nominal
    s = 0.0; n = 0
    for x in f
        isfinite(x) && (s += x; n += 1)
    end
    n > 0 || throw(ArgumentError("read_data: no finite counter readings"))
    return round(s / n; sigdigits=4)
end

function _read_counter(path::AbstractString, tau0::Float64, nominal::Float64, hint::Int,
                       delim::Char, fillgaps::Bool)
    t = Float64[]; v = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, '#')) && continue
        f = _split_fields(s, _delim_for(path, s, delim))
        x = tryparse(Float64, f[end])
        x === nothing && throw(ArgumentError("read_data: non-numeric counter reading \"" * String(s) * "\" in " * String(path)))
        push!(v, _scpi_nan(x))
        if length(f) >= 2
            tt = tryparse(Float64, f[1])
            tt === nothing || push!(t, tt)
        end
    end
    length(v) >= 2 || throw(ArgumentError("read_data: fewer than 2 readings in " * String(path)))
    length(t) == length(v) || empty!(t)
    t0 = _resolve_tau0(tau0, NaN, t, "this counter log (no timestamps; pass tau0 = the gate time)")
    finite = filter(isfinite, v)
    isempty(finite) && throw(ArgumentError("read_data: no finite counter readings in " * String(path)))
    # Time-interval readings (seconds) are phase; anything ≥ 1 is a frequency
    # in Hz, unless the caller says which.
    phase = hint == 0 ? _median_nonan(abs.(finite)) < 1 : hint > 0
    if !phase
        f0 = _nominal(v, nominal)
        v = (v .- f0) ./ f0
    end
    return _Record(_equispace(t, v, fillgaps), t0, phase)
end

# Units in a column name: "(ns)", "[ps]", "_us", " ppb", "Hz"; NaN = Hz.
const _UNITS = [(Regex("(^|[\\s(\\[_/-])" * u * "([\\s)\\]]|\$)"), f)
                for (u, f) in (("ps", 1e-12), ("ns", 1e-9), ("us", 1e-6), ("µs", 1e-6), ("μs", 1e-6),
                               ("ms", 1e-3), ("ppt", 1e-12), ("ppb", 1e-9), ("ppm", 1e-6),
                               ("hz", NaN), ("s", 1.0))]

function _unit_scale(name::AbstractString)
    n = lowercase(name)
    for (re, f) in _UNITS
        occursin(re, n) && return f
    end
    return 1.0
end

const _PHASE_NAMES = r"phase|time ?err|\btie\b|pps|time ?offset|time ?diff|\bx\b"i
const _FREQ_NAMES  = r"freq|\bffe\b|osc|fractional|\by\b"i
const _TIME_NAMES  = r"^(t|mjd|time|timestamp|date|datetime|utc|epoch|elapsed|sec|seconds)\s*([(\[].*)?$"i

# UTC seconds since 2000-01-01 of an ISO 8601 timestamp
# ("2026-11-08T21:41:03.5", a space for the T, then Z, ±hh:mm, ±hhmm or
# nothing); NaN if it is not one.
function _iso_seconds(s::AbstractString)
    m = match(r"^(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2}):(\d{2}(?:\.\d*)?)(Z|([+-])(\d{2}):?(\d{2}))?$", strip(s))
    m === nothing && return NaN
    c = Float64[]
    for i in 1:6
        x = m.captures[i]
        x === nothing && return NaN
        push!(c, parse(Float64, x))
    end
    days = Dates.value(Dates.Date(Int(c[1]), Int(c[2]), Int(c[3])) - Dates.Date(2000, 1, 1))
    sec = days * 86400.0 + c[4] * 3600 + c[5] * 60 + c[6]
    sign, hh, mm = m.captures[8], m.captures[9], m.captures[10]
    (sign === nothing || hh === nothing || mm === nothing) && return sec
    offset = parse(Int, hh) * 3600 + parse(Int, mm) * 60
    return sign == "+" ? sec - offset : sec + offset
end

function _pick_column(names::Vector{String}, column::Int, column_name::String, time::Int)
    if column > 0
        column <= length(names) || throw(ArgumentError("read_data: column=" * string(column) * " exceeds " * string(length(names)) * " columns"))
        return column
    end
    if !isempty(column_name)
        i = findfirst(n -> lowercase(n) == lowercase(column_name), names)
        i === nothing && throw(ArgumentError("read_data: no column named \"" * column_name * "\"; columns: " * join(names, ", ")))
        return i
    end
    for re in (_PHASE_NAMES, _FREQ_NAMES)
        hits = [i for i in eachindex(names) if i != time && occursin(re, names[i])]
        length(hits) == 1 && return hits[1]
    end
    throw(ArgumentError("read_data: cannot tell which column holds the clock; pass column= one of " *
                        join(names, ", ")))
end

function _read_table(path::AbstractString, tau0::Float64, column::Int, column_name::String,
                     nominal::Float64, phase_hint::Int, delim::Char, fillgaps::Bool)
    names = String[]
    rows = Vector{Vector{String}}()
    d::Union{Char,Nothing} = nothing
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, '#')) && continue
        if isempty(names)
            d = _delim_for(path, s, delim)
            names = String[String(strip(f, '"')) for f in _split_fields(s, d)]
            continue
        end
        push!(rows, String[String(strip(f, '"')) for f in _split_fields(s, d)])
    end
    length(rows) >= 2 || throw(ArgumentError("read_data: fewer than 2 rows under the column names in " * String(path)))
    # The time column by name, else the first one holding ISO timestamps.
    time = something(findfirst(n -> occursin(_TIME_NAMES, n), names),
                     findfirst(f -> !isnan(_iso_seconds(f)), rows[1]), 0)
    c = _pick_column(names, column, column_name, time)
    v = Float64[]
    for (k, r) in enumerate(rows)
        length(r) == length(names) || throw(ArgumentError(
            "read_data: row " * string(k) * " has " * string(length(r)) * " fields, expected " * string(length(names))))
        x = tryparse(Float64, r[c])
        push!(v, x === nothing ? NaN : x)
    end
    t = Float64[]
    if time > 0
        for r in rows
            x = tryparse(Float64, r[time])
            push!(t, x === nothing ? _iso_seconds(r[time]) : x)
        end
        any(isnan, t) && empty!(t)
        occursin(r"mjd"i, names[time]) && (t .*= 86400)
    end
    phase = phase_hint == 0 ? !occursin(_FREQ_NAMES, names[c]) || occursin(_PHASE_NAMES, names[c]) :
            phase_hint > 0
    scale = _unit_scale(names[c])
    if isnan(scale)            # absolute frequency in Hz
        f0 = _nominal(v, nominal)
        v = (v .- f0) ./ f0
        phase = false
    elseif scale != 1
        v .*= scale
    end
    t0 = _resolve_tau0(tau0, NaN, t, "this table (no time column)")
    return _Record(_equispace(t, v, fillgaps), t0, phase)
end

# The trim-friendly core of read_data: every argument concrete, NaN / 0 / ""
# meaning "not given", `kind` one of :auto, :phase, :frequency.
function _read_data(path::AbstractString, format::Symbol; tau0::Float64=NaN, column::Int=0,
                    column_name::String="", nominal::Float64=NaN, kind::Symbol=:auto,
                    scaling::Float64=1.0, delim::Char='\0', fillgaps::Bool=false)
    isfile(path) || throw(ArgumentError("read_data: file not found: " * String(path)))
    kind in (:auto, :phase, :frequency) ||
        throw(ArgumentError("read_data: kind must be :auto, :phase or :frequency, got :" * String(kind)))
    fmt = format === :auto ? _detect_format(path, delim) : format
    hint = kind === :phase ? 1 : kind === :frequency ? -1 : 0
    r = if fmt === :timelab
        _read_timelab(path, tau0)
    elseif fmt === :stable32
        _read_stable32(path, tau0, hint, fillgaps)
    elseif fmt === :counter
        _read_counter(path, tau0, nominal, hint, delim, fillgaps)
    elseif fmt === :table
        _read_table(path, tau0, column, column_name, nominal, hint, delim, fillgaps)
    elseif fmt === :columns
        v, t0 = _ingest(path; tau0=isnan(tau0) ? nothing : tau0, time_col=1, value_col=2,
                        header=0, delim=nothing, scaling=1.0, detrend=:none, fillgaps=fillgaps)
        _Record(v, t0, hint >= 0)
    else
        throw(ArgumentError("read_data: unknown format :" * String(fmt) * "; use one of " * join(map(String, INPUT_FORMATS), ", ")))
    end
    scaling == 1 || (r.v .*= scaling)
    phase = hint == 0 ? r.phase : hint > 0
    return _Record(r.v, r.tau0, phase, r.filled)
end

"""
    read_data(path; format=:auto, tau0=nothing, column=nothing, nominal=nothing,
                    kind=:auto, scaling=1.0, delim=nothing, fillgaps=false) → PhaseData or FrequencyData

Read a measurement file in any of the [`INPUT_FORMATS`](@ref), detected with
[`detect_format`](@ref) unless `format` names one. The record type follows
the file: TimeLab logs and counter time-interval readings are phase,
`.frd` files and counter frequency readings are frequency, and a named
column is phase or frequency by its name.

- `tau0`: sample interval in seconds; by default from the file header
  (TimeLab `Sample Interval`, Stable32 `Tau:`) or the median step of the
  time column (seconds, MJD days, or ISO 8601 timestamps).
- `column`: the clock column of a `:table` file, by 1-based index or by
  name. By default the one column whose name says phase / time error /
  PPS, else frequency.
- `nominal`: nominal frequency in Hz for readings in Hz (counter logs, a
  table column named `… (Hz)`); converted to fractional frequency
  `(f − nominal)/nominal`. By default the mean reading rounded to four
  significant figures.
- `kind`: `:phase` or `:frequency` overrides the record type.
- `scaling`: multiplies the samples after any unit conversion from a
  column name (`(ns)`, `[ps]`, `ppb` …).
- `delim`: field delimiter of counter and table files; by default comma
  for `.csv`, tab for `.tsv`, else whitespace.
- `fillgaps`: put samples on an even grid of the time column and fill the
  missing ones (and any `NaN` reading) with [`fillgaps`](@ref).

Counter readings of `9.91E+37` (SCPI's "no reading") become `NaN`.

```julia
p = read_data("trace.tim")                  # TimeLab, τ₀ from the header
y = read_data("53230A.csv"; tau0=1.0)       # counter in Hz → fractional frequency
p = read_data("gpsdo.csv"; column="PPS offset (ns)")
```
"""
function read_data(path::AbstractString; format::Symbol=:auto, tau0::Union{Real,Nothing}=nothing,
                   column::Union{Integer,AbstractString,Nothing}=nothing,
                   nominal::Union{Real,Nothing}=nothing, kind::Symbol=:auto,
                   scaling::Real=1.0, delim::Union{AbstractChar,Nothing}=nothing,
                   fillgaps::Bool=false)
    r = _read_data(path, format;
                   tau0=tau0 === nothing ? NaN : Float64(tau0),
                   column=column isa Integer ? Int(column) : 0,
                   column_name=column isa AbstractString ? String(column) : "",
                   nominal=nominal === nothing ? NaN : Float64(nominal),
                   kind=kind, scaling=Float64(scaling),
                   delim=delim === nothing ? '\0' : Char(delim), fillgaps=fillgaps)
    return r.phase ? PhaseData(r.v, r.tau0; source=String(path)) :
                     FrequencyData(r.v, r.tau0; source=String(path))
end

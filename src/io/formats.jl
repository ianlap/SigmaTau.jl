# Measurement file formats: detection, and readers for the files that are
# not plain numeric columns (TimeLab .tim, Stable32 .phd and .frd, counter
# logs, tables with named columns, BIPM UTC - UTC(k) files, RINEX clock
# files). Plain columns go through read_phase and read_frequency.

"""
    INPUT_FORMATS

The file formats that [`read_data`](@ref) reads, under the names that
[`detect_format`](@ref) returns: `:columns`, `:stable32`, `:timelab`,
`:counter`, `:table`, `:bipm_utc` and `:rinex_clk`.
"""
const INPUT_FORMATS = (:columns, :stable32, :timelab, :counter, :table, :bipm_utc, :rinex_clk)

# What a format reader returns: the samples, τ₀, whether they are phase,
# which samples fillgaps put back (empty when none were) and the start time
# (NaN when unknown).
struct _Record
    v::Vector{Float64}
    tau0::Float64
    phase::Bool
    filled::BitVector
    start::Float64
end
_Record(vf::Tuple{Vector{Float64},BitVector}, tau0::Float64, phase::Bool, start::Float64) =
    _Record(vf[1], tau0, phase, vf[2], start)

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

# The exponent format in which SCPI instruments return a reading:
# "+1.00000000012345E+07".
const _SCPI_REAL = r"^[+-]\d\.\d+E[+-]\d+$"
# The SCPI not-a-number value, which a counter returns for a gate without a
# measurement.
_scpi_nan(v::Float64) = v == 9.91e37 ? NaN : v

# The header of a BIPM UTC - UTC(k) file: "MJD  [UTC-UTC(NIST)]/ns  uA/ns  uB/ns  u/ns".
const _BIPM_UTC_HEAD = r"^\s*MJD\s+\[UTC-UTC\(.+\)\]\s*/\s*ns"i

const _TIMELAB_FIELD = r"^(STR|DBL|S32|BLN|S64|U32) 0x[0-9A-Fa-f]+ \"([^\"]*)\"\s?(.*)$"

function _delim_for(line::AbstractString, delim::Char)
    delim == '\0' || return delim
    occursin(',', line) && return ','
    occursin('\t', line) && return '\t'
    return nothing
end

"""
    detect_format(path; delim=nothing) → Symbol

Name the format of a measurement file, one of [`INPUT_FORMATS`](@ref), from
its extension and its first 200 lines. The tests are made in the order of
the table, and the first that holds decides.

| Format       | Recognized by                                                          |
|:-------------|:-----------------------------------------------------------------------|
| `:timelab`   | a `TIC` block after header fields such as `DBL 0x… "Sample Interval" 0.1` |
| `:stable32`  | the extension `.phd` or `.frd`                                         |
| `:bipm_utc`  | a heading `MJD  [UTC-UTC(k)]/ns` in the first five lines               |
| `:rinex_clk` | a `RINEX VERSION / TYPE` line of type `C` in the first five lines      |
| `:counter`   | one or two fields on every line, the last a reading such as `+1.00000000012345E+07` and the first, if any, a number or an ISO 8601 timestamp |
| `:table`     | a line of column names above rows of numbers                           |
| `:columns`   | anything else                                                          |

A file with a `# Header End` line, which closes the header of a Stable32
export, is `:columns` unless one of the first four tests holds. The counter
and table tests read the first eight lines of data, and split them at
`delim`: by default at commas, or at tabs, when the first of the lines
holds one, and at runs of blanks otherwise.

[`read_data`](@ref) says how each format is read.
"""
detect_format(path::AbstractString; delim::Union{AbstractChar,Nothing}=nothing) =
    _detect_format(path, delim === nothing ? '\0' : Char(delim))

const _RINEX_CLK_HEAD = r"^\s*\d+\.\d+\s+C\b.*RINEX VERSION / TYPE"

function _detect_format(path::AbstractString, delim::Char)
    isfile(path) || throw(ArgumentError("detect_format: file not found: " * String(path)))
    lines = _head_lines(path, 200)
    any(l -> startswith(l, "TIC "), lines) &&
        any(l -> match(_TIMELAB_FIELD, l) !== nothing, lines) && return :timelab
    lowercase(splitext(path)[2]) in (".phd", ".frd") && return :stable32
    any(l -> occursin(_BIPM_UTC_HEAD, l), first(lines, 5)) && return :bipm_utc
    any(l -> occursin(_RINEX_CLK_HEAD, l), first(lines, 5)) && return :rinex_clk
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
    d = _delim_for(data[1], delim)
    # Counter readings, after a timestamp or none.
    all(data) do l
        f = _split_fields(l, d)
        1 <= length(f) <= 2 && occursin(_SCPI_REAL, f[end]) &&
            (length(f) == 1 || _isnumfield(f[1]) || !isnan(_iso_seconds(f[1])))
    end && return :counter
    first_fields = _split_fields(data[1], d)
    if length(data) >= 2 && !all(_isnumfield, first_fields)
        next = _split_fields(data[2], d)
        length(next) == length(first_fields) && any(_isnumfield, next) && return :table
        return :columns
    end
    return :columns
end

# τ₀: the explicit value, else the header's, else the time column's step.
function _resolve_tau0(tau0::Float64, header::Float64, t::Vector{Float64}, what::String)
    isnan(tau0) || return tau0
    isnan(header) || return header
    length(t) >= 2 && return _column_step(t)
    throw(ArgumentError("read_data: " * what * " has no sample interval; pass tau0"))
end

function _read_timelab(path::AbstractString, tau0::Float64)
    fields = Dict{String,String}()
    v = Float64[]
    block = false
    l = _Lines(path)
    try
        while true
            r = _nextline!(l)
            r === nothing && break
            r = _strip(l.buf, r)
            (isempty(r) || l.buf[first(r)] == UInt8(';')) && continue
            if block
                x = _parse_f64(l.buf, r)
                x === nothing && throw(ArgumentError("read_data: non-numeric TimeLab sample \"" *
                                                     _field_string(l.buf, r) * "\" in " * String(path)))
                push!(v, x)
                continue
            end
            line = _field_string(l.buf, r)
            m = match(_TIMELAB_FIELD, line)
            if m !== nothing
                k = m.captures[2]; val = m.captures[3]
                (k === nothing || val === nothing) || (fields[String(k)] = String(strip(val)))
                continue
            end
            name = String(first(split(line)))
            name == "TIC" || throw(ArgumentError(
                "read_data: TimeLab block \"" * name * "\" is not supported (only TIC phase blocks)"))
            block = true
        end
    finally
        close(l)
    end
    block || throw(ArgumentError("read_data: no TIC data block in " * String(path)))
    scale = tryparse(Float64, get(fields, "Scale Factor", "1"))
    scale === nothing || scale == 1 || (v .*= scale)
    si = tryparse(Float64, get(fields, "Sample Interval", ""))
    t0 = _resolve_tau0(tau0, si === nothing ? NaN : si, Float64[], "this TimeLab file")
    return _Record(v, t0, true, falses(0), NaN)
end

function _read_stable32(path::AbstractString, tau0::Float64, hint::Int, fillgaps::Bool)
    nskip, htau, hstart = _scan_comment_header(path)
    # Record type: the header's "Type:" line, else the extension.
    phase = lowercase(splitext(path)[2]) != ".frd"
    for l in _head_lines(path, nskip)
        m = match(r"^Type:\s*(phase|freq)"i, strip(l))
        m === nothing || (phase = lowercase(something(m.captures[1], "phase")) == "phase")
    end
    hint == 0 || (phase = hint > 0)
    (t, v), ncols = _read_cols(path, nskip, nothing, [1, 0])
    length(v) >= 2 || throw(ArgumentError("read_data: fewer than 2 samples in " * String(path)))
    ncols >= 2 || (t = Float64[])
    # A Stable32 time tag is an MJD (days); a large value with sub-day steps.
    start = hstart
    if length(t) >= 2 && t[1] > 30000 && _median_step(t) < 1
        isnan(start) && (start = _mjd_seconds(t[1]))
        t .*= 86400
    end
    t0 = _resolve_tau0(tau0, htau, t, "this Stable32 file")
    return _Record(_gridded(t, v, fillgaps, phase), t0, phase, start)
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

# A time field: seconds (or days) as a number, else an ISO 8601 timestamp;
# NaN if neither. Sets iso[] when it was a timestamp.
function _time_field(buf::Vector{UInt8}, r::UnitRange{Int}, iso::Base.RefValue{Bool})
    x = _parse_f64(buf, r)
    x === nothing || return x
    iso[] = true
    return _iso_seconds(buf, r)
end

function _read_counter(path::AbstractString, tau0::Float64, nominal::Float64, hint::Int,
                       delim::Char, fillgaps::Bool)
    t = Float64[]; v = Float64[]
    iso = Ref(false)
    d = '\0'; first_line = true
    fields = UnitRange{Int}[]
    l = _Lines(path)
    try
        while true
            r = _nextline!(l)
            r === nothing && break
            _blank_or_comment(l.buf, r) && continue
            if first_line
                dd = _delim_for(_field_string(l.buf, _strip(l.buf, r)), delim)
                d = dd === nothing ? '\0' : dd
                first_line = false
            end
            _fields!(fields, l.buf, r, d)
            x = _parse_f64(l.buf, fields[end])
            x === nothing && throw(ArgumentError("read_data: non-numeric counter reading \"" *
                                                 _field_string(l.buf, _strip(l.buf, r)) * "\" in " * String(path)))
            push!(v, _scpi_nan(x))
            length(fields) >= 2 && push!(t, _time_field(l.buf, fields[1], iso))
        end
    finally
        close(l)
    end
    length(v) >= 2 || throw(ArgumentError("read_data: fewer than 2 readings in " * String(path)))
    (length(t) == length(v) && !any(isnan, t)) || empty!(t)
    start = iso[] && !isempty(t) ? t[1] : NaN
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
    return _Record(_gridded(t, v, fillgaps, phase), t0, phase, start)
end

# Units a column name may carry ("(ns)", "[ps]", "_us", " ppb", "ns/s", "Hz"):
# the scale to SI, NaN for absolute Hz, and whether they are time units
# (phase, +1) or frequency units (-1). Ratios of times come first.
const _UNITS = [(Regex("(^|[\\s(\\[_/-])" * u * "([\\s)\\]]|\$)"), f, k)
                for (u, f, k) in (("ps/s", 1e-12, -1), ("ns/s", 1e-9, -1), ("us/s", 1e-6, -1),
                                  ("µs/s", 1e-6, -1), ("μs/s", 1e-6, -1), ("ms/s", 1e-3, -1),
                                  ("s/s", 1.0, -1), ("ps", 1e-12, 1), ("ns", 1e-9, 1), ("us", 1e-6, 1), ("µs", 1e-6, 1),
                                  ("μs", 1e-6, 1), ("ms", 1e-3, 1), ("ppt", 1e-12, -1), ("ppb", 1e-9, -1),
                                  ("ppm", 1e-6, -1), ("hz", NaN, -1), ("s", 1.0, 1))]

# (scale, kind) of a column name's unit; (1, 0) without one.
function _unit(name::AbstractString)
    n = lowercase(name)
    for (re, f, k) in _UNITS
        occursin(re, n) && return f, k
    end
    return 1.0, 0
end

const _PHASE_NAMES = r"phase|time ?err|\btie\b|pps|time ?offset|time ?diff|\bx\b"i
const _FREQ_NAMES  = r"freq|\bffe\b|osc|fractional|ratio|\by\b"i
const _TIME_NAMES  = r"^(t|mjd|time|timestamp|date|datetime|utc|epoch|elapsed|sec|seconds)\s*([(\[].*)?$"i

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
    length(names) == 2 && time > 0 && return 3 - time
    for re in (_PHASE_NAMES, _FREQ_NAMES)
        hits = [i for i in eachindex(names) if i != time && occursin(re, names[i])]
        length(hits) == 1 && return hits[1]
    end
    throw(ArgumentError("read_data: cannot tell which column holds the clock; pass column= one of " *
                        join(names, ", ")))
end

# A table is read in one pass that keeps only its clock and time columns,
# chosen from the names and the first row.
function _read_table(path::AbstractString, tau0::Float64, column::Int, column_name::String,
                     nominal::Float64, phase_hint::Int, delim::Char, fillgaps::Bool)
    names = String[]
    d = '\0'
    c = 0; time = 0
    v = Float64[]; t = Float64[]
    iso = Ref(false)
    fields = UnitRange{Int}[]
    nrow = 0
    l = _Lines(path)
    try
        while true
            r = _nextline!(l)
            r === nothing && break
            _blank_or_comment(l.buf, r) && continue
            if isempty(names)
                dd = _delim_for(_field_string(l.buf, _strip(l.buf, r)), delim)
                d = dd === nothing ? '\0' : dd
                _fields!(fields, l.buf, r, d)
                names = String[_field_string(l.buf, f) for f in fields]
                continue
            end
            _fields!(fields, l.buf, r, d)
            nrow += 1
            length(fields) == length(names) || throw(ArgumentError(
                "read_data: row " * string(nrow) * " has " * string(length(fields)) * " fields, expected " *
                string(length(names))))
            if nrow == 1
                # The time column by name, else the first one holding ISO timestamps.
                time = something(findfirst(n -> occursin(_TIME_NAMES, n), names),
                                 findfirst(f -> !isnan(_iso_seconds(l.buf, f)), fields), 0)
                c = _pick_column(names, column, column_name, time)
            end
            x = _parse_f64(l.buf, fields[c])
            push!(v, x === nothing ? NaN : x)
            time > 0 && push!(t, _time_field(l.buf, fields[time], iso))
        end
    finally
        close(l)
    end
    nrow >= 2 || throw(ArgumentError("read_data: fewer than 2 rows under the column names in " * String(path)))
    any(isnan, t) && empty!(t)
    # On a time axis a row without a reading is a missing reading, which the
    # grid puts back where it belongs, so a column sampled more sparsely than
    # the table keeps its own step.
    if !isempty(t) && any(isnan, v)
        keep = .!isnan.(v)
        count(keep) >= 2 || throw(ArgumentError("read_data: fewer than 2 readings in column \"" * names[c] * "\""))
        t = t[keep]; v = v[keep]
    end
    start = iso[] && !isempty(t) ? t[1] : NaN
    # MJD by name, or by its range unless the name says seconds; else in the
    # time unit the name gives ("time (ms)"), seconds without one or for
    # ISO 8601 timestamps.
    tscale, tkind = isempty(t) ? (1.0, 0) : _unit(names[time])
    if !isempty(t) && (occursin(r"mjd"i, names[time]) ||
                       (tkind <= 0 && !occursin(r"sec|elapsed"i, names[time]) && _is_mjd(t)))
        start = _mjd_seconds(t[1])
        t .*= 86400
    elseif tkind > 0 && tscale != 1 && !iso[]
        t .*= tscale
    end
    scale, unit_kind = _unit(names[c])
    # A unit in the name fixes the record type, else the name's words do; an
    # override must agree with the unit.
    phase = phase_hint != 0 ? phase_hint > 0 : unit_kind != 0 ? unit_kind > 0 :
            !occursin(_FREQ_NAMES, names[c]) || occursin(_PHASE_NAMES, names[c])
    (phase_hint != 0 && unit_kind != 0 && phase_hint != unit_kind) && throw(ArgumentError(
        "read_data: column \"" * names[c] * "\" is in " * (unit_kind > 0 ? "time" : "frequency") *
        " units and cannot be read as " * (phase_hint > 0 ? "phase" : "frequency")))
    if isnan(scale)            # absolute frequency in Hz
        f0 = _nominal(v, nominal)
        v = (v .- f0) ./ f0
        phase = false
    elseif scale != 1
        v .*= scale
    end
    t0 = _resolve_tau0(tau0, NaN, t, "this table (no time column)")
    return _Record(_gridded(t, v, fillgaps, phase), t0, phase, start)
end

# A BIPM UTC - UTC(k) file: rows of an integer MJD and the difference in ns,
# then uncertainties and notes that vary from row to row ("+/-", a Circular T
# tag), so only the first two fields count; the header, blank lines and the
# DOS end-of-file mark are skipped. The record is phase in seconds on the MJD
# axis (five days apart), a missing row a gap.
function _read_bipm_utc(path::AbstractString, tau0::Float64, fillgaps::Bool)
    t = Float64[]; v = Float64[]
    io = open(path, "r")
    try
        for line in eachline(io)
            f = split(line; keepempty=false)
            length(f) >= 2 || continue
            mjd = tryparse(Int, f[1]); x = tryparse(Float64, f[2])
            (mjd === nothing || x === nothing) && continue
            push!(t, Float64(mjd)); push!(v, x * 1e-9)
        end
    finally
        close(io)
    end
    length(t) >= 2 || throw(ArgumentError("read_data: no UTC - UTC(k) rows in " * String(path)))
    start = _mjd_seconds(t[1])
    t .*= 86400
    t0 = _resolve_tau0(tau0, NaN, t, "this BIPM file")
    return _Record(_gridded(t, v, fillgaps, true), t0, true, start)
end

# A RINEX 3 clock file: "AS G08 2025 08 17 00 00 00.000000  2  bias  sigma"
# records after the header, a satellite (AS) or station (AR) clock bias in
# seconds per line (old IGS files write the exponent with D). The chosen clock's
# biases are the phase record, on the MJD axis; an epoch it is missing from is a gap.
function _read_rinex_clk(path::AbstractString, tau0::Float64, column::Int, column_name::String, fillgaps::Bool)
    seen = Set{String}()
    t = Float64[]; v = Float64[]
    target = column_name
    io = open(path, "r")
    try
        body = false
        for line in eachline(io)
            if !body
                body = occursin("END OF HEADER", line)
                continue
            end
            startswith(line, "AS ") || startswith(line, "AR ") || continue
            f = split(line; keepempty=false)
            length(f) >= 10 || continue
            name = String(f[2])
            if !(name in seen)
                push!(seen, name)
                isempty(target) && length(seen) == max(column, 1) && (target = name)
            end
            (isempty(target) ? false : lowercase(name) == lowercase(target)) || continue
            y = tryparse(Int, f[3]); mo = tryparse(Int, f[4]); d = tryparse(Int, f[5])
            h = tryparse(Int, f[6]); mi = tryparse(Int, f[7])
            sec = tryparse(Float64, f[8]); x = tryparse(Float64, replace(f[10], 'D' => 'E', 'd' => 'e'))
            (y === nothing || mo === nothing || d === nothing || h === nothing || mi === nothing ||
             sec === nothing || x === nothing) && continue
            days = Dates.value(Date(y, mo, d) - Date(1858, 11, 17))
            push!(t, days + (3600h + 60mi + sec) / 86400); push!(v, x)
        end
    finally
        close(io)
    end
    if length(t) < 2
        n = string(length(seen))
        throw(ArgumentError("read_data: " * (isempty(seen) ? "no clock records" :
            isempty(target) ? "column=" * string(column) * " exceeds the " * n * " clocks" :
            "no records for clock \"" * target * "\" among the " * n * " clocks") * " in " * String(path)))
    end
    start = _mjd_seconds(t[1])
    t .*= 86400
    t0 = _resolve_tau0(tau0, NaN, t, "this RINEX clock file")
    return _Record(_gridded(t, v, fillgaps, true), t0, true, start)
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
    elseif fmt === :bipm_utc
        _read_bipm_utc(path, tau0, fillgaps)
    elseif fmt === :rinex_clk
        _read_rinex_clk(path, tau0, column, column_name, fillgaps)
    elseif fmt === :columns
        v, t0, start, phase = _ingest(path; tau0=isnan(tau0) ? nothing : tau0, time_col=1, value_col=2,
                                      header=0, delim=delim == '\0' ? nothing : delim, scaling=1.0,
                                      detrend=:none, fillgaps=fillgaps, hint=hint)
        _Record(v, t0, phase, falses(0), start)
    else
        throw(ArgumentError("read_data: unknown format :" * String(fmt) * "; use one of " * join(map(String, INPUT_FORMATS), ", ")))
    end
    scaling == 1 || (r.v .*= scaling)
    phase = hint == 0 ? r.phase : hint > 0
    return _Record(r.v, r.tau0, phase, r.filled, r.start)
end

"""
    read_data(path; format=:auto, tau0=nothing, column=nothing, nominal=nothing,
              kind=:auto, scaling=1.0, delim=nothing,
              fillgaps=false) → PhaseData or FrequencyData

Read a measurement file in one of the [`INPUT_FORMATS`](@ref) and return
its record, phase or frequency as the file has it. The format is found by
[`detect_format`](@ref) unless `format` names it.

# Formats

- `:columns`: numeric columns, read as [`read_phase`](@ref) reads them by
  default. The record is phase, or frequency when the median size of the
  samples is under ``10^{-11}``. Phase samples can be that small, as the
  residuals of a low-noise comparison are: pass `kind=:phase` for them.
- `:stable32`: a Stable32 `.phd` (phase) or `.frd` (frequency) file. A
  `Type:` line in its header overrides the extension. ``τ_0`` is that of
  the `Tau:` header line, else the step of the time tags. With two columns
  or more, the first holds the time tags and the last the samples. Time
  tags above 30000 at steps under one day are read as MJD.
- `:timelab`: a TimeLab `.tim` file. Its `TIC` block holds the phase
  samples in seconds, which are multiplied by the `Scale Factor` of the
  header, and `Sample Interval` is ``τ_0``. A file with another kind of
  block throws an error.
- `:counter`: readings of a counter in the format `+1.00000000012345E+07`,
  one on each line, alone or after a timestamp in seconds or ISO 8601.
  Readings of median size under 1 are taken as time intervals in seconds
  and make a phase record. Larger ones are taken as frequencies in hertz
  and make the frequency record ``(f - f_0)/f_0`` about `nominal`. A
  reading of `9.91E+37`, the SCPI value for not a number, becomes `NaN`.
  Without timestamps `tau0` must be given: it is the gate time.
- `:table`: columns under a line of names, as in a CSV export, described
  below.
- `:bipm_utc`: a `UTC - UTC(k)` file of the BIPM, with one row every five
  days that holds the MJD and the difference in nanoseconds. The record is
  phase in seconds.
- `:rinex_clk`: a RINEX clock file. The record is the clock bias, in
  seconds, of one satellite (`AS` records) or station (`AR` records),
  chosen by `column`. Exponents written with `D` are read. A compressed
  file must be unpacked first.

!!! warning "Gaps in Stable32 files"
    Stable32 marks a missing sample with a value of zero
    ([Stable32 manual, "Gaps & Outliers"](@cite hamilton-2008-stable32-manual)).
    This reader takes a zero as a sample. Set the zeros that mark gaps to
    `NaN` in the record after reading.

# Named columns

In a `:table` file the time column is the first whose name is `t`, `time`,
`timestamp`, `date`, `datetime`, `utc`, `epoch`, `elapsed`, `sec`,
`seconds` or `mjd`, with or without a unit in brackets after it, else the
first column that holds ISO 8601 timestamps. Numbers in it are read as
seconds, or in the unit that its name gives when that is one of the units
of time listed below, as in `time (ms)`: minutes and hours are not among
them. They are read as days when the name contains `mjd`, or when they meet
the MJD conditions of [`read_phase`](@ref) and the name has neither `sec`,
`elapsed` nor one of those units.

The clock column is `column`. By default it is the other column of a file
with a time column and one more, else the only column whose name holds
`phase`, `time err`, `tie`, `pps`, `time offset`, `time diff` or `x`, else
the only one whose name holds `freq`, `ffe`, `osc`, `fractional`, `ratio`
or `y`. An error lists the names when none of these decides.

A unit in the name of the clock column sets the type of the record and
converts the samples. The units of time `ps`, `ns`, `us`, `µs`, `ms` and
`s` make it phase, in seconds. The units `ppt`, `ppb` and `ppm` and the
ratios of times from `ps/s` to `s/s` make it fractional frequency, and so
does `Hz`, about `nominal`. Without a unit the words of the name decide:
the record is frequency when the name holds one of the frequency words and
none of the phase words, and phase otherwise.

A row without a number in the clock column is a missing reading. With a
time column such rows are left out before the grid is fitted, so a column
logged at a longer interval than the table keeps its own interval.

# Timestamps

Where a file has timestamps, they place the samples as in
[`read_phase`](@ref): they give ``τ_0`` when the keyword and the header do
not, and a reading missing from the file leaves a `NaN` sample. ISO 8601
timestamps and MJD values give the record its [`start_time`](@ref).

# Keywords

- `format=:auto`: one of [`INPUT_FORMATS`](@ref), or `:auto` to detect it.
- `tau0`: sample interval in seconds. By default that of the header of the
  file, else that of its timestamps. An error is thrown when the file has
  neither.
- `column`: the clock to read. In a `:table` file, the number of its
  column, counted from 1, or its name. In a `:rinex_clk` file, its position
  in order of first appearance, or its name, such as `"G08"`; by default
  the first clock. Names are matched without regard to case.
- `nominal`: nominal frequency in hertz, for readings in hertz. By default
  the mean reading rounded to four significant figures.
- `kind=:auto`: `:phase` or `:frequency` sets the type of the record
  whatever the file suggests. Counter readings are then read as seconds
  under `:phase` and as hertz under `:frequency`, and a table column whose
  unit is of the other kind throws an error.
- `scaling=1.0`: factor applied to the samples after any conversion of
  units.
- `delim`: character that separates the fields of a `:columns`, `:counter`
  or `:table` file. By default a `:columns` file follows the rule of
  [`read_phase`](@ref). In the other two a comma, else a tab, separates the
  fields when the first line holds one, and runs of blanks otherwise.
- `fillgaps=false`: `true` fills the missing samples with
  [`fillgaps`](@ref), with the reservations stated there: the readings
  that the timestamps show to be absent, and those that are `NaN`.

# Examples

```jldoctest
julia> using SigmaTau

julia> path = joinpath(mktempdir(), "gpsdo.csv");

julia> write(path, "t,PPS offset (ns),temp\\n0,12.5,31\\n1,13.5,31\\n2,11.0,31\\n");

julia> detect_format(path)
:table

julia> p = read_data(path);

julia> p.x ≈ [12.5e-9, 13.5e-9, 11.0e-9]
true
```

```julia
p = read_data("trace.tim")                      # TimeLab, τ₀ from its header
y = read_data("counter.csv"; tau0=1.0)          # readings in hertz at a 1 s gate
p = read_data("utc-nist")                       # BIPM UTC − UTC(NIST)
p = read_data("clocks.clk"; column="G08")       # one satellite clock
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
    return r.phase ? PhaseData(r.v, r.tau0; source=String(path), start=r.start) :
                     FrequencyData(r.v, r.tau0; source=String(path), start=r.start)
end

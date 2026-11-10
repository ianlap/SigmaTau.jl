# Instrument file formats: detection, and readers for the files that are not
# plain numeric columns (TimeLab .tim, Stable32 .phd/.frd, frequency-counter
# logs, logs with named columns). Plain columns go through read_phase /
# read_frequency.

"""
Input formats [`read_data`](@ref) understands, as [`detect_format`](@ref)
names them.
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

# Keysight / Agilent SCPI real format: "+1.00000000012345E+07".
const _SCPI_REAL = r"^[+-]\d\.\d+E[+-]\d+$"
# SCPI's "not a number" reading (9.91E+37), written for a missed gate.
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

Name the layout of a data file from its extension and first lines, one of
[`INPUT_FORMATS`](@ref):

- `:timelab`: a TimeLab `.tim` measurement log (also what the Microchip
  53100A, 5120A/5125A and TimePod write through TimeLab): typed header
  fields such as `DBL 0x… "Sample Interval" 0.1`, then a `TIC n` block of
  phase samples in seconds.
- `:stable32`: a Stable32 `.phd` (phase) or `.frd` (frequency) file. As
  in Stable32, a value of zero is a gap (`NaN`), except the first and last
  point of phase data.
- `:counter`: frequency-counter readings in SCPI real format
  (`+1.00000000012345E+07`, as the Keysight 53230A and 53220A log them),
  one per line, optionally after a timestamp (seconds or ISO 8601).
- `:table`: numeric columns under a line of column names (a CSV export,
  most GPSDO and logging-script output).
- `:bipm_utc`: a BIPM `UTC - UTC(k)` file from the time department's
  `utclab` directory (`MJD  [UTC-UTC(NIST)]/ns  uA/ns …`): one row every five
  days with the difference in nanoseconds, read as phase in seconds.
- `:rinex_clk`: a RINEX 3 clock file (`.clk`, IGS and CODE clock products):
  `AS` satellite and `AR` station clock biases in seconds, many clocks on one
  epoch grid. `column` (1-based, in order of first appearance) or `column_name`
  (`G08`, `ALGO00CAN`) picks the clock, by default the first. Old files'
  `D` exponents (`-.123D-03`) are read too. Unzip `.gz` files first.
- `:columns`: plain numeric columns, read by [`read_phase`](@ref) and
  [`read_frequency`](@ref). `read_data` takes them as phase, unless the
  readings' median size is under 1e-11, which only fractional frequency has.

`delim` is the field delimiter, by default a comma or tab when the first
data line holds one, else whitespace.
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
    # Stable32 writes a gap as a zero; only the end points of phase data may be zero.
    for i in eachindex(v)
        v[i] == 0 && !(phase && (i == 1 || i == length(v))) && (v[i] = NaN)
    end
    # A Stable32 time tag is an MJD (days); a large value with sub-day steps.
    start = hstart
    if length(t) >= 2 && t[1] > 30000 && _median_step(t) < 1
        isnan(start) && (start = _mjd_seconds(t[1]))
        t .*= 86400
    end
    t0 = _resolve_tau0(tau0, htau, t, "this Stable32 file")
    return _Record(_gridded(t, v, fillgaps), t0, phase, start)
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
    return _Record(_gridded(t, v, fillgaps), t0, phase, start)
end

# Units a column name may carry ("(ns)", "[ps]", "_us", " ppb", "Hz"): the
# scale to SI, NaN for absolute Hz, and whether they are time units (phase,
# +1) or frequency units (-1).
const _UNITS = [(Regex("(^|[\\s(\\[_/-])" * u * "([\\s)\\]]|\$)"), f, k)
                for (u, f, k) in (("ps", 1e-12, 1), ("ns", 1e-9, 1), ("us", 1e-6, 1), ("µs", 1e-6, 1),
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
    # MJD by name, or by its range unless the name says seconds.
    if !isempty(t) && (occursin(r"mjd"i, names[time]) ||
                       (_unit(names[time])[2] <= 0 && !occursin(r"sec|elapsed"i, names[time]) && _is_mjd(t)))
        start = _mjd_seconds(t[1])
        t .*= 86400
    end
    phase = phase_hint == 0 ? !occursin(_FREQ_NAMES, names[c]) || occursin(_PHASE_NAMES, names[c]) :
            phase_hint > 0
    scale, unit_kind = _unit(names[c])
    # A unit in the name fixes the record type; an override must agree.
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
    return _Record(_gridded(t, v, fillgaps), t0, phase, start)
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
    return _Record(_gridded(t, v, fillgaps), t0, true, start)
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
    return _Record(_gridded(t, v, fillgaps), t0, true, start)
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
        v, t0, start = _ingest(path; tau0=isnan(tau0) ? nothing : tau0, time_col=1, value_col=2,
                               header=0, delim=nothing, scaling=1.0, detrend=:none, fillgaps=fillgaps)
        # Unnamed readings of median size under 1e-11 are fractional
        # frequency (or a ratio): a time error that small is not a clock.
        f = filter(isfinite, v)
        tiny = !isempty(f) && _median_nonan(abs.(f)) < 1e-11
        _Record(v, t0, hint > 0 || (hint == 0 && !tiny), falses(0), start)
    else
        throw(ArgumentError("read_data: unknown format :" * String(fmt) * "; use one of " * join(map(String, INPUT_FORMATS), ", ")))
    end
    scaling == 1 || (r.v .*= scaling)
    phase = hint == 0 ? r.phase : hint > 0
    return _Record(r.v, r.tau0, phase, r.filled, r.start)
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
  (TimeLab `Sample Interval`, Stable32 `Tau:`) or the grid step of the
  time column (seconds, MJD days, or ISO 8601 timestamps). A column named
  `MJD` is read in days, and so is one whose values lie between MJD 30000
  and 90000 with steps under 0.05 d and more than four decimals, unless its
  name says seconds.
- `column`: the clock column of a `:table` file, by 1-based index or by
  name. By default the only column besides the time column, else the one
  whose name says phase / time error / PPS, else frequency.
- `nominal`: nominal frequency in Hz for readings in Hz (counter logs, a
  table column named `… (Hz)`); converted to fractional frequency
  `(f − nominal)/nominal`. By default the mean reading rounded to four
  significant figures.
- `kind`: `:phase` or `:frequency` overrides the record type.
- `scaling`: multiplies the samples after any unit conversion from a
  column name (`(ns)`, `[ps]`, `ppb` …).
- `delim`: field delimiter of counter and table files; by default a comma
  or tab when the first data line holds one, else whitespace.
- `fillgaps`: fill the missing samples (and any `NaN` reading) with
  [`fillgaps`](@ref).

A time column (seconds, MJD days or ISO 8601 timestamps) places the samples
as in [`read_phase`](@ref): a missing reading leaves a `NaN` sample in its
slot, which the gap-tolerant estimators skip, and a record with none is read
as written. A table row with no reading in the clock column is a missing
reading too, so a column sampled more sparsely than the table keeps its own
step. ISO 8601 and MJD timestamps also set the record's
[`start_time`](@ref). Files are read line by line, keeping only the clock
and time columns.

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
    return r.phase ? PhaseData(r.v, r.tau0; source=String(path), start=r.start) :
                     FrequencyData(r.v, r.tau0; source=String(path), start=r.start)
end

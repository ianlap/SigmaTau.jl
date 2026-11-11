# Readers for files of phase or frequency samples in numeric columns, with
# optional scaling, detrending and gap filling, and the time grid that puts
# timestamped samples in place.

# The delimiter of a file whose first data line is `buf[r]`: tab for a `.tsv`
# whose line holds one, else a comma when the line holds one, else runs of
# blanks (tabs too), so a comma-separated `.txt` and a blank-separated `.csv`
# both read.
function _auto_delim(path::AbstractString, buf::Vector{UInt8}, r::UnitRange{Int})
    lowercase(splitext(path)[2]) == ".tsv" && UInt8('\t') in view(buf, r) && return '\t'
    return UInt8(',') in view(buf, r) ? ',' : '\0'
end

# Columns `cols` of a numeric text file (0 is the last column), read line by
# line so a record of any length costs only those columns. `header` lines are
# skipped, then blank and `#` lines wherever they appear (Stable32's
# "# Header End", notes written mid-record); every row must have the same
# number of fields. The delimiter is `delim`, else the one picked from the
# first data line (see `_auto_delim`). Returns the columns and the field count.
function _read_cols(path::AbstractString, header::Int, delim::Union{Char,Nothing}, cols::Vector{Int})
    isfile(path) || throw(ArgumentError("read: file not found: " * String(path)))
    d = delim === nothing ? '\0' : delim
    out = [Float64[] for _ in cols]
    ncols = -1
    lineno = 0
    fields = UnitRange{Int}[]
    l = _Lines(path)
    try
        while true
            r = _nextline!(l)
            r === nothing && break
            lineno += 1
            (lineno <= header || _blank_or_comment(l.buf, r)) && continue
            ncols < 0 && delim === nothing && (d = _auto_delim(path, l.buf, r))
            _fields!(fields, l.buf, r, d)
            if ncols < 0
                ncols = length(fields)
                for c in cols
                    c <= ncols || throw(ArgumentError("read: column " * string(c) * " exceeds the " *
                                                      string(ncols) * " columns of " * String(path)))
                end
            elseif length(fields) != ncols
                throw(ArgumentError("read: line " * string(lineno) * " of " * String(path) * " has " *
                                    string(length(fields)) * " columns, expected " * string(ncols)))
            end
            for (k, c) in enumerate(cols)
                f = fields[c == 0 ? ncols : c]
                v = _parse_f64(l.buf, f)
                v === nothing && throw(ArgumentError("read: non-numeric entry \"" * _field_string(l.buf, f) *
                                                     "\" on line " * string(lineno) * " of " * String(path)))
                push!(out[k], v)
            end
        end
    finally
        close(l)
    end
    return out, max(ncols, 0)
end

# The field count of the first data line (0 for no data).
function _ncols(path::AbstractString, header::Int, delim::Union{Char,Nothing})
    fields = UnitRange{Int}[]
    lineno = 0
    l = _Lines(path)
    try
        while true
            r = _nextline!(l)
            r === nothing && return 0
            lineno += 1
            (lineno <= header || _blank_or_comment(l.buf, r)) && continue
            _fields!(fields, l.buf, r, delim === nothing ? _auto_delim(path, l.buf, r) : delim)
            return length(fields)
        end
    finally
        close(l)
    end
end

# Median of a NaN-free vector. Hand-rolled instead of `Statistics.median`,
# whose NaN error path prints the whole vector and so blocks
# `juliac --trim` builds. Callers must filter or reject NaN first.
function _median_nonan(v::AbstractVector{<:Real})
    s = sort!(Float64.(v))
    n = length(s)
    n == 0 && throw(ArgumentError("median of an empty collection"))
    return isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

# Median sample spacing of a time column.
function _median_step(t::AbstractVector{<:Real})
    dt = diff(t)
    any(isnan, dt) && throw(ArgumentError("read: time column contains NaN"))
    return _median_nonan(dt)
end

# The even grid of a time column: the slot of each sample and the grid step.
# The step is first the mean of the regular steps (under 1.5 times the
# shortest, so not the gaps), whose sum telescopes, so timestamps rounded to
# a coarse resolution (MJD to 1e-6 day) still average to the true step. The
# step is then the shortest round value of it (1 to 6 significant figures)
# that the scatter of the timestamps cannot tell from their least-squares
# step on the slots it gives (see `_step_tolerance`) and that keeps every
# timestamp within a quarter step of its slot over the whole record, else
# the span over the slots.
# Every step must be within a quarter of a whole number of grid steps, so jittered timestamps keep their order
# and a missing reading leaves an empty slot. No slots (and a NaN step) when
# the column is no time axis: not increasing, or steps that fit no grid.
function _time_grid(t::AbstractVector{<:Real})
    n = length(t)
    none = (Int[], NaN)
    n >= 2 || return none
    for i in 2:n
        t[i] > t[i - 1] || return none
    end
    dt = diff(t)
    lo = minimum(dt)
    s = 0.0; m = 0
    for d in dt
        d < 1.5 * lo && (s += d; m += 1)
    end
    est = s / m
    for digits in 1:6
        step = round(est; sigdigits=digits)
        slots = _slots(dt, step)
        isempty(slots) && continue
        fit, tol = _step_tolerance(t, slots)
        abs(step - fit) <= min(tol, 0.01 * fit) && _on_grid(t, slots, step) && return slots, step
    end
    slots = _slots(dt, est)
    isempty(slots) && return none
    slots[end] > 1 && (est = (t[end] - t[1]) / (slots[end] - 1))
    slots = _slots(dt, est)
    return isempty(slots) ? none : (slots, est)
end

# The least-squares step of the timestamps on their slots, and the most a
# round step may differ from it: the slope error that timestamp errors of
# twice the largest residual can cause at worst, so coarse timestamps admit
# a round step and exact ones only the exact step (to 1e-9 of it).
function _step_tolerance(t::AbstractVector{<:Real}, slots::Vector{Int})
    n = length(slots)
    xb = sum(slots) / n; tb = 0.0
    for i in 1:n
        tb += t[i] - t[1]
    end
    tb /= n
    sxx = 0.0; sxy = 0.0; sax = 0.0
    for i in 1:n
        sxx += (slots[i] - xb)^2
        sax += abs(slots[i] - xb)
        sxy += (slots[i] - xb) * (t[i] - t[1] - tb)
    end
    b = sxy / sxx
    r = 0.0
    for i in 1:n
        r = max(r, abs(t[i] - t[1] - tb - b * (slots[i] - xb)))
    end
    return b, max(2r * sax / sxx, 1e-9 * abs(b))
end

# The slot of each sample on a grid of `step`, each step rounded to a whole
# number of grid steps; empty when a step is off the grid by over a quarter.
function _slots(dt::Vector{Float64}, step::Float64)
    slots = Vector{Int}(undef, length(dt) + 1)
    slots[1] = 1
    for i in eachindex(dt)
        k = dt[i] / step
        k < 2.0^31 || throw(ArgumentError("read: a gap of " * string(k) * " samples in the time column"))
        abs(k - round(k)) <= 0.25 || return Int[]
        slots[i + 1] = slots[i] + round(Int, k)
    end
    return slots
end

# Every timestamp within a quarter step of its slot on the grid from t[1].
function _on_grid(t::AbstractVector{<:Real}, slots::Vector{Int}, step::Float64)
    for i in eachindex(slots)
        abs(t[i] - t[1] - (slots[i] - 1) * step) <= 0.25 * step || return false
    end
    return true
end

# τ₀ from a time column: its grid step, or the median step when it fits no grid.
function _column_step(t::AbstractVector{<:Real})
    step = _time_grid(t)[2]
    return isnan(step) ? _median_step(t) : step
end

# The samples on the grid of their time column, NaN in the slots of missing
# readings (what the gap-tolerant estimators skip); `v` itself when nothing
# is missing. With `fillgaps` the empty slots, and any NaN reading, are
# filled. Returns the samples and the mask of filled ones.
function _gridded(t::Vector{Float64}, v::Vector{Float64}, fillgaps::Bool, phase::Bool)
    s = _time_grid(t)[1]
    g = v
    if !isempty(s) && s[end] != length(v)
        g = fill(NaN, s[end])
        g[s] .= v
    end
    fillgaps || return g, falses(0)
    return _fillgaps(g, phase, :auto), isnan.(g)
end

# Seconds since 2000-01-01T00:00Z, the epoch of `start`, of an MJD.
_mjd_seconds(mjd::Float64) = (mjd - 51544.0) * 86400.0

# A time column in MJD days: values from MJD 30000 to 90000 (1941 to 2105),
# steps under 0.05 d (a sample interval under 72 min), and digits past the
# fourth decimal, which a sub-day MJD needs (1 s is 1.16e-5 d) and seconds
# timestamps in that range rarely carry (43200.1, 43200.2 … at 10 Hz).
function _is_mjd(t::Vector{Float64})
    (length(t) >= 2 && all(isfinite, t) && 30000 < t[1] && t[end] < 90000) || return false
    _median_step(t) < 0.05 || return false
    return any(v -> abs(v * 1e4 - round(v * 1e4)) > 1e-2, t)
end

# Header lines to skip before the data, and the sample interval and start
# time the header states (NaN when it does not):
#
# - a block of leading `#` comment lines, as `save` writes them, with
#   `# tau0: <seconds>` and `# start: <ISO 8601>`;
# - a Stable32 export header: text lines such as `Tau: 1.0` closed by a
#   `# Header End` line within the first 64 lines, all skipped.
function _scan_comment_header(path::AbstractString)
    ncomment = 0
    tau0 = NaN; start = NaN; s32_tau = NaN
    s32_end = 0
    leading = true
    lineno = 0
    # Plain open/close rather than a do-block: the closure would box the
    # reassigned locals.
    io = open(path, "r")
    try
        for line in eachline(io)
            lineno += 1
            if leading && startswith(line, '#')
                ncomment += 1
                m = match(r"^#\s*(tau0|start):\s*(\S+)", line)
                if m !== nothing
                    key = something(m.captures[1], ""); val = something(m.captures[2], "")
                    if key == "tau0"
                        tau0 = something(tryparse(Float64, val), NaN)
                    else
                        start = _iso_seconds(val)
                    end
                end
                # A note such as "# data interval 1.0 s" (allantools' sample files).
                m = match(r"^#\s*(?:data|sample) interval:?\s+(\S+)\s*s\s*$"i, line)
                m === nothing || isnan(tau0) &&
                    (tau0 = something(tryparse(Float64, something(m.captures[1], "")), NaN))
                continue
            end
            leading = false
            lineno > 64 && break
            st = strip(line)
            if st == "# Header End"
                s32_end = lineno
                break
            end
            m = match(r"^Tau:\s*(\S+)", st)
            m === nothing || (s32_tau = something(tryparse(Float64, something(m.captures[1], "")), NaN))
        end
    finally
        close(io)
    end
    s32_end > 0 && return s32_end, (isnan(tau0) ? s32_tau : tau0), start
    return ncomment, tau0, start
end

# Shared by read_phase and read_frequency: the samples after scaling, the
# time grid, gap filling and detrending, the resolved τ₀ and the start time.
function _ingest(path::AbstractString; tau0::Union{Real,Nothing}, time_col::Int, value_col::Int,
                 header::Int, delim::Union{Char,Nothing}, scaling::Real, detrend::Symbol,
                 fillgaps::Bool, hint::Int)
    isfile(path) || throw(ArgumentError("read: file not found: " * String(path)))
    nskip, header_tau0, start = _scan_comment_header(path)
    # A Stable32 header block is itself non-comment header lines, so an
    # explicit `header` that already covers it is not added on top.
    stable32 = nskip > 0 && !startswith(readline(path), '#')
    value_col >= 1 || throw(ArgumentError("read: value_col must be ≥ 1"))
    skip = stable32 ? max(header, nskip) : header + nskip
    # The value column is never its own time axis, and a one-column file read
    # with the default columns is its values alone.
    time_col == value_col && (time_col = 0)
    if time_col == 1 && value_col == 2 && _ncols(path, skip, delim) == 1
        time_col = 0; value_col = 1
    end
    have_time = time_col >= 1
    cols, _ = _read_cols(path, skip, delim, have_time ? [time_col, value_col] : [value_col])
    v = cols[end]
    t = have_time ? cols[1] : Float64[]
    if _is_mjd(t)
        isnan(start) && (start = _mjd_seconds(t[1]))
        t .*= 86400
    end
    scaling == 1 || (v .*= scaling)
    # τ₀: explicit kwarg > `# tau0:` comment header > time column.
    resolved_tau0 = if tau0 !== nothing
        Float64(tau0)
    elseif !isnan(header_tau0)
        header_tau0
    elseif length(t) > 1
        _column_step(t)
    else
        throw(ArgumentError("read: tau0 must be supplied when no time column is present"))
    end
    phase = _unnamed_phase(v, hint)
    v, _ = _gridded(t, v, fillgaps, phase)
    detrend === :none || _detrend_core!(v, detrend)
    return v, resolved_tau0, start, phase
end

# Whether readings are phase: as `hint` says (1 phase, −1 frequency), else by
# size. Unnamed readings of median size under 1e-11 are fractional frequency
# (or a ratio): a time error that small is not a clock.
function _unnamed_phase(v::Vector{Float64}, hint::Int)
    hint == 0 || return hint > 0
    f = filter(isfinite, v)
    return isempty(f) || _median_nonan(abs.(f)) >= 1e-11
end

"""
    read_phase(path; tau0=nothing, time_col=1, value_col=2, header=0, delim=nothing,
               scaling=1.0, detrend=:none, fillgaps=false) → PhaseData

Read a text file of numeric columns as a [`PhaseData`](@ref) record. Each
line holds one sample, in seconds, in column `value_col`, and its time in
column `time_col`. A file with a single column is read as the samples
alone. The record carries `path` as its `source`.

The file is read line by line and only the two columns are kept, so the
memory used is that of the samples. [`read_data`](@ref) reads instrument
files and tables with named columns.

# Header

Lines that begin with `#` at the head of the file are skipped, and three of
them are read:

- `# tau0: 10.0` gives the sample interval in seconds;
- `# start: 2026-11-08T21:41:03Z` gives the UTC time of the first sample
  ([`start_time`](@ref));
- `# data interval 1.0 s`, or `# sample interval: 1.0 s`, gives the sample
  interval when no `# tau0:` line does.

[`save`](@ref) writes the first two, so a saved record is read back without
keywords. A Stable32 export begins with lines of text closed by
`# Header End`. They are skipped, and their `Tau:` line gives the sample
interval. Blank lines and lines that begin with `#` are skipped throughout
the file.

# Time column

The sample interval is the keyword `tau0`, else that of the header, else
the step of the time column. An error is thrown when the file gives none
and the keyword is absent.

The time column also places the samples. Its steps are fitted to an even
grid, and a reading that is missing from the file leaves a `NaN` sample in
its place (see [Timestamps and gaps](@ref)). A file with no reading missing
is read as written.

The grid step starts from the mean of the steps shorter than 1.5 times the
shortest. It is rounded to the fewest significant figures, at most six, for
which two conditions hold: the rounded step differs from the least-squares
step of the timestamps by no more than the uncertainty of that step, which
their scatter about the fit sets, nor by more than 1 %, and every timestamp
lies within a quarter step of its place. When no rounding passes, the step
is the span of the timestamps over the number of steps. Timestamps that
jitter, or that are rounded to a coarse resolution, therefore give a round
``τ_0``, and exact timestamps keep their own step whatever the length of
the record.

A column that is no time axis leaves the samples in file order: one that
does not increase from line to line, or whose steps are not whole multiples
of one step to within a quarter of it. ``τ_0`` is then the median step,
which must be positive, unless the keyword or the header gives it.

A time column is read as MJD when its values lie between 30000 and 90000,
its median step is under 0.05 d and it has digits past the fourth decimal.
It is converted to seconds and gives the record its [`start_time`](@ref).
MJD values at longer steps are read as seconds, so that daily data need
`tau0=86400`. Timestamps in seconds that meet the three conditions are
taken for MJD: pass `tau0` for them, and disregard the start time.

# Keywords

- `tau0`: sample interval in seconds.
- `time_col=1`: column of the timestamps, counted from 1. With `0`, or with
  `time_col == value_col`, the file is read without a time column.
- `value_col=2`: column of the samples.
- `header=0`: number of lines to skip after the leading `#` lines, such as
  a line of column names.
- `delim`: character that separates the fields. By default a tab in a
  `.tsv` file whose first line of data holds one, else a comma when that
  line holds one, else runs of blanks. A field may be enclosed in double
  quotes.
- `scaling=1.0`: factor applied to the samples, `1e-9` for a file in
  nanoseconds.
- `detrend=:none`: trend to remove from the samples, one of the methods of
  [`detrend`](@ref).
- `fillgaps=false`: `true` fills the missing samples with
  [`fillgaps`](@ref), with the reservations stated there: the readings
  that the time column shows to be absent, and those written as `NaN`.

# Examples

```jldoctest
julia> using SigmaTau

julia> path = joinpath(mktempdir(), "phase.txt");

julia> write(path, "0 1.0e-9\\n10 1.2e-9\\n20 1.1e-9\\n40 1.5e-9\\n");

julia> p = read_phase(path);

julia> p.tau0
10.0

julia> print(p.x)
[1.0e-9, 1.2e-9, 1.1e-9, NaN, 1.5e-9]
```

```julia
p = read_phase("phase.csv"; scaling=1e-9)         # samples in nanoseconds
p = read_phase("samples.txt"; tau0=0.1)           # one column, no timestamps
p = read_phase("log.txt"; value_col=4, header=1)  # column 4, below the names
```
"""
function read_phase(path::AbstractString; tau0::Union{Real,Nothing}=nothing, time_col::Int=1,
                    value_col::Int=2, header::Int=0, delim::Union{Char,Nothing}=nothing,
                    scaling::Real=1.0, detrend::Symbol=:none, fillgaps::Bool=false)
    v, t0, start, _ = _ingest(path; tau0=tau0, time_col=time_col, value_col=value_col, header=header,
                              delim=delim, scaling=scaling, detrend=detrend, fillgaps=fillgaps, hint=1)
    return PhaseData(v, t0; source=String(path), start=start)
end

"""
    read_frequency(path; tau0=nothing, time_col=1, value_col=2, header=0, delim=nothing,
                   scaling=1.0, detrend=:none, fillgaps=false) → FrequencyData

Read a text file of numeric columns as a [`FrequencyData`](@ref) record of
fractional-frequency samples. The header, the time column and the keywords
are those of [`read_phase`](@ref).

# Examples

```julia
y = read_frequency("frequency.csv")
y = read_frequency("ffe_ppb.txt"; scaling=1e-9)    # samples in parts per billion
```
"""
function read_frequency(path::AbstractString; tau0::Union{Real,Nothing}=nothing, time_col::Int=1,
                        value_col::Int=2, header::Int=0, delim::Union{Char,Nothing}=nothing,
                        scaling::Real=1.0, detrend::Symbol=:none, fillgaps::Bool=false)
    v, t0, start, _ = _ingest(path; tau0=tau0, time_col=time_col, value_col=value_col, header=header,
                              delim=delim, scaling=scaling, detrend=detrend, fillgaps=fillgaps, hint=-1)
    return FrequencyData(v, t0; source=String(path), start=start)
end

# Readers for phase / frequency data files, with optional scaling, detrending
# and gap filling, and the time grid that puts timestamped samples in place.

# Columns `cols` of a numeric text file (0 is the last column), read line by
# line so a record of any length costs only those columns. `header` lines are
# skipped, then blank and `#` lines wherever they appear (Stable32's
# "# Header End", notes written mid-record); every row must have the same
# number of fields. The delimiter is `delim`, else `,` for `.csv`, tab for
# `.tsv`, else runs of blanks. Returns the columns and the field count.
function _read_cols(path::AbstractString, header::Int, delim::Union{Char,Nothing}, cols::Vector{Int})
    isfile(path) || throw(ArgumentError("read: file not found: " * String(path)))
    ext = lowercase(splitext(path)[2])
    d = delim !== nothing ? delim : ext == ".csv" ? ',' : ext == ".tsv" ? '\t' : '\0'
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
# The step is the median of the regular steps (under 1.5 times the shortest,
# so not the gaps), and every step must be within a quarter of a whole number
# of grid steps, so jittered timestamps keep their order and a missing
# reading leaves an empty slot. No slots (and a NaN step) when the column is
# no time axis: not increasing, or steps that fit no grid.
function _time_grid(t::AbstractVector{<:Real})
    n = length(t)
    none = (Int[], NaN)
    n >= 2 || return none
    for i in 2:n
        t[i] > t[i - 1] || return none
    end
    dt = diff(t)
    lo = minimum(dt)
    step = _median_nonan(filter(d -> d < 1.5 * lo, dt))
    slots = Vector{Int}(undef, n)
    slots[1] = 1
    for i in 2:n
        k = dt[i - 1] / step
        k < 2.0^31 || throw(ArgumentError("read: a gap of " * string(k) * " samples in the time column"))
        abs(k - round(k)) <= 0.25 || return none
        slots[i] = slots[i - 1] + round(Int, k)
    end
    return slots, step
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
function _gridded(t::Vector{Float64}, v::Vector{Float64}, fillgaps::Bool)
    s = _time_grid(t)[1]
    g = v
    if !isempty(s) && s[end] != length(v)
        g = fill(NaN, s[end])
        g[s] .= v
    end
    fillgaps || return g, falses(0)
    f, mask = _howe_fillgaps_core(g)
    return f, BitVector(mask)
end

# Seconds since 2000-01-01T00:00Z, the epoch of `start`, of an MJD.
_mjd_seconds(mjd::Float64) = (mjd - 51544.0) * 86400.0

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
                 fillgaps::Bool)
    isfile(path) || throw(ArgumentError("read: file not found: " * String(path)))
    nskip, header_tau0, start = _scan_comment_header(path)
    # A Stable32 header block is itself non-comment header lines, so an
    # explicit `header` that already covers it is not added on top.
    stable32 = nskip > 0 && !startswith(readline(path), '#')
    value_col >= 1 || throw(ArgumentError("read: value_col must be ≥ 1"))
    have_time = time_col >= 1
    fillgaps && !have_time &&
        throw(ArgumentError("read: fillgaps=true requires a time column (got time_col=0)"))
    cols, _ = _read_cols(path, stable32 ? max(header, nskip) : header + nskip, delim,
                         have_time ? [time_col, value_col] : [value_col])
    v = cols[end]
    t = have_time ? cols[1] : Float64[]
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
    v, _ = _gridded(t, v, fillgaps)
    detrend === :none || _detrend_core!(v, detrend)
    return v, resolved_tau0, start
end

"""
    read_phase(path; tau0=nothing, time_col=1, value_col=2, header=0,
                     delim=nothing, scaling=1.0,
                     detrend=:none, fillgaps=false) → PhaseData

Read a phase-data file and return a [`PhaseData`](@ref) ready for stability
analysis. The returned record carries `path` in its `source` field.

Leading `#` comment lines (the header that [`save`](@ref) writes) are skipped
automatically; a `# tau0: <value>` entry supplies the sample interval when the
`tau0` keyword is omitted and a `# start: <ISO 8601>` entry the record's
[`start_time`](@ref), so files written by `save(path, data)` round-trip with
no keywords at all. Stable32 exports (a text header closed by `# Header End`)
are read the same way, with `Tau:` supplying the sample interval. `#` lines
later in the file are skipped too. The file is read line by line, keeping
only the time and value columns, so multi-day records at high rates load in
the memory of their samples.

The time column places the samples: each step is rounded to a whole number
of median steps, so jittered timestamps stay in order and a missing reading
leaves a `NaN` sample where it belongs, which the gap-tolerant estimators
skip (see [Timestamps and gaps](@ref)). A record with no missing timestamps is read exactly
as written. A column that is not a time axis (not increasing) leaves the
samples in file order.

# Keyword arguments

- `tau0`: sample interval in seconds. When omitted, taken from a
               `# tau0:` comment header if present, else from the time column
               (its median step); required when `time_col=0` and neither
               fallback applies.
- `time_col`: 1-based column index for timestamps (`0` = no time column; the
               samples are then taken in file order).
- `value_col`: 1-based column index for phase samples.
- `header`: number of non-comment header lines to skip (counted after any
               leading `#` comment lines).
- `delim`: field delimiter (`Char`). Auto-detected from extension when
               `nothing` (`,` for `.csv`, tab for `.tsv`, whitespace otherwise).
- `scaling`: multiply phase samples by this factor (e.g. `1e-12` for ps→s).
- `detrend`: one of `:none`, `:mean`, `:endpoint`, `:linear`, `:quadratic`.
               See [`detrend`](@ref).
- `fillgaps`: when `true`, fill the missing samples (and any `NaN` reading)
               with [`fillgaps`](@ref) (Howe's algorithm). Requires a time
               column.

# Examples

```julia
pd = read_phase("phase.csv")                          # 2-col CSV, auto tau0
pd = read_phase("p.tsv"; detrend=:linear)             # remove linear trend
pd = read_phase("p.txt"; time_col=0, tau0=1.0)        # no time column
pd = read_phase("g.csv"; fillgaps=true, scaling=1e-9) # ns→s + Howe fill
```
"""
function read_phase(path::AbstractString; tau0::Union{Real,Nothing}=nothing, time_col::Int=1,
                    value_col::Int=2, header::Int=0, delim::Union{Char,Nothing}=nothing,
                    scaling::Real=1.0, detrend::Symbol=:none, fillgaps::Bool=false)
    v, t0, start = _ingest(path; tau0=tau0, time_col=time_col, value_col=value_col, header=header,
                           delim=delim, scaling=scaling, detrend=detrend, fillgaps=fillgaps)
    return PhaseData(v, t0; source=String(path), start=start)
end

"""
    read_frequency(path; tau0=nothing, time_col=1, value_col=2, header=0,
                         delim=nothing, scaling=1.0,
                         detrend=:none, fillgaps=false) → FrequencyData

Read a fractional-frequency file. Same keyword arguments, comment-header and
time-column handling as [`read_phase`](@ref); only the return type differs.
The returned record carries `path` in its `source` field.
"""
function read_frequency(path::AbstractString; tau0::Union{Real,Nothing}=nothing, time_col::Int=1,
                        value_col::Int=2, header::Int=0, delim::Union{Char,Nothing}=nothing,
                        scaling::Real=1.0, detrend::Symbol=:none, fillgaps::Bool=false)
    v, t0, start = _ingest(path; tau0=tau0, time_col=time_col, value_col=value_col, header=header,
                           delim=delim, scaling=scaling, detrend=detrend, fillgaps=fillgaps)
    return FrequencyData(v, t0; source=String(path), start=start)
end

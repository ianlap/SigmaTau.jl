# Reading the record: a file through SigmaTau's readers, or a synthetic one.

# White FM + random-walk FM, built without FFTs so it matches earlier releases.
# With --clocks, clock --column of the simulated ensemble (the one timescale uses).
function synth_record(o::Options)
    if length(o.models) == 1
        R, q1, q2, q3 = o.models[1]
        tau0 = isnan(o.tau0) ? 1.0 : o.tau0
        z, _ = simulate_clock(ClockModel(; tau0, R, q1, q2, q3), o.synth; rng=Random.Xoshiro(o.seed))
        return PhaseData(z.x, tau0; source="synthetic")
    end
    if o.clocks >= 2 || !isempty(o.models)
        cc = synth_ensemble(o)
        return PhaseData(cc.x[:, o.clock], cc.tau0; source="synthetic")
    end
    rng = Random.Xoshiro(o.seed)
    N = o.synth
    tau0 = isnan(o.tau0) ? 1.0 : o.tau0
    y = 1e-12 .* randn(rng, N) .+ 1e-14 .* cumsum(randn(rng, N))
    return PhaseData(cumsum(y) .* tau0, tau0; source="synthetic")
end

# With `unfilled`, the samples that --fill-gaps put in are taken out again
# (NaN) once the record is read: timescale forms no scale from them.
function load_record(o::Options; unfilled::Bool=false)
    ranged = !(isnan(o.range_start) && isnan(o.range_end))
    data = read_record(o)
    keep = o.remove_drift == "none" && !unfilled && !o.fillgaps ? trues(0) : measured_mask(o, length(data.x))
    if ranged
        idx = range_indices(o, data)
        start = isnan(data.start) ? NaN : data.start + (first(idx) - 1) * data.tau0
        data = PhaseData(data.x[idx], data.tau0; source=data.source, start=start)
        isempty(keep) || (keep = keep[idx])
    end
    count_filled!(o, keep)
    # Phase records are cleaned in frequency space: a glitch is a phase step,
    # i.e. one bad frequency sample, which is gapped, filled and integrated back.
    o.remove_outliers && (data = remove_outliers(data))
    o.remove_drift == "none" || (data = remove_drift!(o, data, keep))
    unfilled && !all(keep) && (data = PhaseData(ifelse.(keep, data.x, NaN), data.tau0; source=data.source, start=data.start))
    return data
end

# Seconds from the record's first sample to the first sample of --range.
range_offset(o::Options, tau0::Float64) =
    isnan(o.range_start) ? 0.0 : max(0, ceil(Int, o.range_start / tau0 - 1e-9)) * tau0

# Phase samples inside --range: sample i is at (i − 1)·τ₀.
function range_indices(o::Options, data::PhaseData)
    N = length(data.x)
    tau0 = Float64(data.tau0)
    i0 = round(Int, range_offset(o, tau0) / tau0) + 1
    i1 = isnan(o.range_end) ? N : min(N, floor(Int, o.range_end / tau0 + 1e-9) + 1)
    i1 - i0 >= 1 || throw(CLIError(EXIT_ANALYSIS, "--range holds fewer than 2 samples of the record"))
    return i0:i1
end

# Which samples were measured, as opposed to put back by --fill-gaps: the
# finite samples of the value (and reference) column read without filling,
# which puts them on the grid of their time column with NaN in its empty
# slots. A filled frequency sample leaves the phase sample after it
# unmeasured. All true without --fill-gaps. Instrument formats report their
# filled samples directly (read_instrument).
function measured_mask(o::Options, N::Int)
    length(o.measured) == N && return copy(o.measured)
    keep = trues(N)
    (o.fillgaps && o.synth <= 0 && input_format(o) === :columns) || return keep
    for c in (o.value_col, o.reference_col)
        c > 0 || continue
        kw = (time_col=o.time_col, value_col=c, header=o.header, scaling=1.0, fillgaps=false)
        v = o.freq ? read_freq(o, kw).y : read_phase_col(o, kw).x
        ok = o.freq ? vcat(trues(1), isfinite.(v)) : isfinite.(v)
        length(ok) == N || return trues(N)
        keep .&= ok
    end
    return keep
end

# Least-squares line (frequency offset) or parabola (offset and a linear
# frequency drift) in the phase, fitted to the finite, measured samples only
# (`keep`), so a record with gaps can still be fitted and filled gaps do not
# weigh on the fit. The fit is removed from the phase and what it removed is
# kept on `o` for the JSON.
function remove_drift!(o::Options, data::PhaseData, keep::BitVector)
    x = Float64[Float64(v) for v in data.x]
    ok = keep .& isfinite.(x)
    quad = o.remove_drift == "quadratic"
    count(ok) >= (quad ? 3 : 2) || throw(CLIError(EXIT_ANALYSIS, "too few samples to fit the drift"))
    c = SigmaTau._index_polyfit(x, ok, quad ? 2 : 1)
    mu = (length(x) - 1) / 2.0
    for i in eachindex(x)
        t = (i - 1) - mu
        x[i] -= c[1] + (c[2] + c[3] * t) * t
    end
    tau0 = Float64(data.tau0)
    # dx/dt at the record's center is the frequency offset; d²x/dt² the drift.
    o.removed_offset = c[2] / tau0
    o.removed_drift = quad ? 2.0 * c[3] / (tau0 * tau0) * 86400.0 : NaN
    return PhaseData(x, data.tau0; source=data.source, start=data.start)
end

# The clock's column, minus the reference column when one is chosen (several
# clocks measured in one file against a common timebase).
function read_record(o::Options)
    o.synth > 0 && return synth_record(o)
    isfile(o.input) || throw(CLIError(EXIT_INPUT, "no such file: " * o.input))
    fmt = input_format(o)
    isempty(o.column) || fmt === :table || fmt === :rinex_clk || throw(usage_error(fmt === :columns ?
        "--column " * o.column * " names a column, and this file has no column names; pick it by number, " *
        "or give --input-format table if its first line names the columns" :
        "--column applies to a table with named columns or a RINEX clock file, not to a " * String(fmt) * " file"))
    fmt === :columns || return read_instrument(o, fmt)
    data = read_column(o, o.value_col)
    o.reference_col > 0 || return data
    o.reference_col == o.value_col && throw(usage_error("the reference column is the clock's column"))
    ref = try
        read_column(o, o.reference_col)
    catch err
        err isa ArgumentError || rethrow()
        throw(CLIError(EXIT_INPUT, "--reference-column " * string(o.reference_col) * " is not a column of this file"))
    end
    length(ref.x) == length(data.x) || throw(CLIError(EXIT_INPUT, "reference and value columns differ in length"))
    return PhaseData(Float64[Float64(data.x[i]) - Float64(ref.x[i]) for i in eachindex(data.x)], data.tau0;
                     source=data.source, start=data.start)
end

# Every keyword is concrete at each call site: a Union{Nothing,Float64} (or
# Char) keyword makes the kwcall dynamic under --trim, so the optional tau0
# and delim get their own branches.
function read_column(o::Options, value_col::Int)
    kw = (time_col=o.time_col, value_col=value_col, header=o.header, scaling=o.scaling, fillgaps=o.fillgaps)
    if o.freq
        fd = read_freq(o, kw)
        return PhaseData(integrated(fd.y, fd.tau0), fd.tau0; source=fd.source, start=fd.start)
    end
    return read_phase_col(o, kw)
end

# The filler was checked with up to 30 % of a record missing (ADEV within a
# few percent of the full record's); past that the results say as much about
# the filler as about the clock, and stderr says so.
const FILL_TESTED = 0.3

share_text(k::Int, n::Int) = 100 * k < n ? "under 1 %" : string(round(Int, 100 * k / n)) * " %"

function count_filled!(o::Options, keep::BitVector)
    o.filled = count(!, keep)
    n = length(keep)
    o.filled > FILL_TESTED * n && isempty(o.note) &&
        (o.note = "--fill-gaps made " * string(o.filled) * " of the " * string(n) * " samples (" * share_text(o.filled, n) *
                  "), more than the 30 % it was tested for: the results describe the filler as much as the clock")
    return o
end

# Same integration as SigmaTau's FrequencyData methods: N frequency samples
# give N + 1 phase points starting at 0, so no sample is lost. A missing
# frequency sample leaves every later phase unknown, so it must be filled.
function integrated(y::Vector{Float64}, tau0::Float64)
    n = count(isnan, y)
    n == 0 || throw(CLIError(EXIT_INPUT, "the frequency record is missing " * string(n) * " samples (" *
                                         share_text(n, length(y)) * " of it); fill them with --fill-gaps"))
    return [0.0; cumsum(y)] .* tau0
end

function read_phase_col(o::Options, kw)
    p = o.input
    if o.delim == '\0'
        return isnan(o.tau0) ? read_phase(p; kw...) : read_phase(p; tau0=o.tau0, kw...)
    end
    return isnan(o.tau0) ? read_phase(p; delim=o.delim, kw...) :
                           read_phase(p; tau0=o.tau0, delim=o.delim, kw...)
end

function read_freq(o::Options, kw)
    p = o.input
    if o.delim == '\0'
        return isnan(o.tau0) ? read_frequency(p; kw...) : read_frequency(p; tau0=o.tau0, kw...)
    end
    return isnan(o.tau0) ? read_frequency(p; delim=o.delim, kw...) :
                           read_frequency(p; tau0=o.tau0, delim=o.delim, kw...)
end

# The reader for this file: the one named by --input-format, else the
# detected one. A table goes through the column reader when its columns are
# given explicitly (as the explorer does after its import dialog).
function input_format(o::Options)
    o.input_format == "auto" || return Symbol(o.input_format)
    fmt = SigmaTau._detect_format(o.input, o.delim)
    fmt === :table && o.columns_set && isempty(o.column) && return :columns
    return fmt
end

# The names on the line above the rows of a table, with that line's number
# (counting every line, as --header does); none for any other file.
function table_head(o::Options)
    (o.synth > 0 || !isfile(o.input)) && return (0, String[])
    fmt = o.input_format == "auto" ? SigmaTau._detect_format(o.input, o.delim) : Symbol(o.input_format)
    fmt === :table || return (0, String[])
    for (i, line) in enumerate(SigmaTau._head_lines(o.input, 200))
        s = strip(line)
        (isempty(s) || startswith(s, '#')) && continue
        return (i, String[strip(f) for f in SigmaTau._split_fields(s, SigmaTau._delim_for(s, o.delim))])
    end
    return (0, String[])
end

table_names(o::Options) = table_head(o)[2]

# The columns of a table, numbered, to pick from.
column_list(names::Vector{String}) = join((string(i, " ", n) for (i, n) in enumerate(names)), ", ")

# timescale and hat read one clock per --column. Without them, a table's
# columns are listed by number, to pick from.
function require_value_cols(o::Options)
    (o.command == "timescale" || o.command == "hat") && o.synth <= 0 && isempty(o.value_cols) && isfile(o.input) || return
    names = table_names(o)
    throw(usage_error(o.command * " needs --column with two or more columns, or --synth" *
                      (isempty(names) ? "" : "; the columns of this file are " * column_list(names))))
end

# --column, --time-column and --reference-column as the readers take them: a
# number as it is, a name looked up in the first line of a table. One clock
# alone in a table or a RINEX clock file goes to their reader, which knows the
# units of a table's columns, by name or by number; with --synth, --column is
# the simulated clock.
function set_columns!(o::Options)
    if o.synth > 0
        isempty(o.columns) && return
        k = length(o.columns) == 1 ? tryparse(Int, o.columns[1]) : nothing
        k === nothing && throw(usage_error("with --synth, --column is the number of one simulated clock"))
        o.clock = k
        return
    end
    several = o.command == "timescale" || o.command == "hat" || o.command == "cross"
    length(o.columns) > 1 && !several &&
        throw(usage_error("`" * o.command * "` reads one clock: give one --column (a list is for timescale, hat and cross)"))
    o.command == "cross" && length(o.columns) > 2 && throw(usage_error("cross takes two --column, one per channel"))
    isempty(o.columns) && isempty(o.time_column) && isempty(o.reference_column) && return
    isfile(o.input) || return             # reading says so
    line, names = table_head(o)
    # A table is read by its own reader, which takes each column's unit from
    # its name; --header asks for its rows as plain columns instead.
    table = !isempty(names) && !o.columns_set
    table && !(isempty(o.time_column) && isempty(o.reference_column)) &&
        throw(usage_error((isempty(o.time_column) ? "--reference-column" : "--time-column") *
                          " applies to plain columns, and this file is a table read by the names of its columns;" *
                          " add --header " * string(line) * " to read its rows as plain columns, in the units of the numbers"))
    if length(o.columns) == 1 && !several && isempty(o.time_column) && isempty(o.reference_column) &&
       tryparse(Int, o.columns[1]) === nothing
        o.column = o.columns[1]           # a table's column or a RINEX clock, by name
        return
    end
    if !isempty(o.time_column)
        o.time_col = column_number(o.time_column, "--time-column", names, 0)
        o.columns_set = true
    end
    isempty(o.reference_column) || (o.reference_col = column_number(o.reference_column, "--reference-column", names, 1))
    cols = Int[column_number(c, "--column", names, 1) for c in o.columns]
    isempty(cols) && return
    o.value_col = cols[1]
    several || return
    length(cols) >= 2 && allunique(cols) || throw(usage_error("--column needs two or more different columns"))
    o.value_cols = cols
    return
end

function column_number(s::String, flag::String, names::Vector{String}, lowest::Int)
    k = tryparse(Int, s)
    if k === nothing
        i = findfirst(==(s), names)
        i === nothing && throw(usage_error(flag * " " * s * ": no column of that name" *
            (isempty(names) ? " (this file has no column names)" : "; the columns are " * column_list(names))))
        k = i
    end
    k >= lowest || throw(usage_error(flag * " must be " * string(lowest) * " or more"))
    return k
end

# TimeLab, Stable32 .phd/.frd, counter and named-column files: one clock per
# file, read by SigmaTau's format readers, which set the record type. With
# a --column number the clock of a table or RINEX clock file is that column.
function read_instrument(o::Options, fmt::Symbol)
    o.reference_col > 0 && throw(usage_error("--reference-column applies to column files only"))
    r = try
        SigmaTau._read_data(o.input, fmt; tau0=o.tau0, column=isempty(o.column) && !isempty(o.columns) ? o.value_col : 0,
                            column_name=o.column, nominal=o.nominal,
                            kind=o.freq ? :frequency : o.phase_given ? :phase : :auto, scaling=o.scaling, delim=o.delim,
                            fillgaps=o.fillgaps, keep_zeros=o.keep_zeros)
    catch e
        e isa ArgumentError || rethrow()
        throw(CLIError(EXIT_INPUT, error_message(e)))
    end
    o.freq = !r.phase
    # A filled frequency sample leaves the phase sample after it unmeasured.
    o.measured = isempty(r.filled) ? falses(0) : r.phase ? .!r.filled : [true; .!r.filled]
    r.phase && return PhaseData(r.v, r.tau0; source=o.input, start=r.start)
    return PhaseData(integrated(r.v, r.tau0), r.tau0; source=o.input, start=r.start)
end

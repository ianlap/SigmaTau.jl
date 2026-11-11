# Reading the record: a file through SigmaTau's readers, or a synthetic one.

detrend_symbol(s::String) = Symbol(s)

# White FM + random-walk FM, built without FFTs so it matches earlier releases.
# With --clocks, clock --clock of the simulated ensemble (the one timescale uses).
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

# With `unfilled`, the samples that --fillgaps put in are taken out again
# (NaN) once the record is read: timescale forms no scale from them.
function load_record(o::Options; unfilled::Bool=false)
    ranged = !(isnan(o.range_start) && isnan(o.range_end))
    # With --range the detrend is fitted to the range alone, after the cut.
    mode = o.detrend
    ranged && (o.detrend = "none")
    data = try
        read_record(o)
    finally
        o.detrend = mode
    end
    keep = o.remove_drift == "none" && !unfilled ? trues(0) : measured_mask(o, length(data.x))
    if ranged
        idx = range_indices(o, data)
        start = isnan(data.start) ? NaN : data.start + (first(idx) - 1) * data.tau0
        data = PhaseData(detrended(o, data.x[idx], data.tau0), data.tau0; source=data.source, start=start)
        isempty(keep) || (keep = keep[idx])
    end
    # Phase records are cleaned in frequency space: a glitch is a phase step,
    # i.e. one bad frequency sample, which is gapped, filled and integrated back.
    o.remove_outliers && (data = remove_outliers(data))
    o.remove_drift == "none" || (data = remove_drift!(o, data, keep))
    unfilled && !all(keep) && (data = PhaseData(ifelse.(keep, data.x, NaN), data.tau0; source=data.source, start=data.start))
    return data
end

# The readers' detrend, on a phase slice: frequency records are detrended in
# frequency, as the readers do, and integrated back.
function detrended(o::Options, x::Vector{Float64}, tau0::Float64)
    mode = detrend_symbol(o.detrend)
    mode === :none && return x
    o.freq || return SigmaTau._detrend_core!(x, mode)
    y = SigmaTau._detrend_core!(diff(x) ./ tau0, mode)
    return [0.0; cumsum(y)] .* tau0
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

# Which samples were measured, as opposed to put back by --fillgaps: the
# finite samples of the value (and reference) column read without filling,
# which puts them on the grid of their time column with NaN in its empty
# slots. A filled frequency sample leaves the phase sample after it
# unmeasured. All true without --fillgaps. Instrument formats report their
# filled samples directly (read_instrument).
function measured_mask(o::Options, N::Int)
    length(o.measured) == N && return copy(o.measured)
    keep = trues(N)
    (o.fillgaps && o.synth <= 0 && input_format(o) === :columns) || return keep
    for c in (o.value_col, o.reference_col)
        c > 0 || continue
        kw = (time_col=o.time_col, value_col=c, header=o.header, scaling=1.0, detrend=:none, fillgaps=false)
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

# The value column, minus the reference column when one is chosen (several
# clocks measured in one file against a common timebase). The readers' detrend
# is linear, so detrending each column and subtracting equals detrending the
# difference.
function read_record(o::Options)
    o.synth > 0 && return synth_record(o)
    isfile(o.input) || throw(CLIError(EXIT_INPUT, "no such file: " * o.input))
    fmt = input_format(o)
    fmt === :columns || return read_instrument(o, fmt)
    data = read_column(o, o.value_col)
    o.reference_col > 0 || return data
    o.reference_col == o.value_col && throw(usage_error("the reference column is the value column"))
    ref = try
        read_column(o, o.reference_col)
    catch err
        err isa ArgumentError || rethrow()
        throw(CLIError(EXIT_INPUT, "--reference-col " * string(o.reference_col) * " is not a column of this file"))
    end
    length(ref.x) == length(data.x) || throw(CLIError(EXIT_INPUT, "reference and value columns differ in length"))
    return PhaseData(Float64[Float64(data.x[i]) - Float64(ref.x[i]) for i in eachindex(data.x)], data.tau0;
                     source=data.source, start=data.start)
end

# Every keyword is concrete at each call site: a Union{Nothing,Float64} (or
# Char) keyword makes the kwcall dynamic under --trim, so the optional tau0
# and delim get their own branches.
function read_column(o::Options, value_col::Int)
    kw = (time_col=o.time_col, value_col=value_col, header=o.header, scaling=o.scaling,
          detrend=detrend_symbol(o.detrend), fillgaps=o.fillgaps)
    if o.freq
        fd = read_freq(o, kw)
        return PhaseData(integrated(fd.y, fd.tau0), fd.tau0; source=fd.source, start=fd.start)
    end
    return read_phase_col(o, kw)
end

# Same integration as SigmaTau's FrequencyData methods: N frequency samples
# give N + 1 phase points starting at 0, so no sample is lost. A missing
# frequency sample leaves every later phase unknown, so it must be filled.
function integrated(y::Vector{Float64}, tau0::Float64)
    n = count(isnan, y)
    n == 0 || throw(CLIError(EXIT_INPUT, "the frequency record is missing " * string(n) *
                                         " samples; fill them with --fillgaps"))
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

# TimeLab, Stable32 .phd/.frd, counter and named-column files: one clock per
# file, read by SigmaTau's format readers, which set the record type. With
# --value-cols (timescale) the clock of a table or RINEX clock file is the
# column being read, by number.
function read_instrument(o::Options, fmt::Symbol)
    o.reference_col > 0 && throw(usage_error("--reference-col applies to column files only"))
    r = try
        SigmaTau._read_data(o.input, fmt; tau0=o.tau0, column=isempty(o.value_cols) ? 0 : o.value_col,
                            column_name=o.column, nominal=o.nominal,
                            kind=o.freq ? :frequency : o.phase_given ? :phase : :auto, scaling=o.scaling, delim=o.delim,
                            fillgaps=o.fillgaps, keep_zeros=o.keep_zeros)
    catch e
        e isa ArgumentError || rethrow()
        throw(CLIError(EXIT_INPUT, cli_wording(error_message(e))))
    end
    o.freq = !r.phase
    # A filled frequency sample leaves the phase sample after it unmeasured.
    o.measured = isempty(r.filled) ? falses(0) : r.phase ? .!r.filled : [true; .!r.filled]
    detrend = detrend_symbol(o.detrend)
    v = r.v
    detrend === :none || SigmaTau._detrend_core!(v, detrend)
    r.phase && return PhaseData(v, r.tau0; source=o.input, start=r.start)
    return PhaseData(integrated(v, r.tau0), r.tau0; source=o.input, start=r.start)
end

# read_data's errors in the words of the command line, where the keyword
# `column` is --column (a name) or --value-cols (numbers).
cli_wording(m::String) =
    replace(replace(replace(m, "pass column= one of" => "pass --column with one of"), "column=" => "column "), "pass tau0" => "pass --tau0")

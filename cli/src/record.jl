# Reading the record: a file through SigmaTau's readers, or a synthetic one.

detrend_symbol(s::String) =
    s == "mean" ? :mean : s == "endpoint" ? :endpoint : s == "linear" ? :linear : :none

# White FM + random-walk FM, built without FFTs so it matches earlier releases.
function synth_record(o::Options)
    rng = Random.Xoshiro(o.seed)
    N = o.synth
    tau0 = isnan(o.tau0) ? 1.0 : o.tau0
    y = 1e-12 .* randn(rng, N) .+ 1e-14 .* cumsum(randn(rng, N))
    return PhaseData(cumsum(y) .* tau0, tau0; source="synthetic")
end

function load_record(o::Options)
    ranged = !(isnan(o.range_start) && isnan(o.range_end))
    # With --range the detrend is fitted to the range alone, after the cut.
    mode = o.detrend
    ranged && (o.detrend = "none")
    data = try
        read_record(o)
    finally
        o.detrend = mode
    end
    keep = o.remove_drift == "none" ? trues(0) : measured_mask(o, length(data.x))
    if ranged
        idx = range_indices(o, data)
        data = PhaseData(detrended(o, data.x[idx], data.tau0), data.tau0; source=data.source)
        isempty(keep) || (keep = keep[idx])
    end
    # Phase records are cleaned in frequency space: a glitch is a phase step,
    # i.e. one bad frequency sample, which is gapped, filled and integrated back.
    o.remove_outliers && (data = remove_outliers(data))
    o.remove_drift == "none" || (data = remove_drift!(o, data, keep))
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
# readers snap the time column to an even grid of its smallest step (see
# SigmaTau's _make_equispaced) and fill the empty slots, so the same snapping
# of the raw time column finds the slots that held data. All true without
# --fillgaps or when the grid does not line up with the record. Instrument
# formats report their filled samples directly (read_instrument).
function measured_mask(o::Options, N::Int)
    length(o.measured) == N && return copy(o.measured)
    keep = trues(N)
    (o.fillgaps && o.synth <= 0 && o.time_col >= 1 && input_format(o) === :columns) || return keep
    t = read_times(o)
    length(t) >= 2 || return keep
    dt = Inf
    for i in 2:length(t)
        d = t[i] - t[i - 1]
        d < dt && (dt = d)
    end
    dt > 0 || return keep
    round(Int, (t[end] - t[1]) / dt) + 1 == N || return keep
    fill!(keep, false)
    for ti in t
        idx = round(Int, (ti - t[1]) / dt) + 1
        1 <= idx <= N && (keep[idx] = true)
    end
    return keep
end

# The time column, read as a value column with the same header and delimiter.
function read_times(o::Options)
    kw = (time_col=0, value_col=o.time_col, header=o.header, scaling=1.0,
          detrend=:none, fillgaps=false, tau0=1.0)
    p = o.input
    pd = o.delim == '\0' ? read_phase(p; kw...) : read_phase(p; delim=o.delim, kw...)
    return Float64[Float64(v) for v in pd.x]
end

det3(a::Float64, b::Float64, c::Float64, d::Float64, e::Float64, f::Float64, g::Float64, h::Float64, k::Float64) =
    a * (e * k - f * h) - b * (d * k - f * g) + c * (d * h - e * g)

# Least-squares fit of a line (frequency offset) or a parabola (offset and a
# linear frequency drift) to the phase, over the finite, measured samples only
# (`keep`), so a record with gaps can still be fitted and filled gaps do not
# weigh on the fit. The fit is removed from the phase and
# what it removed is kept on `o` for the JSON. Time is centred, which zeroes
# the odd power sums and keeps the normal equations well conditioned.
function remove_drift!(o::Options, data::PhaseData, keep::BitVector)
    x = Float64[Float64(v) for v in data.x]
    N = length(x)
    tau0 = Float64(data.tau0)
    mu = (N - 1) / 2.0
    n = 0.0; S1 = 0.0; S2 = 0.0; S3 = 0.0; S4 = 0.0; Sy = 0.0; S1y = 0.0; S2y = 0.0
    for i in 1:N
        xi = x[i]
        (keep[i] && isfinite(xi)) || continue
        t = (i - 1) - mu; t2 = t * t
        n += 1.0; S1 += t; S2 += t2; S3 += t2 * t; S4 += t2 * t2
        Sy += xi; S1y += t * xi; S2y += t2 * xi
    end
    quad = o.remove_drift == "quadratic"
    n >= (quad ? 3 : 2) || throw(CLIError(EXIT_ANALYSIS, "too few samples to fit the drift"))
    c0 = 0.0; c1 = 0.0; c2 = 0.0
    if quad
        # Cramer's rule on [n S1 S2; S1 S2 S3; S2 S3 S4] c = [Sy; S1y; S2y].
        D = det3(n, S1, S2, S1, S2, S3, S2, S3, S4)
        c0 = det3(Sy, S1, S2, S1y, S2, S3, S2y, S3, S4) / D
        c1 = det3(n, Sy, S2, S1, S1y, S3, S2, S2y, S4) / D
        c2 = det3(n, S1, Sy, S1, S2, S1y, S2, S3, S2y) / D
    else
        D = n * S2 - S1 * S1
        c1 = (n * S1y - S1 * Sy) / D
        c0 = (Sy - c1 * S1) / n
    end
    for i in 1:N
        t = (i - 1) - mu
        x[i] -= c0 + c1 * t + c2 * t * t
    end
    # dx/dt at the record's centre is the frequency offset; d²x/dt² the drift.
    o.removed_offset = c1 / tau0
    o.removed_drift = quad ? 2.0 * c2 / (tau0 * tau0) * 86400.0 : NaN
    return PhaseData(x, data.tau0; source=data.source)
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
                     source=data.source)
end

# Every keyword is concrete at each call site: a Union{Nothing,Float64} (or
# Char) keyword makes the kwcall dynamic under --trim, so the optional tau0
# and delim get their own branches.
function read_column(o::Options, value_col::Int)
    kw = (time_col=o.time_col, value_col=value_col, header=o.header, scaling=o.scaling,
          detrend=detrend_symbol(o.detrend), fillgaps=o.fillgaps)
    if o.freq
        fd = read_freq(o, kw)
        # Same integration as SigmaTau's FrequencyData methods: N frequency
        # samples give N + 1 phase points starting at 0, so no sample is lost.
        return PhaseData([0.0; cumsum(fd.y)] .* fd.tau0, fd.tau0; source=fd.source)
    end
    return read_phase_col(o, kw)
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
# file, read by SigmaTau's format readers, which set the record type.
function read_instrument(o::Options, fmt::Symbol)
    o.reference_col > 0 && throw(usage_error("--reference-col applies to column files only"))
    r = SigmaTau._read_data(o.input, fmt; tau0=o.tau0, column_name=o.column, nominal=o.nominal,
                            kind=o.freq ? :frequency : o.phase_given ? :phase : :auto, scaling=o.scaling, delim=o.delim,
                            fillgaps=o.fillgaps)
    o.freq = !r.phase
    # A filled frequency sample leaves the phase sample after it unmeasured.
    isempty(r.filled) || (o.measured = r.phase ? .!r.filled : [true; .!r.filled])
    detrend = detrend_symbol(o.detrend)
    v = r.v
    detrend === :none || SigmaTau._detrend_core!(v, detrend)
    r.phase && return PhaseData(v, r.tau0; source=o.input)
    return PhaseData([0.0; cumsum(v)] .* r.tau0, r.tau0; source=o.input)
end

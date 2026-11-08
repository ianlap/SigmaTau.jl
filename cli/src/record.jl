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
    data = read_record(o)
    # Phase records are cleaned in frequency space: a glitch is a phase step,
    # i.e. one bad frequency sample, which is gapped, filled and integrated back.
    o.remove_outliers && (data = remove_outliers(data))
    return data
end

# The value column, minus the reference column when one is chosen (several
# clocks measured in one file against a common timebase). The readers' detrend
# is linear, so detrending each column and subtracting equals detrending the
# difference.
function read_record(o::Options)
    o.synth > 0 && return synth_record(o)
    isfile(o.input) || throw(CLIError(EXIT_INPUT, "no such file: " * o.input))
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
        return PhaseData(cumsum(fd.y) .* fd.tau0, fd.tau0; source=fd.source)
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

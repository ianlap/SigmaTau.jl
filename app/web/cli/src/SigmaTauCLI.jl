# sigmatau-cli: the analysis core of the web explorer, written to compile
# with `juliac --trim`. It reads one record (a file through SigmaTau's own
# readers, or a synthetic one), runs the requested deviations and prints one
# JSON document on stdout. The JSON contract is in app/web/README.md.
#
#   julia --project=app/web/cli -m SigmaTauCLI --synth 8192 --devs adev,mdev
#   juliac --output-exe sigmatau-cli --bundle build --trim=safe --experimental app/web/cli
#
# Everything here stays type-stable for trimming: options are parsed into a
# concrete struct, estimators are chosen by an if/elseif chain rather than a
# Dict of functions, and JSON is written by hand.
module SigmaTauCLI

using SigmaTau
using Random

# Baked in at compile time; pkgversion needs the source tree at run time.
const SIGMATAU_VERSION = string(pkgversion(SigmaTau))

const ALL_DEVS = ["adev", "mdev", "tdev", "hdev", "mhdev", "htdev",
                  "totdev", "mtotdev", "ttotdev", "htotdev", "mhtotdev",
                  "theo1", "theoh", "pdev", "mtie", "tierms"]
const SIGMA_X_DEVS = ["tdev", "htdev", "ttotdev", "mtie", "tierms"]

struct CLIError <: Exception
    msg::String
end

mutable struct Options
    input::String
    synth::Int
    seed::Int
    kind::String          # "phase" | "frequency"
    tau0::Float64         # NaN = from the file
    time_col::Int
    value_col::Int
    reference_col::Int    # 0 = none; else the value is taken against this column
    header::Int
    scaling::Float64
    detrend::String
    devs::Vector{String}
    grid::String
    ci::Bool
    confidence::Float64
    bias::Bool
    series::Int           # max record points returned for plotting
    fillgaps::Bool        # Howe gap filling (needs a time column; uses FFTW)
    remove_outliers::Bool # replace 5σ outliers and fill (uses FFTW)
    spectrum::Bool        # add S_y(f) and S_x(f) (Welch)
    nperseg::Int          # Welch segment length; 0 = auto
    dynamic::String       # "" = none, "dadev" or "dhdev"
    window::Int           # dynamic window in samples; 0 = N ÷ 8
    step::Int             # dynamic step in samples; 0 = window ÷ 2
end

Options() = Options("", 0, 42, "phase", NaN, 1, 2, 0, 0, 1.0, "none",
                    ["adev", "mdev", "hdev", "mhdev"], "octave", true, 0.683, true, 4000,
                    false, false, false, 0, "", 0, 0)

# ── Argument parsing ──────────────────────────────────────────────────────

function parse_int(s::String, flag::String)
    v = tryparse(Int, s)
    v === nothing && throw(CLIError(flag * " expects an integer, got \"" * s * "\""))
    return v
end

function parse_float(s::String, flag::String)
    v = tryparse(Float64, s)
    v === nothing && throw(CLIError(flag * " expects a number, got \"" * s * "\""))
    return v
end

# The value after the flag at args[i]. A plain function, not a closure over
# `i`: a captured, reassigned local is boxed, which trimming rejects.
function flag_value(args::Vector{String}, i::Int)
    i < length(args) || throw(CLIError(args[i] * " needs a value"))
    return args[i + 1]
end

# "none" (or an empty list) asks for the record and its check only.
function parse_devs(s::String)
    s == "none" && return String[]
    return String[String(d) for d in split(s, ',') if !isempty(d)]
end

function parse_args(args::Vector{String})
    o = Options()
    i = 1
    n = length(args)
    while i <= n
        a = args[i]
        if a == "--input";          o.input = flag_value(args, i); i += 1
        elseif a == "--synth";      o.synth = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--seed";       o.seed = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--kind";       o.kind = flag_value(args, i); i += 1
        elseif a == "--tau0";       o.tau0 = parse_float(flag_value(args, i), a); i += 1
        elseif a == "--time-col";   o.time_col = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--value-col";  o.value_col = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--reference-col"; o.reference_col = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--header";     o.header = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--scaling";    o.scaling = parse_float(flag_value(args, i), a); i += 1
        elseif a == "--detrend";    o.detrend = flag_value(args, i); i += 1
        elseif a == "--devs";       o.devs = parse_devs(flag_value(args, i)); i += 1
        elseif a == "--grid";       o.grid = flag_value(args, i); i += 1
        elseif a == "--confidence"; o.confidence = parse_float(flag_value(args, i), a); i += 1
        elseif a == "--no-ci";      o.ci = false
        elseif a == "--no-bias-correction"; o.bias = false
        elseif a == "--series";     o.series = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--fillgaps";   o.fillgaps = true
        elseif a == "--remove-outliers"; o.remove_outliers = true
        elseif a == "--spectrum";   o.spectrum = true
        elseif a == "--nperseg";    o.nperseg = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--dynamic";    o.dynamic = flag_value(args, i); i += 1
        elseif a == "--window";     o.window = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--step";       o.step = parse_int(flag_value(args, i), a); i += 1
        else throw(CLIError("unknown option " * a))
        end
        i += 1
    end
    isempty(o.input) && o.synth <= 0 && throw(CLIError("give --input PATH or --synth N"))
    o.kind in ("phase", "frequency") || throw(CLIError("--kind must be phase or frequency"))
    o.detrend in ("none", "mean", "endpoint", "linear") ||
        throw(CLIError("--detrend must be none, mean, endpoint or linear"))
    for d in o.devs
        d in ALL_DEVS || throw(CLIError("unknown deviation " * d))
    end
    0 < o.confidence < 1 || throw(CLIError("--confidence must be between 0 and 1"))
    o.dynamic in ("", "dadev", "dhdev") || throw(CLIError("--dynamic must be dadev or dhdev"))
    o.nperseg >= 0 || throw(CLIError("--nperseg must be positive"))
    o.window >= 0 && o.step >= 0 || throw(CLIError("--window and --step must be positive"))
    return o
end

# ── Record ────────────────────────────────────────────────────────────────

function detrend_symbol(s::String)
    s == "mean" && return :mean
    s == "endpoint" && return :endpoint
    s == "linear" && return :linear
    return :none
end

# Synthetic white FM + random-walk FM without FFTs (noise_gen goes through
# FFTW, which does not trim yet).
function synth_record(o::Options)
    rng = Random.Xoshiro(o.seed)
    N = o.synth
    tau0 = isnan(o.tau0) ? 1.0 : o.tau0
    y = 1e-12 .* randn(rng, N) .+ 1e-14 .* cumsum(randn(rng, N))
    return PhaseData(cumsum(y) .* tau0, tau0; source="synthetic")
end

function load_record(o::Options)
    data = read_record(o)
    o.remove_outliers && (data = deglitch(data))
    return data
end

# Remove 5σ outliers in the frequency domain and integrate back. A glitch in a
# phase record is a step, i.e. one bad frequency sample; SigmaTau's
# remove_outliers(::PhaseData) refills the phase sample the step lands on,
# which keeps the step. Going through FrequencyData removes it.
function deglitch(data::PhaseData)
    x = data.x
    N = length(x)
    N < 3 && return data
    t0 = Float64(data.tau0)
    y = Float64[(Float64(x[i + 1]) - Float64(x[i])) / t0 for i in 1:N-1]
    fd = remove_outliers(FrequencyData(y, t0; source=data.source))
    out = Vector{Float64}(undef, N)
    out[1] = Float64(x[1])
    for i in 1:N-1
        out[i + 1] = out[i] + Float64(fd.y[i]) * t0
    end
    return PhaseData(out, t0; source=data.source)
end

# The value column, minus the reference column when one is chosen (several
# clocks measured in one file against a common timebase). The readers' detrend
# is linear, so detrending each column and subtracting equals detrending the
# difference.
function read_record(o::Options)
    o.synth > 0 && return synth_record(o)
    isfile(o.input) || throw(CLIError("no such file: " * o.input))
    data = read_column(o, o.value_col)
    o.reference_col > 0 || return data
    o.reference_col == o.value_col && throw(CLIError("the reference column is the value column"))
    o.reference_col >= 1 || throw(CLIError("--reference-col must be 1 or more"))
    ref = try
        read_column(o, o.reference_col)
    catch err
        err isa ArgumentError || rethrow()
        throw(CLIError("--reference-col " * string(o.reference_col) * " is not a column of this file"))
    end
    length(ref.x) == length(data.x) || throw(CLIError("reference and value columns differ in length"))
    return PhaseData(Float64[Float64(data.x[i]) - Float64(ref.x[i]) for i in eachindex(data.x)], data.tau0;
                     source=data.source)
end

function read_column(o::Options, value_col::Int)
    det = detrend_symbol(o.detrend)
    # Two concrete call sites instead of `tau0 = isnan ? nothing : x`: a
    # Union{Nothing,Float64} keyword makes the kwcall dynamic under --trim.
    if o.kind == "frequency"
        fd = isnan(o.tau0) ?
            read_frequency(o.input; time_col=o.time_col, value_col=value_col,
                           header=o.header, scaling=o.scaling, detrend=det, fillgaps=o.fillgaps) :
            read_frequency(o.input; tau0=o.tau0, time_col=o.time_col, value_col=value_col,
                           header=o.header, scaling=o.scaling, detrend=det, fillgaps=o.fillgaps)
        return PhaseData(cumsum(fd.y) .* fd.tau0, fd.tau0; source=fd.source)
    end
    return isnan(o.tau0) ?
        read_phase(o.input; time_col=o.time_col, value_col=value_col,
                   header=o.header, scaling=o.scaling, detrend=det, fillgaps=o.fillgaps) :
        read_phase(o.input; tau0=o.tau0, time_col=o.time_col, value_col=value_col,
                   header=o.header, scaling=o.scaling, detrend=det, fillgaps=o.fillgaps)
end

# ── Estimators ────────────────────────────────────────────────────────────

function tau_mode(s::String)
    s == "octave" && return Octave
    s == "half-octave" && return HalfOctave
    s == "quarter-octave" && return QuarterOctave
    s == "decade" && return Decade
    s == "half-decade" && return HalfDecade
    s == "all" && return AllTaus
    throw(CLIError("unknown --grid " * s))
end

function run_dev(d::String, data::PhaseData, ms::Vector{Int}, o::Options)
    ci, c, b = o.ci, o.confidence, o.bias
    d == "adev"     && return adev(data, ms; ci=ci, confidence=c)
    d == "mdev"     && return mdev(data, ms; ci=ci, confidence=c)
    d == "tdev"     && return tdev(data, ms; ci=ci, confidence=c)
    d == "hdev"     && return hdev(data, ms; ci=ci, confidence=c)
    d == "mhdev"    && return mhdev(data, ms; ci=ci, confidence=c)
    d == "htdev"    && return htdev(data, ms; ci=ci, confidence=c)
    d == "totdev"   && return totdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    d == "mtotdev"  && return mtotdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    d == "ttotdev"  && return ttotdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    d == "htotdev"  && return htotdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    d == "mhtotdev" && return mhtotdev(data, ms; ci=ci, confidence=c, correct_bias=b)
    d == "theo1"    && return theo1(data, ms; ci=ci, confidence=c, correct_bias=b)
    d == "theoh"    && return theoh(data, ms; ci=ci, confidence=c)
    d == "pdev"     && return pdev(data, ms; ci=ci, confidence=c)
    d == "mtie"     && return mtie(data, ms; ci=ci, confidence=c)
    d == "tierms"   && return tierms(data, ms; ci=ci, confidence=c)
    throw(CLIError("unknown deviation " * d))
end

function kernel_symbol(d::String)
    for (s, k) in (("adev", :adev), ("mdev", :mdev), ("tdev", :tdev), ("hdev", :hdev),
                   ("mhdev", :mhdev), ("htdev", :htdev), ("totdev", :totdev),
                   ("mtotdev", :mtotdev), ("ttotdev", :ttotdev), ("htotdev", :htotdev),
                   ("mhtotdev", :mhtotdev), ("theo1", :theo1), ("theoh", :theoh),
                   ("pdev", :pdev), ("mtie", :mtie), ("tierms", :tierms))
        d == s && return k
    end
    throw(CLIError("unknown deviation " * d))
end

# ── JSON output ───────────────────────────────────────────────────────────

function json_string(io::IO, s::AbstractString)
    print(io, '"')
    for c in s
        if c == '"';      print(io, "\\\"")
        elseif c == '\\'; print(io, "\\\\")
        elseif c == '\n'; print(io, "\\n")
        elseif c == '\r'; print(io, "\\r")
        elseif c == '\t'; print(io, "\\t")
        elseif c < ' ';   print(io, "\\u00"); print(io, string(UInt8(c); base=16, pad=2))
        else              print(io, c)
        end
    end
    print(io, '"')
end

json_number(io::IO, x::Float64) = isfinite(x) ? print(io, x) : print(io, "null")
json_number(io::IO, x::Int) = print(io, x)

function json_array(io::IO, v::AbstractVector{Float64})
    print(io, '[')
    for i in eachindex(v)
        i > firstindex(v) && print(io, ',')
        json_number(io, v[i])
    end
    print(io, ']')
end

function json_array(io::IO, v::AbstractVector{Int})
    print(io, '[')
    for i in eachindex(v)
        i > firstindex(v) && print(io, ',')
        print(io, v[i])
    end
    print(io, ']')
end

function json_array(io::IO, v::AbstractVector{Symbol})
    print(io, '[')
    for i in eachindex(v)
        i > firstindex(v) && print(io, ',')
        json_string(io, String(v[i]))
    end
    print(io, ']')
end

function write_result(io::IO, d::String, r::StabilityResult)
    print(io, "{\"dev\":"); json_string(io, d)
    print(io, ",\"axis\":"); json_string(io, d in SIGMA_X_DEVS ? "sigma_x" : "sigma_y")
    print(io, ",\"tau\":"); json_array(io, r.tau)
    print(io, ",\"sigma\":"); json_array(io, r.dev)
    lo = Float64[b.lo for b in r.ci]
    hi = Float64[b.hi for b in r.ci]
    print(io, ",\"ci_lo\":"); json_array(io, lo)
    print(io, ",\"ci_hi\":"); json_array(io, hi)
    print(io, ",\"edf\":"); json_array(io, r.edf)
    print(io, ",\"noise\":"); json_array(io, r.noise_type)
    print(io, ",\"neff\":"); json_array(io, r.neff)
    print(io, '}')
end

# Evenly strided subset of the phase record, for the Data page plot.
function write_series(io::IO, data::PhaseData, maxpts::Int)
    N = length(data.x)
    step = maxpts <= 0 ? 0 : max(1, cld(N, maxpts))
    idx = step == 0 ? Int[] : collect(1:step:N)
    t = Float64[(i - 1) * data.tau0 for i in idx]
    x = Float64[Float64(data.x[i]) for i in idx]
    print(io, "{\"t\":"); json_array(io, t)
    print(io, ",\"x\":"); json_array(io, x)
    print(io, '}')
end

# Record check for the Data page: gaps, 5σ outliers (SigmaTau's MAD test on
# the phase steps) and the linear frequency drift, per day.
function write_check(io::IO, data::PhaseData)
    x = data.x
    N = length(x)
    gaps = count(isnan, x)
    out = find_outliers(data; nsigma=5)
    # Drift: least-squares slope of y = Δx/τ₀ against time, finite samples only.
    sxy = 0.0; sxx = 0.0; st = 0.0; sy = 0.0; n = 0
    for i in 1:N-1
        yi = (Float64(x[i + 1]) - Float64(x[i])) / data.tau0
        isfinite(yi) || continue
        ti = (i - 1) * data.tau0
        st += ti; sy += yi; sxy += ti * yi; sxx += ti * ti; n += 1
    end
    drift = n > 2 ? (n * sxy - st * sy) / (n * sxx - st * st) * 86400.0 : NaN
    print(io, "{\"samples\":"); print(io, N)
    print(io, ",\"span_s\":"); json_number(io, N * data.tau0)
    print(io, ",\"gaps\":"); print(io, gaps)
    print(io, ",\"outliers\":"); print(io, length(out))
    keep = out[1:min(length(out), 500)]
    print(io, ",\"outlier_t\":"); json_array(io, Float64[(i - 1) * data.tau0 for i in keep])
    print(io, ",\"outlier_x\":"); json_array(io, Float64[Float64(x[i]) for i in keep])
    print(io, ",\"drift_per_day\":"); json_number(io, drift)
    print(io, '}')
end

function write_error(io::IO, msg::String)
    print(io, "{\"error\":"); json_string(io, msg); print(io, "}\n")
end

# ── Entry point ───────────────────────────────────────────────────────────

function analyze(io::IO, o::Options)
    t_start = time_ns()
    data = load_record(o)
    N = length(data.x)
    # Reading and checking a record (--devs none) needs only a few samples;
    # the estimators need at least eight.
    nmin = isempty(o.devs) ? 3 : 8
    N >= nmin || throw(CLIError("the record has only " * string(N) * " samples"))
    mode = tau_mode(o.grid)
    buf = IOBuffer()
    print(buf, "{\"sigmatau\":"); json_string(buf, SIGMATAU_VERSION)
    print(buf, ",\"record\":{\"n\":"); print(buf, N)
    print(buf, ",\"tau0\":"); json_number(buf, Float64(data.tau0))
    print(buf, ",\"kind\":"); json_string(buf, o.kind)
    print(buf, ",\"source\":"); json_string(buf, data.source)
    print(buf, ",\"series\":"); write_series(buf, data, o.series)
    print(buf, ",\"check\":"); write_check(buf, data)
    print(buf, "},\"options\":{\"grid\":"); json_string(buf, o.grid)
    print(buf, ",\"ci\":"); print(buf, o.ci ? "true" : "false")
    print(buf, ",\"confidence\":"); json_number(buf, o.confidence)
    print(buf, ",\"bias_correction\":"); print(buf, o.bias ? "true" : "false")
    print(buf, ",\"fillgaps\":"); print(buf, o.fillgaps ? "true" : "false")
    print(buf, ",\"remove_outliers\":"); print(buf, o.remove_outliers ? "true" : "false")
    print(buf, "},\"results\":[")
    first = true
    for d in o.devs
        ms = tau_values(mode, N, kernel_symbol(d))
        isempty(ms) && continue
        r = run_dev(d, data, ms, o)
        first || print(buf, ',')
        write_result(buf, d, r)
        first = false
    end
    print(buf, ']')
    if o.spectrum
        print(buf, ",\"spectrum\":"); write_spectrum(buf, data, o)
    end
    if !isempty(o.dynamic)
        print(buf, ",\"dynamic\":"); write_dynamic(buf, data, o, mode)
    end
    print(buf, ",\"elapsed_ms\":")
    json_number(buf, (time_ns() - t_start) / 1e6)
    print(buf, "}\n")
    write(io, take!(buf))
    return nothing
end

# ── Spectrum and dynamic deviation ────────────────────────────────────────

# Welch segment: the largest power of two up to N/4 (at least 16), so the
# spectrum reaches low frequencies and still averages four or more segments.
function auto_nperseg(N::Int)
    np = 16
    while 2np <= N ÷ 4
        np *= 2
    end
    return min(np, N)
end

# S_y(f) and S_x(f) on one frequency grid, DC dropped. ℒ(f) is S_x scaled by
# the carrier, so the page computes it.
function write_spectrum(io::IO, data::PhaseData, o::Options)
    N = length(data.x)
    N >= 16 || throw(CLIError("the spectrum needs at least 16 samples"))
    any(isnan, data.x) && throw(CLIError("the spectrum needs a record without gaps; fill them on the Data page"))
    np = o.nperseg > 0 ? min(o.nperseg, N) : auto_nperseg(N)
    sy = Sy(data; nperseg=np)
    sx = Sx(data; nperseg=np)
    n = min(length(sy.freq), length(sx.freq))
    print(io, "{\"nperseg\":"); print(io, np)
    print(io, ",\"f\":"); json_array(io, sx.freq[2:n])
    print(io, ",\"sy\":"); json_array(io, sy.psd[2:n])
    print(io, ",\"sx\":"); json_array(io, sx.psd[2:n])
    print(io, '}')
end

# Sliding-window deviation map: rows are windows (centre times t), columns τ.
function write_dynamic(io::IO, data::PhaseData, o::Options, mode::TauMode)
    N = length(data.x)
    w = o.window > 0 ? min(o.window, N) : max(N ÷ 8, 16)
    st = o.step > 0 ? o.step : max(w ÷ 2, 1)
    hd = o.dynamic == "dhdev"
    ms = tau_values(mode, w, hd ? :hdev : :adev)
    isempty(ms) && throw(CLIError("the dynamic window is too short for any τ"))
    r = hd ? dhdev(data, ms; window=w, step=st) : dadev(data, ms; window=w, step=st)
    print(io, "{\"dev\":"); json_string(io, o.dynamic)
    print(io, ",\"window\":"); print(io, w)
    print(io, ",\"step\":"); print(io, st)
    print(io, ",\"t\":"); json_array(io, r.t)
    print(io, ",\"tau\":"); json_array(io, r.tau)
    print(io, ",\"map\":[")
    for i in 1:size(r.dev, 1)
        i > 1 && print(io, ',')
        json_array(io, Float64[r.dev[i, j] for j in 1:size(r.dev, 2)])
    end
    print(io, "]}")
end

function (@main)(args::Vector{String})::Cint
    io = Core.stdout
    if isempty(args) || args[1] == "--help" || args[1] == "-h"
        println(io, "usage: sigmatau-cli (--input PATH | --synth N) [options]; see app/web/README.md")
        return isempty(args) ? Cint(1) : Cint(0)
    end
    try
        analyze(io, parse_args(args))
        return Cint(0)
    catch err
        msg = err isa CLIError ? err.msg :
              err isa ArgumentError ? (err.msg isa String ? err.msg : "invalid argument") :
              err isa DomainError ? "domain error in the analysis" :
              "analysis failed (" * String(nameof(typeof(err))) * ")"
        write_error(io, msg)
        return Cint(2)
    end
end

end # module

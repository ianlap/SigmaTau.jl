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
    header::Int
    scaling::Float64
    detrend::String
    devs::Vector{String}
    grid::String
    ci::Bool
    confidence::Float64
    bias::Bool
    series::Int           # max record points returned for plotting
end

Options() = Options("", 0, 42, "phase", NaN, 1, 2, 0, 1.0, "none",
                    ["adev", "mdev", "hdev", "mhdev"], "octave", true, 0.683, true, 4000)

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
        elseif a == "--header";     o.header = parse_int(flag_value(args, i), a); i += 1
        elseif a == "--scaling";    o.scaling = parse_float(flag_value(args, i), a); i += 1
        elseif a == "--detrend";    o.detrend = flag_value(args, i); i += 1
        elseif a == "--devs";       o.devs = String[String(s) for s in split(flag_value(args, i), ',') if !isempty(s)]; i += 1
        elseif a == "--grid";       o.grid = flag_value(args, i); i += 1
        elseif a == "--confidence"; o.confidence = parse_float(flag_value(args, i), a); i += 1
        elseif a == "--no-ci";      o.ci = false
        elseif a == "--no-bias-correction"; o.bias = false
        elseif a == "--series";     o.series = parse_int(flag_value(args, i), a); i += 1
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
    o.synth > 0 && return synth_record(o)
    isfile(o.input) || throw(CLIError("no such file: " * o.input))
    det = detrend_symbol(o.detrend)
    # Two concrete call sites instead of `tau0 = isnan ? nothing : x`: a
    # Union{Nothing,Float64} keyword makes the kwcall dynamic under --trim.
    if o.kind == "frequency"
        fd = isnan(o.tau0) ?
            read_frequency(o.input; time_col=o.time_col, value_col=o.value_col,
                           header=o.header, scaling=o.scaling, detrend=det) :
            read_frequency(o.input; tau0=o.tau0, time_col=o.time_col, value_col=o.value_col,
                           header=o.header, scaling=o.scaling, detrend=det)
        return PhaseData(cumsum(fd.y) .* fd.tau0, fd.tau0; source=fd.source)
    end
    return isnan(o.tau0) ?
        read_phase(o.input; time_col=o.time_col, value_col=o.value_col,
                   header=o.header, scaling=o.scaling, detrend=det) :
        read_phase(o.input; tau0=o.tau0, time_col=o.time_col, value_col=o.value_col,
                   header=o.header, scaling=o.scaling, detrend=det)
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

function write_error(io::IO, msg::String)
    print(io, "{\"error\":"); json_string(io, msg); print(io, "}\n")
end

# ── Entry point ───────────────────────────────────────────────────────────

function analyze(io::IO, o::Options)
    t_start = time_ns()
    data = load_record(o)
    N = length(data.x)
    N >= 8 || throw(CLIError("the record has only " * string(N) * " samples"))
    mode = tau_mode(o.grid)
    buf = IOBuffer()
    print(buf, "{\"sigmatau\":"); json_string(buf, SIGMATAU_VERSION)
    print(buf, ",\"record\":{\"n\":"); print(buf, N)
    print(buf, ",\"tau0\":"); json_number(buf, Float64(data.tau0))
    print(buf, ",\"kind\":"); json_string(buf, o.kind)
    print(buf, ",\"source\":"); json_string(buf, data.source)
    print(buf, ",\"series\":"); write_series(buf, data, o.series)
    print(buf, "},\"options\":{\"grid\":"); json_string(buf, o.grid)
    print(buf, ",\"ci\":"); print(buf, o.ci ? "true" : "false")
    print(buf, ",\"confidence\":"); json_number(buf, o.confidence)
    print(buf, ",\"bias_correction\":"); print(buf, o.bias ? "true" : "false")
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
    print(buf, "],\"elapsed_ms\":")
    json_number(buf, (time_ns() - t_start) / 1e6)
    print(buf, "}\n")
    write(io, take!(buf))
    return nothing
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

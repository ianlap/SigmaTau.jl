# Options: one table drives the parser, the help text and the shell
# completions, so the three cannot drift apart.

const COMMANDS = ["analyze", "check", "spectrum", "dynamic", "cross", "drift", "report", "live"]
const COMMAND_HELP = [
    ("analyze",    "FILE", "stability deviations σ(τ) with confidence intervals"),
    ("check",      "FILE", "read a record and summarise it (samples, span, gaps, outliers, drift)"),
    ("spectrum",   "FILE", "Welch power spectral densities S_y(f) and S_x(f)"),
    ("dynamic",    "FILE", "dynamic ADEV or HDEV: σ(τ) in sliding windows"),
    ("cross",      "FILE --with FILE2", "two channels measuring one device: cross Allan deviation and cross-spectrum"),
    ("drift",      "FILE", "fit a frequency drift model: aging per day with its uncertainty"),
    ("live",       "SOURCE", "running σ(τ) from an instrument or the simulator, as JSON lines"),
    ("report",     "FILE", "self-contained HTML report: record, σ(τ) plots and tables (prints to PDF)"),
    ("completion", "SHELL", "print a shell completion script (bash, zsh, fish, powershell)"),
    ("version",    "",     "print the version"),
    ("help",       "[COMMAND]", "print this help, or a command's options"),
]

const SHELLS = ["bash", "zsh", "fish", "powershell"]
const FORMATS = ["csv", "json"]
const INPUT_LAYOUTS = ["auto", "columns", "stable32", "timelab", "counter", "table", "bipm_utc"]
const DETRENDS = ["none", "mean", "endpoint", "linear"]
const DRIFTS = ["none", "linear", "quadratic"]
const GRIDS = ["octave", "half-octave", "quarter-octave", "decade", "half-decade", "all"]
const DYNAMIC_DEVS = ["dadev", "dhdev"]
const DRIFT_MODELS = [String(m) for m in SigmaTau.DRIFT_MODELS]
const MASKS = [[String(m) for m in SigmaTau.TELECOM_MASKS]; "all"]

struct Opt
    long::String
    metavar::String          # "" for a flag
    choices::Vector{String}  # values offered by completion; empty = free
    path::Bool               # value is a file path (completion offers files)
    help::String
    cmds::Vector{String}     # commands that accept it
end

const ALL = filter(!=("live"), COMMANDS)   # the commands that read a FILE
const ANY = COMMANDS
const TABULAR = ["analyze", "check", "spectrum", "dynamic", "cross", "drift"]   # CSV or JSON out (report writes HTML)
const OPTS = [
    # The record
    Opt("--synth", "N", String[], false, "use a synthetic N-sample record (white + random-walk FM) instead of FILE", ALL),
    Opt("--seed", "S", String[], false, "random seed for --synth and the live simulator (default 42)", ANY),
    Opt("--phase", "", String[], false, "FILE holds phase in seconds (default)", ANY),
    Opt("--freq", "", String[], false, "FILE holds fractional frequency", ANY),
    Opt("--input-format", "FMT", INPUT_LAYOUTS, false, "file layout: auto (default), columns, stable32 (.phd/.frd), timelab (.tim), counter (SCPI readings), table (named columns), bipm_utc (BIPM UTC-UTC(k))", ALL),
    Opt("--column", "NAME", String[], false, "clock column of a table file, by name (default: the column named phase / PPS / time error, else frequency)", ALL),
    Opt("--nominal", "HZ", String[], false, "nominal frequency for readings in Hz (default: their mean to 4 significant figures)", ANY),
    Opt("--tau0", "SECONDS", String[], false, "sample interval τ₀ in seconds (default: from the file header or time column)", ANY),
    Opt("--time-col", "N", String[], false, "1-based time column, 0 = none (default 1)", ALL),
    Opt("--value-col", "N", String[], false, "1-based value column (default 2)", ALL),
    Opt("--reference-col", "N", String[], false, "subtract this column from the value column (clocks measured against one reference)", ALL),
    Opt("--header", "N", String[], false, "header lines to skip; '#' lines and a Stable32 header are skipped anyway", ALL),
    Opt("--delim", "CHAR", String[], false, "field delimiter (default: from the extension; '\\t' for tab)", ALL),
    Opt("--scaling", "X", String[], false, "multiply the values by X (1e-9 for nanoseconds)", ALL),
    Opt("--detrend", "MODE", DETRENDS, false, "remove none, mean, endpoint or linear drift (default none)", ALL),
    Opt("--remove-drift", "MODE", DRIFTS, false, "fit and remove from the phase: none (default), linear (frequency offset) or quadratic (offset and linear frequency drift)", ALL),
    Opt("--fillgaps", "", String[], false, "put the record on an even time grid and fill gaps (needs a time column)", ALL),
    Opt("--remove-outliers", "", String[], false, "replace 5σ outliers in frequency and fill them", ALL),
    Opt("--range", "A:B", String[], false, "analyse only seconds A to B from the first sample (either may be empty: 3600: or :86400)", ALL),
    # Analysis
    Opt("--devs", "LIST", String[], false, "comma list of deviations (default adev,mdev,hdev,mhdev); see below", ["analyze", "report", "live"]),
    Opt("--grid", "GRID", GRIDS, false, "τ grid: octave (default), half-octave, quarter-octave, decade, half-decade, all", ["analyze", "dynamic", "cross", "report", "live"]),
    Opt("--m", "LIST", String[], false, "explicit averaging factors m (τ = m·τ₀), e.g. 1,2,4,10; replaces --grid", ["analyze", "cross", "report", "live"]),
    Opt("--tau-min", "SECONDS", String[], false, "leave out averaging times below SECONDS", ["analyze", "report"]),
    Opt("--tau-max", "SECONDS", String[], false, "leave out averaging times above SECONDS", ["analyze", "report"]),
    Opt("--confidence", "P", String[], false, "confidence level of the intervals (default 0.683)", ["analyze", "report"]),
    Opt("--no-ci", "", String[], false, "skip confidence intervals, EDF and noise identification", ["analyze", "report"]),
    Opt("--no-bias-correction", "", String[], false, "leave the total estimators and Thêo1 uncorrected", ["analyze", "report"]),
    Opt("--spectrum", "", String[], false, "also return the spectrum (JSON only)", ["analyze"]),
    Opt("--fit", "", String[], false, "also return the power-law noise model fitted to MHDEV (JSON only)", ["analyze"]),
    Opt("--mask", "NAME", MASKS, false, "also check MTIE and TDEV against an ITU-T wander mask, or all of them (JSON only): G811, G813_opt1, G8262_opt1, G823_SSU ..., all", ["analyze"]),
    Opt("--nperseg", "N", String[], false, "Welch segment length (default: largest power of 2 ≤ N/4)", ["analyze", "spectrum", "cross"]),
    Opt("--with", "FILE2", String[], true, "the second channel, read with the same options as FILE", ["cross"]),
    Opt("--dev", "NAME", DYNAMIC_DEVS, false, "dadev (default) or dhdev", ["dynamic"]),
    Opt("--window", "N", String[], false, "window length in samples (default: largest power of 2 ≤ N/8)", ["dynamic"]),
    Opt("--step", "N", String[], false, "window step in samples (default window/2)", ["dynamic"]),
    Opt("--model", "MODEL", DRIFT_MODELS, false, "linear (default), quadratic, log (a·ln(1+b·t)+c) or diffusion (a+b·√(t+c))", ["drift"]),
    # Live
    Opt("--hz", "", String[], false, "readings are counter frequencies in Hz, against --nominal (default: the first reading to 4 significant figures)", ["live"]),
    Opt("--query", "CMD", String[], false, "SCPI query sent before every reading, e.g. READ? (default: READ? for tcp://, none otherwise)", ["live"]),
    Opt("--every", "N", String[], false, "write a line every N readings (default 1)", ["live"]),
    Opt("--count", "N", String[], false, "stop after N readings (default: when the source ends)", ["live"]),
    Opt("--fast", "", String[], false, "run the simulator flat out instead of one reading per τ₀", ["live"]),
    # Output
    Opt("--series", "N", String[], false, "include up to N record points in the JSON (default 0)", TABULAR),
    Opt("--format", "FMT", FORMATS, false, "csv (default) or json", TABULAR),
    Opt("--out", "PATH", String[], true, "write to PATH instead of stdout", ALL),
    Opt("--quiet", "", String[], false, "no informational messages on stderr", ANY),
    Opt("--help", "", String[], false, "show this command's options", ANY),
]

mutable struct Options
    command::String
    input::String
    synth::Int
    seed::Int
    freq::Bool
    tau0::Float64          # NaN = from the file
    time_col::Int
    value_col::Int
    reference_col::Int     # 0 = none
    header::Int
    delim::Char            # '\0' = from the extension
    scaling::Float64
    detrend::String
    fillgaps::Bool
    remove_outliers::Bool
    devs::Vector{Symbol}
    grid::String
    m::Vector{Int}         # empty = use the grid
    ci::Bool
    confidence::Float64
    bias::Bool
    spectrum::Bool
    nperseg::Int           # 0 = library default
    dev::String
    window::Int            # 0 = default
    step::Int              # 0 = window ÷ 2
    series::Int
    format::String
    out::String
    quiet::Bool
    help::Bool
    remove_drift::String
    removed_offset::Float64   # fitted and removed by --remove-drift; NaN = none
    removed_drift::Float64    # per day
    input_format::String
    column::String            # "" = by name
    nominal::Float64          # NaN = from the readings
    columns_set::Bool         # --time-col, --value-col or --header given
    measured::BitVector       # phase samples the format readers did not fill; empty = all
    with::String              # cross: the second channel's file
    phase_given::Bool         # --phase given (else instrument readers pick the kind)
    model::String
    hz::Bool
    query::String             # "\0" = the default for the source
    every::Int
    count::Int                # 0 = until the source ends
    fast::Bool
    range_start::Float64      # --range, seconds from the first sample; NaN = open
    range_end::Float64
    tau_min::Float64          # NaN = no limit
    tau_max::Float64
    fit::Bool                 # --fit: power-law model in the JSON
    mask::String              # --mask: ITU-T mask checked in the JSON; "" = none
end

Options() = Options("", "", 0, 42, false, NaN, 1, 2, 0, 0, '\0', 1.0, "none", false, false,
                    collect(Symbol, DEFAULT_DEVIATIONS), "octave", Int[], true,
                    Float64(DEFAULT_CONFIDENCE), true, false, 0, "dadev", 0, 0, 0, "csv", "", false, false,
                    "none", NaN, NaN, "auto", "", NaN, false, falses(0), "", false, "linear",
                    false, " ", 1, 0, false, NaN, NaN, NaN, NaN, false, "")

quantity(o::Options) = o.freq ? "frequency" : "phase"

function find_opt(name::String)
    for op in OPTS
        op.long == name && return op
    end
    return nothing
end

function parse_int(s::String, flag::String)
    v = tryparse(Int, s)
    v === nothing && throw(usage_error(flag * " expects an integer, got \"" * s * "\""))
    return v
end

function parse_float(s::String, flag::String)
    v = tryparse(Float64, s)
    (v === nothing || !isfinite(v)) && throw(usage_error(flag * " expects a finite number, got \"" * s * "\""))
    return v
end

# "A:B" in seconds; an empty side is open.
function parse_range!(o::Options, s::String)
    k = findfirst(==(':'), s)
    k === nothing && throw(usage_error("--range expects START:END in seconds, got \"" * s * "\""))
    a, b = s[1:prevind(s, k)], s[nextind(s, k):end]
    o.range_start = isempty(a) ? NaN : parse_float(String(a), "--range")
    o.range_end = isempty(b) ? NaN : parse_float(String(b), "--range")
    !isnan(o.range_start) && !isnan(o.range_end) && o.range_end <= o.range_start &&
        throw(usage_error("--range end must be after its start"))
    return o
end

function check_choice(s::String, flag::String, choices::Vector{String})
    s in choices || throw(usage_error(flag * " must be one of " * join(choices, ", ") * ", got \"" * s * "\""))
    return s
end

function estimator_names()
    return String[String(e.name) for e in estimators()]
end

function parse_devs(s::String)
    names = estimator_names()
    devs = Symbol[]
    for tok in split(s, ',')
        t = String(strip(tok))
        isempty(t) && continue
        if t in DYNAMIC_DEVS
            throw(usage_error(t * " is a dynamic deviation; use `sigmatau dynamic --dev " * t * "`"))
        end
        t in names || throw(usage_error("unknown deviation \"" * t * "\"; valid: " * join(names, ", ")))
        push!(devs, Symbol(t))
    end
    isempty(devs) && throw(usage_error("--devs: no deviations given"))
    return devs
end

function parse_m(s::String)
    ms = Int[]
    for tok in split(s, ',')
        t = String(strip(tok))
        isempty(t) && continue
        v = tryparse(Int, t)
        (v === nothing || v < 1) && throw(usage_error("--m expects positive integers, got \"" * t * "\""))
        push!(ms, v)
    end
    isempty(ms) && throw(usage_error("--m: no averaging factors given"))
    return sort!(unique!(ms))
end

function parse_delim(s::String)
    s == "\\t" && return '\t'
    length(s) == 1 || throw(usage_error("--delim expects one character, got \"" * s * "\""))
    return s[1]
end

function set_option!(o::Options, name::String, v::String)
    if name == "--synth";              o.synth = parse_int(v, name)
    elseif name == "--seed";           o.seed = parse_int(v, name)
    elseif name == "--tau0";           o.tau0 = parse_float(v, name)
    elseif name == "--input-format";   o.input_format = check_choice(v, name, INPUT_LAYOUTS)
    elseif name == "--column";         o.column = v
    elseif name == "--nominal";        o.nominal = parse_float(v, name)
    elseif name == "--time-col";       o.time_col = parse_int(v, name); o.columns_set = true
    elseif name == "--value-col";      o.value_col = parse_int(v, name); o.columns_set = true
    elseif name == "--reference-col";  o.reference_col = parse_int(v, name)
    elseif name == "--header";         o.header = parse_int(v, name); o.columns_set = true
    elseif name == "--delim";          o.delim = parse_delim(v)
    elseif name == "--scaling";        o.scaling = parse_float(v, name)
    elseif name == "--detrend";        o.detrend = check_choice(v, name, DETRENDS)
    elseif name == "--remove-drift";   o.remove_drift = check_choice(v, name, DRIFTS)
    elseif name == "--devs";           o.devs = parse_devs(v)
    elseif name == "--grid";           tau_mode(v); o.grid = v
    elseif name == "--m";              o.m = parse_m(v)
    elseif name == "--confidence";     o.confidence = parse_float(v, name)
    elseif name == "--range";          parse_range!(o, v)
    elseif name == "--tau-min";        o.tau_min = parse_float(v, name)
    elseif name == "--tau-max";        o.tau_max = parse_float(v, name)
    elseif name == "--nperseg";        o.nperseg = parse_int(v, name)
    elseif name == "--dev";            o.dev = check_choice(v, name, DYNAMIC_DEVS)
    elseif name == "--window";         o.window = parse_int(v, name)
    elseif name == "--step";           o.step = parse_int(v, name)
    elseif name == "--model";          o.model = check_choice(v, name, DRIFT_MODELS)
    elseif name == "--series";         o.series = parse_int(v, name)
    elseif name == "--query";          o.query = v
    elseif name == "--every";          o.every = parse_int(v, name)
    elseif name == "--count";          o.count = parse_int(v, name)
    elseif name == "--format";         o.format = check_choice(v, name, FORMATS)
    elseif name == "--out";            o.out = v
    elseif name == "--with";           o.with = v
    elseif name == "--mask";           o.mask = check_choice(v, name, MASKS)
    else throw(usage_error("option " * name * " takes no value"))
    end
    return o
end

function set_flag!(o::Options, name::String)
    if name == "--phase";                  o.freq = false; o.phase_given = true
    elseif name == "--freq";               o.freq = true
    elseif name == "--fillgaps";           o.fillgaps = true
    elseif name == "--remove-outliers";    o.remove_outliers = true
    elseif name == "--no-ci";              o.ci = false
    elseif name == "--no-bias-correction"; o.bias = false
    elseif name == "--spectrum";           o.spectrum = true
    elseif name == "--hz";                 o.hz = true
    elseif name == "--fast";               o.fast = true
    elseif name == "--fit";                o.fit = true
    elseif name == "--quiet";              o.quiet = true
    elseif name == "--help";               o.help = true
    end
    return o
end

function parse_args!(o::Options, args::Vector{String})
    i = 2
    n = length(args)
    while i <= n
        a = args[i]
        if a == "-h"
            o.help = true
        elseif startswith(a, "--")
            eq = findfirst(==('='), a)
            name = eq === nothing ? a : a[1:prevind(a, eq)]
            op = find_opt(name)
            op === nothing && throw(usage_error("unknown option " * name))
            o.command in op.cmds || throw(usage_error(name * " does not apply to `" * o.command * "`"))
            if isempty(op.metavar)
                eq === nothing || throw(usage_error("option " * name * " takes no value"))
                set_flag!(o, name)
            elseif eq !== nothing
                set_option!(o, name, a[nextind(a, eq):end])
            else
                i < n || throw(usage_error("option " * name * " needs a value"))
                set_option!(o, name, args[i + 1])
                i += 1
            end
        elseif startswith(a, "-") && a != "-"
            throw(usage_error("unknown option " * a))
        else
            isempty(o.input) || throw(usage_error("unexpected argument \"" * a * "\" (one FILE per call)"))
            o.input = a
        end
        i += 1
    end
    o.help && return o
    ("--phase" in args && "--freq" in args) && throw(usage_error("give only one of --phase and --freq"))
    if o.command == "live"
        isempty(o.input) && throw(usage_error("live: give a SOURCE: sim, tcp://HOST:PORT or a device path"))
        o.hz && o.freq && throw(usage_error("give only one of --freq and --hz"))
        o.every >= 1 || throw(usage_error("--every must be 1 or more"))
        o.count >= 0 || throw(usage_error("--count must be 0 or more"))
    elseif isempty(o.input) && o.synth <= 0
        throw(usage_error(o.command * ": give a FILE or --synth N (see `sigmatau help " * o.command * "`)"))
    end
    !isempty(o.input) && o.synth > 0 && throw(usage_error("give either a FILE or --synth, not both"))
    o.command == "cross" && (isempty(o.with) || o.synth > 0) &&
        throw(usage_error("cross: give FILE and the second channel with --with FILE2"))
    0 < o.confidence < 1 || throw(usage_error("--confidence must be between 0 and 1"))
    o.nperseg >= 0 || throw(usage_error("--nperseg must be positive"))
    o.window >= 0 && o.step >= 0 || throw(usage_error("--window and --step must be positive"))
    o.series >= 0 || throw(usage_error("--series must be 0 or more"))
    o.reference_col >= 0 || throw(usage_error("--reference-col must be 1 or more"))
    o.time_col >= 0 || throw(usage_error("--time-col must be 0 (no time column) or more"))
    o.value_col >= 1 || throw(usage_error("--value-col must be 1 or more"))
    o.header >= 0 || throw(usage_error("--header must be 0 or more"))
    isnan(o.tau0) || o.tau0 > 0 || throw(usage_error("--tau0 must be positive"))
    isnan(o.nominal) || o.nominal > 0 || throw(usage_error("--nominal must be positive"))
    o.scaling != 0 || throw(usage_error("--scaling must not be 0"))
    o.spectrum && o.format != "json" &&
        throw(usage_error("--spectrum adds to the JSON output; use --format json, or `sigmatau spectrum`"))
    o.fit && o.format != "json" && throw(usage_error("--fit adds to the JSON output; use --format json"))
    isempty(o.mask) || o.format == "json" || throw(usage_error("--mask adds to the JSON output; use --format json"))
    return o
end

# ── Help ──────────────────────────────────────────────────────────────────

# `s` then spaces up to `width` columns. Base.rpad goes through repeat(::Char),
# which trimming cannot resolve.
function padded(io::IO, s::String, width::Int)
    print(io, s)
    for _ in 1:(width - length(s))
        print(io, ' ')
    end
end

function print_help(io::IO)
    print(io, "sigmatau ", SIGMATAU_VERSION, ": frequency stability analysis with SigmaTau.jl\n\n")
    print(io, "Usage: sigmatau COMMAND [FILE] [options]\n\nCommands:\n")
    for (c, arg, h) in COMMAND_HELP
        print(io, "  "); padded(io, c * (isempty(arg) ? "" : " " * arg), 22); print(io, h, "\n")
    end
    print(io, "\nRun `sigmatau help COMMAND` for a command's options.\n")
    print(io, "Output is CSV by default; --format json gives the versioned JSON schema.\n")
    print(io, "Exit codes: 0 ok, 1 usage, 2 input, 3 analysis, 4 output.\n\n")
    print(io, "Examples:\n")
    print(io, "  sigmatau analyze clock.dat --tau0 1 --devs adev,mdev,tdev\n")
    print(io, "  sigmatau analyze freq.txt --freq --m 1,10,100 --format json --out result.json\n")
    print(io, "  sigmatau check clock.dat --time-col 1 --value-col 3 --fillgaps\n")
    print(io, "  sigmatau analyze trace.tim                      # TimeLab, τ₀ from the header\n")
    print(io, "  sigmatau spectrum --synth 65536\n")
    print(io, "  sigmatau dynamic clock.dat --dev dhdev --window 4096\n")
    print(io, "  sigmatau drift ocxo.dat --freq --model log\n")
    print(io, "  sigmatau live tcp://192.168.1.20:5025 --hz --nominal 10e6 --tau0 1\n")
    print(io, "  sigmatau report clock.dat --devs adev,tdev --out report.html\n")
end

function print_command_help(io::IO, cmd::String)
    if cmd == "completion"
        print(io, "Usage: sigmatau completion SHELL\n\nShells: ", join(SHELLS, ", "), "\n")
        return
    elseif !(cmd in COMMANDS)
        print(io, "sigmatau: no help for \"", cmd, "\"; commands: ", join(COMMANDS, ", "), "\n")
        return
    end
    for (c, arg, h) in COMMAND_HELP
        c == cmd && print(io, "Usage: sigmatau ", c, " ", arg, " [options]\n\n", uppercasefirst(h), ".\n")
    end
    if cmd == "live"
        print(io, "SOURCE is sim (a simulated clock), tcp://HOST:PORT (an LXI instrument's SCPI socket,\n",
              "usually port 5025) or a serial device path set up with stty. Each reading is the first\n",
              "number on a line: phase in seconds by default. Every line of output is one JSON object.\n\nOptions:\n")
    else
        print(io, "FILE is a whitespace, comma or tab separated table (time, value by default).\n\nOptions:\n")
    end
    sigs = String[]
    for op in OPTS
        cmd in op.cmds && push!(sigs, isempty(op.metavar) ? op.long : op.long * " " * op.metavar)
    end
    w = maximum(length, sigs) + 2
    k = 0
    for op in OPTS
        cmd in op.cmds || continue
        k += 1
        print(io, "  "); padded(io, sigs[k], w); print(io, op.help, "\n")
    end
    if cmd == "live"
        print(io, "\nDeviations: adev, mdev, hdev, mhdev\n")
    elseif cmd == "analyze" || cmd == "report"
        print(io, "\nDeviations: ", join(estimator_names(), ", "), "\n")
    end
end

# Options: one table drives the parser, the help text and the shell
# completions, so the three cannot drift apart.

const COMMANDS = ["analyze", "check", "spectrum", "dynamic", "cross", "drift", "holdover", "steer", "timescale", "hat", "model", "phasenoise", "report", "live"]
const COMMAND_HELP = [
    ("analyze",    "FILE", "stability deviations σ(τ) with confidence intervals"),
    ("check",      "FILE", "read a record and summarize it (samples, span, gaps, outliers, drift and the statistics of its samples), with their histogram and autocorrelation on request"),
    ("spectrum",   "FILE", "Welch power spectral densities S_y(f) and S_x(f)"),
    ("dynamic",    "FILE", "dynamic ADEV or HDEV: σ(τ) in sliding windows"),
    ("cross",      "FILE --with FILE2", "two channels measuring one device: cross Allan deviation and cross-spectrum"),
    ("drift",      "FILE", "fit a frequency drift model: aging per day with its uncertainty"),
    ("holdover",   "FILE", "holdover test: the time error after a training stretch against holdover time, and the time to reach a limit"),
    ("steer",      "FILE", "steer the clock of a record to its reference in replay: σ(τ) free-running, steered and as the clock model predicts"),
    ("timescale",  "FILE --value-cols LIST", "form a timescale from several clocks read against one reference: weights and stability"),
    ("hat",        "FILE --value-cols LIST", "the stability of each of three or more clocks measured against one another: N-cornered hat or Groslambert covariance"),
    ("model",      "--synth N --noise-model R:Q1:Q2:Q3", "the σ(τ) that a clock noise model predicts (MDEV, ADEV and MHDEV), without a simulation"),
    ("phasenoise", "TABLE --carrier HZ", "the σ(τ) that a table of the phase noise ℒ(f), or of another spectrum, gives over its band"),
    ("live",       "SOURCE", "running σ(τ) from an instrument or the simulator, as JSON lines"),
    ("report",     "FILE", "self-contained HTML report: record, σ(τ) plots and tables (prints to PDF)"),
    ("completion", "SHELL", "print a shell completion script (bash, zsh, fish, powershell)"),
    ("version",    "",     "print the version"),
    ("help",       "[COMMAND]", "print this help, or a command's options"),
]

const SHELLS = ["bash", "zsh", "fish", "powershell"]
const FORMATS = ["csv", "json"]
const INPUT_LAYOUTS = ["auto", "columns", "stable32", "timelab", "counter", "table", "bipm_utc", "rinex_clk"]
const DETRENDS = ["none", "mean", "endpoint", "linear"]
const DRIFTS = ["none", "linear", "quadratic"]
const BIASES = ["exact", "published", "stable32", "none"]
const NOISES = ["auto", "fit", "WHPM", "FLPM", "WHFM", "FLFM", "RWFM", "FWFM", "RRFM"]
const GRIDS = ["octave", "half-octave", "quarter-octave", "decade", "half-decade", "all"]
const DYNAMIC_DEVS = ["dadev", "dhdev"]
const STEPS = ["best", "real"]
const DRIFT_MODELS = [String(m) for m in SigmaTau.DRIFT_MODELS]
const DRIFT_METHODS = [String(m) for m in SigmaTau.DRIFT_METHODS]
const TIMESCALES = ["at1", "at2", "ewfa", "greenhall", "jones_tryon", "stein"]
const NOISE_NAMES = ["WHPM", "FLPM", "WHFM", "FLFM", "RWFM", "FWFM", "RRFM"]   # α = 2 … −4
const HAT_DEVS = [String(d) for d in SigmaTau._HAT_DEVIATIONS]
const MASKS = [[String(m) for m in SigmaTau.TELECOM_MASKS]; "all"]
const QUANTITIES = ["phase", "frequency"]
const DENSITIES = ["L", "Sphi", "Sx", "Sy"]
const PHASENOISE_DEVS = [:adev, :mdev, :tdev, :hdev, :mhdev]

struct Opt
    long::String
    metavar::String          # "" for a flag
    choices::Vector{String}  # values offered by completion; empty = free
    path::Bool               # value is a file path (completion offers files)
    help::String
    cmds::Vector{String}     # commands that accept it
end

const ANY = COMMANDS
const BATCH = filter(!=("live"), COMMANDS)    # one result per call, to stdout or --out
const CLOCK = filter(!=("phasenoise"), COMMANDS)   # a record or a clock model, not a spectrum table
const READS = filter(!=("model"), CLOCK)      # the commands that read a record
const ALL = filter(!=("live"), READS)         # … from a FILE or --synth
const TABULAR = filter(!=("report"), ALL)     # CSV or JSON out (report writes HTML; model is always JSON)
const OPTS = [
    # The record
    Opt("--synth", "N", String[], false, "use a synthetic N-sample record (white + random-walk FM) instead of FILE; for model, the length of the record whose averaging times it takes", filter(!=("live"), CLOCK)),
    Opt("--seed", "S", String[], false, "random seed for --synth and the live simulator (default 42)", READS),
    Opt("--clocks", "K", String[], false, "with --synth: simulate K clocks (2 to 64) against a perfect reference, each noisier than the last by a factor of 2; the analysis reads clock --clock, timescale uses all K (default 4)", ALL),
    Opt("--noise-model", "R:Q1:Q2:Q3", String[], false, "with --synth: one clock's noise, white PM variance R (s²) and diffusion q1, q2, q3 (white FM, random-walk FM, random-run FM); repeat per clock, which sets the number of clocks. `model` predicts its σ(τ)", filter(!=("live"), CLOCK)),
    Opt("--clock", "I", String[], false, "with --synth --clocks: the simulated clock to read, 1 to K (default 1)", ALL),
    Opt("--phase", "", String[], false, "the values are phase in seconds (default)", READS),
    Opt("--freq", "", String[], false, "the values are fractional frequency", READS),
    Opt("--input-format", "FMT", INPUT_LAYOUTS, false, "file layout: auto (default), columns, stable32 (.phd/.frd), timelab (.tim), counter (SCPI readings), table (named columns), bipm_utc (BIPM UTC-UTC(k)), rinex_clk (RINEX clock file; pick the clock with --column)", ALL),
    Opt("--column", "NAME", String[], false, "clock column of a table file, or the clock of a RINEX clock file (G08, ALGO00CAN), by name (default: the column named phase / PPS / time error, else frequency)", ALL),
    Opt("--nominal", "HZ", String[], false, "nominal frequency for readings in Hz (default: their mean, or for live the first reading, to 4 significant figures)", READS),
    Opt("--tau0", "SECONDS", String[], false, "sample interval τ₀ in seconds (default: from the file header or time column; 1 for --synth, model and live)", CLOCK),
    Opt("--time-col", "N", String[], false, "1-based time column, 0 = none (default 1)", ALL),
    Opt("--value-col", "N", String[], false, "1-based value column (default 2)", [ALL; "phasenoise"]),
    Opt("--reference-col", "N", String[], false, "subtract this column from the value column (clocks measured against one reference)", ALL),
    Opt("--header", "N", String[], false, "header lines to skip; '#' lines and a Stable32 header are skipped anyway", [ALL; "phasenoise"]),
    Opt("--delim", "CHAR", String[], false, "field delimiter (default: a comma where the first data line has one, a tab in a table or .tsv file, else blanks; '\\t' for tab)", [ALL; "phasenoise"]),
    Opt("--scaling", "X", String[], false, "multiply the values by X (1e-9 for nanoseconds)", ALL),
    Opt("--detrend", "MODE", DETRENDS, false, "remove none, mean, endpoint or linear drift (default none)", ALL),
    Opt("--remove-drift", "MODE", DRIFTS, false, "fit and remove from the phase: none (default), linear (frequency offset) or quadratic (offset and linear frequency drift)", ALL),
    Opt("--fillgaps", "", String[], false, "put the record on an even time grid and fill gaps and NaN readings", ALL),
    Opt("--keep-zeros", "", String[], false, "read the zeros of a Stable32 .phd/.frd file as readings, not gaps (quantized data)", ALL),
    Opt("--remove-outliers", "", String[], false, "replace 5σ outliers in frequency and fill them", ALL),
    Opt("--range", "A:B", String[], false, "analyze only seconds A to B from the first sample (either may be empty: 3600: or :86400)", ALL),
    # Analysis
    Opt("--devs", "LIST", String[], false, "comma list of deviations (default adev,mdev,hdev,mhdev; for hat, adev); see below", ["analyze", "hat", "report", "live", "phasenoise"]),
    Opt("--grid", "GRID", GRIDS, false, "τ grid: octave (default), half-octave, quarter-octave, decade, half-decade, all", ["analyze", "dynamic", "cross", "steer", "hat", "report", "live"]),
    Opt("--m", "LIST", String[], false, "explicit averaging factors m (τ = m·τ₀), e.g. 1,2,4,10; replaces --grid. For check --acf, the one averaging factor of the autocorrelation (default 1)", ["analyze", "check", "cross", "steer", "hat", "report", "live"]),
    Opt("--tau-min", "SECONDS", String[], false, "leave out averaging times below SECONDS", ["analyze", "report"]),
    Opt("--tau-max", "SECONDS", String[], false, "leave out averaging times above SECONDS", ["analyze", "report"]),
    Opt("--confidence", "P", String[], false, "confidence level of the intervals (default 0.683); for holdover, the share of trials the quantile covers (default 0.95)", ["analyze", "hat", "report", "holdover"]),
    Opt("--no-ci", "", String[], false, "skip confidence intervals, EDF and noise identification", ["analyze", "report"]),
    Opt("--bias", "MODE", BIASES, false, "bias correction of the total estimators and Thêo1: exact (default), published (the published tables only), stable32 (what Stable32 applies to each) or none", ["analyze", "report"]),
    Opt("--no-bias-correction", "", String[], false, "the same as --bias none", ["analyze", "report"]),
    Opt("--noise", "TYPE", NOISES, false, "noise type behind the intervals and bias factors: auto (identified at each τ, the default), fit (the strongest term of the power-law model fitted to MHDEV) or one of WHPM, FLPM, WHFM, FLFM, RWFM, FWFM, RRFM for every τ", ["analyze", "report"]),
    Opt("--spectrum", "", String[], false, "also return the spectrum (JSON only)", ["analyze"]),
    Opt("--fit", "", String[], false, "also return the power-law noise model fitted to MHDEV (JSON only); with holdover, the error that model predicts and each noise type's share", ["analyze", "holdover"]),
    Opt("--component", "NOISE:A:B", String[], false, "build the --fit model by hand: fit NOISE (WHPM, FLPM, WHFM, FLFM, RWFM, FWFM, RRFM or an α) on MHDEV between A and B seconds (either may be empty), after the components before it; repeat in order", ["analyze"]),
    Opt("--mask", "NAME", MASKS, false, "check MTIE and TDEV against an ITU-T wander mask, or all of them: G811, G8262_opt1, G8272_PRTC_A, G8272_PRTC_B, G823_SSU ..., all; JSON gains the limits and the checks, CSV a verdict on stderr", ["analyze"]),
    Opt("--nperseg", "N", String[], false, "Welch segment length (default: largest power of 2 ≤ N/4)", ["analyze", "spectrum", "cross"]),
    Opt("--with", "FILE2", String[], true, "the second channel, read with the same options as FILE", ["cross"]),
    Opt("--dev", "NAME", DYNAMIC_DEVS, false, "dadev (default) or dhdev", ["dynamic"]),
    Opt("--window", "N", String[], false, "window length in samples (default: largest power of 2 ≤ N/8)", ["dynamic"]),
    Opt("--step", "N", String[], false, "window step in samples (default window/2)", ["dynamic"]),
    Opt("--value-cols", "LIST", String[], false, "1-based value columns, one clock each, all read against the same reference (with --synth: a simulated ensemble, see --clocks)", ["timescale", "hat"]),
    Opt("--pairs", "", String[], false, "the three --value-cols are the differences A−B, B−C and C−A, each measured on its own, not clocks read against one reference", ["hat"]),
    Opt("--no-reference", "", String[], false, "leave the reference out: the hat of the clocks' differences from one another alone (three or more clocks)", ["hat"]),
    Opt("--covariance", "", String[], false, "the Groslambert covariance of the difference records in place of the hat: ADEV, with the noise of each --pairs record's instrument kept out", ["hat"]),
    Opt("--klts", "", String[], false, "--covariance with the interval of Lantz et al. (2019), for three clocks, given also where a clock is not resolved", ["hat"]),
    Opt("--algorithm", "ALG", TIMESCALES, false, "at1 (default, inverse-variance weights), at2 (at1 with a Kalman frequency update), ewfa (equal weights), or greenhall / jones_tryon / stein (Kalman scales on clock models fitted from the data; need 3+ clocks)", ["timescale"]),
    Opt("--pivot", "K", String[], false, "clock whose reading gives the greenhall, jones_tryon or stein scale, 1 to the number of clocks (default 1)", ["timescale"]),
    Opt("--steps", "MODE", STEPS, false, "at2 only: detect clock frequency steps; best recomputes the scale from the step, real keeps the values already issued and corrects forward", ["timescale"]),
    Opt("--flicker", "", String[], false, "greenhall and jones_tryon only: keep each clock's flicker FM and flicker PM as Markov states where the fit finds them, instead of folding flicker PM into the reading noise", ["timescale"]),
    Opt("--wmax", "W", String[], false, "largest weight one clock may take, 0 to 1 (default 1, no cap)", ["timescale"]),
    Opt("--train", "SECONDS", String[], false, "training stretch the frequency is learned on (default: a quarter of the record); for steer, the stretch the clock model is fitted to, with results from its end on", ["holdover", "steer"]),
    Opt("--grow", "K", String[], false, "in place of --train: a stretch K times each holdover time (2, or 3 with --fit-drift, is the test behind TDEV and HTDEV)", ["holdover"]),
    Opt("--tau", "LIST", String[], false, "holdover times in seconds, e.g. 60,3600,14400 (default: the octave multiples of τ₀ up to the training stretch); for phasenoise, the averaging times (default: 1, 2 and 5 in each decade from 1/f_h to 0.1/f₁, for a band from f₁ to f_h)", ["holdover", "phasenoise"]),
    Opt("--fit-drift", "", String[], false, "fit a parabola to the stretch, so the frequency drift is learned and carried on (default: a line); for steer, a drift state in the clock model, steered out", ["holdover", "steer"]),
    Opt("--no-anchor", "", String[], false, "start from the fit's value at the end of the stretch, not from the last reading", ["holdover"]),
    Opt("--kalman", "", String[], false, "predict with a Kalman filter in place of the least-squares fit; its clock model is fitted to the first training stretch alone, and `claimed` is the error the filter itself expects", ["holdover"]),
    Opt("--limit", "SECONDS", String[], false, "also give the holdover time at which the error reaches SECONDS (1e-6 for a microsecond)", ["holdover"]),
    Opt("--time-constant", "SECONDS", String[], false, "time constant of a critically damped loop: gains by pole placement", ["steer"]),
    Opt("--integral", "", String[], false, "with --time-constant: a third gain, on the running sum of the phase", ["steer"]),
    Opt("--cost", "P:F:U", String[], false, "in place of --time-constant: the gains that minimize P·phase² + F·frequency² + U·steer² (LQG), e.g. 1e-3:1e-3:1e12", ["steer"]),
    Opt("--interval", "SECONDS", String[], false, "seconds between steers, a whole multiple of τ₀ (default τ₀)", ["steer"]),
    Opt("--quantity", "Q", QUANTITIES, false, "phase or frequency: the samples the statistics, --histogram and --acf describe (default: those of the record)", ["check"]),
    Opt("--histogram", "", String[], false, "the histogram of the samples; with CSV, its bins in place of the summary", ["check"]),
    Opt("--bins", "N", String[], false, "number of histogram bins (default: Scott's rule, at most 1000)", ["check"]),
    Opt("--acf", "", String[], false, "the autocorrelation of the samples at averaging factor --m; with CSV, its lags in place of the summary", ["check"]),
    Opt("--lags", "K", String[], false, "largest autocorrelation lag (default: half the series, at most 1000)", ["check"]),
    Opt("--model", "MODEL", DRIFT_MODELS, false, "linear (default), quadratic, log (a·ln(1+b·t)+c) or diffusion (a+b·√(t+c))", ["drift"]),
    Opt("--method", "METHOD", DRIFT_METHODS, false, "estimator of the linear model: frequency_fit (default, a least-squares line through the frequency), phase_fit (a least-squares parabola through the phase), bisection (the three-point fit), second_difference (mean second difference of the phase) or four_point", ["drift"]),
    Opt("--carrier", "HZ", String[], false, "carrier frequency ν₀ in Hz, which ℒ(f) and S_φ(f) need", ["phasenoise"]),
    Opt("--density", "S", DENSITIES, false, "the spectral density that the table holds: L (default, ℒ(f) in dBc/Hz), Sphi (S_φ(f) in rad²/Hz), Sx (S_x(f) in s²/Hz) or Sy (S_y(f) in 1/Hz)", ["phasenoise"]),
    Opt("--band", "A:B", String[], false, "integrate from A to B hertz only, inside the table (default: the whole table)", ["phasenoise"]),
    # Live
    Opt("--hz", "", String[], false, "readings are counter frequencies in Hz, against --nominal", ["live"]),
    Opt("--query", "CMD", String[], false, "SCPI query sent before every reading, e.g. READ? (default: READ? for tcp://, none otherwise)", ["live"]),
    Opt("--every", "N", String[], false, "write a line every N readings (default 1)", ["live"]),
    Opt("--count", "N", String[], false, "stop after N readings (default: when the source ends)", ["live"]),
    Opt("--fast", "", String[], false, "run the simulator flat out instead of one reading per τ₀", ["live"]),
    # Output
    Opt("--series", "N", String[], false, "include up to N record points in the JSON (default 0)", TABULAR),
    Opt("--format", "FMT", FORMATS, false, "csv (default) or json; for model, json alone", [TABULAR; "model"; "phasenoise"]),
    Opt("--out", "PATH", String[], true, "write to PATH instead of stdout", BATCH),
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
    keep_zeros::Bool
    remove_outliers::Bool
    devs::Vector{Symbol}
    grid::String
    m::Vector{Int}         # empty = use the grid
    ci::Bool
    confidence::Float64
    bias::String
    noise::String
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
    method::String            # drift: estimator of the linear model
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
    value_cols::Vector{Int}   # timescale: one value column per clock
    algorithm::String
    wmax::Float64
    components::Vector{Tuple{Int,Float64,Float64}}   # --component: α, τ from, τ to (NaN = open)
    clocks::Int               # --clocks: simulated clocks; 0 = not given
    clock::Int                # --clock: which simulated clock to read
    models::Vector{NTuple{4,Float64}}   # --noise-model: R, q1, q2, q3 per simulated clock
    mask::String              # --mask: ITU-T mask checked in the JSON; "" = none
    pivot::Int                # --pivot: clock the Kalman timescales read at (1-based)
    steps::String             # --steps: at2 frequency-step detection, "" = off
    train::Float64            # holdover: training stretch in seconds; NaN = a quarter of the record
    taus::Vector{Float64}     # holdover: holdover times in seconds; empty = the octave grid
    fit_drift::Bool           # holdover: parabola instead of a line
    anchor::Bool              # holdover: start from the last reading
    limit::Float64            # holdover: NaN = none
    grow::Int                 # holdover: the stretch as a multiple of each holdover time; 0 = --train
    kalman::Bool              # holdover: a Kalman filter as the predictor
    flicker::Bool             # --flicker: flicker states in the Kalman scales' clock models
    time_constant::Float64    # steer: NaN = not given
    integral::Bool            # steer: an integral gain
    cost::Vector{Float64}     # steer: --cost weights on phase, frequency and steer; empty = not given
    interval::Float64         # steer: seconds between steers; NaN = τ₀
    quantity::String          # check: "phase" or "frequency"; "" = the record's
    histogram::Bool           # check: --histogram
    bins::Int                 # check: histogram bins; 0 = Scott's rule
    acf::Bool                 # check: --acf
    lags::Int                 # check: largest ACF lag; -1 = the default
    pairs::Bool               # hat: the columns are the pairwise differences
    reference::Bool           # hat: the reference counts as a clock
    covariance::Bool          # hat: the Groslambert covariance in place of the hat
    klts::Bool                # hat: the covariance with the KLTS interval
    carrier::Float64          # phasenoise: carrier frequency in Hz; NaN = not given
    density::String           # phasenoise: the spectral density that the table holds
    band::Vector{Float64}     # phasenoise: --band limits in Hz; empty = the whole table
end

Options() = Options("", "", 0, 42, false, NaN, 1, 2, 0, 0, '\0', 1.0, "none", false, false, false,
                    collect(Symbol, DEFAULT_DEVIATIONS), "octave", Int[], true,
                    Float64(DEFAULT_CONFIDENCE), "exact", "auto", false, 0, "dadev", 0, 0, 0, "csv", "", false, false,
                    "none", NaN, NaN, "auto", "", NaN, false, falses(0), "", false, "linear", "frequency_fit",
                    false, "\0", 1, 0, false, NaN, NaN, NaN, NaN, false, Int[], "at1", 1.0,
                    Tuple{Int,Float64,Float64}[], 0, 1, NTuple{4,Float64}[], "", 1, "",
                    NaN, Float64[], false, true, NaN, 0, false, false,
                    NaN, false, Float64[], NaN, "", false, 0, false, -1, false, true, false, false, NaN, "L", Float64[])

quantity(o::Options) = o.freq ? "frequency" : "phase"

# --bias as the library's `correct_bias`.
correct_bias(o::Options) = o.bias == "exact" ? true : o.bias == "none" ? false : Symbol(o.bias)

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

# "R:Q1:Q2:Q3", all finite and not negative.
function parse_noise_model!(o::Options, s::String)
    parts = split(s, ':')
    length(parts) == 4 || throw(usage_error("--noise-model expects R:Q1:Q2:Q3, got \"" * s * "\""))
    v = Float64[parse_float(String(p), "--noise-model") for p in parts]
    all(>=(0), v) || throw(usage_error("--noise-model levels must not be negative"))
    push!(o.models, (v[1], v[2], v[3], v[4]))
end

# "P:F:U": weights on phase, frequency and steer.
function parse_cost(s::String)
    parts = split(s, ':')
    length(parts) == 3 || throw(usage_error("--cost expects P:F:U, got \"" * s * "\""))
    v = Float64[parse_float(String(p), "--cost") for p in parts]
    v[1] > 0 && v[2] >= 0 && v[3] > 0 || throw(usage_error("--cost: the phase and steer weights must be positive and the frequency weight not negative"))
    return v
end

# "NOISE:A:B", A and B in seconds; an empty side is open.
function parse_component!(o::Options, s::String)
    parts = split(s, ':')
    length(parts) == 3 || throw(usage_error("--component expects NOISE:A:B, got \"" * s * "\""))
    name = String(parts[1])
    α = tryparse(Int, name)
    if α === nothing
        k = findfirst(==(uppercase(name)), NOISE_NAMES)
        k === nothing && throw(usage_error("--component: unknown noise \"" * name * "\"; use " * join(NOISE_NAMES, ", ") * " or an α"))
        α = 3 - k
    end
    a = isempty(parts[2]) ? NaN : parse_float(String(parts[2]), "--component")
    b = isempty(parts[3]) ? NaN : parse_float(String(parts[3]), "--component")
    push!(o.components, (α, a, b))
    o.fit = true
end

# "A:B" in hertz, both given.
function parse_band(s::String)
    parts = split(s, ':')
    length(parts) == 2 || throw(usage_error("--band expects A:B in hertz, got \"" * s * "\""))
    v = Float64[parse_float(String(p), "--band") for p in parts]
    0 < v[1] < v[2] || throw(usage_error("--band expects two rising frequencies above 0 Hz"))
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

# "2,3,5": value columns in the order given.
function parse_cols(s::String)
    cols = Int[]
    for tok in split(s, ',')
        v = tryparse(Int, strip(tok))
        (v === nothing || v < 1) && throw(usage_error("--value-cols expects column numbers, got \"" * s * "\""))
        push!(cols, v)
    end
    length(cols) >= 2 && allunique(cols) || throw(usage_error("--value-cols needs two or more different columns"))
    return cols
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

# "60,3600,14400": holdover times in seconds, rising.
function parse_taus(s::String)
    taus = Float64[]
    for tok in split(s, ',')
        v = parse_float(String(strip(tok)), "--tau")
        v > 0 || throw(usage_error("--tau expects positive times in seconds"))
        push!(taus, v)
    end
    return sort!(unique!(taus))
end

function parse_delim(s::String)
    s == "\\t" && return '\t'
    length(s) == 1 || throw(usage_error("--delim expects one character, got \"" * s * "\""))
    return s[1]
end

function set_option!(o::Options, name::String, v::String)
    if name == "--synth";              o.synth = parse_int(v, name)
    elseif name == "--seed";           o.seed = parse_int(v, name)
    elseif name == "--clocks";         o.clocks = parse_int(v, name)
    elseif name == "--clock";          o.clock = parse_int(v, name)
    elseif name == "--noise-model";    parse_noise_model!(o, v)
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
    elseif name == "--bias";           o.bias = check_choice(v, name, BIASES)
    elseif name == "--noise";          o.noise = check_choice(v, name, NOISES)
    elseif name == "--devs";           o.devs = parse_devs(v)
    elseif name == "--grid";           tau_mode(v); o.grid = v
    elseif name == "--m";              o.m = parse_m(v)
    elseif name == "--confidence";     o.confidence = parse_float(v, name)
    elseif name == "--range";          parse_range!(o, v)
    elseif name == "--component";      parse_component!(o, v)
    elseif name == "--tau-min";        o.tau_min = parse_float(v, name)
    elseif name == "--tau-max";        o.tau_max = parse_float(v, name)
    elseif name == "--nperseg";        o.nperseg = parse_int(v, name)
    elseif name == "--dev";            o.dev = check_choice(v, name, DYNAMIC_DEVS)
    elseif name == "--window";         o.window = parse_int(v, name)
    elseif name == "--step";           o.step = parse_int(v, name)
    elseif name == "--model";          o.model = check_choice(v, name, DRIFT_MODELS)
    elseif name == "--quantity";       o.quantity = check_choice(v, name, QUANTITIES)
    elseif name == "--bins";           o.bins = parse_int(v, name); o.bins >= 1 || throw(usage_error("--bins must be 1 or more"))
    elseif name == "--lags";           o.lags = parse_int(v, name); o.lags >= 0 || throw(usage_error("--lags must be 0 or more"))
    elseif name == "--method";         o.method = check_choice(v, name, DRIFT_METHODS)
    elseif name == "--train";          o.train = parse_float(v, name)
    elseif name == "--grow";           o.grow = parse_int(v, name)
    elseif name == "--tau";            o.taus = parse_taus(v)
    elseif name == "--limit";          o.limit = parse_float(v, name)
    elseif name == "--time-constant";  o.time_constant = parse_float(v, name)
    elseif name == "--cost";           o.cost = parse_cost(v)
    elseif name == "--interval";       o.interval = parse_float(v, name)
    elseif name == "--carrier";        o.carrier = parse_float(v, name)
    elseif name == "--density";        o.density = check_choice(v, name, DENSITIES)
    elseif name == "--band";           o.band = parse_band(v)
    elseif name == "--value-cols";     o.value_cols = parse_cols(v)
    elseif name == "--algorithm";      o.algorithm = check_choice(v, name, TIMESCALES)
    elseif name == "--wmax";           o.wmax = parse_float(v, name)
    elseif name == "--pivot";          o.pivot = parse_int(v, name)
    elseif name == "--steps";          o.steps = check_choice(v, name, STEPS)
    elseif name == "--series";         o.series = parse_int(v, name)
    elseif name == "--query";          o.query = v
    elseif name == "--every";          o.every = parse_int(v, name)
    elseif name == "--count";          o.count = parse_int(v, name)
    elseif name == "--format"
        # Refused here, so that a later error on the line is still answered in JSON.
        o.command == "model" && v != "json" && throw(usage_error("model writes JSON only"))
        o.format = check_choice(v, name, FORMATS)
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
    elseif name == "--keep-zeros";         o.keep_zeros = true
    elseif name == "--remove-outliers";    o.remove_outliers = true
    elseif name == "--no-ci";              o.ci = false
    elseif name == "--no-bias-correction"; o.bias = "none"
    elseif name == "--spectrum";           o.spectrum = true
    elseif name == "--hz";                 o.hz = true
    elseif name == "--fast";               o.fast = true
    elseif name == "--fit";                o.fit = true
    elseif name == "--histogram";          o.histogram = true
    elseif name == "--acf";                o.acf = true
    elseif name == "--fit-drift";          o.fit_drift = true
    elseif name == "--no-anchor";          o.anchor = false
    elseif name == "--kalman";             o.kalman = true
    elseif name == "--integral";           o.integral = true
    elseif name == "--flicker";            o.flicker = true
    elseif name == "--pairs";              o.pairs = true
    elseif name == "--no-reference";       o.reference = false
    elseif name == "--covariance";         o.covariance = true
    elseif name == "--klts";               o.klts = true
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
    elseif o.command == "model"
        isempty(o.input) && o.synth > 0 && !isempty(o.models) ||
            throw(usage_error("model: give --synth N and a clock with --noise-model R:Q1:Q2:Q3, and no FILE"))
    elseif o.command == "phasenoise"
        isempty(o.input) && throw(usage_error("phasenoise: give the TABLE of the spectrum (see `sigmatau help phasenoise`)"))
        o.density in ("L", "Sphi") && isnan(o.carrier) && throw(usage_error("phasenoise: ℒ(f) and S_φ(f) need the carrier frequency, --carrier HZ"))
        isnan(o.carrier) || o.carrier > 0 || throw(usage_error("--carrier must be positive"))
        all(in(PHASENOISE_DEVS), o.devs) || throw(usage_error("phasenoise converts to " * join(String[String(d) for d in PHASENOISE_DEVS], ", ")))
    elseif isempty(o.input) && o.synth <= 0
        throw(usage_error(o.command * ": give a FILE or --synth N (see `sigmatau help " * o.command * "`)"))
    end
    !isempty(o.input) && o.synth > 0 && throw(usage_error("give either a FILE or --synth, not both"))
    if o.synth <= 0
        (o.clocks != 0 || o.clock != 1 || !isempty(o.models)) && throw(usage_error("--clocks, --clock and --noise-model go with --synth"))
    else
        nm = length(o.models)
        nm <= 64 || throw(usage_error("--noise-model: at most 64 clocks"))
        nm > 0 && o.clocks != 0 && o.clocks != nm && throw(usage_error("--clocks is " * string(o.clocks) * " but " * string(nm) * " --noise-model were given"))
        nm == 0 && o.clocks != 0 && !(2 <= o.clocks <= 64) && throw(usage_error("--clocks must be between 2 and 64"))
        1 <= o.clock <= max(o.clocks, nm, 1) || throw(usage_error("--clock must be between 1 and the number of clocks"))
        o.clocks == 0 && nm == 0 && o.clock != 1 && throw(usage_error("--clock needs --clocks"))
    end
    o.command == "cross" && (isempty(o.with) || o.synth > 0) &&
        throw(usage_error("cross: give FILE and the second channel with --with FILE2"))
    0 < o.confidence < 1 || throw(usage_error("--confidence must be between 0 and 1"))
    o.method == "frequency_fit" || o.model == "linear" ||
        throw(usage_error("--method " * o.method * " estimates a linear drift and goes with --model linear"))
    (isnan(o.train) || o.train > 0) && (isnan(o.limit) || o.limit > 0) ||
        throw(usage_error("--train and --limit must be positive"))
    o.grow == 0 || isnan(o.train) && o.grow > 1 + o.fit_drift ||
        throw(usage_error("--grow goes without --train and is 2 or more, 3 or more with --fit-drift"))
    o.kalman && (o.fit || !o.anchor || o.grow > 0) &&
        throw(usage_error("--kalman goes with none of --fit, --no-anchor and --grow"))
    if o.command == "steer"
        isnan(o.time_constant) == isempty(o.cost) && throw(usage_error("steer: give one of --time-constant SECONDS and --cost P:F:U"))
        (isnan(o.time_constant) || o.time_constant > 0) && (isnan(o.interval) || o.interval > 0) ||
            throw(usage_error("--time-constant and --interval must be positive"))
        o.integral && isnan(o.time_constant) && throw(usage_error("--integral goes with --time-constant"))
    end
    if o.command == "timescale"
        0 < o.wmax <= 1 || throw(usage_error("--wmax must be above 0 and at most 1"))
        isempty(o.steps) || o.algorithm == "at2" || throw(usage_error("--steps needs --algorithm at2"))
    end
    o.nperseg >= 0 || throw(usage_error("--nperseg must be positive"))
    o.window >= 0 && o.step >= 0 || throw(usage_error("--window and --step must be positive"))
    o.series >= 0 || throw(usage_error("--series must be 0 or more"))
    o.reference_col >= 0 || throw(usage_error("--reference-col must be 1 or more"))
    o.time_col >= 0 || throw(usage_error("--time-col must be 0 (no time column) or more"))
    o.value_col >= 1 || throw(usage_error("--value-col must be 1 or more"))
    o.header >= 0 || throw(usage_error("--header must be 0 or more"))
    if !isempty(o.value_cols)
        isempty(o.column) || throw(usage_error("--column names one clock; give the clocks of a timescale with --value-cols"))
        o.value_col = o.value_cols[1]
    end
    isnan(o.tau0) || o.tau0 > 0 || throw(usage_error("--tau0 must be positive"))
    isnan(o.nominal) || o.nominal > 0 || throw(usage_error("--nominal must be positive"))
    o.scaling != 0 || throw(usage_error("--scaling must not be 0"))
    o.spectrum && o.format != "json" &&
        throw(usage_error("--spectrum adds to the JSON output; use --format json, or `sigmatau spectrum`"))
    if o.command == "check"
        o.bins == 0 || o.histogram || throw(usage_error("--bins goes with --histogram"))
        (o.lags < 0 && isempty(o.m)) || o.acf || throw(usage_error("--lags and --m go with --acf"))
        length(o.m) > 1 && throw(usage_error("check --acf takes one averaging factor --m"))
        o.histogram && o.acf && o.format != "json" &&
            throw(usage_error("--histogram and --acf each print a table; ask for one, or use --format json"))
    end
    o.fit && o.format != "json" && throw(usage_error("--fit and --component add to the JSON output; use --format json"))
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
        sig = isempty(arg) ? c : c * " " * arg
        print(io, "  ")
        # A signature that leaves no gap before the description takes a line of its own.
        length(sig) <= 20 ? padded(io, sig, 22) : (print(io, sig, "\n  "); padded(io, "", 22))
        print(io, h, "\n")
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
    print(io, "  sigmatau holdover gpsdo.dat --train 86400 --tau 3600,14400,86400 --limit 1e-6\n")
    print(io, "  sigmatau steer maser.dat --time-constant 86400 --format json\n")
    print(io, "  sigmatau hat clocks.dat --value-cols 2,3 --devs adev,hdev\n")
    print(io, "  sigmatau phasenoise ocxo.csv --carrier 10e6 --devs adev,mdev\n")
    print(io, "  sigmatau live tcp://192.168.1.20:5025 --hz --nominal 10e6 --tau0 1\n")
    print(io, "  sigmatau report clock.dat --devs adev,tdev --out report.html\n")
end

function print_command_help(io::IO, cmd::String)
    if cmd == "completion"
        print(io, "Usage: sigmatau completion SHELL\n\nShells: ", join(SHELLS, ", "), "\n")
        return
    end
    for (c, arg, h) in COMMAND_HELP
        c == cmd && print(io, "Usage: sigmatau ", c, " ", arg, " [options]\n\n", uppercasefirst(h), ".\n")
    end
    if cmd == "live"
        print(io, "SOURCE is sim (a simulated clock), tcp://HOST:PORT (an LXI instrument's SCPI socket,\n",
              "usually port 5025) or a serial device path set up with stty. Each reading is the first\n",
              "number on a line: phase in seconds by default. Every line of output is one JSON object.\n\nOptions:\n")
    elseif cmd == "model"
        print(io, "The averaging times are the octave multiples of τ₀ that a record of N samples reaches.\n",
              "The output is JSON.\n\nOptions:\n")
    elseif cmd == "phasenoise"
        print(io, "TABLE is a whitespace, comma or tab separated table: the Fourier frequency in Hz, then the\n",
              "spectrum. Between its points the spectrum is a straight line on a log-log plot. Each result is\n",
              "the deviation of the noise inside the band of the table. `below` and `above` are what the\n",
              "spectrum would add if it went on below and above the band on the slope of its lowest and\n",
              "highest decade, as fractions of the variance found: blank, or null in JSON, where that has no\n",
              "limit, as above the band for ADEV and HDEV of white or flicker PM.\n\nOptions:\n")
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
    elseif cmd == "hat"
        print(io, "\nDeviations: ", join(HAT_DEVS, ", "), "\n")
    elseif cmd == "phasenoise"
        print(io, "\nDeviations: ", join(String[String(d) for d in PHASENOISE_DEVS], ", "), "\n")
    end
end

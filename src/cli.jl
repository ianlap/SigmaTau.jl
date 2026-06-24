# cli.jl — non-interactive command-line interface for SigmaTau.
#
# Public entry point: `SigmaTau.cli_main(args)::Cint` (intentionally NOT exported;
# reachable as `SigmaTau.cli_main`). `bin/sigmatau.jl` is a three-line launcher
# that wraps it. All output flows through injected `out`/`err` IO streams so the
# CLI is fully testable in-process and PackageCompiler-friendly. `cli_main` never
# calls `exit` — it always returns a `Cint` exit code; only the launcher exits.
#
# This is a thin wrapper over the existing public API (`read_phase` /
# `read_frequency` / `stability`); it touches no numerical kernel. There is NO
# interactive / TUI / REPL behavior — output is plain stdout, errors go to stderr.

# ──────────────────────────────────────────────────────────────────────────
# ── Centralized metadata (single source of truth) ─────────────────────────
#
# These tables drive (a) argument parsing/validation, (b) `--help` rendering,
# and (c) all four shell-completion generators, so the command/option/value
# vocabularies can never drift between them.

"""
Deviation symbols that `analyze` (via [`stability`](@ref)) can compute — exactly
the `_dispatch_dev` set in `suite.jl`. Used for `--devs` validation and the
"valid deviations" error list.
"""
const SUITE_DEVS = (:adev, :mdev, :tdev, :hdev, :mhdev, :htdev,
                    :totdev, :mtotdev, :ttotdev, :htotdev, :mhtotdev,
                    :mtie, :tierms, :pdev, :theo1, :theoh)

# Dynamic (time-resolved) deviations. Offered in the shell-completion vocabulary
# (the user asked for them there) but REJECTED by `analyze`: they return a
# two-axis σ_y(t,τ) `DynamicStabilityResult` — visualized as a waterfall or a
# 3-D surface of σ vs τ over time, NOT a heatmap — which does not fit the per-τ
# CSV/JSON table schema. `analyze --devs ...,dadev` exits 1 with a pointer to the
# dynamic Julia API.
const DYNAMIC_DEVS = (:dadev, :dhdev)

# Full vocabulary offered for `--devs` completion (static first, dynamic last).
const COMPLETION_DEVS = (SUITE_DEVS..., DYNAMIC_DEVS...)
const COMPLETION_DEVS_STR = Tuple(string(d) for d in COMPLETION_DEVS)

# `--tau` grid tokens → TauMode enum values.
const TAU_MODE_MAP = Dict(
    "octave"        => Octave,
    "decade"        => Decade,
    "alltaus"       => AllTaus,
    "halfoctave"    => HalfOctave,
    "quarteroctave" => QuarterOctave,
    "halfdecade"    => HalfDecade,
)
const TAU_CHOICES    = ("octave", "decade", "alltaus", "halfoctave", "quarteroctave", "halfdecade")
const TAU_COMPLETE   = ("octave", "decade")          # terse completion set (per spec)
const FORMAT_CHOICES = ("csv", "json")
const BOOL_CHOICES   = ("true", "false")
const DETREND_CHOICES = ("none", "mean", "endpoint", "linear")

"""Descriptor for one `analyze` option, consumed by the parser, help, and completions."""
struct CliOpt
    long::String                    # e.g. "--devs"
    takes_arg::Bool                 # true for value options, false for flags
    is_path::Bool                   # true when the value is a file/dir (drives file completion)
    choices::Tuple{Vararg{String}}  # () = free value; non-empty = enumerated candidates
    metavar::String                 # "LIST"/"N"/"FLOAT"/"PATH"; "" for flags
    help::String
end

"""Top-level sub-command descriptor."""
struct CliCmd
    name::String
    help::String
end

# Every `analyze` option. `--header/--delim/--scaling/--detrend` map 1:1 to the
# existing `read_phase`/`read_frequency` keyword arguments.
const ANALYZE_OPTS = (
    CliOpt("--phase",      false, false, (),                 "",      "treat FILE as phase data (default)"),
    CliOpt("--freq",       false, false, (),                 "",      "treat FILE as fractional-frequency data"),
    CliOpt("--tau0",       true,  false, (),                 "FLOAT", "base sample interval τ₀ in seconds"),
    CliOpt("--time-col",   true,  false, (),                 "N",     "1-based time column (0 = none)"),
    CliOpt("--value-col",  true,  false, (),                 "N",     "1-based value column"),
    CliOpt("--devs",       true,  false, COMPLETION_DEVS_STR, "LIST", "comma list of deviations (default adev,mdev,hdev,mhdev)"),
    CliOpt("--tau",        true,  false, TAU_CHOICES,        "GRID",  "tau-grid spacing (default octave; exclusive with --taus)"),
    CliOpt("--taus",       true,  false, (),                 "LIST",  "explicit comma list of integer averaging factors m (τ = m·τ₀)"),
    CliOpt("--ci",         true,  false, BOOL_CHOICES,       "BOOL",  "compute confidence intervals (default true)"),
    CliOpt("--confidence", true,  false, (),                 "FLOAT", "CI confidence level in (0,1) (default 0.683)"),
    CliOpt("--format",     true,  false, FORMAT_CHOICES,     "FMT",   "output format: csv or json (default csv)"),
    CliOpt("--out",        true,  true,  (),                 "PATH",  "write output to PATH (default: stdout)"),
    CliOpt("--header",     true,  false, (),                 "N",     "non-comment header lines to skip"),
    CliOpt("--delim",      true,  false, (),                 "CHAR",  "field delimiter (default: auto from extension)"),
    CliOpt("--scaling",    true,  false, (),                 "FLOAT", "multiply samples by this factor"),
    CliOpt("--detrend",    true,  false, DETREND_CHOICES,    "MODE",  "detrend before analysis: none|mean|endpoint|linear"),
    CliOpt("--quiet",      false, false, (),                 "",      "suppress informational messages"),
    CliOpt("--debug",      false, false, (),                 "",      "print stack traces on error"),
)
const _OPT_BY_LONG = Dict(o.long => o for o in ANALYZE_OPTS)

const COMMANDS = (
    CliCmd("analyze",    "compute stability deviations from a data file"),
    CliCmd("completion", "print a shell completion script (bash|zsh|fish|powershell)"),
    CliCmd("version",    "print the SigmaTau version"),
    CliCmd("help",       "print this help message"),
)
const SHELLS = ("bash", "zsh", "fish", "powershell")

# Space-joined vocabularies for the completion-script generators.
const _SUBCMDS_SP      = join((c.name for c in COMMANDS), " ")
const _ANALYZE_OPTS_SP = join((o.long for o in ANALYZE_OPTS), " ")
const _DEVS_SP         = join(COMPLETION_DEVS_STR, " ")
const _SHELLS_SP       = join(SHELLS, " ")
const _TAU_SP          = join(TAU_COMPLETE, " ")
const _FMT_SP          = join(FORMAT_CHOICES, " ")
const _BOOL_SP         = join(BOOL_CHOICES, " ")
const _DETREND_SP      = join(DETREND_CHOICES, " ")

# ──────────────────────────────────────────────────────────────────────────
# ── Errors ────────────────────────────────────────────────────────────────

"""
Internal CLI exception carrying the process exit code. `showvalid` requests the
"valid deviations" list be appended to the message (for unknown-deviation errors).
"""
struct CliError <: Exception
    code::Cint
    msg::String
    showvalid::Bool
end
CliError(code::Integer, msg::AbstractString; showvalid::Bool=false) =
    CliError(Cint(code), String(msg), showvalid)

# Human-readable message for an arbitrary exception (unwrap ArgumentError's text).
_emsg(e) = e isa ArgumentError ? e.msg : sprint(showerror, e)

# A downstream consumer closing the pipe early (e.g. `sigmatau analyze … | head`)
# surfaces as an EPIPE IOError. That is the normal Unix idiom, not a failure, so
# we exit cleanly rather than printing a spurious "output error".
_is_broken_pipe(e) = e isa Base.IOError &&
    (e.code == Base.UV_EPIPE || occursin("EPIPE", e.msg) || occursin("broken pipe", e.msg))

# Print a failure and return its exit code. A `CliError` keeps its own code;
# any other exception is wrapped with `defaultcode`/`prefix` (the per-stage
# mapping). Stack traces print only under `--debug`.
function _fail(e, err::IO, debug::Bool, bt, defaultcode::Integer, prefix::AbstractString)
    ce = e isa CliError ? e : CliError(defaultcode, string(prefix, _emsg(e)))
    cause = e isa CliError ? nothing : e
    println(err, "sigmatau: ", ce.msg)
    ce.showvalid && println(err, "  valid deviations: ", join(string.(SUITE_DEVS), ", "))
    if debug
        if cause !== nothing
            showerror(err, cause)
            println(err)
        end
        for frame in stacktrace(bt)
            println(err, "  ", frame)
        end
    end
    return ce.code
end

# ──────────────────────────────────────────────────────────────────────────
# ── Argument value parsers (all throw CliError(1, …) on bad input) ─────────

_split_eq(a::AbstractString) =
    (j = findfirst('=', a)) === nothing ? (a, nothing) : (a[1:j-1], a[j+1:end])

function _parse_bool(name, v)
    v == "true"  && return true
    v == "false" && return false
    throw(CliError(1, "$name expects true or false, got '$v'"))
end

function _parse_float(name, v)
    x = tryparse(Float64, v)
    x === nothing && throw(CliError(1, "$name expects a number, got '$v'"))
    return x
end

function _parse_int(name, v)
    x = tryparse(Int, v)
    x === nothing && throw(CliError(1, "$name expects an integer, got '$v'"))
    return x
end

function _parse_devs(v)
    syms = Symbol[]
    for tok in split(v, ',')
        t = strip(tok)
        isempty(t) && continue
        s = Symbol(t)
        if s in SUITE_DEVS
            push!(syms, s)
        elseif s in DYNAMIC_DEVS
            throw(CliError(1,
                "'$t' is a dynamic (time-resolved) deviation σ_y(t,τ) and cannot be " *
                "computed by `analyze`; use the SigmaTau `$t` Julia API (it returns a " *
                "DynamicStabilityResult — a waterfall / 3-D surface of σ vs τ over time)."))
        else
            throw(CliError(1, "unknown deviation '$t'"; showvalid=true))
        end
    end
    isempty(syms) && throw(CliError(1, "--devs: no deviations given"))
    return Tuple(syms)
end

function _parse_taus(v)
    ms = Int[]
    for tok in split(v, ',')
        t = strip(tok)
        isempty(t) && continue
        m = tryparse(Int, t)
        (m === nothing || m < 1) &&
            throw(CliError(1, "--taus expects positive integer averaging factors, got '$t'"))
        push!(ms, m)
    end
    isempty(ms) && throw(CliError(1, "--taus: no averaging factors given"))
    return ms
end

function _parse_tau_mode(v)
    key = lowercase(v)
    haskey(TAU_MODE_MAP, key) ||
        throw(CliError(1, "--tau: unknown grid '$v' (valid: $(join(TAU_CHOICES, ", ")))"))
    return TAU_MODE_MAP[key]
end

function _parse_detrend(v)
    v in DETREND_CHOICES ||
        throw(CliError(1, "--detrend: unknown mode '$v' (valid: $(join(DETREND_CHOICES, ", ")))"))
    return Symbol(v)
end

function _parse_delim(v)
    s = v == "\\t" ? "\t" : v       # accept the literal escape "\t" for a tab
    length(s) == 1 || throw(CliError(1, "--delim expects a single character, got '$v'"))
    return only(s)
end

# ──────────────────────────────────────────────────────────────────────────
# ── `analyze` argument parsing ─────────────────────────────────────────────

# Returns `(help=true,)` for `analyze --help`, else a fully-resolved settings
# NamedTuple. Throws `CliError(1, …)` for any usage error.
function _parse_analyze(args)
    file = nothing
    phase_set = false
    freq_set  = false
    raw   = Dict{String,String}()
    flags = Set{String}()

    i, n = 1, length(args)
    while i <= n
        a = args[i]
        if a == "-h" || a == "--help"
            return (help = true,)
        elseif a == "--phase"
            phase_set = true; i += 1
        elseif a == "--freq"
            freq_set = true; i += 1
        elseif a == "--quiet" || a == "--debug"
            push!(flags, a); i += 1
        elseif startswith(a, "--")
            name, inlineval = _split_eq(a)
            o = get(_OPT_BY_LONG, name, nothing)
            o === nothing && throw(CliError(1, "unknown option: $name"))
            o.takes_arg || throw(CliError(1, "option $name takes no value"))
            if inlineval !== nothing
                raw[name] = inlineval; i += 1
            else
                i + 1 <= n || throw(CliError(1, "option $name requires a value"))
                raw[name] = args[i + 1]; i += 2
            end
        elseif a != "-" && startswith(a, "-")
            throw(CliError(1, "unknown option: $a"))
        else
            file === nothing || throw(CliError(1, "unexpected extra argument: $a"))
            file = a; i += 1
        end
    end

    file === nothing &&
        throw(CliError(1, "analyze: missing required FILE argument (see `sigmatau analyze --help`)"))
    phase_set && freq_set && throw(CliError(1, "specify only one of --phase / --freq"))
    haskey(raw, "--tau") && haskey(raw, "--taus") &&
        throw(CliError(1, "--tau and --taus are mutually exclusive"))

    devs = haskey(raw, "--devs") ? _parse_devs(raw["--devs"]) : DEFAULT_DEVIATIONS
    taus = haskey(raw, "--taus") ? _parse_taus(raw["--taus"]) :
           haskey(raw, "--tau")  ? _parse_tau_mode(raw["--tau"]) : Octave
    ci = haskey(raw, "--ci") ? _parse_bool("--ci", raw["--ci"]) : true
    confidence = if haskey(raw, "--confidence")
        c = _parse_float("--confidence", raw["--confidence"])
        (0 < c < 1) || throw(CliError(1, "--confidence must be in (0,1), got $c"))
        c
    else
        DEFAULT_CONFIDENCE
    end
    format = get(raw, "--format", "csv")
    format in FORMAT_CHOICES ||
        throw(CliError(1, "--format: unknown format '$format' (valid: $(join(FORMAT_CHOICES, ", ")))"))

    readkw = (
        tau0      = haskey(raw, "--tau0")      ? _parse_float("--tau0", raw["--tau0"])         : nothing,
        time_col  = haskey(raw, "--time-col")  ? _parse_int("--time-col", raw["--time-col"])   : 1,
        value_col = haskey(raw, "--value-col") ? _parse_int("--value-col", raw["--value-col"]) : 2,
        header    = haskey(raw, "--header")    ? _parse_int("--header", raw["--header"])        : 0,
        delim     = haskey(raw, "--delim")     ? _parse_delim(raw["--delim"])                   : nothing,
        scaling   = haskey(raw, "--scaling")   ? _parse_float("--scaling", raw["--scaling"])    : 1.0,
        detrend   = haskey(raw, "--detrend")   ? _parse_detrend(raw["--detrend"])               : :none,
    )

    return (help = false, file = file, kind = freq_set ? :frequency : :phase,
            devs = devs, taus = taus, ci = ci, confidence = confidence,
            format = format, out = get(raw, "--out", nothing),
            quiet = ("--quiet" in flags), debug = ("--debug" in flags),
            readkw = readkw)
end

# ──────────────────────────────────────────────────────────────────────────
# ── Output writers — pure `(io, suite, source)` functions ──────────────────

# RFC-4180-quote a field that contains a comma, quote, or newline.
function _csv_field(s::AbstractString)
    any(c -> c == ',' || c == '"' || c == '\n' || c == '\r', s) || return s
    return string('"', replace(s, '"' => "\"\""), '"')
end

# Numeric cell: blank for nothing / NaN / Inf (never the literal "NaN"/"nothing").
_csv_num(x) = x === nothing ? "" : (x isa Real && !isfinite(x)) ? "" : string(x)

const _CSV_HEADER = "source,quantity,deviation,tau,dev,ci_lower,ci_upper,edf,neff,noise_type"

function _write_csv(io::IO, suite::StabilitySuite, source::AbstractString)
    println(io, _CSV_HEADER)
    q   = string(suite.data_kind)
    src = _csv_field(source)
    for r in suite.results
        dev_sym = string(r.deviation_type)
        has_ci  = !isempty(r.ci)
        has_edf = !isempty(r.edf)
        has_nt  = !isempty(r.noise_type)
        for i in eachindex(r.tau)
            ci_lo = has_ci  ? r.ci[i].lo : nothing
            ci_hi = has_ci  ? r.ci[i].hi : nothing
            edf   = has_edf ? r.edf[i]   : nothing
            nt    = has_nt  ? string(r.noise_type[i]) : ""
            println(io, src, ',', q, ',', dev_sym, ',',
                    _csv_num(r.tau[i]), ',', _csv_num(r.dev[i]), ',',
                    _csv_num(ci_lo), ',', _csv_num(ci_hi), ',', _csv_num(edf), ',',
                    string(r.neff[i]), ',', nt)
        end
    end
    return nothing
end

# Map nothing / NaN / Inf to `nothing` so JSON.print emits `null`.
_json_num(x) = (x === nothing || (x isa Real && !isfinite(x))) ? nothing : x

function _build_json(suite::StabilitySuite, source::AbstractString)
    analysis = OrderedDict{String,Any}(
        "quantity"   => string(suite.data_kind),
        "tau0"       => suite.tau0,
        "n"          => suite.n,
        "ci"         => !isnan(suite.confidence),
        "confidence" => isnan(suite.confidence) ? nothing : suite.confidence,
        "tau_mode"   => string(suite.tau_mode),
        "deviations" => [string(r.deviation_type) for r in suite.results],
    )
    results = Any[]
    for r in suite.results
        has_ci  = !isempty(r.ci)
        has_edf = !isempty(r.edf)
        has_nt  = !isempty(r.noise_type)
        points = Any[]
        for i in eachindex(r.tau)
            push!(points, OrderedDict{String,Any}(
                "tau"        => _json_num(r.tau[i]),
                "dev"        => _json_num(r.dev[i]),
                "ci_lower"   => has_ci  ? _json_num(r.ci[i].lo) : nothing,
                "ci_upper"   => has_ci  ? _json_num(r.ci[i].hi) : nothing,
                "edf"        => has_edf ? _json_num(r.edf[i])   : nothing,
                "neff"       => r.neff[i],
                "noise_type" => has_nt  ? string(r.noise_type[i]) : nothing,
            ))
        end
        push!(results, OrderedDict{String,Any}(
            "deviation" => string(r.deviation_type), "points" => points))
    end
    return OrderedDict{String,Any}(
        "schema_version" => 1,
        "source"         => source,
        "analysis"       => analysis,
        "results"        => results,
    )
end

function _write_json(io::IO, suite::StabilitySuite, source::AbstractString)
    JSON.print(io, _build_json(suite, source), 2)
    println(io)
    return nothing
end

function _write_output(settings, suite::StabilitySuite, out::IO)
    writer = settings.format == "json" ? _write_json : _write_csv
    if settings.out === nothing
        writer(out, suite, settings.file)
    else
        open(settings.out, "w") do io
            writer(io, suite, settings.file)
        end
    end
    return nothing
end

# ──────────────────────────────────────────────────────────────────────────
# ── Help ───────────────────────────────────────────────────────────────────

_pkg_version() = string(pkgversion(@__MODULE__))

_print_version(out::IO) = println(out, "sigmatau ", _pkg_version())

function _print_analyze_options(out::IO)
    println(out, "Options for 'analyze' (sigmatau analyze FILE [options]):")
    sigs = String[o.metavar == "" ? o.long : string(o.long, " ", o.metavar) for o in ANALYZE_OPTS]
    w = maximum(length, sigs)
    for (o, sig) in zip(ANALYZE_OPTS, sigs)
        println(out, "  ", rpad(sig, w), "  ", o.help)
    end
    println(out, "  Deviations (--devs): ", join(string.(SUITE_DEVS), ", "))
end

function _print_help(out::IO)
    println(out, "sigmatau — SigmaTau.jl clock-stability analysis (command-line interface)")
    println(out)
    println(out, "Usage:")
    println(out, "  sigmatau <command> [options]")
    println(out)
    println(out, "Commands:")
    for c in COMMANDS
        println(out, "  ", rpad(c.name, 12), c.help)
    end
    println(out)
    _print_analyze_options(out)
    println(out)
    println(out, "Exit codes: 0 ok · 1 usage · 2 input/read · 3 analysis · 4 output/write")
    println(out)
    println(out, "Examples:")
    println(out, "  sigmatau analyze phase.csv --devs adev,mdev --tau octave --format json")
    println(out, "  sigmatau analyze freq.txt --freq --tau0 1 --taus 1,2,4,8 --out out.csv")
    println(out, "  sigmatau completion bash > /usr/local/etc/bash_completion.d/sigmatau")
end

function _print_analyze_help(out::IO)
    println(out, "Usage: sigmatau analyze FILE [options]")
    println(out)
    println(out, "Compute stability deviations from a phase or frequency data file.")
    println(out, "By default column 1 is time and column 2 is the value (override with")
    println(out, "--time-col / --value-col). --taus values are integer averaging factors m")
    println(out, "(τ = m·τ₀); --tau selects a named geometric grid instead.")
    println(out)
    _print_analyze_options(out)
    println(out)
    println(out, "Exit codes: 0 ok · 1 usage · 2 input/read · 3 analysis · 4 output/write")
end

# ──────────────────────────────────────────────────────────────────────────
# ── Shell completion generators ────────────────────────────────────────────
#
# Each returns a STATIC script string (no file reads, no analysis at completion
# time), with vocabularies interpolated from the metadata tables so they cannot
# drift. Native file-path completion is preserved for the FILE positional and
# path options (`--out`, plus future `--config`/`--input`/`--reference`).

# Substitute @@TOKEN@@ placeholders (raw strings keep shell `$` literal).
_fill(tpl) = replace(tpl,
    "@@DEVS@@"    => _DEVS_SP,
    "@@TAU@@"     => _TAU_SP,
    "@@FMT@@"     => _FMT_SP,
    "@@BOOL@@"    => _BOOL_SP,
    "@@DETREND@@" => _DETREND_SP,
    "@@SUBCMDS@@" => _SUBCMDS_SP,
    "@@SHELLS@@"  => _SHELLS_SP,
    "@@OPTS@@"    => _ANALYZE_OPTS_SP,
)

function _completion_bash()
    tpl = raw"""
# sigmatau bash completion.
# Install:  sigmatau completion bash > /usr/local/etc/bash_completion.d/sigmatau
#      or:  source <(sigmatau completion bash)
_sigmatau() {
    local cur prev cword
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    cword=$COMP_CWORD

    # Value completion keyed on the preceding option.
    case "$prev" in
        --devs)    COMPREPLY=( $(compgen -W "@@DEVS@@" -- "${cur##*,}") ); return ;;
        --tau)     COMPREPLY=( $(compgen -W "@@TAU@@" -- "$cur") ); return ;;
        --format)  COMPREPLY=( $(compgen -W "@@FMT@@" -- "$cur") ); return ;;
        --ci)      COMPREPLY=( $(compgen -W "@@BOOL@@" -- "$cur") ); return ;;
        --detrend) COMPREPLY=( $(compgen -W "@@DETREND@@" -- "$cur") ); return ;;
        --out|--config|--input|--reference)
                   COMPREPLY=(); return ;;   # empty -> `-o default` does file/dir completion
        --tau0|--time-col|--value-col|--taus|--confidence|--header|--delim|--scaling)
                   COMPREPLY=(); return ;;   # free value; default filename completion is harmless
    esac

    if [[ $cword -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "@@SUBCMDS@@" -- "$cur") ); return
    fi

    if [[ "${COMP_WORDS[1]}" == "completion" && $cword -eq 2 ]]; then
        COMPREPLY=( $(compgen -W "@@SHELLS@@" -- "$cur") ); return
    fi

    if [[ "${COMP_WORDS[1]}" == "analyze" ]]; then
        if [[ "$cur" == -* ]]; then
            COMPREPLY=( $(compgen -W "@@OPTS@@" -- "$cur") )
        else
            COMPREPLY=()        # FILE positional -> native file/dir completion via `-o default`
        fi
        return
    fi

    COMPREPLY=()
}
complete -o bashdefault -o default -F _sigmatau sigmatau
"""
    return _fill(tpl)
end

function _completion_zsh()
    tpl = raw"""
#compdef sigmatau
# sigmatau zsh completion.
# Install:  sigmatau completion zsh > "${fpath[1]}/_sigmatau"   (then restart shell)
_sigmatau() {
  local curcontext="$curcontext" state line
  typeset -A opt_args
  _arguments -C '1: :->cmd' '*:: :->args'
  case $state in
    cmd)
      _values 'command' @@SUBCMDS@@
      ;;
    args)
      case $words[1] in
        completion)
          _values 'shell' @@SHELLS@@
          ;;
        analyze)
          _arguments -s \
            '--phase[treat FILE as phase data (default)]' \
            '--freq[treat FILE as fractional-frequency data]' \
            '--tau0[base sample interval]:tau0:' \
            '--time-col[1-based time column]:n:' \
            '--value-col[1-based value column]:n:' \
            '--devs[comma list of deviations]:devs:_values -s , deviation @@DEVS@@' \
            '--tau[tau-grid spacing]:grid:(@@TAU@@)' \
            '--taus[explicit averaging factors]:taus:' \
            '--ci[confidence intervals]:bool:(@@BOOL@@)' \
            '--confidence[CI confidence level]:level:' \
            '--format[output format]:fmt:(@@FMT@@)' \
            '--out[write output to PATH]:path:_files' \
            '--header[header lines to skip]:n:' \
            '--delim[field delimiter]:char:' \
            '--scaling[scale factor]:factor:' \
            '--detrend[detrend mode]:mode:(@@DETREND@@)' \
            '--quiet[suppress informational messages]' \
            '--debug[show stack traces on error]' \
            '*:file:_files'
          ;;
      esac
      ;;
  esac
}
_sigmatau "$@"
"""
    return _fill(tpl)
end

function _completion_fish()
    tpl = raw"""
# sigmatau fish completion.
# Install:  sigmatau completion fish > ~/.config/fish/completions/sigmatau.fish

# Subcommands (-f: the first word is a command, not a file).
complete -c sigmatau -f -n '__fish_use_subcommand' -a 'analyze'    -d 'compute stability deviations'
complete -c sigmatau -f -n '__fish_use_subcommand' -a 'completion' -d 'print a shell completion script'
complete -c sigmatau -f -n '__fish_use_subcommand' -a 'version'    -d 'print version'
complete -c sigmatau -f -n '__fish_use_subcommand' -a 'help'       -d 'print help'

complete -c sigmatau -f -n '__fish_seen_subcommand_from completion' -a '@@SHELLS@@'

# analyze flags (no value)
complete -c sigmatau -n '__fish_seen_subcommand_from analyze' -l phase -d 'phase data (default)'
complete -c sigmatau -n '__fish_seen_subcommand_from analyze' -l freq  -d 'frequency data'
complete -c sigmatau -n '__fish_seen_subcommand_from analyze' -l quiet -d 'suppress informational messages'
complete -c sigmatau -n '__fish_seen_subcommand_from analyze' -l debug -d 'show stack traces on error'

# analyze enumerated-value options (-x: require a value, do not offer files)
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l devs    -a '@@DEVS@@'    -d 'comma list of deviations'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l tau     -a '@@TAU@@'     -d 'tau-grid spacing'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l format  -a '@@FMT@@'     -d 'output format'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l ci      -a '@@BOOL@@'    -d 'confidence intervals'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l detrend -a '@@DETREND@@' -d 'detrend mode'

# analyze free-value options (-x: require a value, do not offer files)
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l tau0       -d 'base sample interval'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l time-col   -d '1-based time column'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l value-col  -d '1-based value column'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l taus       -d 'explicit averaging factors'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l confidence -d 'CI confidence level'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l header     -d 'header lines to skip'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l delim      -d 'field delimiter'
complete -c sigmatau -x -n '__fish_seen_subcommand_from analyze' -l scaling    -d 'scale factor'

# analyze path option: require a value AND force file completion.
complete -c sigmatau -r -F -n '__fish_seen_subcommand_from analyze' -l out -d 'write output to PATH'

# Positional FILE: the analyze option rules above are not marked -f, so fish
# offers files and directories by default for any non-option argument.
"""
    return _fill(tpl)
end

# PowerShell arrays use comma-separated quoted literals: 'a','b','c'.
_ps_arr(items) = join(("'" * string(x) * "'" for x in items), ",")

function _completion_powershell()
    tpl = raw"""
# sigmatau PowerShell completion.
# Install (add to $PROFILE):  sigmatau completion powershell | Out-String | Invoke-Expression
Register-ArgumentCompleter -Native -CommandName sigmatau -ScriptBlock {
    param($wordToComplete, $commandAst, $cursorPosition)

    $elements = @($commandAst.CommandElements | ForEach-Object { $_.ToString() })
    $sub  = if ($elements.Count -ge 2) { $elements[1] } else { '' }
    $prev = if ($elements.Count -ge 1) { $elements[-1] } else { '' }
    if ($wordToComplete) { $prev = if ($elements.Count -ge 2) { $elements[-2] } else { '' } }

    $subcommands = @@SUBCMDS_C@@
    $shells      = @@SHELLS_C@@
    $opts        = @@OPTS_C@@
    $devs        = @@DEVS_C@@

    function script:_res($items) {
        $items | Where-Object { $_ -like "$wordToComplete*" } | ForEach-Object {
            [System.Management.Automation.CompletionResult]::new($_, $_, 'ParameterValue', $_)
        }
    }

    switch ($prev) {
        '--devs'    { return _res $devs }
        '--tau'     { return _res @(@@TAU_C@@) }
        '--format'  { return _res @(@@FMT_C@@) }
        '--ci'      { return _res @(@@BOOL_C@@) }
        '--detrend' { return _res @(@@DETREND_C@@) }
        '--out'     { return }   # return nothing -> native path completion
    }

    if ($elements.Count -le 2 -and -not $wordToComplete.StartsWith('-')) {
        if ($elements.Count -le 1 -or $elements[1] -eq $wordToComplete) { return _res $subcommands }
    }
    if ($sub -eq 'completion') { return _res $shells }
    if ($sub -eq 'analyze') {
        if ($wordToComplete.StartsWith('-')) { return _res $opts }
        return   # FILE positional -> native path completion
    }
    return
}
"""
    return replace(tpl,
        "@@SUBCMDS_C@@" => _ps_arr(c.name for c in COMMANDS),
        "@@SHELLS_C@@"  => _ps_arr(SHELLS),
        "@@OPTS_C@@"    => _ps_arr(o.long for o in ANALYZE_OPTS),
        "@@DEVS_C@@"    => _ps_arr(COMPLETION_DEVS_STR),
        "@@TAU_C@@"     => _ps_arr(TAU_COMPLETE),
        "@@FMT_C@@"     => _ps_arr(FORMAT_CHOICES),
        "@@BOOL_C@@"    => _ps_arr(BOOL_CHOICES),
        "@@DETREND_C@@" => _ps_arr(DETREND_CHOICES),
    )
end

function _run_completion(args, out::IO, err::IO)
    if isempty(args)
        println(err, "sigmatau: completion requires a shell (valid: ", join(SHELLS, ", "), ")")
        return Cint(1)
    end
    if args[1] == "-h" || args[1] == "--help"
        println(out, "Usage: sigmatau completion <shell>")
        println(out, "Shells: ", join(SHELLS, ", "))
        return Cint(0)
    end
    sh = lowercase(args[1])
    script =
        sh == "bash" ? _completion_bash() :
        sh == "zsh"  ? _completion_zsh() :
        sh == "fish" ? _completion_fish() :
        (sh == "powershell" || sh == "pwsh") ? _completion_powershell() : nothing
    if script === nothing
        println(err, "sigmatau: unknown shell '", args[1], "' (valid: ", join(SHELLS, ", "), ")")
        return Cint(1)
    end
    print(out, script)
    endswith(script, "\n") || println(out)
    return Cint(0)
end

# ──────────────────────────────────────────────────────────────────────────
# ── `analyze` command ──────────────────────────────────────────────────────

function _run_analyze(args, out::IO, err::IO)
    debug = ("--debug" in args)             # honored even if parsing fails
    settings = try
        _parse_analyze(args)
    catch e
        return _fail(e, err, debug, catch_backtrace(), 1, "")
    end
    get(settings, :help, false) && (_print_analyze_help(out); return Cint(0))
    debug = settings.debug

    data = try
        settings.kind === :phase ?
            read_phase(settings.file; settings.readkw...) :
            read_frequency(settings.file; settings.readkw...)
    catch e
        return _fail(e, err, debug, catch_backtrace(), 2, "input error: ")
    end

    suite = try
        stability(data; devs = settings.devs, taus = settings.taus,
                  ci = settings.ci, confidence = settings.confidence)
    catch e
        return _fail(e, err, debug, catch_backtrace(), 3, "analysis error: ")
    end

    try
        _write_output(settings, suite, out)
    catch e
        _is_broken_pipe(e) && return Cint(0)   # consumer closed the pipe (e.g. | head)
        return _fail(e, err, debug, catch_backtrace(), 4, "output error: ")
    end

    if settings.out !== nothing && !settings.quiet
        println(err, "sigmatau: wrote ", length(suite.results), " deviation(s) to ", settings.out)
    end
    return Cint(0)
end

# ──────────────────────────────────────────────────────────────────────────
# ── Top-level entry point ──────────────────────────────────────────────────

"""
    SigmaTau.cli_main(args=ARGS; out=stdout, err=stderr) -> Cint

Non-interactive command-line entry point for SigmaTau. Returns a process exit
code (`0` ok, `1` usage error, `2` input/read error, `3` analysis error,
`4` output/write error) and never calls `exit`, so it is safe to call from tests
and from a PackageCompiler app. `out`/`err` are injectable for testing.

Commands: `analyze FILE [options]`, `completion SHELL`, `version`, `help`.
See `cli_main(["help"])` for the full option list. The thin launcher
`bin/sigmatau.jl` wraps this as `exit(SigmaTau.cli_main(ARGS))`.
"""
function cli_main(args::AbstractVector{<:AbstractString}=ARGS;
                  out::IO=stdout, err::IO=stderr)::Cint
    isempty(args) && (_print_help(out); return Cint(0))
    cmd  = args[1]
    rest = view(args, 2:length(args))
    if cmd in ("-h", "--help", "help")
        _print_help(out); return Cint(0)
    elseif cmd in ("-V", "--version", "version")
        _print_version(out); return Cint(0)
    elseif cmd == "completion"
        return _run_completion(rest, out, err)
    elseif cmd == "analyze"
        return _run_analyze(rest, out, err)
    elseif startswith(cmd, "-")
        println(err, "sigmatau: unknown option: ", cmd)
        println(err, "Run 'sigmatau --help' for usage.")
        return Cint(1)
    else
        println(err, "sigmatau: unknown command: ", cmd)
        println(err, "Run 'sigmatau --help' for usage.")
        return Cint(1)
    end
end

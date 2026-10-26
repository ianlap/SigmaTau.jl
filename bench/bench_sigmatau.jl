# Time SigmaTau's deviation kernels on a deterministic synthetic record, or on
# a phase record read from a file, and write the samples to disk so
# bench_allantools.py can time the same samples.
#
#   julia --project=. -t auto bench/bench_sigmatau.jl [N] [m_max]
#   julia --project=. -t auto bench/bench_sigmatau.jl --input <path> [--tau0 <s>] [--tag <name>] [m_max]
#
# Synthetic (default): N = 25_000 white-PM phase samples (seed 20260601) at
# τ₀ = 1 s, octave grid up to m_max = 512. Each kernel is called once to
# compile, then timed as the minimum of `reps` one-shot calls over the full τ
# grid (ci=false and raw totals, so the comparison is kernel against kernel;
# allantools returns no CI either).
#
# --input: the last whitespace-separated column of each non-empty, non-comment
# line is taken as phase (seconds or any unit; timing does not depend on it),
# sampled at --tau0 (default 1 s). Each kernel is compiled, then timed with a
# single call. Results are keyed by --tag (default "measured"), not the path.
#
# Both modes also time the default call (ci=true, bias-corrected totals) once
# per kernel, after a compiling call, as the end-to-end figure.
#
# Writes bench/results/input_<key>.txt and bench/results/sigmatau_<key>.json,
# where <key> is N for the synthetic record and the tag for --input.

using SigmaTau
using Random
using Printf

function parse_args(args)
    opts = Dict{String,String}()
    pos = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--input", "--tau0", "--tag")
            i < length(args) || error("$a needs a value")
            opts[a[3:end]] = args[i+1]
            i += 2
        elseif startswith(a, "--")
            error("unknown option $a")
        else
            push!(pos, a)
            i += 1
        end
    end
    return opts, pos
end

function read_last_column(path)
    x = Float64[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, '#') || startswith(s, '%')) && continue
        push!(x, parse(Float64, last(split(s))))
    end
    return x
end

const OPTS, POS = parse_args(ARGS)
const FROM_FILE = haskey(OPTS, "input")
const OUT = joinpath(@__DIR__, "results")
mkpath(OUT)

if FROM_FILE
    x     = read_last_column(OPTS["input"])
    tau0  = parse(Float64, get(OPTS, "tau0", "1.0"))
    key   = get(OPTS, "tag", "measured")
    M_MAX = length(POS) >= 1 ? parse(Int, POS[1]) : 512
    source = "file"
else
    N0    = length(POS) >= 1 ? parse(Int, POS[1]) : 25_000
    M_MAX = length(POS) >= 2 ? parse(Int, POS[2]) : 512
    tau0  = parse(Float64, get(OPTS, "tau0", "1.0"))
    x     = randn(Xoshiro(20260601), N0) .* 1e-9     # white PM, seconds
    key   = string(N0)
    source = "synthetic white PM, seed 20260601"
end
const N = length(x)

m_values = [2^k for k in 0:floor(Int, log2(M_MAX))]
open(joinpath(OUT, "input_$key.txt"), "w") do io
    for v in x
        @printf(io, "%.17e\n", v)
    end
end
pd = PhaseData(x, tau0)

kernels = (
    (:adev,    adev),    (:mdev,    mdev),    (:hdev,    hdev),    (:tdev, tdev),
    (:totdev,  totdev),  (:mtotdev, mtotdev), (:htotdev, htotdev),
)
# Totals are bias-corrected by default; the raw kernel is what allantools computes.
call(f, pd, m) = f in (totdev, mtotdev, htotdev) ? f(pd, m; ci=false, correct_bias=false) :
                                                   f(pd, m; ci=false)

function best_of(g, reps)
    best = Inf
    for _ in 1:reps
        t = @elapsed g()
        best = min(best, t)
    end
    return best
end

cpu = Sys.cpu_info()[1].model
println("SigmaTau $(pkgversion(SigmaTau)), Julia $VERSION, $(Threads.nthreads()) threads, ",
        "$cpu ($(Sys.CPU_THREADS) logical CPUs)")
println("$source, N = $N, τ₀ = $tau0 s, m ≤ $M_MAX")
results  = Dict{String,Float64}()
defaults = Dict{String,Float64}()
reps_used = Dict{String,Int}()
for (name, f) in kernels
    call(f, pd, m_values)                                   # compile
    reps = FROM_FILE ? 1 : (name in (:mtotdev, :htotdev) ? 3 : 20)
    t = best_of(() -> call(f, pd, m_values), reps)
    results[string(name)] = t
    reps_used[string(name)] = reps
    f(pd, m_values)                                         # compile the default path
    d = @elapsed f(pd, m_values)                            # default call: ci=true, bias-corrected totals
    defaults[string(name)] = d
    @printf("  %-8s kernel %10.3f ms (min of %2d)   default call %10.3f ms\n",
            name, 1e3 * t, reps, 1e3 * d)
end

function write_seconds(io, label, d, last)
    println(io, "  \"$label\": {")
    for (i, (k, _)) in enumerate(kernels)
        @printf(io, "    \"%s\": %.6e%s\n", k, d[string(k)], i < length(kernels) ? "," : "")
    end
    println(io, last ? "  }" : "  },")
end

open(joinpath(OUT, "sigmatau_$key.json"), "w") do io
    println(io, "{")
    println(io, "  \"library\": \"SigmaTau.jl $(pkgversion(SigmaTau))\",")
    println(io, "  \"julia\": \"$VERSION\", \"threads\": $(Threads.nthreads()),")
    println(io, "  \"cpu\": \"$cpu\", \"logical_cpus\": $(Sys.CPU_THREADS),")
    println(io, "  \"source\": \"$source\",")
    println(io, "  \"N\": $N, \"tau0\": $tau0, \"m_values\": $(m_values),")
    println(io, "  \"reps\": {", join(["\"$k\": $(reps_used[string(k)])" for (k, _) in kernels], ", "), "},")
    write_seconds(io, "seconds", results, false)
    write_seconds(io, "seconds_default_call", defaults, true)
    println(io, "}")
end
println("wrote ", joinpath(OUT, "sigmatau_$key.json"))

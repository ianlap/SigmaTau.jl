# Time SigmaTau's deviation kernels on a deterministic synthetic record and
# write the record to disk so bench_allantools.py can time the same samples.
#
#   julia --project=. -t auto bench/bench_sigmatau.jl [N] [m_max]
#
# Defaults: N = 25_000 white-PM phase samples at τ₀ = 1 s, octave grid up to
# m_max = 512. Each kernel is called once to compile, then timed as the
# minimum of `reps` one-shot calls over the full τ grid (ci=false, so the
# comparison is kernel against kernel; allantools returns no CI either).
# Writes bench/results/input_<N>.txt and bench/results/sigmatau_<N>.json.

using SigmaTau
using Random
using Printf

const N     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 25_000
const M_MAX = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 512
const OUT   = joinpath(@__DIR__, "results")
mkpath(OUT)

m_values = [2^k for k in 0:floor(Int, log2(M_MAX))]
x = randn(Xoshiro(20260601), N) .* 1e-9          # white PM, seconds
open(joinpath(OUT, "input_$N.txt"), "w") do io
    for v in x
        @printf(io, "%.17e\n", v)
    end
end
pd = PhaseData(x, 1.0)

kernels = (
    (:adev,    adev),    (:mdev,    mdev),    (:hdev,    hdev),    (:tdev, tdev),
    (:totdev,  totdev),  (:mtotdev, mtotdev), (:htotdev, htotdev),
)
# Totals are bias-corrected by default; the raw kernel is what allantools computes.
call(f, pd, m) = f in (totdev, mtotdev, htotdev) ? f(pd, m; ci=false, correct_bias=false) :
                                                   f(pd, m; ci=false)

function best_of(f, pd, m, reps)
    best = Inf
    for _ in 1:reps
        t = @elapsed call(f, pd, m)
        best = min(best, t)
    end
    return best
end

println("SigmaTau $(pkgversion(SigmaTau)), Julia $VERSION, $(Threads.nthreads()) threads, N = $N, m ≤ $M_MAX")
results = Dict{String,Float64}()
for (name, f) in kernels
    call(f, pd, m_values)                                   # compile
    reps = name in (:mtotdev, :htotdev) ? 3 : 20
    t = best_of(f, pd, m_values, reps)
    results[string(name)] = t
    @printf("  %-8s %10.3f ms\n", name, 1e3 * t)
end

open(joinpath(OUT, "sigmatau_$N.json"), "w") do io
    println(io, "{")
    println(io, "  \"library\": \"SigmaTau.jl $(pkgversion(SigmaTau))\",")
    println(io, "  \"julia\": \"$VERSION\", \"threads\": $(Threads.nthreads()),")
    println(io, "  \"N\": $N, \"m_values\": $(m_values),")
    println(io, "  \"seconds\": {")
    ks = collect(keys(results))
    for (i, k) in enumerate(ks)
        @printf(io, "    \"%s\": %.6e%s\n", k, results[k], i < length(ks) ? "," : "")
    end
    println(io, "  }\n}")
end
println("wrote ", joinpath(OUT, "sigmatau_$N.json"))

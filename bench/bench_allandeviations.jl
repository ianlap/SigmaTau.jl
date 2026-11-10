# Time AllanDeviations.jl on the record written by bench_sigmatau.jl.
#
#   julia --project=tools/validation/allandeviations bench/bench_allandeviations.jl [N]
#   julia --project=tools/validation/allandeviations bench/bench_allandeviations.jl --tag <name>
#
# Reads bench/results/input_<key>.txt and the τ grid from sigmatau_<key>.json
# (<key> is N, default 25000, or the tag), compiles each kernel with one call,
# then times it as the minimum of 20 calls (one call for a tagged record) and
# writes bench/results/allandeviations_<key>.json. AllanDeviations.jl has no
# modified total or Hadamard total. Single-threaded, as the package is.
using AllanDeviations, Printf

const OUT = joinpath(@__DIR__, "results")
const TAGGED = length(ARGS) >= 2 && ARGS[1] == "--tag"
const KEY = TAGGED ? ARGS[2] : (isempty(ARGS) ? "25000" : ARGS[1])

sig = read(joinpath(OUT, "sigmatau_$KEY.json"), String)
grab(k) = match(Regex("\"$k\": ([^,\\n]+)"), sig).captures[1]
tau0 = parse(Float64, grab("tau0"))
ms = parse.(Int, split(strip(match(r"\"m_values\": \[([^\]]*)\]", sig).captures[1]), ','))
x = parse.(Float64, readlines(joinpath(OUT, "input_$KEY.txt")))
taus = ms .* tau0
rate = 1 / tau0

# SigmaTau's adev and hdev are the overlapping estimators, AllanDeviations' default.
kernels = (("adev", allandev), ("mdev", mallandev), ("hdev", hadamarddev),
           ("tdev", timedev), ("totdev", totaldev))
reps = TAGGED ? 1 : 20
println("AllanDeviations ", pkgversion(AllanDeviations), ", Julia ", VERSION, ", N = ", length(x))
secs = Dict{String,Float64}()
for (name, f) in kernels
    f(x, rate; taus)                                    # compile
    secs[name] = minimum(@elapsed(f(x, rate; taus)) for _ in 1:reps)
    @printf("  %-8s %10.3f ms (min of %d)\n", name, 1e3 * secs[name], reps)
end

open(joinpath(OUT, "allandeviations_$KEY.json"), "w") do io
    println(io, "{")
    println(io, "  \"library\": \"AllanDeviations.jl $(pkgversion(AllanDeviations))\", \"julia\": \"$VERSION\",")
    println(io, "  \"N\": $(length(x)), \"tau0\": $tau0, \"m_values\": $ms, \"reps\": $reps,")
    println(io, "  \"seconds\": {")
    for (i, (name, _)) in enumerate(kernels)
        @printf(io, "    \"%s\": %.6e%s\n", name, secs[name], i < length(kernels) ? "," : "")
    end
    println(io, "  }")
    println(io, "}")
end
println("wrote ", joinpath(OUT, "allandeviations_$KEY.json"))

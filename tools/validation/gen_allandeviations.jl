# Regenerate the AllanDeviations.jl reference table for the cross-check suite.
#
#     julia --project=tools/validation/allandeviations -e 'using Pkg; Pkg.instantiate()'
#     julia --project=tools/validation/allandeviations tools/validation/gen_allandeviations.jl [DIR]
#
# Runs every estimator AllanDeviations.jl has on the committed power-law
# records (noise/alpha_*.txt, N = 4096) at octave averaging factors and writes
# DIR/devs.csv (default test/fixtures/validation/allandeviations/) in the
# layout of the allantools table, with allantools' estimator names, so the
# same SigmaTau calls check both.
using AllanDeviations, Printf

const FIX = normpath(joinpath(@__DIR__, "..", "..", "test", "fixtures", "validation"))
const OUT = isempty(ARGS) ? joinpath(FIX, "allandeviations") : ARGS[1]

const ESTIMATORS = (
    ("adev",  (x, t) -> allandev(x, 1.0; taus=t, overlapping=false)),
    ("oadev", (x, t) -> allandev(x, 1.0; taus=t)),
    ("mdev",  (x, t) -> mallandev(x, 1.0; taus=t)),
    ("tdev",  (x, t) -> timedev(x, 1.0; taus=t)),
    ("hdev",  (x, t) -> hadamarddev(x, 1.0; taus=t, overlapping=false)),
    ("ohdev", (x, t) -> hadamarddev(x, 1.0; taus=t)),
    ("totdev", (x, t) -> totaldev(x, 1.0; taus=t)),
    ("mtie",  (x, t) -> mtie(x, 1.0; taus=t)),
)

readrec(path) = [parse(Float64, l) for l in eachline(path) if !isempty(strip(l)) && !startswith(l, '#')]

mkpath(OUT)
open(joinpath(OUT, "devs.csv"), "w") do io
    println(io, "dataset,estimator,m,tau,n,dev")
    for a in (2, 1, 0, -1, -2, -3, -4)
        name = a > 0 ? "alpha_p$a" : a == 0 ? "alpha_0" : "alpha_m$(-a)"
        x = readrec(joinpath(FIX, "noise", name * ".txt"))
        taus = [2.0^k for k in 0:floor(Int, log2(length(x) ÷ 4))]
        for (est, f) in ESTIMATORS
            r = f(x, taus)
            for (t, n, d) in zip(r.tau, r.count, r.deviation)
                @printf(io, "alpha_%+d,%s,%d,%.17e,%d,%.17e\n", a, est, round(Int, t), t, n, d)
            end
        end
    end
end
open(joinpath(OUT, "VERSION"), "w") do io
    println(io, "AllanDeviations ", pkgversion(AllanDeviations))
end

# make_noise_fixtures.jl — writes the power-law phase records used by the
# cross-check suite (test/validation/). Run once; the records are committed,
# so allantools and SigmaTau read the same numbers.
#
#   julia --project=. tools/validation/make_noise_fixtures.jl
#
# One record per integer α from +2 (white PM) to −4 (random-run FM),
# N = 4096 samples at τ0 = 1 s, Kasdin–Walter synthesis with fixed seeds.
using SigmaTau, FFTW, Random, Printf

const N = 4096
const OUT = joinpath(@__DIR__, "..", "..", "test", "fixtures", "validation", "noise")
mkpath(OUT)

name(α) = α > 0 ? "alpha_p$(α).txt" : (α == 0 ? "alpha_0.txt" : "alpha_m$(-α).txt")

for α in 2:-1:-4
    rng = Xoshiro(20261108 + 17 * (α + 4))
    x = SigmaTau._gen_powerlaw_phase(α, N; rng=rng)
    x ./= maximum(abs, x)                # unit-scale so %.12e keeps full meaning
    open(joinpath(OUT, name(α)), "w") do io
        println(io, "# SigmaTau cross-check fixture: phase, tau0 = 1 s, N = $N, alpha = $α")
        println(io, "# Kasdin-Walter synthesis, Xoshiro seed $(20261108 + 17 * (α + 4))")
        for v in x
            @printf(io, "%.12e\n", v)
        end
    end
end
println("wrote ", readdir(OUT))

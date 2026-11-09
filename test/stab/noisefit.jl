# noisefit.jl — the interactive fit recovers noise_gen's levels one component
# at a time, and undo takes them back.

using Test
using SigmaTau
using Random

@testset "interactive noise fit" begin
    h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30)
    p = noise_gen(PhaseData, 65536, 1.0; h=h, rng=MersenneTwister(1))
    r = mhdev(p)                                     # τ = 1, 2, 4, …, 8192 s
    f = noise_fit(r)
    add_component!(f, :WHPM, 1:2)
    add_component!(f, 0, 6:9)
    add_component!(f, :RWFM, 12:14)
    @test [c.alpha for c in f.components] == [2, 0, -2]
    @test f.components[1].h ≈ 1e-22 rtol=0.25
    @test f.components[2].h ≈ 1e-24 rtol=0.1
    @test f.components[3].h ≈ 1e-30 rtol=0.5
    @test powerlaw_sigma(f, r.tau; deviation=:mhdev) ≈ r.dev rtol=0.1
    @test noise_residual(f) ≈ r.dev .^ 2 .- powerlaw_sigma(f, r.tau; deviation=:mhdev) .^ 2
    @test clock_model(f).q1 == f.components[2].h / 2

    # A component is fitted to what the others leave, and undo restores it.
    before = noise_residual(f)
    undo!(f)
    @test length(f.components) == 2 && noise_residual(f)[14] > before[14]
    @test noise_fit(mhtotdev(p)).response == :mhdev
    a = add_component!(noise_fit(adev(p; ci=false)), :WHPM, 1:2)   # no intervals: equal weights
    @test add_component!(a, 0, 4:7).components[2].h ≈ 1e-24 rtol=0.1

    @test_throws ArgumentError add_component!(f, 0, 6:9)          # already in
    @test_throws ArgumentError add_component!(noise_fit(adev(p)), -4, 1:3)
    @test_throws ArgumentError add_component!(f, -2, 1:2)         # residual ≤ 0 there
    @test_throws ArgumentError noise_fit(tdev(p))
    @test_throws ArgumentError undo!(noise_fit(r))
end

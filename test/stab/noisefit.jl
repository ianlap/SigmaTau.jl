# noisefit.jl — the interactive fit recovers noise_gen's levels whatever the
# order the components are added in, and undo takes them back.

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
    # White PM is fitted where white FM adds a quarter to σ²: it reads 25 %
    # high unless it is refitted once white FM is in the model.
    @test f.components[1].h ≈ 1e-22 rtol=0.03
    @test f.components[2].h ≈ 1e-24 rtol=0.03
    @test f.components[3].h ≈ 1e-30 rtol=0.5
    g = add_component!(add_component!(add_component!(noise_fit(r), :RWFM, 12:14), 0, 6:9), :WHPM, 1:2)
    @test [c.h for c in g.components] ≈ [c.h for c in reverse(f.components)]
    @test powerlaw_sigma(f, r.tau; deviation=:mhdev) ≈ r.dev rtol=0.1
    @test noise_residual(f) ≈ r.dev .^ 2 .- powerlaw_sigma(f, r.tau; deviation=:mhdev) .^ 2
    @test clock_model(f).q1 == f.components[2].h / 2

    # Undo drops the last component and refits the rest.
    before = noise_residual(f)
    undo!(f)
    @test length(f.components) == 2 && noise_residual(f)[14] > before[14]
    @test [c.h for c in f.components] ≈ [c.h for c in add_component!(add_component!(noise_fit(r), :WHPM, 1:2), 0, 6:9).components]
    @test noise_fit(mhtotdev(p)).response == :mhdev
    a = add_component!(noise_fit(adev(p; ci=false)), :WHPM, 1:2)   # no intervals: equal weights
    @test add_component!(a, 0, 4:7).components[2].h ≈ 1e-24 rtol=0.1

    @test_throws ArgumentError add_component!(f, 0, 6:9)          # already in
    @test_throws ArgumentError add_component!(noise_fit(adev(p)), -4, 1:3)
    over = add_component!(noise_fit(r), 0, 12:14)                 # white FM fitted on random walk
    @test_throws ArgumentError add_component!(over, 2, 6:9)       # residual ≤ 0 there
    @test_throws ArgumentError noise_fit(tdev(p))
    @test_throws ArgumentError undo!(noise_fit(r))
end

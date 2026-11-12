# noisefit.jl — the interactive fit recovers noise_gen's levels whatever the
# order the components are added in, and undo takes them back.

using Test
using SigmaTau
using Random

@testset "interactive noise fit" begin
    h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30)
    # On the model's own σ(τ) the levels come back exactly, whatever the order.
    # Fitted once in sequence, white PM reads 25 % high: white FM adds a
    # quarter to σ² where it is fitted.
    τ = 2.0 .^ (0:13)
    exact = StabilityResult(:mhdev, τ, powerlaw_sigma(h, τ; deviation=:mhdev, tau0=1.0), Symbol[],
                            SigmaTau._CIBound[], Float64[], zeros(Int, 14), 1.0, 65536, NaN)
    e = add_component!(add_component!(add_component!(noise_fit(exact), :WHPM, 1:2), 0, 6:9), :RWFM, 12:14)
    @test [c.h for c in e.components] ≈ [1e-22, 1e-24, 1e-30] rtol=1e-6
    g = add_component!(add_component!(add_component!(noise_fit(exact), :RWFM, 12:14), 0, 6:9), :WHPM, 1:2)
    @test [c.h for c in g.components] ≈ [1e-30, 1e-24, 1e-22] rtol=1e-6

    p = noise_gen(PhaseData, 65536, 1.0; h=h, rng=MersenneTwister(1))
    r = mhdev(p)                                     # τ = 1, 2, 4, …, 8192 s
    f = noise_fit(r)
    add_component!(f, :WHPM, 1:2)
    add_component!(f, 0, 6:9)
    add_component!(f, :RWFM, 12:14)
    @test [c.alpha for c in f.components] == [2, 0, -2]
    @test powerlaw_sigma(f, r.tau; deviation=:mhdev) ≈ r.dev rtol=0.1
    @test noise_residual(f) ≈ r.dev .^ 2 .- powerlaw_sigma(f, r.tau; deviation=:mhdev) .^ 2
    @test clock_model(f).q1 == f.components[2].h / 2
    @test_throws "not 2.0 s" clock_model(f; tau0=2.0)
    # The model covers the τ the fit saw: flicker states end at its longest τ.
    fl = NoiseFit(r, :mhdev, [(alpha=0, h=1e-24, range=1:3), (alpha=-1, h=1e-26, range=10:14)])
    hl = Dict(0 => 1e-24, -1 => 1e-26)
    @test clock_model(fl).fm_band == clock_model(hl; tau0=1.0, tau_max=8192.0).fm_band != clock_model(hl; tau0=1.0).fm_band

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
    @test_throws "unknown noise type :PINK; use an integer α or one of WHPM," add_component!(noise_fit(adev(p)), :PINK, 1:3)
    over = add_component!(noise_fit(r), 0, 12:14)                 # white FM fitted on random walk
    @test_throws ArgumentError add_component!(over, 2, 6:9)       # residual ≤ 0 there
    @test_throws ArgumentError noise_fit(tdev(p))
    @test_throws ArgumentError undo!(noise_fit(r))
end

# test/stab/metadata.jl — the record metadata a StabilityResult carries.
using Test, Random, SigmaTau

@testset "StabilityResult carries tau0, n and confidence" begin
    x = cumsum(randn(MersenneTwister(5), 1024)) .* 1e-9
    pd = PhaseData(x, 2.0)
    r = adev(pd; confidence = 0.95)
    @test (r.tau0, r.n, r.confidence) == (2.0, 1024, 0.95)
    @test isnan(adev(pd; ci = false).confidence)
    @test isnan(mtie(pd).confidence)            # no interval model
    # 1023 frequency samples integrate to 1024 phase points (x[1] = 0).
    @test adev(FrequencyData(diff(x), 2.0)).n == 1024
    @test stability(pd)[:mdev].tau0 == 2.0
    acc = StreamingStability(:adev, 2.0, [1, 2])
    append!(acc, x[1:100])
    @test (snapshot(acc).tau0, snapshot(acc).n) == (2.0, 100)
    legacy = StabilityResult(:adev, [1.0], [1.0], Symbol[], SigmaTau._CIBound[], Float64[], [1])
    @test isnan(legacy.tau0) && legacy.n == 0 && isnan(legacy.confidence)
    @test psd_y(pd; taper = :hamming).taper === :hamming
end

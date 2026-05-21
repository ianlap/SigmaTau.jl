using Test
using SigmaTau

@testset "Shared types" begin
    @testset "PhaseData" begin
        p = PhaseData([1.0, 2.0, 3.0], 1.0)
        @test p.x == [1.0, 2.0, 3.0]
        @test p.tau0 == 1.0
        @test p isa AbstractTimingData
    end

    @testset "FrequencyData" begin
        f = FrequencyData([0.1, 0.2], 0.5)
        @test f.y == [0.1, 0.2]
        @test f.tau0 == 0.5
        @test f isa AbstractTimingData
    end

    @testset "StabilityResult fields" begin
        r = StabilityResult(:adev, [1.0], [0.5], Symbol[], Float64[], Float64[], Float64[])
        @test r.deviation_type === :adev
        @test r.tau == [1.0]
        @test r.dev == [0.5]
        @test isempty(r.noise_type)
        @test isempty(r.edf)
    end
end

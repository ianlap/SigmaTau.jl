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

    @testset "constructor validation" begin
        # tau0 must be positive
        @test_throws ArgumentError PhaseData([1.0, 2.0], 0.0)
        @test_throws ArgumentError PhaseData([1.0, 2.0], -1.0)
        @test_throws ArgumentError FrequencyData([0.1, 0.2], 0.0)
        @test_throws ArgumentError FrequencyData([0.1, 0.2], -2.0)
        # need at least 2 samples
        @test_throws ArgumentError PhaseData(Float64[], 1.0)
        @test_throws ArgumentError PhaseData([1.0], 1.0)
        @test_throws ArgumentError FrequencyData(Float64[], 1.0)
        @test_throws ArgumentError FrequencyData([0.1], 1.0)
        # integer tau0 is accepted and stored as Float64
        p = PhaseData([1.0, 2.0], 1)
        @test p.tau0 === 1.0
        f = FrequencyData([0.1, 0.2], 2)
        @test f.tau0 === 2.0
        # Float32 samples still construct (parameterization preserved)
        p32 = PhaseData(Float32[1, 2, 3], 1.0)
        @test p32 isa PhaseData{Float32}
        @test p32.tau0 === 1.0
    end
end

using Test
using SigmaTau

@testset "Shared types" begin
    @testset "PhaseData" begin
        p = PhaseData([1.0, 2.0, 3.0], 1.0)
        @test p.x == [1.0, 2.0, 3.0]
        @test p.tau0 == 1.0
        @test p isa AbstractTimingData
        # tau0 defaults to 1.0
        @test PhaseData([1.0, 2.0, 3.0]).tau0 == 1.0
    end

    @testset "FrequencyData" begin
        f = FrequencyData([0.1, 0.2], 0.5)
        @test f.y == [0.1, 0.2]
        @test f.tau0 == 0.5
        @test f isa AbstractTimingData
        # tau0 defaults to 1.0
        @test FrequencyData([0.1, 0.2]).tau0 == 1.0
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

    @testset "show summaries" begin
        # Compact one-line summaries, not full-array dumps.
        r = StabilityResult(:adev, [1.0, 2.0, 4.0], [3e-10, 2e-10, 1e-10],
                            Symbol[], Float64[], Float64[], Float64[])
        s = sprint(show, r)
        @test occursin("StabilityResult", s)
        @test occursin("adev", s)
        @test occursin("3 pts", s)
        @test occursin("no CI", s)
        # CI-populated result reports "with CI".
        rci = StabilityResult(:adev, [1.0], [3e-10], [:wpm], [2e-10], [4e-10], [10.0])
        @test occursin("with CI", sprint(show, rci))
        # Zero-point result does not error.
        r0 = StabilityResult(:adev, Float64[], Float64[], Symbol[], Float64[], Float64[], Float64[])
        @test occursin("0 pts", sprint(show, r0))

        sp = SpectralResult(:Sy, [0.1, 0.2], [1.0, 0.5], :per_Hz, 256, 128, :hann)
        @test occursin("SpectralResult", sprint(show, sp))
        @test occursin("Sy", sprint(show, sp))
        @test occursin("2 bins", sprint(show, sp))
    end
end

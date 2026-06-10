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
        CIB = @NamedTuple{lo::Float64, hi::Float64}
        r = StabilityResult(:adev, [1.0], [0.5], Symbol[], CIB[], Float64[], [10])
        @test r.deviation_type === :adev
        @test r.tau == [1.0]
        @test r.dev == [0.5]
        @test isempty(r.noise_type)
        @test isempty(r.ci)
        @test isempty(r.edf)
        @test r.neff == [10]
        # Accessors are empty-in empty-out and typed Vector{Float64}.
        @test ci_lower(r) isa Vector{Float64} && isempty(ci_lower(r))
        @test ci_upper(r) isa Vector{Float64} && isempty(ci_upper(r))
        # Tuple field access and accessor agreement on a CI-carrying result.
        rci = StabilityResult(:adev, [1.0, 2.0], [0.5, 0.4], [:WHFM, :WHFM],
                              CIB[(lo=0.4, hi=0.6), (lo=0.3, hi=0.5)],
                              [10.0, 5.0], [10, 8])
        @test rci.ci[2].lo == 0.3
        @test rci.ci[2].hi == 0.5
        @test ci_lower(rci) == [0.4, 0.3]
        @test ci_upper(rci) == [0.6, 0.5]
        @test rci.dev .- ci_lower(rci) ≈ [0.1, 0.1]
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
        CIB = @NamedTuple{lo::Float64, hi::Float64}
        # Compact one-line summaries, not full-array dumps.
        r = StabilityResult(:adev, [1.0, 2.0, 4.0], [3e-10, 2e-10, 1e-10],
                            Symbol[], CIB[], Float64[], [98, 96, 92])
        s = sprint(show, r)
        @test occursin("StabilityResult", s)
        @test occursin("adev", s)
        @test occursin("3 pts", s)
        @test occursin("no CI", s)
        # CI-populated result reports "with CI".
        rci = StabilityResult(:adev, [1.0], [3e-10], [:wpm],
                              CIB[(lo=2e-10, hi=4e-10)], [10.0], [98])
        @test occursin("with CI", sprint(show, rci))
        # Zero-point result does not error.
        r0 = StabilityResult(:adev, Float64[], Float64[], Symbol[], CIB[], Float64[], Int[])
        @test occursin("0 pts", sprint(show, r0))

        sp = SpectralResult(:Sy, [0.1, 0.2], [1.0, 0.5], :per_Hz, 256, 128, :hann)
        @test occursin("SpectralResult", sprint(show, sp))
        @test occursin("Sy", sprint(show, sp))
        @test occursin("2 bins", sprint(show, sp))
    end

    @testset "text/plain table" begin
        CIB = @NamedTuple{lo::Float64, hi::Float64}
        # Full REPL display: aligned table with %.4g rounding (display only —
        # the stored vectors keep full precision).
        rci = StabilityResult(:adev, [1.0, 2.0], [3.14159265e-10, 2.71828182e-10],
                              [:WHFM, :WHFM],
                              CIB[(lo=2.5e-10, hi=4.2e-10), (lo=2.2e-10, hi=3.5e-10)],
                              [10.123456, 5.654321], [98, 96])
        s = sprint(show, MIME"text/plain"(), rci)
        for col in ("tau", "dev", "neff", "ci_lo", "ci_hi", "edf", "noise")
            @test occursin(col, s)
        end
        @test occursin("3.142e-10", s)          # rounded for display
        @test !occursin("3.14159265e-10", s)    # not full precision
        @test occursin("10.12", s)
        @test occursin("WHFM", s)
        @test occursin("98", s)                 # neff column carries the counts
        @test rci.dev[1] == 3.14159265e-10      # underlying value untouched

        # No CI ⇒ only the tau/dev/neff columns appear.
        rno = StabilityResult(:adev, [1.0], [0.5], Symbol[], CIB[], Float64[], [10])
        sno = sprint(show, MIME"text/plain"(), rno)
        @test !occursin("ci_lo", sno)
        @test occursin("dev", sno)
        @test occursin("neff", sno)

        # Long results truncate to 11 rows plus an ellipsis line.
        long = StabilityResult(:adev, collect(1.0:20.0), fill(1e-10, 20),
                               Symbol[], CIB[], Float64[], fill(50, 20))
        sl = sprint(show, MIME"text/plain"(), long)
        @test occursin("… 20 rows total", sl)
        @test count(==('\n'), sl) == 13         # header + 11 rows + ellipsis
    end
end

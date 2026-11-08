using Test
using SigmaTau

const CIB = @NamedTuple{lo::Float64, hi::Float64}

@testset "Shared types" begin
    @testset "data records" begin
        p = PhaseData([1.0, 2.0, 3.0], 2.0)
        f = FrequencyData([0.1, 0.2], 0.5)
        @test (p.x, p.tau0, p.source) == ([1.0, 2.0, 3.0], 2.0, "user")
        @test (f.y, f.tau0, f.source) == ([0.1, 0.2], 0.5, "user")
        @test PhaseData([1.0, 2.0]).tau0 === FrequencyData([0.1, 0.2], 1).tau0 === 1.0
        @test PhaseData([1.0, 2.0]; source="a.dat").source == "a.dat"
        @test PhaseData(Float32[1, 2, 3], 1.0) isa PhaseData{Float32}
        for T in (PhaseData, FrequencyData)
            @test_throws ArgumentError T([1.0, 2.0], 0.0)
            @test_throws ArgumentError T([1.0, 2.0], -1.0)
            @test_throws ArgumentError T([1.0], 1.0)
        end
        # One-line show; the source appears only when it is not "user".
        @test sprint(show, p) == "PhaseData{Float64}(3 samples, τ₀=2 s)"
        @test !occursin("source", sprint(show, f))
        @test occursin("source=\"freq.csv\"",
                       sprint(show, FrequencyData([0.1, 0.2]; source="freq.csv")))
    end

    @testset "StabilityResult" begin
        r = StabilityResult(:adev, [1.0, 2.0], [0.5, 0.4], [:WHFM, :WHFM],
                            CIB[(lo=0.4, hi=0.6), (lo=0.3, hi=0.5)], [10.0, 5.0], [10, 8])
        @test (ci_lower(r), ci_upper(r)) == ([0.4, 0.3], [0.6, 0.5])
        r0 = StabilityResult(:adev, [1.0], [0.5], Symbol[], CIB[], Float64[], [10])
        @test ci_lower(r0) isa Vector{Float64} && isempty(ci_upper(r0))

        @test sprint(show, r) == "StabilityResult(:adev, 2 pts, τ∈[1, 2] s, with CI)"
        @test occursin("no CI", sprint(show, r0))
        empty = StabilityResult(:adev, Float64[], Float64[], Symbol[], CIB[], Float64[], Int[])
        @test occursin("0 pts", sprint(show, empty))
        sp = SpectralResult(:Sy, [0.1, 0.2], [1.0, 0.5], :per_Hz, 256, 128, :hann)
        @test all(occursin(sprint(show, sp)), ("SpectralResult", "Sy", "2 bins"))
    end

    @testset "text/plain table" begin
        # Aligned table rounded to %.4g for display; stored values keep full precision.
        r = StabilityResult(:adev, [1.0, 2.0], [3.14159265e-10, 2.71828182e-10],
                            [:WHFM, :WHFM],
                            CIB[(lo=2.5e-10, hi=4.2e-10), (lo=2.2e-10, hi=3.5e-10)],
                            [10.123456, 5.654321], [98, 96])
        s = sprint(show, MIME"text/plain"(), r)
        @test all(occursin(s), ("tau", "dev", "neff", "ci_lo", "ci_hi", "edf", "noise",
                                "3.142e-10", "10.12", "WHFM", "98"))
        @test !occursin("3.14159265e-10", s)

        sno = sprint(show, MIME"text/plain"(),
                     StabilityResult(:adev, [1.0], [0.5], Symbol[], CIB[], Float64[], [10]))
        @test occursin("neff", sno) && !occursin("ci_lo", sno)

        long = StabilityResult(:adev, collect(1.0:20.0), fill(1e-10, 20),
                               Symbol[], CIB[], Float64[], fill(50, 20))
        sl = sprint(show, MIME"text/plain"(), long)
        @test occursin("… 20 rows total", sl)
        @test count(==('\n'), sl) == 13         # header + 11 rows + ellipsis
    end
end

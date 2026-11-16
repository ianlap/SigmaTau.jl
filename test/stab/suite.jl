# test/stab/suite.jl — stability() compute-all + StabilitySuite.
# Standalone-runnable and included from test/stab/runtests.jl.

using Test
using Random
using SigmaTau

@testset "stability() compute-all + StabilitySuite" begin
    Random.seed!(11)
    pd = PhaseData(cumsum(randn(1024)) .* 1e-9, 1.0)

    @testset "ordered results, indexing, keys, metadata" begin
        s = stability(pd; devs=(:adev, :mdev, :hdev), m=Octave, ci=false)
        @test s isa StabilitySuite
        @test length(s) == 3
        @test keys(s) == [:adev, :mdev, :hdev]        # request order preserved
        @test haskey(s, :mdev)
        @test !haskey(s, :pdev)
        @test s[:adev].deviation_type === :adev
        @test s[2].deviation_type === :mdev           # positional indexing
        @test_throws KeyError s[:pdev]
        @test [r.deviation_type for r in s] == [:adev, :mdev, :hdev]
        @test s.data_kind === :phase
        @test s.n == 1024
        @test s.tau0 == 1.0
        @test s.tau_mode === :Octave
        @test isnan(s.confidence)                      # ci=false ⇒ NaN
    end

    @testset "matches individual calls on the per-deviation grid" begin
        s = stability(pd; devs=(:adev, :hdev, :mtie), m=Octave, ci=false)
        @test s[:adev].dev == adev(pd, tau_values(Octave, 1024, :adev); ci=false).dev
        @test s[:adev].tau == adev(pd, tau_values(Octave, 1024, :adev); ci=false).tau
        @test s[:hdev].dev == hdev(pd, tau_values(Octave, 1024, :hdev); ci=false).dev
        @test s[:mtie].dev == mtie(pd, tau_values(Octave, 1024, :mtie)).dev
    end

    @testset "ci / confidence propagate" begin
        s = stability(pd; devs=(:adev,), ci=true, confidence=0.9)
        @test !isempty(s[:adev].ci)
        @test s.confidence == 0.9
        s0 = stability(pd; devs=(:adev,), ci=false)
        @test isempty(s0[:adev].ci)
        # neff is populated either way.
        @test s0[:adev].neff == s[:adev].neff
        @test all(>(0), s0[:adev].neff)
    end

    @testset "kwarg filtering: correct_bias reaches total family, not adev" begin
        # RWFM phase so the TOTVAR bias correction (B = 1 − 0.75·τ/T) actually
        # bites, making correct_bias=false vs =true numerically distinguishable.
        Random.seed!(2026)
        rw = PhaseData(cumsum(cumsum(randn(1024))) .* 1e-12, 1.0)
        ms = tau_values(Octave, 1024, :totdev)
        s = stability(rw; devs=(:adev, :totdev), m=Octave, ci=false, correct_bias=false)
        @test s[:adev].dev   == adev(rw, tau_values(Octave, 1024, :adev); ci=false).dev
        @test s[:totdev].dev == totdev(rw, ms; ci=false, correct_bias=false).dev
        # correct_bias=false differs from the corrected default — confirms the
        # kwarg actually reached the total family (and adev silently ignores it).
        @test s[:totdev].dev != totdev(rw, ms; ci=false, correct_bias=true).dev
    end

    @testset "default devs, explicit Vector{Int}, and TauMode taus" begin
        sd = stability(pd; ci=false)
        @test keys(sd) == collect(DEFAULT_DEVIATIONS)   # (:adev, :mdev, :hdev, :mhdev)
        se = stability(pd; devs=(:adev,), m=[1, 2, 4, 8], ci=false)
        @test se[:adev].tau == [1.0, 2.0, 4.0, 8.0]
        @test se.tau_mode === :explicit
        sde = stability(pd; devs=(:adev,), m=Decade, ci=false)
        @test sde.tau_mode === :Decade
    end

    @testset "FrequencyData entry point" begin
        Random.seed!(12)
        fd = FrequencyData(randn(512) .* 1e-9, 1.0)
        s = stability(fd; devs=(:adev, :hdev), m=Octave, ci=false)
        @test s.data_kind === :frequency
        @test s.n == 513   # 512 frequency samples integrate to 513 phase points
        @test s[:adev].dev == adev(fd, tau_values(Octave, 513, :adev); ci=false).dev
        # The grid comes from the N + 1 phase points: 1025 samples reach m = 512.
        fd2 = FrequencyData(randn(1025) .* 1e-9, 1.0)
        @test 512 in round.(Int, adev(fd2; ci=false).tau)
    end

    @testset "unknown deviation errors" begin
        @test_throws ArgumentError stability(pd; devs=(:nope,), ci=false)
    end

    @testset "show summary" begin
        s = stability(pd; devs=(:adev, :hdev), m=Octave, ci=false)
        compact = sprint(show, s)
        @test occursin("StabilitySuite", compact)
        @test occursin("2 devs", compact)
        @test occursin("no CI", compact)
        # Multi-line MIME view lists each deviation as its own row.
        pretty = sprint(show, MIME("text/plain"), s)
        @test occursin("adev", pretty)
        @test occursin("hdev", pretty)
        @test occursin("\n", pretty)
        # One row per τ with the deviations side by side, a dash where an
        # estimator has no value at that τ.
        @test occursin(r"tau +adev +hdev", pretty)
        @test occursin(r"\n +1 +\d\.\d{3}e-\d+ +\d\.\d{3}e-\d+", pretty)
        gappy = StabilitySuite([StabilityResult(:adev, [1.0, 2.0, 4.0], [3.0, 2.0, 1.0], Symbol[], SigmaTau._CIBound[], Float64[], [9, 8, 7]),
                                StabilityResult(:tdev, [1.0, 2.0], [5.0, 6.0], Symbol[], SigmaTau._CIBound[], Float64[], [9, 8])],
                               :phase, 1.0, 10, NaN, :explicit)
        @test occursin(r"\n +4 +1 +-$", sprint(show, MIME("text/plain"), gappy))
        # Past 12 rows the table keeps the first six and the last five.
        taus = 2.0 .^ (0:15)
        long = StabilitySuite([StabilityResult(:adev, taus, 1 ./ taus, Symbol[], SigmaTau._CIBound[], Float64[], fill(9, 16))],
                              :phase, 1.0, 100, NaN, :explicit)
        shown = sprint(show, MIME("text/plain"), long)
        @test occursin("⋮  5 more rows", shown) && occursin(r"\n +3\.277e\+04 +", shown)
    end
end

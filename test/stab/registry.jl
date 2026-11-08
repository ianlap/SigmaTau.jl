# test/stab/registry.jl — estimators(), deviation(), tau_mode(), summarize()
# and the shared spectrum / dynamic-window defaults.
using Test, Random, SigmaTau

@testset "Registry and shared defaults" begin
    rng = MersenneTwister(11)
    x = cumsum(randn(rng, 2048)) .* 1e-9
    pd = PhaseData(x, 1.0)
    fd = FrequencyData(diff(x), 1.0)

    @testset "estimators()" begin
        e = estimators()
        names = [r.name for r in e]
        @test length(e) == 16
        @test allunique(names)
        @test all(r -> r.family in (:allan, :total, :theo, :time), e)
        @test all(r -> r.axis in (:sigma_y, :sigma_x), e)
        @test Set(r.name for r in e if r.axis === :sigma_x) ==
              Set([:tdev, :htdev, :ttotdev, :mtie, :tierms])
        @test Set(r.name for r in e if r.correct_bias) ==
              Set([:totdev, :mtotdev, :ttotdev, :htotdev, :mhtotdev, :theo1])
        @test all(s -> s in names, DEFAULT_DEVIATIONS)
        @test all(r -> !isempty(r.label) && !isempty(r.description), e)
        e[1] = (; e[1]..., label = "x")        # a copy: the registry is unchanged
        @test estimators()[1].label == "ADEV"
    end

    @testset "deviation() matches the estimator functions" begin
        fns = Dict(:adev => adev, :mdev => mdev, :tdev => tdev, :hdev => hdev,
                   :mhdev => mhdev, :htdev => htdev, :totdev => totdev,
                   :mtotdev => mtotdev, :ttotdev => ttotdev, :htotdev => htotdev,
                   :mhtotdev => mhtotdev, :theo1 => theo1, :theoh => theoh,
                   :pdev => pdev, :mtie => mtie, :tierms => tierms)
        small = PhaseData(x[1:512], 1.0)
        for r in estimators()
            f = fns[r.name]
            a = deviation(r.name, small)
            b = f(small)
            @test a.deviation_type === r.name
            @test a.tau == b.tau
            @test isequal(a.dev, b.dev)
            @test isequal(a.ci, b.ci)
            c = deviation(r.name, small, HalfOctave; ci=false)
            @test isequal(c.dev, f(small, HalfOctave; ci=false).dev)
        end
        @test deviation(:adev, pd, 1:4; confidence=0.95).ci == adev(pd, 1:4; confidence=0.95).ci
        @test deviation(:totdev, pd, [1, 2]; correct_bias=false).dev ==
              totdev(pd, [1, 2]; correct_bias=false).dev
        @test deviation(:adev, fd, [1, 2]).dev == adev(fd, [1, 2]).dev
        @test_throws ArgumentError deviation(:nope, pd)
        @test_throws ArgumentError deviation(:nope, pd, [1, 2])
        @test_throws ArgumentError deviation(:dadev, pd)
    end

    @testset "stability runs through deviation()" begin
        s = stability(pd; devs = ["adev", :totdev], correct_bias = false)
        @test keys(s) == [:adev, :totdev]
        @test s[:totdev].dev == totdev(pd; correct_bias = false).dev
        @test_throws ArgumentError stability(pd; devs = (:nope,))
    end

    @testset "tau_mode()" begin
        for (s, m) in (("octave", Octave), ("half-octave", HalfOctave),
                       ("halfoctave", HalfOctave), ("Half_Octave", HalfOctave),
                       ("quarter-octave", QuarterOctave), ("quarteroctave", QuarterOctave),
                       ("decade", Decade), ("half-decade", HalfDecade),
                       ("halfdecade", HalfDecade), ("all", AllTaus),
                       ("alltaus", AllTaus), (" AllTaus ", AllTaus))
            @test tau_mode(s) === m
        end
        for m in instances(TauMode)
            @test tau_mode(string(m)) === m
        end
        @test_throws ArgumentError tau_mode("weekly")
    end

    @testset "summarize()" begin
        t = 0:999
        q = PhaseData(0.5e-12 .* t .^ 2, 2.0)       # y slope = 1e-12/τ0² per s
        s = summarize(q)
        @test s.n == 1000
        @test s.tau0 == 2.0
        @test s.span == 2000.0
        @test s.gaps == 0
        @test s.outliers == Int[]
        @test s.drift_per_day ≈ 0.25e-12 * 86400 rtol = 1e-9
        y = randn(MersenneTwister(3), 500) .* 1e-12
        y[100] = 1e-9
        y[200] = NaN
        f = summarize(FrequencyData(y, 1.0))
        @test f.gaps == 1
        @test f.outliers == find_outliers(FrequencyData(y, 1.0); nsigma = 5)
        @test 100 in f.outliers
        @test isfinite(f.drift_per_day)
        @test isnan(summarize(PhaseData([0.0, 1.0, 2.0])).drift_per_day)
    end

    @testset "dynamic defaults" begin
        @test SigmaTau._default_window(2048) == 256
        @test SigmaTau._default_window(100) == 16
        @test SigmaTau._default_window(10) == 10
        r = dadev(pd)
        @test r.window == 256
        @test r.tau == adev(PhaseData(x[1:256], 1.0)).tau
        @test size(r.dev, 1) == length(1:128:(2048 - 256 + 1))
        @test dhdev(pd).window == 256
        @test dadev(pd, [1, 2]).window == 256
        @test dadev(fd).window == SigmaTau._default_window(length(fd))
        @test dadev(pd, Decade; window = 512).window == 512
    end
end

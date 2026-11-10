# test/stab/api_inputs.jl — the estimators and records accept any integer
# vector of averaging factors and any real vector of samples.
using Test, Random, SigmaTau

@testset "API input types" begin
    rng = MersenneTwister(7)
    x = cumsum(randn(rng, 1024)) .* 1e-9
    pd = PhaseData(x, 1.0)
    fd = FrequencyData(diff(x), 1.0)
    ms = [1, 2, 4, 8]

    @testset "records from any real vector" begin
        @test PhaseData(view(x, 1:100)).x == x[1:100]
        @test PhaseData(1:10, 2.0).x == collect(1.0:10.0)
        @test PhaseData(Float32[1, 2, 3]) isa PhaseData{Float32}
        @test FrequencyData(view(diff(x), 1:50), 0.5).y == diff(x)[1:50]
        @test FrequencyData([1, 2, 3]).y == [1.0, 2.0, 3.0]
        @test PhaseData(1:10; source="r").source == "r"
        @test length(pd) == 1024
        @test length(fd) == 1023
        @test deviation(:adev, fd).tau == adev(fd).tau      # both on the 1024 phase samples
    end

    @testset "m as a range or other integer vector" begin
        fns = (adev, mdev, tdev, hdev, mhdev, htdev, totdev, mtotdev, ttotdev,
               htotdev, mhtotdev, mtie, tierms, pdev, theoh)
        forms = (2 .^ (0:3), Int32.(ms), view(ms, :), UInt8.(ms))
        @test isempty([(f, typeof(m)) for f in fns, m in forms
                       if !isequal(f(pd, m).dev, f(pd, ms).dev)])
        @test isempty([f for f in fns if f(fd, 2 .^ (0:3)).tau != f(fd, ms).tau])
        @test theo1(pd, 10:10:40).tau == theo1(pd, [10, 20, 30, 40]).tau
        @test dadev(pd, 1:4; window=256).tau == dadev(pd, [1, 2, 3, 4]; window=256).tau
        @test dhdev(pd, 1:4; window=256).tau == dhdev(pd, [1, 2, 3, 4]; window=256).tau
    end

    @testset "m below 1 and gaps an estimator cannot skip are errors" begin
        rejects(f, data, m) = try f(data, m); false catch e; e isa ArgumentError end
        fns = (adev, mdev, tdev, hdev, mhdev, htdev, totdev, mtotdev, ttotdev,
               htotdev, mhtotdev, mtie, tierms, pdev, theo1, theoh, dadev, dhdev)
        @test all(f -> rejects(f, pd, [0, 1]) && rejects(f, fd, [-1]), fns)
        gapped = PhaseData([x[1:99]; NaN; x[101:end]], 1.0)
        @test all(f -> rejects(f, gapped, ms), (totdev, mtotdev, ttotdev, htotdev, mhtotdev,
                                                pdev, theo1, theoh, dadev, dhdev))
    end

    @testset "identify_noise on a record" begin
        ref = identify_noise(x, ms)
        @test identify_noise(pd, ms) == ref
        @test identify_noise(pd, 2 .^ (0:3)) == ref
        @test identify_noise(pd, ms; dmax=3) == identify_noise(x, ms; dmax=3)
        @test identify_noise(fd, ms) == identify_noise(cumsum([0.0; diff(x)]), ms)
        @test length(identify_noise(pd, Octave)) == length(adev(pd, Octave).tau)
    end

    @testset "stability forwards confidence and rejects dynamic devs" begin
        # mtie / tierms have no interval model; confidence is passed through
        # for uniformity and must not error.
        s = stability(pd; devs=[:adev, :mtie, :tierms], confidence=0.95)
        @test s[:adev].ci == adev(pd; confidence=0.95).ci
        @test s[:tierms].dev == tierms(pd).dev
        err = try stability(pd; devs=[:dadev]); nothing catch e; e end
        @test err isa ArgumentError
        @test occursin("dynamic", sprint(showerror, err))
        @test stability(pd; devs=[:adev], m=1:4)[:adev].tau == adev(pd, 1:4).tau
    end
end

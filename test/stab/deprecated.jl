# test/stab/deprecated.jl — names replaced in 0.6 still work, with a warning,
# and the StabilityResult record metadata added in 0.6.
using Test, Random, SigmaTau

@testset "0.6 renames and deprecations" begin
    rng = MersenneTwister(5)
    x = cumsum(randn(rng, 1024)) .* 1e-9
    pd = PhaseData(x, 2.0)

    @testset "StabilityResult carries tau0, n and confidence" begin
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
    end

    @testset "taper and identify_noise remove_drift" begin
        @test psd_y(pd; taper = :hamming).taper === :hamming
        @test psd_y(pd).window === :hann                 # old field name still reads
        ms = [1, 2, 4]
        @test identify_noise(x, ms; remove_drift = false) ==
              (@test_deprecated identify_noise(x, ms; detrend = false))
    end

    @testset "deprecated names forward to the new ones" begin
        @test (@test_deprecated Sy(pd; nperseg = 128)).psd == psd_y(pd; nperseg = 128).psd
        @test (@test_deprecated Sx(pd; window = :hamming)).psd == psd_x(pd; taper = :hamming).psd
        @test (@test_deprecated L(pd; f_carrier = 1e7)).psd == ssb_phase_noise(pd; f_carrier = 1e7).psd
        @test (@test_deprecated stability(pd; devs = (:adev,), taus = 1:3))[:adev].tau ==
              adev(pd, 1:3).tau
        rs = [adev(PhaseData(x .+ k .* 1e-10 .* randn(rng, 1024), 2.0), Octave; ci = false) for k in 1:3]
        @test (@test_deprecated nch(rs...))[1].dev == ncornered_hat(rs...)[1].dev
        mktempdir() do d
            r = adev(pd)
            p = joinpath(d, "r.tsv")
            @test_deprecated save_result(p, r)
            @test (@test_deprecated load_result(p)).dev == r.dev
            s = stability(pd; devs = (:adev,))
            q = joinpath(d, "s.tsv")
            @test_deprecated save_suite(q, s)
            @test (@test_deprecated load_suite(q))[:adev].dev == s[:adev].dev
        end
    end
end

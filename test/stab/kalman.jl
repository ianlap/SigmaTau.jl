# kalman.jl — the single-clock filter and smoother against simulated truth,
# holdover in closed form, gaps, the innovation check, and flicker states.

using Test
using SigmaTau
using Random
using Statistics

@testset "Kalman filter" begin
    m = ClockModel(; tau0=1.0, R=1e-20, q1=1e-22, q2=1e-28)
    p, X = simulate_clock(m, 5000; rng=MersenneTwister(3))
    kf = kalman_filter(p, m)
    ks = kalman_smooth(kf)
    late = 500:5000
    rms(e) = sqrt(mean(abs2, e))

    @testset "errors match the covariance" begin
        for (r, s) in ((kf, 1), (kf, 2), (ks, 1), (ks, 2))
            @test rms(r.x[late, s] .- X[late, s]) ≈ sqrt(mean(r.P[s, s, late])) rtol=0.25   # frequency errors are slow, so few independent
        end
        @test mean(ks.P[1, 1, late]) < mean(kf.P[1, 1, late])
        @test ks.smoothed && ks.innovation == kf.innovation
        @test_throws ArgumentError kalman_smooth(ks)
    end

    @testset "innovation check" begin
        @test innovation_check(kf).ok
        wrong = innovation_check(kalman_filter(p, ClockModel(; tau0=1.0, R=1e-20, q1=1e-20, q2=1e-28)))
        @test !wrong.ok && wrong.nis < wrong.nis_bounds[1]
    end

    @testset "holdover and gaps" begin
        P = kf.P[:, :, end]
        τ = 100.0
        @test holdover_error(kf, [τ])[1] ≈ sqrt(P[1, 1] + 2τ * P[1, 2] + τ^2 * P[2, 2] + 1e-22τ + 1e-28τ^3 / 3)
        g = copy(p.x); g[2000:2100] .= NaN
        kg = kalman_filter(PhaseData(g, 1.0), m)
        @test all(isnan, kg.innovation[2000:2100]) && kg.P[1, 1, 2100] > kf.P[1, 1, 2100]
        @test kg.P[1, 1, 3000] ≈ kf.P[1, 1, 3000] rtol=1e-3       # the filter recovers
    end

    @testset "noiseless readings and drift" begin
        exact = kalman_filter(p, ClockModel(; tau0=1.0, q1=1e-22, q2=1e-28))
        @test exact.x[:, 1] ≈ p.x && maximum(exact.P[1, 1, :]) < 1e-12 * 1e-20
        m3 = ClockModel(; tau0=10.0, R=1e-20, q1=1e-22, q2=1e-28, q3=1e-40)
        p3, X3 = simulate_clock(m3, 3000; rng=MersenneTwister(4))
        k3 = kalman_filter(p3, m3)
        @test size(k3.x, 2) == 3
        @test rms(k3.x[2000:end, 3] .- X3[2000:end, 3]) < sqrt(mean(k3.P[3, 3, 2000:end]))
    end

    @testset "prediction error and tuning" begin
        line = PhaseData(1e-9 .+ 2e-12 .* (0:999), 1.0)
        @test prediction_error(kalman_filter(line, m), line, [1, 10]; burnin=100) ≈ [0, 0] atol=1e-18
        start = ClockModel(; tau0=1.0, R=1e-20, q1=1e-21, q2=1e-29)
        t = tune_q(p, start; horizons=[10, 100], points=5)
        @test t.model.q1 / t.scale ≈ 1e-22 rtol=1e-9       # on the grid, one decade down
        @test t.scale ≈ 1 rtol=0.1
        @test size(t.cost) == (5, 5) && t.edge == any(in((1, 5)), Tuple(argmin(t.cost)))
        @test minimum(t.cost) == t.cost[argmin(t.cost)] <= t.cost[3, 3]
        # Predictions do not see the size of the levels: from a start a
        # hundred times too large, the scale brings the model back.
        big = tune_q(p, ClockModel(; tau0=1.0, R=1e-18, q1=1e-19, q2=1e-27); horizons=[10, 100], points=5)
        @test [big.model.R, big.model.q1, big.model.q2] ≈ [t.model.R, t.model.q1, t.model.q2] rtol=1e-3
        @test_throws ArgumentError tune_q(p, ClockModel(; tau0=1.0, q1=1e-22))
    end

    @testset "flicker FM states" begin
        # A flicker floor of 1.5e-13 between white FM and random-walk FM, from
        # noise_gen's spectral generator. The model with Markov states
        # predicts better than the white-only fit and owns up to its error.
        h = Dict(2 => 8π^2 * 1e-26, 0 => 2e-24, -1 => 1.6e-26, -2 => 1.7e-31)
        pf = noise_gen(PhaseData, 65536, 1.0; h=h, rng=MersenneTwister(2001))
        mf = clock_model(h; tau0=1.0, flicker=true)
        kf = kalman_filter(pf, mf)
        kw = kalman_filter(pf, clock_model(powerlaw_fit(adev(pf); alphas=[2, 0, -2])))
        @test size(kf.x, 2) == 7
        @test innovation_check(kf).nis ≈ 1 rtol=0.05
        # At 100 and 1000 s the record holds enough independent spans to
        # test on; at 10 000 s it holds six.
        err = prediction_error(kf, pf, [100, 1000]; burnin=8192)
        @test err[2] < 0.9 * prediction_error(kw, pf, [1000]; burnin=8192)[1]
        @test holdover_error(kf, [100.0, 1000.0]) ≈ err rtol=0.3
        @test tune_q(pf, mf; points=3).model.fm_band == mf.fm_band
        # A record of the model itself has flicker FM's ADEV, and the filter
        # follows its frequency, Markov states included.
        ms = ClockModel(; tau0=1.0, flicker_fm=1e-26)
        ps, Xs = simulate_clock(ms, 2^17; rng=MersenneTwister(3))
        @test size(Xs, 2) == 9
        @test adev(ps, [1, 16, 256]; ci=false).dev ≈ powerlaw_sigma(Dict(-1 => 1e-26), [1.0, 16.0, 256.0]; tau0=1.0) rtol=0.1
        @test cor(vec(sum(kalman_filter(ps, ms).x[:, 2:9]; dims=2)), vec(sum(Xs[:, 2:9]; dims=2))) > 0.85
    end

    @testset "flicker PM states" begin
        # Flicker PM on white FM and random-walk FM. With Markov phase states
        # the innovations have the claimed variance, where the fold into R
        # overstates it, and the next reading is predicted better.
        h = Dict(1 => 2.2e-22, 0 => 2e-24, -2 => 1.7e-31)
        pp = noise_gen(PhaseData, 65536, 1.0; h=h, rng=MersenneTwister(5))
        mp = clock_model(h; tau0=1.0, flicker=true)
        kp, k0 = kalman_filter(pp, mp), kalman_filter(pp, clock_model(h; tau0=1.0))
        @test (size(kp.x, 2), size(k0.x, 2)) == (6, 2)
        @test innovation_check(kp).nis ≈ 1 rtol=0.05
        @test innovation_check(k0).nis < 0.7
        @test prediction_error(kp, pp, [1]; burnin=4096)[1] < 0.95 * prediction_error(k0, pp, [1]; burnin=4096)[1]
        # A record of the model itself: flicker PM's ADEV, white innovations,
        # and a phase estimate as good as it claims, the reading's flicker
        # kept out of it.
        ps, Xs = simulate_clock(mp, 2^16; rng=MersenneTwister(5))
        ks = kalman_filter(ps, mp)
        @test adev(ps, [1, 8, 64]; ci=false).dev ≈ powerlaw_sigma(h, [1.0, 8.0, 64.0]; tau0=1.0) rtol=0.05
        @test innovation_check(ks; confidence=0.9999).ok
        e = rms(ks.x[1000:end, 1] .- Xs[1000:end, 1])
        @test e ≈ sqrt(mean(ks.P[1, 1, 1000:end])) rtol=0.15
        @test e < 0.9 * rms(ps.x[1000:end] .- Xs[1000:end, 1])
    end

    @test_throws ArgumentError kalman_filter(PhaseData(p.x, 2.0), m)
    @test_throws ArgumentError kalman_filter(p, m; x0=[0.0, 0.0, 0.0], P0=zeros(3, 3))
    @test_throws "no noise" kalman_filter(p, ClockModel(; tau0=1.0))
end

# kalman.jl — the single-clock filter and smoother against simulated truth,
# holdover in closed form, gaps, and the innovation check.

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

    @test_throws ArgumentError kalman_filter(PhaseData(p.x, 2.0), m)
    @test_throws ArgumentError kalman_filter(p, m; x0=[0.0, 0.0, 0.0], P0=zeros(3, 3))
end

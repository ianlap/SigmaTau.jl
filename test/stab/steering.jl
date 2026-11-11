# steering.jl — pole placement against the formulas of Matsakis (2019, 2020),
# cost-minimizing gains against those printed by Koppang & Leland (1999), the
# replayed loop against the closed-loop recursion, and the predicted steady
# state against a simulated clock.

using Test
using SigmaTau
using LinearAlgebra
using Random
using Statistics

@testset "Steering" begin
    m = ClockModel(; tau0=1.0, R=1e-20, q1=1e-22, q2=1e-26)
    g = steering_gains(30.0, 1.0)

    @testset "pole placement" begin
        r = exp(-1 / 30)
        @test g.phase ≈ (1 - r)^2 && g.frequency ≈ 1 - r^2 && g.integral == 0     # eqs 3.10, 3.11
        @test steering_poles(steering_gains(1e5, 60.0)).time_constants ≈ [1e5, 1e5] rtol=1e-6
        # Off the critical line: the poles of eq. 3.9, their decay (3.12) and ringing (3.13).
        τ, g₁, g₂ = 2.0, 0.02, 0.3
        p = steering_poles(SteeringGains(g₁, g₂, 0, τ))
        mid = (2 - τ * g₁ - g₂) / 2
        @test sort(p.poles; by=imag) ≈ [mid - sqrt(complex(mid^2 - (1 - g₂))), mid + sqrt(complex(mid^2 - (1 - g₂)))]
        @test p.time_constants ≈ fill(-τ / log(sqrt(1 - g₂)), 2)
        @test p.frequencies ≈ fill(atan(sqrt(4(1 - g₂) / (2 - τ * g₁ - g₂)^2 - 1)) / (2π * τ), 2)
        @test steering_gains(p.poles, τ).phase ≈ g₁ && steering_gains(p.poles, τ).frequency ≈ g₂
        # The edge of stability, 4 − τg₁ − 2g₂ = 0 (eq. 3.17).
        @test steering_poles(SteeringGains(0.99, 1.5, 0, 1.0)).stable
        @test !steering_poles(SteeringGains(1.01, 1.5, 0, 1.0)).stable
        @test !steering_poles(SteeringGains(-1e-3, 0.5, 0, 1.0)).stable
        # Three gains: the poles are the eigenvalues of Φ − BG on (Σp, p, f),
        # the roots of the cubic of Matsakis (2020), eq. 3.10.
        z = [0.9, 0.8 + 0.1im, 0.8 - 0.1im]
        k = steering_gains(z, τ)
        A = [1 1 τ; 0 1 τ; 0 0 1] .- [τ, τ, 1] * [k.integral k.phase k.frequency]
        @test sort(eigvals(A); by=imag) ≈ sort(z; by=imag)
        @test sort(steering_poles(k).poles; by=imag) ≈ sort(z; by=imag)
        J, P, D = τ * k.integral, τ * k.phase, k.frequency
        @test all(abs(x^3 + (J + P + D - 3) * x^2 + (3 - P - 2D) * x + D - 1) < 1e-12 for x in z)
        c = steering_gains(30.0, 1.0; integral=true)
        @test [c.integral, c.phase, c.frequency] ≈ [1 - 3r + 3r^2 - r^3, 1 - 3r^2 + 2r^3, 1 - r^3]   # eqs 3.15–3.17
        @test steering_poles(c).time_constants ≈ fill(30.0, 3) rtol=1e-3
        @test_throws ArgumentError steering_gains([0.9, 0.8 + 0.1im], 1.0)
        @test_throws ArgumentError steering_gains([0.9], 1.0)
        @test_throws ArgumentError steering_gains(-1.0, 1.0)
    end

    @testset "quadratic cost" begin
        k(τ; kw...) = (l = lqg_gains(τ; kw...); (l.phase, l.frequency))
        near(a, b; rtol) = all(isapprox.(a, b; rtol))
        # The three designs of Koppang & Leland (1999), §IV, a steer each hour, as printed.
        @test near(k(3600.0; phase=1e-3, frequency=1e-3, steer=1e5), (6.503e-5, 0.5771); rtol=5e-4)
        @test near(k(3600.0; phase=1e-3, frequency=1e-3, steer=1e12), (3.138e-8, 0.01498); rtol=5e-4)
        @test near(k(3600.0; phase=1e-4, frequency=1e9, steer=1e12), (9.8376e-9, 0.03221); rtol=5e-4)
        # Free steers take the whole error out at once; costly ones tend to the
        # gains of the continuous double integrator, √ρ and √(2√ρ + ρ_f).
        @test near(k(60.0; phase=1.0, steer=1e-30), (1 / 60, 1.0); rtol=1e-9)
        @test near(k(1.0; phase=1e-16, steer=1.0), (1e-8, sqrt(2e-8)); rtol=1e-3)
        @test near(k(1.0; phase=1e-16, frequency=1e-6, steer=1.0), (1e-8, sqrt(2e-8 + 1e-6)); rtol=1e-3)
        l = lqg_gains(3600.0; phase=1e-3, steer=1e12)
        @test l.integral == 0 && l.interval == 3600 && steering_poles(l).stable
        # The cost that the predicted steady state gives is least at these gains.
        w = (phase=1.0, frequency=1e3, steer=1e6)
        J(g) = (q = steer(m, g); w.phase * q.phase^2 + w.frequency * q.frequency^2 + w.steer * q.steers^2)
        l = lqg_gains(1.0; w...)
        @test all(J(SteeringGains(a * l.phase, b * l.frequency, 0, 1.0)) > J(l) for (a, b) in ((1.05, 1), (0.95, 1), (1, 1.05), (1, 0.95)))
        @test_throws ArgumentError lqg_gains(1.0; phase=0.0, steer=1.0)
    end

    p, X = simulate_clock(m, 100_000; rng=MersenneTwister(5))
    s = steer(p, m, g)
    rms(v) = sqrt(mean(abs2, v))

    @testset "replay" begin
        # The steered estimate obeys x̂⁺ = (Φ − BG) x̂ + K ν⁺ (eqs A.11, A.12);
        # K ν⁺ is what the free-running filter adds to its own prediction.
        kf = kalman_filter(p, m)
        Φ = [1 1; 0 1.0]
        A = Φ .- [1, 1] * [g.phase g.frequency]
        x̂ = [s.phase s.frequency]'
        @test x̂[:, 2:end] .- A * x̂[:, 1:end-1] ≈ kf.x'[:, 2:end] .- Φ * kf.x'[:, 1:end-1] rtol=1e-9
        @test cumsum(s.steers)[1:end-1] ≈ s.correction[2:end]
        @test s.steered.x .- p.x ≈ cumsum(s.correction)
        # A steer every tenth sample, with the gains for that interval.
        t = steer(p, m, steering_gains(30.0, 10.0))
        @test findall(!iszero, t.steers) == 10:10:100_000
        @test rms(t.steered.x) ≈ rms(s.steered.x) rtol=0.02
        # Under a drift the model does not hold the phase stands off. A drift
        # state steers it out; an integral gain centres the estimate.
        d = PhaseData(p.x .+ 1e-11 / 2 .* (0:99_999) .^ 2, 1.0)
        off = mean(steer(d, m, g).steered.x[1000:end])
        md = ClockModel(; tau0=1.0, R=m.R, q1=m.q1, q2=m.q2, drift=true)
        @test off > 1e-8 && abs(mean(steer(d, md, g).steered.x[1000:end])) < 1e-3 * off
        @test abs(mean(steer(d, m, steering_gains(30.0, 1.0; integral=true)).phase[1000:end])) < 1e-3 * off
        x = copy(p.x); x[500:520] .= NaN
        @test count(isnan, steer(PhaseData(x, 1.0), m, g).steered.x) == 21
        @test_throws ArgumentError steer(p, m, steering_gains(30.0, 1.5))
    end

    @testset "steady state" begin
        q = steer(m, g; tau=[1.0, 16.0, 256.0, 16384.0])
        e = X[1000:end, 1] .+ (s.steered.x.-p.x)[1000:end]          # the clock itself, steered
        @test rms(e) ≈ q.phase rtol=0.1
        @test rms(X[1000:end, 2] .+ s.correction[1000:end]) ≈ q.frequency rtol=0.1
        @test rms(s.steers[1000:end]) ≈ q.steers rtol=0.05
        a = adev(PhaseData(e, 1.0), [1, 16, 256]; ci=false).dev
        @test a[1:2] ≈ q.adev[1:2] rtol=0.03
        @test a[3] ≈ q.adev[3] rtol=0.15                             # a few hundred loop times in the record
        @test q.adev[4] ≈ sqrt(3) * q.phase / 16384 rtol=1e-3        # far past the loop: uncorrelated phase
        # The settled filter is the one a long run reaches.
        f = SigmaTau._steady_filter(m)
        @test (I - f.K * f.h') * f.P ≈ kalman_filter(p, m).P[:, :, end] rtol=1e-6
        # Flicker FM and PM as Markov states, with an integral gain.
        fm = ClockModel(; tau0=1.0, R=1e-24, q1=1e-24, flicker_fm=1e-24, flicker_pm=1e-22)
        c = steering_gains(30.0, 1.0; integral=true)
        pf, Xf = simulate_clock(fm, 100_000; rng=MersenneTwister(6))
        sf = steer(pf, fm, c)
        @test rms(Xf[1000:end, 1] .+ (sf.steered.x.-pf.x)[1000:end]) ≈ steer(fm, c).phase rtol=0.1
        @test rms(sf.steers[1000:end]) ≈ steer(fm, c).steers rtol=0.05
        @test steer(m, g).tau == 2.0 .^ (0:9)
        @test_throws ArgumentError steer(m, steering_gains(30.0, 10.0))
        @test_throws ArgumentError steer(m, SteeringGains(1.01, 1.5, 0, 1.0))
        @test_throws ArgumentError steer(ClockModel(; tau0=1.0, q1=1e-22, q3=1e-40), g)
    end
end

# ensemble_validation.jl — the Kalman smoother and the timescales against independent references.

using Test
using SigmaTau
using Random
using Statistics
using LinearAlgebra

@testset "ensemble validation" begin
    @testset "smoother = exact batch posterior, gap included" begin
        m = ClockModel(; tau0=1.0, R=1e-20, q1=1e-22, q2=1e-28)
        p, _ = simulate_clock(m, 60; rng=MersenneTwister(21))
        z = copy(p.x); z[20:30] .= NaN; z[1:3] .= NaN
        x0, P0 = zeros(2), Diagonal([1e-14, 1e-20])
        ks = kalman_smooth(kalman_filter(PhaseData(z, m.tau0), m; x0, P0))
        # Information form of the Gaussian over the whole state path, in the prior's units.
        Φ, Q = SigmaTau.state_transition(m, 1.0), SigmaTau.process_noise(m, 1.0)
        sc = sqrt.(diag(P0)); T = Diagonal(sc)
        Φs, Qi = T \ Φ * T, inv(T \ Q / T)
        N = length(z); Λ = zeros(2N, 2N); η = zeros(2N); ix(k) = 2k-1:2k
        Λ[ix(1), ix(1)] += inv(T \ P0 / T)
        for k in 2:N
            a, b = ix(k-1), ix(k)
            Λ[a, a] += Φs' * Qi * Φs; Λ[b, b] += Qi; Λ[a, b] -= Φs' * Qi; Λ[b, a] -= Qi * Φs
        end
        for k in findall(isfinite, z)
            Λ[2k-1, 2k-1] += sc[1]^2 / m.R; η[2k-1] += sc[1] * z[k] / m.R
        end
        μ = inv(Symmetric(Λ)) * η
        xb = reduce(hcat, (sc .* μ[ix(k)] for k in 1:N))'
        @test maximum(abs.(ks.x .- xb) ./ sqrt.(diag(inv(Symmetric(Λ)))[1:2:end] .* sc[1]^2)) < 1e-6
    end

    @testset "timescale :rkt = Greenhall's optimal weights" begin
        ms = [ClockModel(; tau0=1.0, q1=q1, q2=q2) for (q1, q2) in ((1e-24, 1e-31), (4e-24, 1e-30), (1e-23, 1e-29))]
        cc, _ = simulate_ensemble(ms, 2000; reference=0, rng=MersenneTwister(2))
        rk = timescale(cc; algorithm=:rkt, models=ms)
        # Unreduced Jones-Tryon filter on the differences to clock 1 (Greenhall 2003, Thm 1 and 2).
        k = 3; Φ = Matrix{Float64}(I, 2k, 2k); Q = zeros(2k, 2k); X = zeros(2k); P = zeros(2k, 2k)
        for (i, m) in enumerate(ms)
            j = [i, k+i]; x0, P0 = SigmaTau._kalman_start(cc.x[:, i], 1.0, 2)
            Φ[j, j] = SigmaTau.state_transition(m, 1.0); Q[j, j] = SigmaTau.process_noise(m, 1.0)
            X[j] = x0; P[j, j] = P0
        end
        H = zeros(k-1, 2k); H[:, 1] .= -1; H[1, 2] = 1; H[2, 3] = 1
        Pyy = zeros(k, k); dy = 0.0; dw = 0.0
        for n in 1:2000
            if n > 1; X = Φ*X; P = Φ*P*Φ' + Q; end
            K = P*H' / (H*P*H'); X += K*(cc.x[n, 2:k] .- cc.x[n, 1] .- H*X); P = Symmetric(P - K*H*P)
            dy = max(dy, maximum(abs.(rk.y[n, :] .- X[k+1:2k]) ./ sqrt.(diag(P)[k+1:2k])))
            if n > 3
                λ = (Pyy + Diagonal([m.q1 + m.q2/3 for m in ms])) \ ones(k)
                dw = max(dw, maximum(abs.(rk.weights[n, :] .- λ ./ sum(λ))))
            end
            Pyy = Matrix(P[k+1:2k, k+1:2k])
        end
        @test dy < 1e-6      # the reduction leaves the frequency estimates alone
        @test dw < 1e-6      # and its weights minimise the variance of the scale step
    end

    @testset "EWFA does not step when a clock drops out" begin
        ms = [ClockModel(; tau0=1.0, q1=q, q2=1e-34) for q in (1e-26, 2e-26, 3e-26, 5e-26)]
        cc, _ = simulate_ensemble(ms, 6000; reference=0, rng=MersenneTwister(5))
        x = cc.x .+ (0:5999) .* [0, 3e-11, -2e-11, 1e-11]'      # distinct frequency offsets
        xg = copy(x); xg[2000:4000, 3] .= NaN
        u = timescale_phase(timescale(ClockComparisons(x, 1.0); algorithm=:ewfa)).x
        ug = timescale_phase(timescale(ClockComparisons(xg, 1.0); algorithm=:ewfa)).x
        @test abs(ug[2000] - ug[1999] - (u[2000] - u[1999])) < 1e-12
    end

    @testset "AT1 with frequency offsets and short or gappy records" begin
        ms = [ClockModel(; tau0=1.0, q1=1e-24) for _ in 1:4]
        cc, _ = simulate_ensemble(ms, 4000; reference=0, rng=MersenneTwister(2))
        x = copy(cc.x); x[:, 2] .+= 1e-10 .* (0:3999)             # 100× the one-second noise
        @test_broken all(timescale(ClockComparisons(x, 1.0)).weights[end, :] .> 0)   # all clocks rejected at sample 2
        @test_broken timescale(ClockComparisons(cc.x[1:12, :], 1.0)) isa Any          # BoundsError under 16 samples
        xn = copy(cc.x); xn[3, 2] = NaN
        @test_broken timescale(ClockComparisons(xn, 1.0)) isa Any                     # NaN in the first 10 %
    end
end

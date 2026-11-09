# timescale.jl — EWFA and AT1 on simulated ensembles with a known truth.

using Test
using SigmaTau
using Random
using Statistics

@testset "timescales" begin
    q = [1.0, 4.0, 9.0, 16.0] .* 1e-24                    # white FM levels
    ms = [ClockModel(; tau0=1.0, q1=qi) for qi in q]
    cc, _ = simulate_ensemble(ms, 50000; reference=0, rng=MersenneTwister(1))
    adev1(r) = adev(timescale_phase(r), [1]; ci=false).dev[1]
    late = 25000:50000

    ew = timescale(cc; algorithm=:ewfa)
    @test all(ew.weights .≈ 0.25)
    @test adev1(ew) ≈ sqrt(sum(q)) / 4 rtol=0.05            # σ² of an equal mean

    at = timescale(cc)
    wopt = (1 ./ q) ./ sum(1 ./ q)
    @test vec(mean(at.weights[late, :]; dims=1)) ≈ wopt rtol=0.05
    @test adev1(at) ≈ sqrt(sum(wopt .^ 2 .* q)) rtol=0.05   # the inverse-variance optimum
    @test adev1(at) < sqrt(q[1])                            # better than the best clock
    @test all(sum(at.weights; dims=2) .≈ 1)
    @test timescale_phase(at).x ≈ at.x[:, 3] .- cc.x[:, 3]   # the same through any clock
    # A short record where the Groslambert start fails for the best clock.
    rw = [ClockModel(; tau0=1.0, q1=1e-24 * 4.0^k, q2=1e-30 * 4.0^k) for k in 0:3]
    short, _ = simulate_ensemble(rw, 4096; reference=0, rng=Xoshiro(42))
    @test all(<(0), diff(vec(mean(timescale(short).weights[2048:end, :]; dims=1))))

    @testset "caps, gaps and outliers" begin
        capped = timescale(cc; wmax=0.4)
        @test maximum(capped.weights[late, :]) <= 0.4 + 1e-12
        x = copy(cc.x); x[30000:30100, 1] .= NaN; x[40000, 2] += 1e-9
        r = timescale(ClockComparisons(x, 1.0))
        @test all(iszero, r.weights[30000:30100, 1]) && all(isnan, r.error[30000:30100, 1])
        @test r.weights[40000, 2] == 0                      # a 1 ns jump is left out
        @test abs(timescale_phase(r).x[40000] - timescale_phase(at).x[40000]) < 1e-11
    end

    @testset "reduced Kalman" begin
        rk = timescale(cc; algorithm=:rkt, models=ms)
        @test vec(mean(rk.weights[late, :]; dims=1)) ≈ wopt rtol=0.01   # Greenhall 2003: ∝ 1/q1
        @test adev1(rk) ≈ sqrt(sum(wopt .^ 2 .* q)) rtol=0.05
        # Dissimilar clocks with white PM: quieter than every clock out to 100 s.
        mix = [ClockModel(; tau0=1.0, R=1e-24, q1=isodd(k) ? 1e-22 : 1e-24, q2=isodd(k) ? 1e-32 : 1e-28) for k in 1:4]
        cm, _ = simulate_ensemble(mix, 50000; reference=0, rng=Xoshiro(3))
        τs = [1, 10, 100]
        best = minimum(reduce(hcat, (adev(cm[i], τs; ci=false).dev for i in 1:4)); dims=2)
        @test all(adev(timescale_phase(timescale(cm; algorithm=:rkt, models=mix)), τs; ci=false).dev .< vec(best))
        x = copy(cc.x); x[30000:30100, 1] .= NaN
        g = timescale_phase(timescale(ClockComparisons(x, 1.0); algorithm=:rkt, models=ms))
        @test maximum(abs, g.x[29000:31000] .- timescale_phase(rk).x[29000:31000]) < 5e-11
        @test_throws ArgumentError timescale(cc; algorithm=:rkt)
    end

    @test_throws ArgumentError timescale(cc; algorithm=:at3)
    @test_throws ArgumentError timescale(cc; wmax=0.2)
end

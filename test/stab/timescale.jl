# timescale.jl — EWFA, AT1 and the reduced Kalman timescale on simulated ensembles with a known truth.

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

    wopt = (1 ./ q) ./ sum(1 ./ q)
    wlate(r) = vec(mean(r.weights[late, :]; dims=1))
    at = timescale(cc)                                      # the paper's bias term K (eqs 5-6)
    @test wlate(at) ≈ wopt rtol=0.1
    @test adev1(at) ≈ sqrt(sum(wopt .^ 2 .* q)) rtol=0.05   # the inverse-variance optimum
    @test adev1(at) < sqrt(q[1])                            # better than the best clock
    @test all(sum(at.weights; dims=2) .≈ 1)
    @test timescale_phase(at).x ≈ at.x[:, 3] .- cc.x[:, 3]   # the same through any clock
    ex = timescale(cc; bias=:exact)                         # exact for independent errors
    @test wlate(ex) ≈ wopt rtol=0.02
    # The paper's constants: m from τ_min (eq. 10), per clock.
    @test SigmaTau._filter_constant(1.0, 1.0) ≈ (sqrt(5 / 3) - 1) / 2
    @test SigmaTau._filter_constant(100.0, 1.0) > SigmaTau._filter_constant(10.0, 1.0)
    @test adev1(timescale(cc; tau_min=[50, 200, 200, 200])) ≈ adev1(at) rtol=0.05
    # A short record where the Groslambert-style start is poor for the best clock.
    rw = [ClockModel(; tau0=1.0, q1=1e-24 * 4.0^k, q2=1e-30 * 4.0^k) for k in 0:3]
    short, _ = simulate_ensemble(rw, 4096; reference=0, rng=Xoshiro(42))
    @test all(<(0), diff(vec(mean(timescale(short).weights[2048:end, :]; dims=1))))

    @testset "caps, gaps and outliers" begin
        capped = timescale(cc; wmax=0.4)
        @test maximum(capped.weights[late, :]) <= 0.4 + 1e-12
        x = copy(cc.x); x[30000:30100, 1] .= NaN
        gap = timescale(ClockComparisons(x, 1.0))
        @test all(iszero, gap.weights[30000:30100, 1]) && all(isnan, gap.error[30000:30100, 1])
        clean = timescale(ClockComparisons(x, 1.0); outlier=4)
        x[40000, 2] += 1e-9
        r = timescale(ClockComparisons(x, 1.0); outlier=4)
        @test r.weights[40000, 2] == 0                      # a 1 ns jump is left out
        @test abs(timescale_phase(r).x[40000] - timescale_phase(clean).x[40000]) < 1e-11
        @test timescale(ClockComparisons(x, 1.0)).weights[40000, 2] > 0   # not by default
    end

    @testset "offsets, gaps near the start and short records (AT1)" begin
        # Frequency offsets 1000 times the clocks' noise: weights unchanged,
        # also with the outlier test on, which never rejects every clock.
        off = ClockComparisons(cc.x .+ (0:size(cc.x, 1)-1) .* [1e-9 -2e-9 3e-9 0.5e-9], 1.0)
        @test wlate(timescale(off)) ≈ wlate(at) rtol=0.1
        @test wlate(timescale(off; outlier=4)) ≈ wlate(at) rtol=0.1
        # τ₀ = 1000 s keeps its weights.
        slow = timescale(ClockComparisons(cc.x[1:5000, :], 1000.0))
        @test vec(mean(slow.weights[2500:end, :]; dims=1)) ≈ wopt rtol=0.15
        x = copy(cc.x); x[1:5, 2] .= NaN; x[7, 3] = NaN
        @test wlate(timescale(ClockComparisons(x, 1.0))) ≈ wlate(at) rtol=0.1
        for n in (3, 10, 15)
            r = timescale(ClockComparisons(cc.x[1:n, :], 1.0))
            @test all(isfinite, r.x) && all(sum(r.weights; dims=2) .≈ 1)
        end
    end

    @testset "offsets, late starters and short records (EWFA)" begin
        # Frequency offsets of 1 ns/s per clock change nothing, a clock that
        # starts late joins at its first reading, and a short record runs.
        off = cc.x .+ (0:size(cc.x, 1)-1) .* [1e-9 -2e-9 3e-9 0.5e-9]
        @test all(timescale(ClockComparisons(off, 1.0); algorithm=:ewfa).weights[end, :] .≈ 0.25)
        x = copy(cc.x); x[1:5, 2] .= NaN
        joined = timescale(ClockComparisons(x, 1.0); algorithm=:ewfa)
        @test all(isfinite, timescale_phase(joined).x) && all(iszero, joined.weights[1:5, 2])
        @test joined.weights[end, 2] ≈ 0.25
        @test abs(adev1(joined) - adev1(ew)) < 0.05 * adev1(ew)
        @test size(timescale(ClockComparisons(cc.x[1:3, :], 1.0); algorithm=:ewfa).x, 1) == 3
    end

    @testset "Greenhall Kalman" begin
        rk = timescale(cc; algorithm=:greenhall, models=ms)
        @test vec(mean(rk.weights[late, :]; dims=1)) ≈ wopt rtol=0.01   # Greenhall 2003: ∝ 1/q1
        @test adev1(rk) ≈ sqrt(sum(wopt .^ 2 .* q)) rtol=0.05
        # Dissimilar clocks with white PM: quieter than every clock out to 100 s.
        mix = [ClockModel(; tau0=1.0, R=1e-24, q1=isodd(k) ? 1e-22 : 1e-24, q2=isodd(k) ? 1e-32 : 1e-28) for k in 1:4]
        cm, _ = simulate_ensemble(mix, 50000; reference=0, rng=Xoshiro(3))
        τs = [1, 10, 100]
        best = minimum(reduce(hcat, (adev(cm[i], τs; ci=false).dev for i in 1:4)); dims=2)
        @test all(adev(timescale_phase(timescale(cm; algorithm=:greenhall, models=mix)), τs; ci=false).dev .< vec(best))
        x = copy(cc.x); x[30000:30100, 1] .= NaN
        g = timescale_phase(timescale(ClockComparisons(x, 1.0); algorithm=:greenhall, models=ms))
        @test maximum(abs, g.x[29000:31000] .- timescale_phase(rk).x[29000:31000]) < 5e-11
        @test_throws ArgumentError timescale(cc; algorithm=:greenhall)
        # Noiseless readings: every pivot gives the same scale; with white PM
        # they differ by the pivots' own reading noise.
        @test timescale_phase(timescale(cc; algorithm=:greenhall, models=ms, pivot=3)).x ≈ timescale_phase(rk).x atol=1e-12
        @test_throws ArgumentError timescale(cc; algorithm=:greenhall, models=ms, pivot=5)
        # Jones–Tryon: the same filter unreduced; the scale is a clock's reading
        # less its phase estimate, equal at every clock for noiseless readings.
        jt = timescale(cc; algorithm=:jones_tryon, models=ms)
        @test jt.algorithm === :jones_tryon && all(isnan, jt.weights)
        @test maximum(abs, jt.y .- rk.y) < 1e-3 * std(rk.y[late, 1])   # frequency estimates agree (Theorem 1)
        @test timescale_phase(timescale(cc; algorithm=:jones_tryon, models=ms, pivot=3)).x ≈ timescale_phase(jt).x atol=1e-12
        @test adev1(rk) < adev1(jt) < 1.6 * adev1(rk)             # the reduction quiets the scale at short τ
        a1 = timescale_phase(timescale(cm; algorithm=:greenhall, models=mix)).x
        a2 = timescale_phase(timescale(cm; algorithm=:greenhall, models=mix, pivot=2)).x
        @test 1e-13 < std(a1 .- a2) < 3e-12                    # about the readings' noise, √R = 1e-12 s
    end

    @test_throws ArgumentError timescale(cc; algorithm=:at3)
    @test_throws ArgumentError timescale(cc; wmax=0.2)
end

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

    @testset "EWFA is the plain mean of the phases" begin
        # Kettering 2023, section 3.1: the mean of the clocks' frequency steps
        # steers the scale, which without gaps is the mean of their phases.
        m = timescale_phase(ew).x
        @test m .- m[1] ≈ -(vec(sum(cc.x; dims=2)) ./ size(cc.x, 2) .- sum(cc.x[1, :]) / size(cc.x, 2)) atol=1e-15
        @test_throws MethodError timescale(cc; algorithm=:ewfa, tau_y=10.0)
    end

    @testset "AT2" begin
        # Same weights as AT1 (the time update is AT1's); frequency by the
        # Kalman update, which in steady state is AT1's filter (eq 17).
        a2 = timescale(cc; algorithm=:at2)
        @test wlate(a2) ≈ wlate(at) rtol=0.05
        @test all(sum(a2.weights; dims=2) .≈ 1)
        τs = [1, 10, 100, 1000]
        @test adev(timescale_phase(a2), τs; ci=false).dev ≈ adev(timescale_phase(at), τs; ci=false).dev rtol=0.1
        # Units of eqs 13-18: from m(m+1) = σα²/(σ_β²τ) the steady state of the
        # recursion has gain P̂/(σα²+P̂) = 1/(1+m), the AT1 constant of eq 17.
        for (m, σα²) in ((0.5, 4e-26), (30.0, 1e-24), (4000.0, 3e-23))
            s = σα² / (m * (m + 1))
            P = SigmaTau._at2_steady(σα², s)
            P̂ = P + s
            @test σα² * P̂ / (σα² + P̂) ≈ P rtol=1e-12            # a fixed point of (15)
            @test σα² / P̂ ≈ m rtol=1e-9                           # eq 17 with τ-scaled σ_β
        end
        # Random-walk FM clocks: the same steady state, so the same stability.
        rw2 = timescale(short; algorithm=:at2)
        @test adev(timescale_phase(rw2), [1, 100]; ci=false).dev ≈ adev(timescale_phase(timescale(short)), [1, 100]; ci=false).dev rtol=0.15
        x = copy(cc.x); x[30000:30100, 1] .= NaN; x[1:5, 2] .= NaN
        g = timescale(ClockComparisons(x, 1.0); algorithm=:at2)
        @test all(isfinite, timescale_phase(g).x) && all(iszero, g.weights[30000:30100, 1])
        @test all(isfinite, timescale(ClockComparisons(cc.x[1:10, :], 1.0); algorithm=:at2).x)
    end

    @testset "AT2 frequency-step detector" begin
        sms = [ClockModel(; tau0=1.0, q1=q, q2=q * 1e-8) for q in (1e-22, 2e-22, 4e-22, 8e-22)]
        sc, _ = simulate_ensemble(sms, 8000; reference=0, rng=MersenneTwister(5))
        sx = copy(sc.x); sx[4001:end, 2] .+= (1:4000) .* 5e-11     # clock 2 steps in frequency
        stepped = ClockComparisons(sx, 1.0)
        off = timescale(stepped; algorithm=:at2)
        on = timescale(stepped; algorithm=:at2, steps=true)
        # Caught at once, the clock kept out of the scale for τ_min and the scale
        # rerun from the step: much steadier at long τ than letting the filter adapt.
        @test all(iszero, on.weights[4001:4020, 2]) && off.weights[4005, 2] > 0.1
        @test adev(timescale_phase(on), [1000]; ci=false).dev[1] < 0.6 * adev(timescale_phase(off), [1000]; ci=false).dev[1]
        @test timescale(stepped; algorithm=:at2, steps=false).x == off.x
        # :real keeps what was issued: identical up to the catch, and the scale
        # stays continuous where the rerun corrects it.
        re = timescale(stepped; algorithm=:at2, steps=:real)
        @test re.weights[1:4000, :] == on.weights[1:4000, :] && re.x[1:4000, :] == off.x[1:4000, :]
        @test all(iszero, re.weights[4010:4020, 2]) && all(isfinite, timescale_phase(re).x)
        @test maximum(abs, diff(timescale_phase(re).x)) < 20 * maximum(abs, diff(timescale_phase(off).x))
        @test adev(timescale_phase(re), [1000]; ci=false).dev[1] < adev(timescale_phase(off), [1000]; ci=false).dev[1]
        @test timescale(stepped; algorithm=:at2, steps=true).x == on.x
        @test_throws ArgumentError timescale(stepped; algorithm=:at2, steps=:fast)
        # Eq 20 reads the same however the units are chosen: a scaled phase
        # (all variances ∝ c²) and a scaled τ₀ (frequencies ∝ 1/a) change nothing.
        for (c, a) in ((1e3, 1.0), (1.0, 10.0))
            u = timescale(ClockComparisons(sx .* c, a); algorithm=:at2, steps=true)
            @test u.weights ≈ on.weights rtol=1e-6
        end
        @test_throws ArgumentError timescale(stepped; steps=true)
        # At the default 5σ a step-free ensemble triggers nothing (at the paper's
        # 4σ the noisier clocks do: 53 false steps in 24 such runs).
        free = timescale(sc; algorithm=:at2)
        @test timescale(sc; algorithm=:at2, steps=:best).weights == free.weights
        @test timescale(sc; algorithm=:at2, steps=:real).weights == free.weights
    end

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
        # The implicit mean of the true phases: the readout at R = 0, quieter than it with white PM.
        cq, tq = simulate_ensemble(ms, 5000; reference=0, rng=MersenneTwister(9))
        rq = timescale(cq; algorithm=:greenhall, models=ms)
        dsum = PhaseData(timescale_phase(rq).x .+ implicit_mean(rq, tq).x, 1.0)    # equal up to sign and a frequency offset
        @test adev(dsum, [1, 10, 100]; ci=false).dev[1] < 1e-3 * adev(timescale_phase(rq), [1]; ci=false).dev[1]
        cn, tn = simulate_ensemble(mix, 5000; reference=0, rng=MersenneTwister(9))
        rn = timescale(cn; algorithm=:greenhall, models=mix)
        @test adev(implicit_mean(rn, tn), [1]; ci=false).dev[1] < 0.8 * adev(timescale_phase(rn), [1]; ci=false).dev[1]
        @test_throws ArgumentError implicit_mean(timescale(cc), tq)
        a1 = timescale_phase(timescale(cm; algorithm=:greenhall, models=mix)).x
        a2 = timescale_phase(timescale(cm; algorithm=:greenhall, models=mix, pivot=2)).x
        @test 1e-13 < std(a1 .- a2) < 3e-12                    # about the readings' noise, √R = 1e-12 s
    end

    @testset "Stein Kalman (KAS-2)" begin
        # Stein 2003 eqs 4.1–7.3: weights ∝ 1/variance of the matching noise (Levine 2012 §XVI).
        st = timescale(cc; algorithm=:stein, models=ms)
        @test st.algorithm === :stein
        @test vec(mean(st.weights[late, :]; dims=1)) ≈ wopt rtol=1e-12
        @test adev1(st) ≈ sqrt(sum(wopt .^ 2 .* q)) rtol=0.05
        # The same optimum as Greenhall's reduced Kalman scale, to a few percent at every τ.
        τs = [1, 10, 100]
        gr = timescale(cc; algorithm=:greenhall, models=ms)
        @test adev(timescale_phase(st), τs; ci=false).dev ≈ adev(timescale_phase(gr), τs; ci=false).dev rtol=0.05
        # Clocks differing in white FM and random-walk FM: quieter than each clock, as in the paper's figs 2–4.
        mix = [ClockModel(; tau0=1.0, q1=isodd(i) ? 1e-22 : 1e-24, q2=isodd(i) ? 1e-32 : 1e-28) for i in 1:4]
        cm, _ = simulate_ensemble(mix, 20000; reference=0, rng=Xoshiro(3))
        best = minimum(reduce(hcat, (adev(cm[i], τs; ci=false).dev for i in 1:4)); dims=2)
        sm = timescale(cm; algorithm=:stein, models=mix)
        @test all(adev(timescale_phase(sm), τs; ci=false).dev .< vec(best))
        @test vec(sm.weights[end, :]) ≈ [1, 100, 1, 100] ./ 202 rtol=1e-12
        # The pivot is a choice of reference: another pivot gives the same weights and a scale as stable.
        s3 = timescale(cc; algorithm=:stein, models=ms, pivot=3)
        @test s3.weights[late, :] ≈ st.weights[late, :]
        @test adev(timescale_phase(s3), τs; ci=false).dev ≈ adev(timescale_phase(st), τs; ci=false).dev rtol=0.02
        # wmax cuts the weights as in KAS-2; a gap in one clock is carried through.
        sc = timescale(cc; algorithm=:stein, models=ms, wmax=0.3)
        @test maximum(sc.weights[late, :]) ≈ 0.3 && all(sum(sc.weights; dims=2) .≈ 1)
        xg = copy(cc.x); xg[30000:30100, 2] .= NaN
        sg = timescale_phase(timescale(ClockComparisons(xg, 1.0); algorithm=:stein, models=ms))
        @test maximum(abs, sg.x[29000:31000] .- timescale_phase(st).x[29000:31000]) < 5e-11
        @test_throws ArgumentError timescale(cc; algorithm=:stein)
        @test_throws ArgumentError timescale(cc; algorithm=:stein, models=ms, pivot=9)
        @test_throws ArgumentError timescale(cc; algorithm=:stein, models=ms, wmax=0.2)
    end

    @test_throws ArgumentError timescale(cc; algorithm=:at3)
    @test_throws ArgumentError timescale(cc; wmax=0.2)
end

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
    @test at.xhat === at.x                                  # no filter state apart from the scale's own phases
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
        # Out of the scale the clock is still estimated (Weiss & Weissert 1991, §4),
        # so its frequency holds the step 20 samples on (0.96 ± 0.08 of it over 200 seeds).
        @test on.y[4020, 2] - on.y[3990, 2] ≈ 5e-11 rtol=0.45
        # :real keeps what was issued: identical up to the catch, and the scale
        # stays continuous where the rerun corrects it.
        re = timescale(stepped; algorithm=:at2, steps=:real)
        @test re.weights[1:4000, :] == on.weights[1:4000, :] && re.x[1:4000, :] == off.x[1:4000, :]
        @test all(iszero, re.weights[4010:4020, 2]) && all(isfinite, timescale_phase(re).x)
        @test re.y[4020, 2] - re.y[3990, 2] ≈ 5e-11 rtol=0.45
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
        # Readings at every fourth sample: the interval of the equations is the
        # time between readings (Weiss & Weissert 1991, §2.8), so the scale is
        # that of the readings alone, the caught step included.
        x4 = fill(NaN, 32000, 4); x4[1:4:end, :] = sx
        for kw in ((; algorithm=:ewfa), (; algorithm=:at1), (; algorithm=:at2, steps=:real))
            @test timescale(ClockComparisons(x4, 0.25); kw...).x[1:4:end, :] ≈ timescale(stepped; kw...).x rtol=1e-9
        end
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

    @testset "a clock back from a gap joins anew" begin
        # No paper gives a rule. The clock has no weight at its first reading
        # after the gap and the scale does not step there.
        x = cc.x[1:8000, :]; x[3001:5000, 1] .= NaN
        cg = ClockComparisons(x, 1.0)
        function rejoin(; kw...)                            # step at the return over the scale's rms step
            r = timescale(cg; kw...)
            e = timescale_phase(r).x
            @test r.weights[5001, 1] == 0 && r.weights[5002, 1] > 0
            abs(e[5001] - e[5000]) / sqrt(mean(abs2, diff(e[1:3000])))
        end
        @test rejoin(; algorithm=:ewfa) < 8
        @test rejoin(; algorithm=:at1) < 8
        @test rejoin(; algorithm=:at2, steps=true) < 8
        @test rejoin(; algorithm=:stein, models=ms, pivot=2) < 8
    end

    @testset "a clock that joins late starts from its own readings" begin
        # Weiss & Weissert 1991, §2.7: ε² starts from the clock's own
        # characterization. The best clock, with a rate of 10⁻¹¹, joins at
        # sample 20001: it has its weight at once (1.02 of the optimum, 0.99 to
        # 1.06 over 40 seeds) and the scale keeps its rate (within 4×10⁻¹⁴).
        x = copy(cc.x); x[:, 1] .+= (0:49999) .* 1e-11; x[1:20000, 1] .= NaN
        lj = ClockComparisons(x, 1.0)
        @test mean(timescale(lj).weights[20002:22000, 1]) ≈ wopt[1] rtol=0.1
        e = timescale_phase(timescale(lj; algorithm=:ewfa)).x
        @test abs((e[end] - e[20001]) / 29999 - (e[20000] - e[1]) / 19999) < 2e-13
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
        # Two noiseless clocks and one with reading noise: the weights stay ∝ 1/q1,
        # not the gain's row mean (1/3 each for the noisy one).
        mx = [ClockModel(; tau0=1.0, R=r, q1=q, q2=1e-30) for (r, q) in zip([0.0, 0.0, 1e-25], [1e-24, 4e-24, 2.5e-25])]
        cx, _ = simulate_ensemble(mx, 4000; reference=0, rng=Xoshiro(5))
        @test vec(mean(timescale(cx; algorithm=:greenhall, models=mx).weights[2000:end, :]; dims=1)) ≈ [0.190, 0.048, 0.762] atol=0.01
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
        # Greenhall 2006 §7: clocks 1 and 2 read without noise, clock 3 with noise of
        # variance 10. The implicit mean is Σ aᵢ(xᵢ − x̂ᵢ) (3.13). The corrected clocks
        # stay about it with the variances printed in §9.1, and its ADEV is the theory
        # of §9.2 (the paper's Fig. 3), evaluated here at τ = 1 and 8.
        g7 = [ClockModel(; tau0=1.0, R, q1, q2) for (R, q1, q2) in ((0.0, 4.0, 0.37), (0.0, 25.0, 0.031), (10.0, 4.0, 0.37))]
        c7, t7 = simulate_ensemble(g7, 20000; reference=0, rng=MersenneTwister(9))
        r7 = timescale(c7; algorithm=:greenhall, models=g7)
        iem = implicit_mean(r7, t7)
        @test vec(mean(abs2, t7[:, 1, :] .- r7.xhat .- iem.x; dims=1)) ≈ [1.101, 1.101, 2.31] rtol=0.08
        @test adev(iem, [1, 8]; ci=false).dev ≈ [1.405, 0.844] rtol=0.05
        @test r7.xhat[:, 1] ≈ r7.x[:, 1]                       # the pivot is read at its own estimate
        # Without reading noise the corrected clocks agree: the readout is the mean, reversed.
        g0 = [ClockModel(; tau0=1.0, q1=m.q1, q2=m.q2) for m in g7]
        c0, t0 = simulate_ensemble(g0, 5000; reference=0, rng=MersenneTwister(9))
        r0 = timescale(c0; algorithm=:greenhall, models=g0)
        @test implicit_mean(r0, t0).x ≈ -timescale_phase(r0).x
        @test r0.xhat ≈ r0.x
        # A clock without a reading has no weight (§10): the readout stays the mean
        # through the gap, and the scale does not step when the clock returns.
        xg = copy(c0.x); xg[1001:3000, 3] .= NaN
        rg = timescale(ClockComparisons(xg, 1.0); algorithm=:greenhall, models=g0)
        eg = timescale_phase(rg).x
        @test all(iszero, rg.weights[1001:3000, 3])
        @test implicit_mean(rg, t0).x ≈ -eg
        @test abs(eg[3001] - eg[3000]) < 12 * sqrt(mean(abs2, diff(eg[1:1000])))
        @test_throws ArgumentError implicit_mean(timescale(cc), t0)
        a1 = timescale_phase(timescale(cm; algorithm=:greenhall, models=mix)).x
        a2 = timescale_phase(timescale(cm; algorithm=:greenhall, models=mix, pivot=2)).x
        @test 1e-13 < std(a1 .- a2) < 3e-12                    # about the readings' noise, √R = 1e-12 s
    end

    @testset "Stein Kalman (KAS-2)" begin
        # Stein 2003 eqs 4.1–7.3: weights ∝ 1/variance of the matching noise (Levine 2012 §XVI).
        st = timescale(cc; algorithm=:stein, models=ms)
        @test st.algorithm === :stein
        @test st.xhat[2:end, :] ≈ st.x[2:end, :]              # noiseless readings: the estimate is the readout, past the start
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
        @test_throws "no noise" timescale(cc; algorithm=:stein, models=[ClockModel(; tau0=1.0), ClockModel(; tau0=1.0), ms[3]])
    end

    @testset "flicker states in the Kalman scales" begin
        # Flicker FM: the clocks' frequencies follow the truth through the
        # Markov states, the reduction leaves the unreduced filter's
        # frequencies as they are, and Stein's three-state clocks refuse.
        mf = [ClockModel(; tau0=1.0, q1=1e-26, flicker_fm=1e-26, tau_max=1e3),
              ClockModel(; tau0=1.0, q1=1e-26, flicker_fm=4e-26, tau_max=1e3), ClockModel(; tau0=1.0, q1=1e-24)]
        cf, tf = simulate_ensemble(mf, 8192; reference=0, rng=MersenneTwister(11))
        gf = timescale(cf; algorithm=:greenhall, models=mf)
        rest = 2000:8192
        @test cor(gf.y[rest, 1] .- gf.y[rest, 2], tf[rest, 2, 1] .- tf[rest, 2, 2]) > 0.75
        @test maximum(abs, timescale(cf; algorithm=:jones_tryon, models=mf).y .- gf.y) < 1e-6 * std(gf.y[rest, 1])
        @test_throws ArgumentError timescale(cf; algorithm=:stein, models=mf)
        # Flicker PM on two of three equal clocks: the weights stay equal,
        # and past the shortest τ the scale is quieter than with that
        # flicker folded into R.
        hp = [0.0, 2.2e-22, 2.2e-21]
        mp = [ClockModel(; tau0=1.0, q1=1e-24, q2=1e-30, flicker_pm=h) for h in hp]
        fold = [clock_model(Dict(1 => h, 0 => 2e-24, -2 => 1e-30 / 2π^2); tau0=1.0) for h in hp]
        cp, _ = simulate_ensemble(mp, 16384; reference=0, rng=MersenneTwister(31))
        gp = timescale(cp; algorithm=:greenhall, models=mp)
        scale(r) = adev(timescale_phase(r), [16, 32, 64]; ci=false).dev
        @test gp.weights[end, :] ≈ fill(1 / 3, 3) rtol=0.01
        @test all(scale(gp) .< 0.95 .* scale(timescale(cp; algorithm=:greenhall, models=fold)))
        @test maximum(abs, timescale(cp; algorithm=:jones_tryon, models=mp).y .- gp.y) < 1e-6 * std(gp.y[rest, 1])
        # Two noiseless clocks and one whose readings carry flicker PM states: those
        # readings are not noiseless, so only the first two share a weight (∝ 1/q1).
        mq = [ClockModel(; tau0=1.0, q1=q, q2=1e-30, flicker_pm=i == 3 ? 2.2e-20 : 0.0) for (i, q) in enumerate([1e-24, 4e-24, 2.5e-25])]
        cq, _ = simulate_ensemble(mq, 8192; reference=0, rng=Xoshiro(5))
        wq = vec(mean(timescale(cq; algorithm=:greenhall, models=mq).weights[4000:end, :]; dims=1))
        @test wq ≈ [0.407, 0.106, 0.488] atol=0.02
        @test wq[1] / wq[2] ≈ 4 rtol=0.05
    end

    @test_throws ArgumentError timescale(cc; algorithm=:at3)
    @test_throws ArgumentError timescale(cc; wmax=0.2)
end

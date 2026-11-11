# holdover.jl — the measured holdover test against the closed-form error
# variances on a simulated clock, an exact parabola with gaps, the time to a
# limit, and the stretch that grows with τ against TDEV and HTDEV.

using Test
using SigmaTau
using Random

@testset "Holdover test" begin
    m = ClockModel(; tau0=1.0, q1=1e-22, q2=1e-26)
    p, _ = simulate_clock(m, 100_000; rng=MersenneTwister(7))
    T, τ = 256.0, [64.0, 256.0, 1024.0]
    r = τ ./ T
    h = holdover(p; train=T, tau=τ)

    @testset "measured against predicted" begin
        @test h.rms ≈ holdover(m; train=T, tau=τ).sigma rtol=0.1
        @test h.rms ≈ sqrt.(m.q1 .* τ .* (1 .+ 6r ./ 5) .+ m.q2 .* (τ .^ 3 ./ 3 .+ 13T .* τ .^ 2 ./ 35)) rtol=0.1
        @test h.quantile ≈ 1.96 .* h.rms rtol=0.1          # Gaussian errors
        @test all(h.worst .> h.quantile) && h.n == 100_000 - 256 + 1 .- Int.(τ)
        @test holdover(p; train=T, tau=τ, anchor=false, drift=true).rms ≈
              holdover(m; train=T, tau=τ, anchor=false, drift=true).sigma rtol=0.1
        # The noises no clock model holds.
        hα = Dict(2 => 1e-20, 1 => 1e-21, -1 => 1e-24)
        f = noise_gen(PhaseData, 100_000, 1.0; h=hα, rng=MersenneTwister(2))
        @test holdover(f; train=T, tau=τ).rms ≈ holdover(hα; tau0=1.0, train=T, tau=τ).sigma rtol=0.1
    end

    @testset "predicted against the closed forms" begin
        Tm, t = 2048.0, [128.0, 2048.0, 8192.0]
        r = t ./ Tm
        ln = log.(r ./ (1 .+ r))
        unit = Dict(0 => 1.0, -1 => 1.0, -2 => 1.0)
        # Delporte et al. (2000) with k_α₋₂ = h_α/4π²: eqs 17–19 for a line and
        # 11–13 for a parabola, on white, flicker and random-walk FM.
        line = holdover(unit; tau0=1.0, train=Tm, tau=t, anchor=false)
        @test all(isapprox.(line.budget, hcat(
            Tm / 15 .* (9r .^ 2 .+ 9r .+ 1),
            Tm^2 / 12 .* (12r .^ 4 .+ 24r .^ 3 .+ 20r .^ 2 .+ 8r .+ 1 .+ 2 .* log1p.(r) .* (6r .^ 2 .+ 6r .+ 1) .+
                          2r .^ 3 .* ln .* (6r .^ 2 .+ 15r .+ 8)),
            2π^2 * Tm^3 / 105 .* (35r .^ 3 .+ 39r .^ 2 .+ 11r .+ 1)); rtol=4e-3))
        par = holdover(unit; tau0=1.0, train=Tm, tau=t, anchor=false, drift=true)
        @test all(isapprox.(par.budget, hcat(
            3Tm / 70 .* (50r .^ 4 .+ 100r .^ 3 .+ 69r .^ 2 .+ 19r .+ 1),
            Tm^2 / 32 .* (192r .^ 6 .+ 576r .^ 5 .+ 692r .^ 4 .+ 424r .^ 3 .+ 136r .^ 2 .+ 20r .+ 1 .+
                          96r .^ 3 .* ln .* (2r .^ 4 .+ 7r .^ 3 .+ 9r .^ 2 .+ 5r .+ 1)),
            π^2 * Tm^3 / 630 .* (450r .^ 4 .+ 690r .^ 3 .+ 303r .^ 2 .+ 42r .+ 2)); rtol=4e-3))
        # From the last reading: the clock's own noise over τ plus the learned frequency's.
        cm = ClockModel(; tau0=1.0, R=1.0, q1=1.0, q2=1.0)
        a = holdover(cm; train=Tm, tau=t)
        @test a.alphas == [2, 0, -2] && a.sigma .^ 2 ≈ vec(sum(a.budget; dims=2))
        @test all(isapprox.(a.budget, hcat(2 .+ 12 .* r ./ Tm .* (1 .+ r), t .* (1 .+ 6r ./ 5),
                                            t .^ 3 ./ 3 .+ 13Tm .* t .^ 2 ./ 35); rtol=2e-3))
        @test a.quantile ≈ 1.96 .* a.sigma rtol=1e-4
        @test holdover(cm; train=Tm, tau=t, limit=a.quantile[2]).time_to_limit ≈ t[2]
        @test_throws ArgumentError holdover(Dict(-3 => 1.0); tau0=1.0, train=Tm)
        @test_throws ArgumentError holdover(unit; train=Tm)
    end

    @testset "Kalman filter as the predictor" begin
        k = holdover(p, m; train=2000, tau=τ)
        @test k.rms ≈ k.claimed rtol=0.05                  # the filter's claim holds on its own model
        @test all(k.rms .< h.rms)
        @test k.rms ≈ holdover(p; train=180, tau=τ).rms rtol=0.1    # the line at its best stretch, √(42q1/13q2)
        @test k.claimed .^ 2 ≈ holdover_error(kalman_filter(p, m), τ) .^ 2 .+ m.R rtol=1e-6
        @test k.n == 100_000 - 2000 + 1 .- Int.(τ) && !k.anchor
        low = holdover(p, ClockModel(; tau0=1.0, q1=m.q1 / 100, q2=m.q2 / 100); train=2000, tau=τ)
        @test all(low.rms .> 5 .* low.claimed)             # a model with too little noise claims too little
        g = copy(p.x); g[[5000, 5001, 70_000]] .= NaN
        @test holdover(PhaseData(g, 1.0), m; train=2000, tau=τ).n == k.n .- 6          # less the starts and the readings at a gap
        @test holdover(p, m).train == 25_000 && isempty(h.claimed)
        # Flicker FM and PM as Markov states: the claim holds there too.
        fm = ClockModel(; tau0=1.0, R=1e-24, q1=1e-24, flicker_fm=1e-24, flicker_pm=1e-22)
        f = holdover(simulate_clock(fm, 100_000; rng=MersenneTwister(1))[1], fm; train=2000, tau=τ)
        @test f.rms ≈ f.claimed rtol=0.05
        # A flicker FM clock, its model fitted to the first stretch: the flicker has to stay in.
        ff = noise_gen(PhaseData, 8192, 1.0; h=Dict(-1 => 3e-25), rng=MersenneTwister(5))
        fk = holdover(ff, clock_model(powerlaw_fit(PhaseData(ff.x[1:2048], 1.0))); train=2048, tau=[1, 4])
        @test fk.rms ≈ fk.claimed rtol=0.15
        @test_throws ArgumentError holdover(p, m; train=2000, tau=[1e6])
        @test_throws ArgumentError holdover(p, ClockModel(; tau0=1.0))      # no noise: nothing to weigh a reading with
    end

    @testset "drift and gaps" begin
        D, n = 1e-15, 100
        x = 1e-6 .+ 3e-10 .* (0:999) .+ D / 2 .* (0:999) .^ 2
        t = [10.0, 100.0]
        # A line's frequency is that of the stretch's middle, (n − 1)/2 samples old.
        @test holdover(PhaseData(x, 1.0); train=n, tau=t).mean ≈ D / 2 .* t .* (t .+ (n - 1)) rtol=1e-9
        x[[150, 151, 400]] .= NaN
        g = holdover(PhaseData(x, 1.0); train=n, tau=t, drift=true)
        @test maximum(g.worst) < 1e-20 && g.n == [885, 797]    # less the starts and readings at a gap
    end

    @testset "time to a limit and defaults" begin
        at(limit) = holdover(p; train=T, tau=τ, limit).time_to_limit
        @test at(sqrt(h.quantile[2] * h.quantile[3])) ≈ sqrt(τ[2] * τ[3])
        @test (at(h.quantile[1] / 2), at(2h.quantile[3])) == (0.0, Inf) && isnan(h.time_to_limit)
        d = holdover(p)
        @test d.train == 25_000 && d.tau == 2.0 .^ (0:14) && d.anchor && !d.drift
        @test occursin("trained on 2.5e+04 s", sprint(show, MIME"text/plain"(), d))
    end

    @testset "a stretch that grows with tau" begin
        t = [1.0, 64.0, 1024.0]
        g = holdover(p; grow=2, tau=t)
        @test g.rms[2:3] ≈ holdover(m; grow=2, tau=t[2:3]).sigma rtol=0.1
        # At τ₀ the error is the second difference TVAR is defined on; with a
        # parabola and three times τ, the third difference of HTVAR.
        @test g.rms[1] ≈ sqrt(6) * tdev(p, [1]; ci=false).dev[1]
        @test holdover(p; grow=3, tau=t, drift=true).rms[1] ≈ sqrt(20) * htdev(p, [1]; ci=false).dev[1]
        # Above τ₀, a constant over √6 TDEV: √(8/5) on white FM, √(452/231) on random-walk FM.
        σ(α) = holdover(Dict(α => 1.0); tau0=1.0, grow=2, tau=[4096.0]).sigma[1]
        @test σ(0) ≈ sqrt(8 / 5 * 4096 / 2) rtol=1e-3
        @test σ(-2) ≈ sqrt(452 / 231 * 11π^2 * 4096.0^3 / 10) rtol=1e-3
        d = holdover(p; grow=2)
        @test isnan(d.train) && d.grow == 2 && d.tau == 2.0 .^ (0:14) && occursin("trained on 2 τ", sprint(show, d))
        @test_throws ArgumentError holdover(p; grow=2, train=T)
        @test_throws ArgumentError holdover(p; grow=2, drift=true)
        @test_throws ArgumentError holdover(m; grow=2)
    end

    @test_throws ArgumentError holdover(p; train=0.5)
    @test_throws ArgumentError holdover(p; train=T, tau=[1e6])
    @test_throws ArgumentError holdover(p; confidence=1)
end

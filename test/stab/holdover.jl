# holdover.jl — the measured holdover test against the closed-form error
# variances on a simulated clock, an exact parabola with gaps, and the time
# to a limit.

using Test
using SigmaTau
using Random

@testset "Holdover test" begin
    m = ClockModel(; tau0=1.0, q1=1e-22, q2=1e-26)
    p, _ = simulate_clock(m, 100_000; rng=MersenneTwister(7))
    T, τ = 256.0, [64.0, 256.0, 1024.0]
    r = τ ./ T
    h = holdover(p; train=T, tau=τ)

    @testset "white and random-walk FM variances" begin
        @test h.rms ≈ sqrt.(m.q1 .* τ .* (1 .+ 6r ./ 5) .+ m.q2 .* (τ .^ 3 ./ 3 .+ 13T .* τ .^ 2 ./ 35)) rtol=0.1
        @test h.quantile ≈ 1.96 .* h.rms rtol=0.1          # Gaussian errors
        @test all(h.worst .> h.quantile) && h.n == 100_000 - 256 + 1 .- Int.(τ)
        # Delporte et al. (2000) with k₋₂ = q1/2π² and k₋₄ = q2/8π⁴: eqs 17 and
        # 19 for a line, 11 and 13 for a parabola.
        line = holdover(p; train=T, tau=τ, anchor=false)
        @test line.rms ≈ sqrt.(2m.q1 * T / 15 .* (9r .^ 2 .+ 9r .+ 1) .+
                               m.q2 * T^3 / 105 .* (35r .^ 3 .+ 39r .^ 2 .+ 11r .+ 1)) rtol=0.1
        par = holdover(p; train=T, tau=τ, anchor=false, drift=true)
        @test par.rms ≈ sqrt.(3m.q1 * T / 35 .* (50r .^ 4 .+ 100r .^ 3 .+ 69r .^ 2 .+ 19r .+ 1) .+
                              m.q2 * T^3 / 1260 .* (450r .^ 4 .+ 690r .^ 3 .+ 303r .^ 2 .+ 42r .+ 2)) rtol=0.1
        # Their eq. 18, a line on flicker FM, with k₋₃ = h₋₁/4π².
        f = noise_gen(PhaseData, 100_000, 1.0; h=Dict(-1 => 1e-24), rng=MersenneTwister(2))
        @test holdover(f; train=T, tau=τ, anchor=false).rms ≈ sqrt.(1e-24 * T^2 / 12 .*
            (12r .^ 4 .+ 24r .^ 3 .+ 20r .^ 2 .+ 8r .+ 1 .+ 2 .* log1p.(r) .* (6r .^ 2 .+ 6r .+ 1) .+
             2r .^ 3 .* log.(r ./ (1 .+ r)) .* (6r .^ 2 .+ 15r .+ 8))) rtol=0.1
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

    @test_throws ArgumentError holdover(p; train=0.5)
    @test_throws ArgumentError holdover(p; train=T, tau=[1e6])
    @test_throws ArgumentError holdover(p; confidence=1)
end

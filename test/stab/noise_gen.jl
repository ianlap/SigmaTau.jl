# test/stab/noise_gen.jl — the calibrated power-law noise generator.

using Test
using Random
using SigmaTau

@testset "noise_gen — calibrated power-law noise generator" begin
    @testset "single-α σ_y(τ=τ₀) target hit exactly per realization" begin
        Random.seed!(42)
        N    = 4096
        tau0 = 1.0
        for (α, σ) in ((-2, 1e-13), (-1, 5e-13), (0, 1e-12),
                       (1, 2e-12), (2, 3e-12))
            p = noise_gen(PhaseData, N, tau0; sigma1 = Dict(α => σ))
            @test p isa PhaseData
            @test length(p.x) == N
            σ_emp = _adev_core(p.x, [1], tau0)[1]
            @test isapprox(σ_emp, σ; rtol = 1e-10)
        end
    end

    @testset "FrequencyData path matches PhaseData via cumsum" begin
        Random.seed!(7)
        f = noise_gen(FrequencyData, 1024, 1.0; sigma1 = Dict(0 => 1e-12))
        Random.seed!(7)
        p = noise_gen(PhaseData, 1024, 1.0; sigma1 = Dict(0 => 1e-12))
        @test f isa FrequencyData
        @test length(f.y) == 1024
        @test cumsum(f.y) .* f.tau0 ≈ p.x
    end

    @testset "composite mixture: independent components sum in quadrature" begin
        # Two comparable-amplitude noises (WFM + FFM) at the same σ_y(τ₀).
        # Single-α calibration pegs each component to exactly σ at τ=τ₀;
        # being statistically independent, their sum has variance ≈ 2σ²
        # at τ=τ₀ (the cross-term is zero in expectation; sample
        # covariance fluctuates around zero at scale 1/√N).
        Random.seed!(99)
        σ = 1e-12
        N = 8192
        p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(0 => σ, -1 => σ))
        @test p isa PhaseData
        @test length(p.x) == N
        σ_emp = _adev_core(p.x, [1], 1.0)[1]
        # Allow ±20% — the independent-components cross-term has a
        # sample-size-limited spread; well within tolerance at N=8192.
        @test isapprox(σ_emp, sqrt(2) * σ; rtol = 0.2)
    end

    @testset "h-mode: the generated spectrum carries the requested coefficient" begin
        using Statistics: var, median
        # An `h` component is scaled analytically on the Kasdin–Walter
        # filter's low-frequency PSD, 2τ₀(2πτ₀f)^α per unit innovation.
        # Exact discrete checks for the three integer-order processes
        # (sample variances, so a few percent of spread at N = 65536):
        τ0 = 1.0; h = 1e-24; N = 65536
        y_rw = noise_gen(FrequencyData, N, τ0; h = Dict(-2 => h), rng = Xoshiro(3)).y
        @test isapprox(var(diff(y_rw)), 2π^2 * τ0 * h; rtol = 0.05)      # random-walk innovations
        x_wp = noise_gen(PhaseData, N, τ0; h = Dict(2 => h), rng = Xoshiro(3)).x
        @test isapprox(var(x_wp), h / (8π^2 * τ0); rtol = 0.05)           # white phase, S_x = h/(4π²)
        y_wf = noise_gen(FrequencyData, N, τ0; h = Dict(0 => h), rng = Xoshiro(3)).y
        @test isapprox(var(y_wf), h / (2τ0); rtol = 0.05)                 # white frequency
        # Welch estimate of S_y(f)/(h f^α) well below Nyquist, every α.
        for α in (2, 1, 0, -1, -2)
            p  = noise_gen(FrequencyData, N, τ0; h = Dict(α => h), rng = Xoshiro(7))
            sy = psd_y(p)
            keep = (sy.freq .> 0.003 / τ0) .& (sy.freq .< 0.03 / τ0)
            ratio = median(sy.psd[keep] ./ (h .* sy.freq[keep] .^ α))
            @test isapprox(ratio, 1.0; rtol = 0.15)
        end
    end

    @testset "h-mode RWFM at τ₀ ≠ 1: discrete random walk, σ²(τ₀) = π² h_{-2} τ₀" begin
        # For the sampled random walk the innovation variance is q = 2π²τ₀h
        # and ADEV²(τ₀) = q/2 = π²τ₀h, not the continuous-time (2π²/3)hτ₀
        # (the latter holds for τ ≫ τ₀). Sample estimate at N = 65536.
        using Statistics: var
        tau0 = 0.5
        h_m2 = 1e-26
        p = noise_gen(FrequencyData, 65536, tau0; h = Dict(-2 => h_m2), rng = Xoshiro(13))
        @test isapprox(var(diff(p.y)), 2π^2 * tau0 * h_m2; rtol = 0.05)
        x = cumsum(p.y) .* tau0
        @test isapprox(_adev_core(x, [1], tau0)[1], sqrt(π^2 * tau0 * h_m2); rtol = 0.05)
    end

    @testset "determinism: seed → same output" begin
        Random.seed!(123)
        a = noise_gen(PhaseData, 512, 1.0; sigma1 = Dict(0 => 1e-12))
        Random.seed!(123)
        b = noise_gen(PhaseData, 512, 1.0; sigma1 = Dict(0 => 1e-12))
        @test a.x == b.x
    end

    @testset "rng kwarg: independent reproducible streams" begin
        # Equal-seed streams reproduce; different seeds differ; the rng
        # path does not disturb the global RNG (so parallel MC realizations
        # can each carry their own stream).
        a = noise_gen(PhaseData, 512, 1.0; sigma1 = Dict(0 => 1e-12), rng = Xoshiro(123))
        b = noise_gen(PhaseData, 512, 1.0; sigma1 = Dict(0 => 1e-12), rng = Xoshiro(123))
        c = noise_gen(PhaseData, 512, 1.0; sigma1 = Dict(0 => 1e-12), rng = Xoshiro(999))
        @test a.x == b.x
        @test a.x != c.x
        y1 = SigmaTau._gen_powerlaw_y(0.0, 256; rng = Xoshiro(7))
        y2 = SigmaTau._gen_powerlaw_y(0.0, 256; rng = Xoshiro(7))
        @test y1 == y2
    end

    @testset "argument validation" begin
        @test_throws ArgumentError noise_gen(PhaseData, 4096, 1.0)
        @test_throws ArgumentError noise_gen(PhaseData, 4096, 1.0;
            sigma1 = Dict(0 => 1e-12), h = Dict(0 => 1e-24))
        @test_throws ArgumentError noise_gen(PhaseData, 3, 1.0;
            sigma1 = Dict(0 => 1e-12))
        @test_throws ArgumentError noise_gen(PhaseData, 256, -1.0;
            sigma1 = Dict(0 => 1e-12))
        @test_throws ArgumentError noise_gen(PhaseData, 256, 1.0;
            sigma1 = Dict(3 => 1e-12))
        @test_throws ArgumentError noise_gen(PhaseData, 256, 1.0;
            sigma1 = Dict(0 => -1e-12))
        @test_throws ArgumentError noise_gen(PhaseData, 256, 1.0;
            h = Dict(0 => -1e-24))
    end

    @testset "zero amplitude → zero-vector contribution" begin
        Random.seed!(0)
        p = noise_gen(PhaseData, 256, 1.0; sigma1 = Dict(0 => 0.0))
        @test all(iszero, p.x)
    end
end

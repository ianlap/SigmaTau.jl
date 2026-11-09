# powerlaw.jl — the power-law model against SP1065's closed forms, and fits
# that recover the h_α of noise_gen records.

using Test
using SigmaTau
using Random

@testset "power-law model" begin
    @testset "σ(τ) matches the closed forms" begin
        σ²(dev, α, τ) = powerlaw_sigma(Dict(α => 1.0), [τ]; deviation=dev, tau0=1.0)[1]^2
        fh = 0.5
        @test σ²(:adev, 2, 1.0)    ≈ 3fh / (4π^2)              rtol=1e-12   # white PM, every m
        @test σ²(:adev, 2, 37.0)   ≈ 3fh / (4π^2 * 37^2)       rtol=1e-12
        @test σ²(:adev, 0, 1e3)    ≈ 1 / 2e3                   rtol=1e-6
        @test σ²(:adev, -1, 1e3)   ≈ 2log(2)                   rtol=1e-5
        @test σ²(:adev, -2, 1e3)   ≈ 2π^2 / 3 * 1e3            rtol=1e-5
        @test σ²(:mdev, 0, 1e3)    ≈ 1 / 4e3                   rtol=1e-5
        @test σ²(:mdev, -2, 1e3)   ≈ 11π^2 / 20 * 1e3          rtol=1e-5
        @test σ²(:hdev, 0, 1e3)    ≈ 1 / 2e3                   rtol=1e-6
        @test σ²(:hdev, 0, 1e6)    ≈ 1 / 2e6                   rtol=1e-6    # past the resolved lobes
        @test σ²(:mdev, -1, 1.0)   ≈ σ²(:adev, -1, 1.0)        rtol=1e-12   # m = 1: modified = plain
        @test σ²(:mhdev, -4, 1.0)  ≈ σ²(:hdev, -4, 1.0)        rtol=1e-12
    end

    @testset "spectra" begin
        h = Dict(:WHPM => 1e-22, 0 => 1e-24, -2 => 1e-30)
        f = [1e-4, 1e-2, 1.0, 1e2]
        Sy = 1e-22 .* f .^ 2 .+ 1e-24 .+ 1e-30 ./ f .^ 2
        @test powerlaw_psd(h, f) ≈ Sy
        @test powerlaw_psd(h, f; quantity=:Sx) ≈ Sy ./ (2π .* f) .^ 2
        L = powerlaw_psd(h, f; quantity=:L, f_carrier=10e6)
        @test L ≈ 10 .* log10.(0.5 .* (2π * 10e6)^2 .* Sy ./ (2π .* f) .^ 2)
        fit = powerlaw_fit(f, L; quantity=:L, f_carrier=10e6, alphas=[2, 0, -2])
        @test fit.h ≈ [1e-22, 1e-24, 1e-30] rtol=1e-8
    end

    @testset "fits recover noise_gen's h_α" begin
        h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30)
        p = noise_gen(PhaseData, 65536, 1.0; h=h, rng=MersenneTwister(1))
        for f in (powerlaw_fit(p), powerlaw_fit(p; deviation=:adev), powerlaw_fit(psd_y(p; nperseg=4096)))
            @test f.alphas == [2, 1, 0, -1, -2]
            @test f.h[1] ≈ 1e-22 rtol=0.05
            @test f.h[3] ≈ 1e-24 rtol=0.1
            @test f.h[5] ≈ 1e-30 rtol=0.5
            @test f.h[2] * 0.5 < 0.05 * f.h[1] * 0.5^2     # spurious flicker PM is small at f_h
        end
        f = powerlaw_fit(p)
        @test powerlaw_sigma(f, [1.0, 100.0]) ≈ adev(p, [1, 100]; ci=false).dev rtol=0.02
        @test powerlaw_fit(FrequencyData(diff(p.x), 1.0); deviation=:hdev, alphas=0:-1:-4).alphas == [0, -1, -2, -3, -4]
    end

    @testset "bad input" begin
        p = noise_gen(PhaseData, 1024, 1.0; h=Dict(0 => 1e-24), rng=MersenneTwister(2))
        @test_throws ArgumentError powerlaw_fit(p; deviation=:totdev)
        @test_throws ArgumentError powerlaw_fit(p; deviation=:mdev, alphas=[0, -3])
        @test_throws ArgumentError powerlaw_sigma(Dict(0 => 1e-24), [1.5]; tau0=1.0)
        @test_throws ArgumentError powerlaw_sigma(Dict(0 => 1e-24), [1.0])
        @test_throws ArgumentError powerlaw_psd(Dict(0 => 1e-24), [1.0]; quantity=:L)
    end
end

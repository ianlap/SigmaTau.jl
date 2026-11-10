# clockmodel.jl — the clock model's Φ and Q in closed form, its mapping to
# the h_α, a model recovered from a noise_gen record's power-law fit, and the
# Markov states of the flicker FM and flicker PM.

using Test
using SigmaTau
using Random
using LinearAlgebra

@testset "clock model" begin
    m3 = ClockModel(; tau0=1.0, R=1e-20, q1=3e-24, q2=2e-30, q3=5e-40)
    m2 = ClockModel(; tau0=1.0, q1=3e-24, q2=2e-30)

    @testset "Φ and Q" begin
        @test m3.drift && !m2.drift
        @test state_transition(m3, 2.0) == [1 2 2; 0 1 2; 0 0 1]
        @test state_transition(m2, 2.0) == [1 2; 0 1]
        Q = process_noise(m3, 2.0)
        @test Q ≈ [3e-24*2 + 2e-30*8/3 + 5e-40*32/20   2e-30*2 + 5e-40*16/8   5e-40*8/6
                   2e-30*2 + 5e-40*16/8                 2e-30*2 + 5e-40*8/3    5e-40*2
                   5e-40*8/6                            5e-40*2                5e-40*2]
        @test process_noise(m2, 2.0) ≈ process_noise(ClockModel(; tau0=1.0, q1=3e-24, q2=2e-30, drift=true), 2.0)[1:2, 1:2]
        for m in (m2, m3)    # Q composes over steps
            Φ = state_transition(m, 3.0)
            @test process_noise(m, 5.0) ≈ Φ * process_noise(m, 2.0) * Φ' + process_noise(m, 3.0)
        end
    end

    @testset "h_α mapping" begin
        h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30, -4 => 1e-40)
        m = clock_model(h; tau0=0.1)
        @test (m.R, m.q1, m.q2, m.q3) == (1e-22 / (8π^2 * 0.1), 0.5e-24, 2π^2 * 1e-30, 8π^4 * 1e-40)
        @test clock_model(Dict(:WHFM => 1e-24); tau0=1.0).q1 == 0.5e-24
        @test powerlaw_sigma(m, [0.1, 10.0]; deviation=:hdev) ≈ powerlaw_sigma(h, [0.1, 10.0]; deviation=:hdev, tau0=0.1)
        @test powerlaw_psd(m2, [1e-3, 1.0]) ≈ 2 * 3e-24 .+ 2e-30 / (2π^2) ./ [1e-3, 1.0] .^ 2
        # Flicker PM folds into R: the white PM with the PM's Allan variance
        # at the last lag where the PM still leads the FM.
        pm, fm = Dict(2 => 1e-22, 1 => 2e-22), Dict(0 => 2e-24)
        k = SigmaTau._pm_handover(merge(pm, fm), 1.0)
        σ(h, m) = powerlaw_sigma(h, [m]; tau0=1.0)[1]
        @test σ(pm, k) > σ(fm, k) && σ(pm, k + 1) < σ(fm, k + 1)
        mf = clock_model(merge(pm, fm); tau0=1.0)
        @test mf.R ≈ σ(pm, k)^2 * k^2 / 3 && mf.q1 == 1e-24
        @test mf.R > σ(pm, 1)^2 / 3 > 1e-22 / 8π^2
        @test clock_model(Dict(1 => 1e-26, 0 => 2e-24); tau0=1.0).R ≈ σ(Dict(1 => 1e-26), 1)^2 / 3
        run = merge(pm, Dict(-4 => 1e-36))
        kh = SigmaTau._pm_handover(run, 1.0)
        σh(h, m) = powerlaw_sigma(h, [m]; tau0=1.0, deviation=:hdev)[1]
        @test 1 < kh < 2^16 && σh(pm, kh) > σh(Dict(-4 => 1e-36), kh)
        @test clock_model(run; tau0=1.0).R ≈ σ(pm, kh)^2 * kh^2 / 3
    end

    @testset "from a record's fit" begin
        h = Dict(2 => 1e-22, 0 => 1e-24, -2 => 1e-30)
        p = noise_gen(PhaseData, 65536, 1.0; h=h, rng=MersenneTwister(1))
        m = clock_model(powerlaw_fit(p; alphas=[2, 0, -2]))
        @test !m.drift
        @test m.R  ≈ 1e-22 / 8π^2 rtol=0.05
        @test m.q1 ≈ 0.5e-24      rtol=0.1
        @test m.q2 ≈ 2π^2 * 1e-30 rtol=0.5
        @test powerlaw_sigma(m, [1.0, 100.0]) ≈ adev(p, [1, 100]; ci=false).dev rtol=0.05
    end

    @testset "flicker states" begin
        hf, hp = 1.6e-26, 2.2e-22
        m = ClockModel(; tau0=1.0, R=1e-24, q1=1e-24, q2=3.3e-30, q3=1e-42, flicker_fm=hf, flicker_pm=hp)
        # States only where each flicker holds a quarter of the Allan variance,
        # a decade apart; none for a flicker that never does.
        Tf, Tp = SigmaTau._fm_taus(m), SigmaTau._pm_taus(m)
        @test Tf ≈ 56.2 .* 10.0 .^ (0:3) rtol=0.01
        @test Tp == [1.0, 10.0, 100.0]
        @test size(state_transition(m, 1.0)) == (10, 10)
        @test SigmaTau._nstates(ClockModel(; tau0=1.0, q1=1e-24, flicker_fm=1e-30, tau_max=1e3)) == 2
        # Φ and Q are those of dx = (y + Σm) dt, dmⱼ = −mⱼ/Tⱼ dt + √(2σ²/Tⱼ) dW
        # with σ² = h₋₁ ln 10, and the same for the phase states with
        # σ² = h₁ ln 10/4π², by Van Loan's matrix exponential.
        F, W = zeros(10, 10), zeros(10)
        F[1, 2] = F[2, 3] = 1
        W[1:3] = [m.q1 + SigmaTau._fm_q1(m), m.q2, m.q3]
        for (j, T) in enumerate(Tf)
            F[1, 3+j], F[3+j, 3+j], W[3+j] = 1, -1 / T, 2hf * log(10) / T
        end
        for (j, T) in enumerate(Tp)
            F[7+j, 7+j], W[7+j] = -1 / T, 2hp * log(10) / (4π^2 * T)
        end
        E = exp([-F Diagonal(W); zeros(10, 10) F'] .* 7.0)
        Φ = E[11:20, 11:20]'
        @test state_transition(m, 7.0) ≈ Φ rtol=1e-10
        @test process_noise(m, 7.0) ≈ Φ * E[1:10, 11:20] rtol=1e-10
        Φ3 = state_transition(m, 3.0)
        @test process_noise(m, 5.0) ≈ Φ3 * process_noise(m, 2.0) * Φ3' + process_noise(m, 3.0)
        # One frequency state in the paper's own form: Davis, Greenhall &
        # Stacey (2005) eqs 2.7–2.9 with R = 1/T and σ_M² = 2σ²/T.
        one = ClockModel(; tau0=1.0, flicker_fm=hf, fm_band=(40, 40))
        x, σM² = 7 / 40, 2hf * log(10) / 40
        Q = process_noise(one, 7.0)
        @test state_transition(one, 7.0)[1, 3] ≈ (1 - exp(-x)) * 40
        @test Q[1, 1] - 7 * SigmaTau._fm_q1(one) ≈ σM² * 7^3 * (-3 / 2 + x + 2exp(-x) - exp(-2x) / 2) / x^3
        @test Q[1, 3] ≈ σM² * 7^2 * (1 / 2 - exp(-x) + exp(-2x) / 2) / x^2
        @test Q[3, 3] ≈ σM² * 7 * (1 - exp(-2x)) / 2x
        # A slow state loses no digits: its phase variance is (2/3) σ² δ³/T to first order.
        slow = ClockModel(; tau0=1.0, flicker_fm=hf, fm_band=(1e6, 1e6))
        @test process_noise(slow, 1.0)[1, 1] - SigmaTau._fm_q1(slow) ≈ 2 / 3 * hf * log(10) * 1e-6 rtol=1e-6
        # The states alone are flicker: the Allan variance of their sum plus
        # the white top-up, Σ σ² g(τ/Tⱼ)/(τ/Tⱼ)² + q1/τ for the frequency
        # states (the paper's eq. 2.5) and Σ σ² (1 − a)(3 − a)/τ² + 3R/τ²
        # for the phase states, against the power law; exact at τ₀.
        τ = [1.0, 4.0, 16.0, 256.0, 4096.0]
        fm = ClockModel(; tau0=1.0, flicker_fm=hf, fm_band=(2, 2e5))
        avf(τ) = sum(hf * log(10) * SigmaTau._markov_g(τ / T) / (τ / T)^2 for T in SigmaTau._fm_taus(fm)) + SigmaTau._fm_q1(fm) / τ
        @test sqrt.(avf.(τ)) ≈ powerlaw_sigma(fm, τ) rtol=0.03
        @test avf(1.0) ≈ 2hf rtol=1e-12
        pm = ClockModel(; tau0=1.0, flicker_pm=hp, pm_band=(1, 1e4))
        avp(τ) = (sum(SigmaTau._pm_var(pm) * (1 - exp(-τ / T)) * (3 - exp(-τ / T)) for T in SigmaTau._pm_taus(pm)) +
                  3 * SigmaTau._reading_var(pm)) / τ^2
        @test sqrt.(avp.(τ)) ≈ powerlaw_sigma(pm, τ) rtol=0.03
        @test avp(1.0) ≈ powerlaw_sigma(pm, [1.0])[1]^2 rtol=1e-12
        # From a fit: kept only on request.
        hs = Dict(1 => hp, 0 => 2e-24, -1 => hf, -2 => 1e-31)
        off, on = clock_model(hs; tau0=1.0), clock_model(hs; tau0=1.0, flicker=true)
        @test (off.flicker_fm, off.flicker_pm) == (0, 0) && (on.flicker_fm, on.flicker_pm) == (hf, hp)
        @test clock_model(hs; tau0=1.0, flicker=true, fm_band=(10, 1000)).fm_band == (10.0, 1000.0)
        @test_throws ArgumentError ClockModel(; tau0=1.0, flicker_fm=hf, fm_band=(0.5, 100.0))
    end

    @testset "bad input" begin
        @test_throws ArgumentError ClockModel(; tau0=0.0)
        @test_throws ArgumentError ClockModel(; tau0=1.0, q1=-1.0)
        @test_throws ArgumentError ClockModel(; tau0=1.0, q3=1.0, drift=false)
        @test_throws ArgumentError clock_model(powerlaw_fit([1e-3, 1e-2, 1e-1], [1.0, 1.0, 1.0]; alphas=[0]))
    end
end

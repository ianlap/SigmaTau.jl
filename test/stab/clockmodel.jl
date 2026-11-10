# clockmodel.jl — the clock model's Φ and Q in closed form, its mapping to
# the h_α, and a model recovered from a noise_gen record's power-law fit.

using Test
using SigmaTau
using Random

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

    @testset "bad input" begin
        @test_throws ArgumentError ClockModel(; tau0=0.0)
        @test_throws ArgumentError ClockModel(; tau0=1.0, q1=-1.0)
        @test_throws ArgumentError ClockModel(; tau0=1.0, q3=1.0, drift=false)
        @test_throws ArgumentError clock_model(powerlaw_fit([1e-3, 1e-2, 1e-1], [1.0, 1.0, 1.0]; alphas=[0]))
    end
end

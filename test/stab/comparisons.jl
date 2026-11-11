# comparisons.jl — the multi-clock record and the ensemble simulator.

using Test
using SigmaTau
using Random
using Statistics

@testset "clock comparisons" begin
    a = PhaseData([0.0, 1.0, 2.0], 1.0); b = PhaseData([0.0, 0.5, NaN], 1.0)
    c = ClockComparisons([a, b]; names=["Cs", "Rb"], reference="UTC(lab)")
    @test size(c) == (3, 2) && c.names == ["Cs", "Rb"]
    @test c[2].x[1:2] == [0.0, 0.5] && isnan(c[2].x[3])
    @test c[1, 2].x[2] == 0.5 && c[1, 2].source == "Cs − Rb"
    @test_throws ArgumentError ClockComparisons([a, PhaseData([0.0, 1.0], 1.0)])
    @test_throws ArgumentError ClockComparisons([a, PhaseData(a.x, 2.0)])
    @test_throws ArgumentError ClockComparisons(ones(3, 2), 1.0; names=["x", "x"])
    @test_throws ArgumentError ClockComparisons(ones(3, 1), 1.0)

    # Three clocks against clock 1: its column is zero, a difference of two
    # readings is the difference of the true phases plus both channels' noise,
    # and each clock's ADEV comes out at its model's level.
    ms = [ClockModel(; tau0=1.0, q1=q) for q in (1e-24, 4e-24, 9e-24)]
    ms[2] = ClockModel(; tau0=1.0, R=1e-22, q1=4e-24)
    cc, truth = simulate_ensemble(ms, 20000; rng=MersenneTwister(7))
    @test cc.reference == "clock 1" && all(iszero, cc.x[:, 1])
    @test std(cc[2, 3].x .- (truth[:, 1, 2] .- truth[:, 1, 3])) ≈ 1e-11 rtol=0.05
    # ADEV² = Σq1/τ + 3R/τ² (white FM plus the white PM of the readings).
    for (i, j) in ((2, 1), (3, 1))
        @test adev(cc[i, j], [100]; ci=false).dev[1]^2 ≈ (ms[i].q1 + ms[j].q1) / 100 + 3ms[i].R / 100^2 rtol=0.2
    end
    # Models fitted from the comparisons alone: the white FM level comes back
    # (R and the random walk are only weakly set by one record), and the
    # Kalman timescale runs on them.
    fit = clock_models(cc)
    @test length(fit) == 3 && all(m -> m.tau0 == 1.0, fit)
    @test [m.q1 for m in fit] ≈ [m.q1 for m in ms] rtol=0.3
    @test timescale(cc; algorithm=:greenhall, models=fit) isa TimescaleResult
    # Missing readings: the terms they touch are left out of the covariance.
    # A clock unread for a tenth of the record and one that joins late (within
    # 0.13 over 20 seeds); every clock read at every tenth sample, fitted on
    # the octaves of ten (within 0.44); with three clocks, ten readings of one
    # leave no clock enough terms.
    x = copy(cc.x); x[5001:7000, 2] .= NaN; x[1:3000, 3] .= NaN
    @test [m.q1 for m in clock_models(ClockComparisons(x, 1.0))] ≈ [m.q1 for m in ms] rtol=0.3
    x = copy(cc.x); x[[t % 10 != 1 for t in 1:20000], :] .= NaN
    @test [m.q1 for m in clock_models(ClockComparisons(x, 1.0))] ≈ [m.q1 for m in ms] rtol=0.5
    x = copy(cc.x); x[1:19990, 3] .= NaN
    @test_throws "clock 1 has 1 usable τ points" clock_models(ClockComparisons(x, 1.0))
    # Quiet clocks beside clocks 300 times noisier: the fit weighs the couples
    # by their noise, so the quiet ones are read against each other (within
    # 0.19 over 60 seeds; with equal weights, 0 to 2500 times the truth).
    mix = [ClockModel(; tau0=1.0, q1=q, q2=1e-6q) for q in (1e-27, 4e-27, 2e-27, 1e-22, 2e-22, 4e-22)]
    cm, _ = simulate_ensemble(mix, 20000; reference=0, rng=Xoshiro(2))
    @test all(isapprox(f.q1, m.q1; rtol=0.3) for (f, m) in zip(clock_models(cm), mix))
    # Flicker PM on a reading is fitted as PM and does not pass for white FM.
    x = copy(cc.x)
    x[:, 3] .+= noise_gen(PhaseData, 20000, 1.0; h=Dict(1 => 2e-21), rng=MersenneTwister(8)).x
    noisy = ClockComparisons(x, 1.0)
    @test clock_models(noisy)[3].q1 ≈ ms[3].q1 rtol=0.3
    @test clock_models(noisy; alphas=[2, 0, -2])[3].q1 > 2ms[3].q1
    @test_throws ArgumentError clock_models(ClockComparisons(ones(10, 2), 1.0))
    # With flicker: the levels come back and each clock carries its states;
    # the reference's readings are exact, so it shows no flicker PM.
    fl = [ClockModel(; tau0=1.0, q1=1e-24 * k, q2=1e-31, flicker_fm=1.6e-26 * k, flicker_pm=2.2e-22) for k in 1:3]
    cf, _ = simulate_ensemble(fl, 65536; rng=MersenneTwister(4))
    ff = clock_models(cf; flicker=true)
    @test [m.flicker_fm for m in ff] ≈ [m.flicker_fm for m in fl] rtol=0.5
    @test [m.flicker_pm for m in ff] ≈ [0, 2.2e-22, 2.2e-22] rtol=0.3
    @test all(m -> m.fm_band[1] > 1, ff) && ff[2].pm_band[1] == 1
    ideal, _ = simulate_ensemble(ms[[1, 3]], 100; reference=0)
    @test ideal.reference == "ideal reference" && !iszero(ideal.x[:, 1])
    @test_throws ArgumentError simulate_ensemble(ms[1:1], 10)
    @test_throws ArgumentError simulate_ensemble([ms[1], ClockModel(; tau0=2.0)], 10)
end

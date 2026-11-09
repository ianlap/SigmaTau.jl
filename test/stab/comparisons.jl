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
    ideal, _ = simulate_ensemble(ms[[1, 3]], 100; reference=0)
    @test ideal.reference == "ideal reference" && !iszero(ideal.x[:, 1])
    @test_throws ArgumentError simulate_ensemble(ms[1:1], 10)
    @test_throws ArgumentError simulate_ensemble([ms[1], ClockModel(; tau0=2.0)], 10)
end

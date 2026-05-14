# umbrella_smoke.jl — top-level `using SigmaTau` re-export smoke test.
#
# Verifies that every public symbol the package promises lands in the caller's
# namespace under a bare `using SigmaTau`, that both submodules are reachable
# without explicit qualification, and that the FrequencyData dispatches all
# round-trip through `_freq_to_phase` correctly.

using Test
using LinearAlgebra
using SigmaTau

@testset "Submodules accessible" begin
    @test SigmaTau.Stab isa Module
    @test SigmaTau.Est  isa Module
    # Internal kernels reachable through the dotted path even though only the
    # high-level wrappers are re-exported.
    @test isdefined(SigmaTau.Stab, :_adev_core)
    @test isdefined(SigmaTau.Est, :state_transition)
end

@testset "Stab public API re-exported" begin
    # Must be accessible bare (not via SigmaTau.Stab.adev).
    for sym in (:adev, :mdev, :tdev, :hdev, :mhdev, :htdev,
                :totdev, :mtotdev, :ttotdev, :htotdev, :mhtotdev,
                :mtie, :pdev,
                :identify_noise, :calculate_edf, :confidence_intervals,
                :bias_correction, :DEFAULT_CONFIDENCE,
                :save_result, :load_result,
                :_adev_core, :_mdev_core, :_tdev_core,
                :_hdev_core, :_mhdev_core,
                :_totdev_core, :_mtotdev_core, :_htotdev_core, :_mhtotdev_core,
                :_mtie_core, :_pdev_core,
                :ldev)        # deprecated alias, still re-exported
        @test isdefined(@__MODULE__, sym)
    end
end

@testset "Est public API re-exported" begin
    for sym in (:AbstractClockModel, :TwoStateClock, :ThreeStateClock, :RelativisticClock,
                :nstates, :state_transition, :process_noise,
                :measurement_matrix, :measurement_noise,
                :AbstractEstimator, :StandardKalmanFilter,
                :UDFactorizedFilter, :KuramotoOscillator,
                :predict!, :update!, :prop!,
                :PIDController, :step!, :steer_to_correction)
        @test isdefined(@__MODULE__, sym)
    end
end

@testset "Shared types re-exported" begin
    for sym in (:AbstractTimingData, :PhaseData, :FrequencyData, :StabilityResult)
        @test isdefined(@__MODULE__, sym)
    end
end

@testset "End-to-end Stab call works under bare using" begin
    # Quadratic phase: ADEV(τ) = √2 at m=1 by hand (see test/stab `ADEV Core`).
    p = PhaseData(collect(1.0:100.0) .^ 2, 1.0)
    r = adev(p, [1, 2]; calc_ci=false)
    @test r.deviation_type === :adev
    @test r.dev[1] ≈ sqrt(2.0)
    @test length(r.tau) == 2
end

@testset "End-to-end Est call works under bare using" begin
    m  = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=0.0)
    est = StandardKalmanFilter([0.0, 0.0, 0.0], Matrix(1e-18 * I(3)))
    update!(est, m, 1e-9)
    prop!(est, m, 5.0)
    @test est.k == 1
    @test length(est.x) == 3
    @test size(est.P) == (3, 3)
end

@testset "FrequencyData dispatch on every deviation" begin
    using Random
    Random.seed!(20260509)
    tau0 = 1.0
    y = randn(200) .* 1e-9
    fd = FrequencyData(y, tau0)
    pd_eq = PhaseData(cumsum(y) .* tau0, tau0)
    ms = [1, 2, 4]

    for f in (adev, mdev, tdev, hdev, mhdev, htdev,
              totdev, mtotdev, htotdev, mhtotdev,
              mtie, pdev)
        a = f(fd,    ms; calc_ci=false)
        b = f(pd_eq, ms; calc_ci=false)
        @test a.dev ≈ b.dev
    end
end

@testset "ldev is a forwarding alias for htdev" begin
    pd = PhaseData(cumsum(randn(64)), 1.0)
    a = htdev(pd, [1, 2]; calc_ci=false)
    b = ldev(pd, [1, 2]; calc_ci=false)
    @test a.dev == b.dev
    @test a.deviation_type == b.deviation_type == :htdev
end

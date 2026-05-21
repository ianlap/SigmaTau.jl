# umbrella_smoke.jl — top-level `using SigmaTau` re-export smoke test.
#
# Verifies that every public symbol the package promises lands in the caller's
# namespace under a bare `using SigmaTau`, and that the FrequencyData
# dispatches all round-trip through `_freq_to_phase` correctly. As of v0.3.0
# the `Stab` submodule is flattened into the umbrella and the `Est` surface
# has moved to ClockEnsemble.jl; both assertions are pinned negative below.

using Test
using SigmaTau

@testset "Flattened — no Stab/Est submodules" begin
    # Stab is no longer a submodule (flattened into the umbrella).
    @test !isdefined(SigmaTau, :Stab)
    # Est moved to ClockEnsemble.jl.
    @test !isdefined(SigmaTau, :Est)
    # Internal kernels still reachable through `SigmaTau.` directly.
    @test isdefined(SigmaTau, :_adev_core)
end

@testset "Public API re-exported" begin
    # Must be accessible bare under `using SigmaTau`.
    for sym in (:adev, :mdev, :tdev, :hdev, :mhdev, :htdev,
                :totdev, :mtotdev, :ttotdev, :htotdev, :mhtotdev,
                :mtie, :pdev,
                :identify_noise, :calculate_edf, :confidence_intervals,
                :bias_correction, :DEFAULT_CONFIDENCE,
                :save_result, :load_result,
                :read_phase, :read_frequency,
                :detrend, :fillgaps,
                :noise_gen,
                :_adev_core, :_mdev_core, :_tdev_core,
                :_hdev_core, :_mhdev_core,
                :_totdev_core, :_mtotdev_core, :_htotdev_core, :_mhtotdev_core,
                :_mtie_core, :_pdev_core,
                :ldev)        # deprecated alias, still re-exported
        @test isdefined(@__MODULE__, sym)
    end
end

@testset "Shared types re-exported" begin
    for sym in (:AbstractTimingData, :PhaseData, :FrequencyData, :StabilityResult)
        @test isdefined(@__MODULE__, sym)
    end
end

@testset "End-to-end deviation call works under bare using" begin
    # Quadratic phase: ADEV(τ) = √2 at m=1 by hand (see test/stab `ADEV Core`).
    p = PhaseData(collect(1.0:100.0) .^ 2, 1.0)
    r = adev(p, [1, 2]; calc_ci=false)
    @test r.deviation_type === :adev
    @test r.dev[1] ≈ sqrt(2.0)
    @test length(r.tau) == 2
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

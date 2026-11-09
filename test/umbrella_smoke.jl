# umbrella_smoke.jl — the exported surface under a bare `using SigmaTau`, and
# FrequencyData dispatch on every estimator.

using Test
using Random
using SigmaTau

@testset "exported surface" begin
    public = Set([
        :AbstractTimingData, :PhaseData, :FrequencyData, :StabilityResult,
        :StabilitySuite, :SpectralResult, :DynamicStabilityResult,
        :ci_lower, :ci_upper, :DEFAULT_CONFIDENCE, :DEFAULT_DEVIATIONS,
        :TauMode, :AllTaus, :Octave, :HalfOctave, :QuarterOctave, :Decade,
        :HalfDecade, :tau_values, :tau_mode,
        :adev, :mdev, :tdev, :hdev, :mhdev, :htdev, :totdev, :mtotdev,
        :ttotdev, :htotdev, :mhtotdev, :mtie, :tierms, :pdev, :theo1, :theoh,
        :dadev, :dhdev, :stability, :deviation, :estimators, :summarize,
        :ncornered_hat, :StreamingStability, :snapshot, :nsamples,
        :identify_noise, :noise_id_nnls, :calculate_edf, :confidence_intervals, :bias_correction,
        :LiveStability, :SimulatedSource, :InstrumentSource, :read_reading,
        :psd_y, :psd_x, :ssb_phase_noise, :noise_gen,
        :cross_psd_y, :cross_psd_x, :cross_adev, :CrossSpectralResult, :CrossStabilityResult,
        :read_phase, :read_frequency, :read_data, :detect_format, :INPUT_FORMATS,
        :save, :load,
        :detrend, :fillgaps, :find_outliers, :remove_outliers,
        :drift_fit, :DriftFit])
    deprecated = Set([:Sy, :Sx, :L, :nch, :save_result, :save_suite,
                      :load_result, :load_suite])
    @test setdiff(Set(names(SigmaTau)), [:SigmaTau]) == public ∪ deprecated
    # The kernels stay internal but reachable as SigmaTau._adev_core etc.
    @test isdefined(SigmaTau, :_adev_core)
end

@testset "FrequencyData dispatch on every estimator" begin
    y = randn(MersenneTwister(20260509), 200) .* 1e-9
    fd = FrequencyData(y, 1.0)
    pd = PhaseData([0.0; cumsum(y)], 1.0)
    fns = (adev, mdev, tdev, hdev, mhdev, htdev, totdev, mtotdev, ttotdev,
           htotdev, mhtotdev, mtie, tierms, pdev)
    @test isempty([f for f in fns if !(f(fd, [1, 2, 4]; ci=false).dev ≈
                                       f(pd, [1, 2, 4]; ci=false).dev)])
end

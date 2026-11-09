module SigmaTau

using DocStringExtensions
using Statistics
using Distributions
using FFTW
using LinearAlgebra: I, Diagonal, Symmetric, diag, eigen
using Dates
using Printf

# ── Shared types ────────────────────────────────────────────────────────
include("types.jl")

# ── IO: file readers, detrend, gap fill, result round-trip ──────────────
include("io/results.jl")
include("io/detrend.jl")
include("io/fillgaps.jl")
include("io/read.jl")
include("io/formats.jl")
include("io/save.jl")
include("io/outliers.jl")

"""
Default confidence level (0.683, one sigma, the Stable32 convention) for
every deviation with a χ² confidence model. Override per call with
`confidence=0.95`. `mtie` and `tierms` have no interval model and ignore it.
"""
const DEFAULT_CONFIDENCE = 0.683

# ── Clock-stability analysis ────────────────────────────────────────────
# Dependency order: kernels (raw _*_core) → noise ID / EDF stats → grids and
# conversion helpers → spectral → public deviation API → compute-all suite.
include("kernels.jl")
include("noise.jl")
include("edf.jl")
include("grids.jl")
include("spectral.jl")
include("gaps.jl")
include("cross.jl")
include("drift.jl")
include("powerlaw.jl")
include("clockmodel.jl")
include("kalman.jl")
include("deviations.jl")
include("registry.jl")
include("suite.jl")
include("ensemble.jl")

# ── Real-time streaming accumulators and live sources ──────────────────
include("streaming.jl")
include("live.jl")

# ── Names replaced in 0.6 (deprecated, exported until 0.7) ──────────────
include("deprecated.jl")

# ── Flat exports ────────────────────────────────────────────────────────
export AbstractTimingData, PhaseData, FrequencyData, StabilityResult, StabilitySuite
export SpectralResult, DynamicStabilityResult
export ci_lower, ci_upper

export save, load
export read_phase, read_frequency, read_data, detect_format, INPUT_FORMATS
export detrend, fillgaps
export find_outliers, remove_outliers

export DEFAULT_CONFIDENCE

export TauMode, AllTaus, Octave, HalfOctave, QuarterOctave, Decade, HalfDecade
export tau_values

# The deviation kernels (`_*_core`) are not exported. They are reachable as
# `SigmaTau._adev_core` etc.; the leading underscore marks them unsupported.

export identify_noise, noise_id_nnls, calculate_edf, confidence_intervals, bias_correction

export adev, mdev, tdev
export hdev, mhdev, htdev
export totdev, mtotdev, ttotdev, htotdev, mhtotdev
export mtie, tierms, pdev
export dadev, dhdev
export ncornered_hat, groslambert
export drift_fit, DriftFit
export powerlaw_fit, PowerLawFit, powerlaw_sigma, powerlaw_psd
export ClockModel, clock_model, state_transition, process_noise
export KalmanResult, kalman_filter, kalman_smooth, holdover_error, innovation_check, simulate_clock
export prediction_error, tune_q
export theo1, theoh
export stability, DEFAULT_DEVIATIONS
export estimators, deviation, tau_mode, summarize
export StreamingStability, snapshot, nsamples
export LiveStability, SimulatedSource, InstrumentSource, read_reading

export noise_gen
export psd_y, psd_x, ssb_phase_noise
export cross_psd_y, cross_psd_x, cross_adev, CrossSpectralResult, CrossStabilityResult

# Plot recipes for `StabilityResult` live in the `SigmaTauRecipesBaseExt`
# package extension and load automatically when `RecipesBase` (or `Plots`) is.

# ── Precompile workload ─────────────────────────────────────────────────
# Run every public entry point once on tiny deterministic records during
# package precompilation, so the first real call in a fresh session is fast
# instead of paying full JIT compilation. Values are irrelevant; the records
# are just long enough for every kernel to produce valid windows at m ≤ 2.
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    _pc_x = cumsum(sinpi.((1:128) ./ 7.3)) ./ 1.0e9   # deterministic "phase"
    @compile_workload begin
        pd = PhaseData(_pc_x, 1.0)
        fd = FrequencyData(diff(_pc_x), 1.0)
        ms = [1, 2]
        for f in (adev, mdev, tdev, hdev, mhdev, htdev,
                  totdev, mtotdev, ttotdev, htotdev, mhtotdev, mtie, tierms, pdev)
            f(pd, ms)
            f(fd, ms)
        end
        theo1(pd, [2, 4])
        theo1(fd, [2, 4])
        theoh(pd)
        stability(pd)
        identify_noise(_pc_x, ms)
        detrend(pd)
        detrend(fd; method=:quadratic)
        fillgaps(pd)
        noise_gen(PhaseData, 32, 1.0; sigma1=Dict(0 => 1e-12))
        psd_y(fd)
        psd_x(pd)
        ssb_phase_noise(pd; f_carrier=10e6)
        r = adev(pd, ms)
        show(IOBuffer(), MIME"text/plain"(), r)
        _pc_path = joinpath(mktempdir(), "r.tsv")
        save(_pc_path, r)
        load(_pc_path)
        _pc_acc = StreamingStability(:adev, 1.0, ms)
        append!(_pc_acc, _pc_x[1:32])
        snapshot(_pc_acc)
        nsamples(_pc_acc)
        drift_fit(fd)
        _pc_fit = powerlaw_fit(pd; alphas=[2, 0, -2])
        powerlaw_sigma(_pc_fit, [1.0, 2.0])
        cross_adev(pd, pd, ms)
        cross_psd_y(pd, pd)
        groslambert(pd, pd, pd, ms)
        _pc_live = LiveStability([:adev], 1.0, ms)
        push!(_pc_live, 0.0)
        snapshot(_pc_live)
    end
end

end # module SigmaTau

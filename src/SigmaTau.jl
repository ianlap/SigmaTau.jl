module SigmaTau

using DocStringExtensions
using Statistics
using Distributions
using StaticArrays
using DelimitedFiles
using FFTW
using Dates
using Printf

# ── Shared types ────────────────────────────────────────────────────────
include("types.jl")

# ── IO: file readers, detrend, gap fill, result round-trip ──────────────
# IO functions return top-level types (PhaseData / FrequencyData /
# StabilityResult). DelimitedFiles, FFTW, and Statistics (for `median`)
# are imported above so included files don't need their own directives.
include("io/results.jl")
include("io/detrend.jl")
include("io/fillgaps.jl")
include("io/read.jl")
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
include("deviations.jl")
include("suite.jl")

# ── Real-time streaming accumulators (Dobrogowski–Kasznia scheme) ───────
include("streaming.jl")

# ── Flat exports ────────────────────────────────────────────────────────
export AbstractTimingData, PhaseData, FrequencyData, StabilityResult, StabilitySuite
export SpectralResult, DynamicStabilityResult
export ci_lower, ci_upper

export save, save_result, load_result, save_suite, load_suite
export read_phase, read_frequency
export detrend, fillgaps
export find_outliers, remove_outliers

export DEFAULT_CONFIDENCE

export TauMode, AllTaus, Octave, HalfOctave, QuarterOctave, Decade, HalfDecade
export tau_values

# The deviation kernels (`_*_core`) are not exported. They are reachable as
# `SigmaTau._adev_core` etc.; the leading underscore marks them unsupported.

export identify_noise, calculate_edf, confidence_intervals, bias_correction

export adev, mdev, tdev
export hdev, mhdev, htdev
export totdev, mtotdev, ttotdev, htotdev, mhtotdev
export mtie, tierms, pdev
export dadev, dhdev
export nch
export theo1, theoh
export stability, DEFAULT_DEVIATIONS
export StreamingStability, snapshot, nsamples

export noise_gen
export Sy, Sx, L

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
        Sy(fd)
        Sx(pd)
        L(pd; f_carrier=10e6)
        r = adev(pd, ms)
        show(IOBuffer(), MIME"text/plain"(), r)
        _pc_path = joinpath(mktempdir(), "r.tsv")
        save_result(_pc_path, r)
        load_result(_pc_path)
        _pc_acc = StreamingStability(:adev, 1.0, ms)
        append!(_pc_acc, _pc_x[1:32])
        snapshot(_pc_acc)
        nsamples(_pc_acc)
    end
end

end # module SigmaTau

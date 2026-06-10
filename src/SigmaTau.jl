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
Package-wide default confidence factor used by every public deviation API
with a χ²/EDF confidence model (`adev`, `mdev`, `tdev`, `hdev`, `mhdev`,
`htdev`, `totdev`, `mtotdev`, `ttotdev`, `htotdev`, `mhtotdev`, `pdev`,
`theo1`, `theoh`) when `confidence` is not supplied. (`mtie` and `tierms`
have no published EDF model and ignore it.)

Set to 0.683 (1-sigma) — the time-and-frequency stability convention used
by Stable32, allantools' published error bars, and the
Greenhall–Riley uncertainty papers. Override per call by passing
`confidence=0.95` (or any other level) explicitly.
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

# ── Flat exports ────────────────────────────────────────────────────────
export AbstractTimingData, PhaseData, FrequencyData, StabilityResult, StabilitySuite
export SpectralResult
export ci_lower, ci_upper

export save, save_result, load_result, save_suite, load_suite
export read_phase, read_frequency
export detrend, fillgaps
export find_outliers, remove_outliers

export DEFAULT_CONFIDENCE

export TauMode, AllTaus, Octave, HalfOctave, QuarterOctave, Decade, HalfDecade
export tau_values

# Internal deviation kernels (`_*_core`) are intentionally NOT exported — they
# are reachable as `SigmaTau._adev_core` etc. for power users, but the leading
# underscore marks them unsupported and they are kept out of the bare namespace.

export identify_noise, calculate_edf, confidence_intervals, bias_correction

export adev, mdev, tdev
export hdev, mhdev, htdev
export totdev, mtotdev, ttotdev, htotdev, mhtotdev
export mtie, tierms, pdev
export nch
export theo1, theoh
export stability, DEFAULT_DEVIATIONS

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
    end
end

end # module SigmaTau

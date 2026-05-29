module SigmaTau

using DocStringExtensions
using Statistics
using Distributions
using StaticArrays
using DelimitedFiles
using FFTW
using Dates

# ── Shared types ────────────────────────────────────────────────────────
include("types/abstract.jl")
include("types/phase_data.jl")
include("types/frequency_data.jl")
include("types/stability_result.jl")
include("types/stability_suite.jl")

# ── IO: file readers, detrend, gap fill, result round-trip ──────────────
# IO functions return top-level types (PhaseData / FrequencyData /
# StabilityResult). DelimitedFiles, FFTW, and Statistics (for `median`)
# are imported above so included files don't need their own directives.
include("io/results.jl")
include("io/detrend.jl")
include("io/fillgaps.jl")
include("io/read.jl")

"""
Package-wide default confidence factor used by every public deviation API
(`adev`, `mdev`, `hdev`, `tdev`, `mhdev`, `htdev`, `totdev`, `mtotdev`,
`htotdev`, `mhtotdev`) when `confidence` is not supplied.

Set to 0.683 (1-sigma) — the time-and-frequency stability convention used
by Stable32, AllanLab, allantools' published error bars, and the
Greenhall–Riley uncertainty papers. Override per call by passing
`confidence=0.95` (or any other level) explicitly.
"""
const DEFAULT_CONFIDENCE = 0.683

# ── Stab: clock-stability analysis ──────────────────────────────────────
include("stab/core/allan.jl")
include("stab/core/hadamard.jl")
include("stab/core/total.jl")
include("stab/core/mtie.jl")
include("stab/core/pdev.jl")

include("stab/noise/lag1.jl")
include("stab/noise/synth.jl")
include("stab/noise/gen.jl")
include("stab/stats/edf.jl")

include("stab/utils.jl")
include("stab/taus.jl")

include("stab/api/allan.jl")
include("stab/api/hadamard.jl")
include("stab/api/total.jl")
include("stab/api/mtie.jl")
include("stab/api/pdev.jl")
include("stab/api/suite.jl")

# ── Flat exports ────────────────────────────────────────────────────────
export AbstractTimingData, PhaseData, FrequencyData, StabilityResult, StabilitySuite

export save_result, load_result, save_suite, load_suite
export read_phase, read_frequency
export detrend, fillgaps

export DEFAULT_CONFIDENCE
export TOTDEV_DETREND_DEFAULT, TOTAL_FAMILY_DETREND_DEFAULT

export TauMode, AllTaus, Octave, HalfOctave, QuarterOctave, Decade, HalfDecade
export tau_values

# Internal deviation kernels (`_*_core`) are intentionally NOT exported — they
# are reachable as `SigmaTau._adev_core` etc. for power users, but the leading
# underscore marks them unsupported and they are kept out of the bare namespace.

export identify_noise, calculate_edf, confidence_intervals, bias_correction

export adev, mdev, tdev
export hdev, mhdev, htdev
export ldev   # deprecated alias for htdev — remove in a future release
export totdev, mtotdev, ttotdev, htotdev, mhtotdev
export mtie, pdev
export stability, DEFAULT_DEVIATIONS

export noise_gen

# Plot recipes for `StabilityResult` live in the `SigmaTauRecipesBaseExt`
# package extension and load automatically when `RecipesBase` (or `Plots`) is.

end # module SigmaTau

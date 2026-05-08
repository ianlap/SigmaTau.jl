module SigmaTauStability

using Statistics
using Distributions
using DocStringExtensions
using SigmaTauBase

"""
Package-wide default confidence factor used by every public deviation API
(`adev`, `mdev`, `hdev`, `tdev`, `mhdev`, `ldev`, `totdev`, `mtotdev`,
`htotdev`, `mhtotdev`) when `confidence` is not supplied.

Set to 0.683 (1-sigma) — the time-and-frequency stability convention used
by Stable32, AllanLab, allantools' published error bars, and the
Greenhall–Riley uncertainty papers. Override per call by passing
`confidence=0.95` (or any other level) explicitly.
"""
const DEFAULT_CONFIDENCE = 0.683

include("core/allan.jl")
include("core/hadamard.jl")
include("core/total.jl")

include("noise/lag1.jl")
include("noise/synth.jl")
include("stats/edf.jl")

include("utils.jl")

include("api/allan.jl")
include("api/hadamard.jl")
include("api/total.jl")

export _adev_core, _mdev_core, _tdev_core
export _hdev_core, _mhdev_core
export _totdev_core, _mtotdev_core, _htotdev_core, _mhtotdev_core

export identify_noise, calculate_edf, confidence_intervals, bias_correction
export DEFAULT_CONFIDENCE

export adev, mdev, tdev
export hdev, mhdev, ldev
export totdev, mtotdev, htotdev, mhtotdev

end # module

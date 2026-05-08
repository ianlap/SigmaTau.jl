module SigmaTauStability

using Statistics
using Distributions
using DocStringExtensions
using SigmaTauBase

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

export adev, mdev, tdev
export hdev, mhdev, ldev
export totdev, mtotdev, htotdev, mhtotdev

end # module

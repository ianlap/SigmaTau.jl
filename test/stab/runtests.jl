using Test
using Random
using SigmaTau
using SigmaTau: NEFF_RELIABLE
using SigmaTau: _adev_core, _mdev_core, _tdev_core, _hdev_core, _mhdev_core,
                _totdev_core, _mtotdev_core, _htotdev_core, _mhtotdev_core,
                _mtie_core, _pdev_core

include("reference_kernels.jl")
const LK = ReferenceKernels

include("kernels.jl")
include("stable32gen.jl")
include("noise_types.jl")
include("noise_id.jl")
include("defaults.jl")
include("taus.jl")
include("suite.jl")
include("spectral.jl")
include("noise_gen.jl")
include("cross.jl")
include("mhtotdev_mc.jl")
include("pdev_edf_mc.jl")
include("theo.jl")
include("tierms.jl")
include("gaps.jl")
include("drift.jl")
include("powerlaw.jl")
include("clockmodel.jl")
include("kalman.jl")
include("ensemble.jl")
include("dynamic.jl")
include("streaming.jl")
include("live.jl")
include("api_inputs.jl")
include("registry.jl")
include("deprecated.jl")
include("nnls_noise_id.jl")

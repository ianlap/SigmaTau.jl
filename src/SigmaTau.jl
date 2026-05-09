module SigmaTau

using DocStringExtensions
using Reexport

# ── Shared types ────────────────────────────────────────────────────────
include("types/abstract.jl")
include("types/phase_data.jl")
include("types/frequency_data.jl")
include("types/stability_result.jl")

export AbstractTimingData, PhaseData, FrequencyData, StabilityResult

# ── Stab: clock-stability analysis ──────────────────────────────────────
module Stab
    using ..SigmaTau: AbstractTimingData, PhaseData, FrequencyData,
                      StabilityResult
    using Statistics
    using Distributions
    using DocStringExtensions

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

    include("stab/core/allan.jl")
    include("stab/core/hadamard.jl")
    include("stab/core/total.jl")
    include("stab/core/mtie.jl")
    include("stab/core/pdev.jl")

    include("stab/noise/lag1.jl")
    include("stab/noise/synth.jl")
    include("stab/stats/edf.jl")

    include("stab/utils.jl")

    include("stab/api/allan.jl")
    include("stab/api/hadamard.jl")
    include("stab/api/total.jl")
    include("stab/api/mtie.jl")
    include("stab/api/pdev.jl")

    export _adev_core, _mdev_core, _tdev_core
    export _hdev_core, _mhdev_core
    export _totdev_core, _mtotdev_core, _htotdev_core, _mhtotdev_core
    export _mtie_core, _pdev_core

    export identify_noise, calculate_edf, confidence_intervals, bias_correction
    export DEFAULT_CONFIDENCE

    export adev, mdev, tdev
    export hdev, mhdev, htdev
    export ldev   # deprecated alias for htdev — remove in a future release
    export totdev, mtotdev, htotdev, mhtotdev
    export mtie, pdev
end

# ── Est: clock estimation ───────────────────────────────────────────────
module Est
    using ..SigmaTau: AbstractTimingData, PhaseData, FrequencyData,
                      StabilityResult
    using LinearAlgebra
    using StaticArrays
    using DocStringExtensions

    include("est/models/clocks.jl")
    include("est/estimators/filters.jl")

    export AbstractClockModel, TwoStateClock, ThreeStateClock, RelativisticClock
    export nstates, state_transition, process_noise, measurement_matrix, measurement_noise
    export AbstractEstimator, StandardKalmanFilter, UDFactorizedFilter, KuramotoOscillator
    export predict!, update!
    export PIDController, step!, steer_to_correction
end

# ── Flatten submodule exports onto the umbrella ─────────────────────────
@reexport using .Stab
@reexport using .Est

# Make submodules themselves accessible without qualification.
export Stab, Est

# Plot recipes for `StabilityResult` live in the `SigmaTauRecipesBaseExt`
# package extension and load automatically when `RecipesBase` (or `Plots`) is.

end # module SigmaTau

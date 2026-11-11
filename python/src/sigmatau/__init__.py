"""sigmatau: analysis of clock and oscillator stability, the Python package of SigmaTau.jl.

Each function calls the Julia function of the same name in a Julia runtime
inside the Python process, through juliacall, so the numbers are those of the
Julia package. The names and keywords are the same, and the documentation at
https://ianlap.github.io/SigmaTau.jl gives the equations, the sources cited in
these docstrings and what differs between the two languages.
"""

__version__ = "0.6.0"

from ._julia import julia_version, sigmatau_jl_version
from .cross import cross_adev, cross_psd_x, cross_psd_y
from .deviations import (
    DEFAULT_CONFIDENCE,
    adev,
    deviation,
    estimators,
    hdev,
    htdev,
    htotdev,
    mdev,
    mhdev,
    mhtotdev,
    mtie,
    mtotdev,
    pdev,
    tdev,
    theo1,
    theoh,
    tierms,
    totdev,
    ttotdev,
)
from .drift import DRIFT_MODELS, drift_fit
from .dynamic import dadev, dhdev
from .edf import bias_correction, calculate_edf, confidence_intervals
from .ensemble import (
    ClockComparisons,
    ClockModel,
    HoldoverPrediction,
    HoldoverResult,
    KalmanResult,
    NoiseFit,
    SteeringGains,
    SteeringPrediction,
    SteeringResult,
    TimescaleResult,
    clock_model,
    clock_models,
    holdover,
    holdover_error,
    implicit_mean,
    innovation_check,
    kalman_filter,
    kalman_smooth,
    lqg_gains,
    noise_fit,
    prediction_error,
    process_noise,
    simulate_clock,
    simulate_ensemble,
    state_transition,
    steer,
    steering_gains,
    steering_poles,
    timescale,
    timescale_phase,
    tune_q,
)
from .grids import (
    AllTaus,
    Decade,
    HalfDecade,
    HalfOctave,
    Octave,
    QuarterOctave,
    TauMode,
    tau_values,
)
from .io import (
    detect_format,
    detrend,
    fillgaps,
    find_outliers,
    load,
    read_data,
    read_frequency,
    read_phase,
    remove_outliers,
    save,
    summarize,
)
from .masks import TELECOM_MASKS, MaskCheck, mask_check, telecom_mask
from .noise import identify_noise, noise_gen, noise_id_nnls
from .powerlaw import powerlaw_fit, powerlaw_psd, powerlaw_sigma
from .spectral import integrated_jitter, psd_x, psd_y, ssb_phase_noise
from .streaming import StreamingStability
from .suite import DEFAULT_DEVIATIONS, groslambert, ncornered_hat, stability
from .types import (
    CrossSpectralResult,
    CrossStabilityResult,
    DriftFit,
    DynamicStabilityResult,
    FrequencyData,
    PhaseData,
    PowerLawFit,
    SpectralResult,
    StabilityResult,
    StabilitySuite,
    start_time,
)

__all__ = [
    # data + result types
    "PhaseData",
    "FrequencyData",
    "StabilityResult",
    "StabilitySuite",
    "SpectralResult",
    "DynamicStabilityResult",
    "DriftFit",
    "PowerLawFit",
    "CrossSpectralResult",
    "CrossStabilityResult",
    # tau grids
    "TauMode",
    "AllTaus",
    "Octave",
    "HalfOctave",
    "QuarterOctave",
    "Decade",
    "HalfDecade",
    "tau_values",
    # static deviations
    "deviation",
    "estimators",
    "adev",
    "mdev",
    "tdev",
    "hdev",
    "mhdev",
    "htdev",
    "totdev",
    "mtotdev",
    "ttotdev",
    "htotdev",
    "mhtotdev",
    "mtie",
    "tierms",
    "pdev",
    "theo1",
    "theoh",
    "DEFAULT_CONFIDENCE",
    # telecom masks
    "TELECOM_MASKS",
    "MaskCheck",
    "telecom_mask",
    "mask_check",
    # dynamic deviations
    "dadev",
    "dhdev",
    # compute-all suite + clock ensembles
    "stability",
    "StreamingStability",
    "DEFAULT_DEVIATIONS",
    "ncornered_hat",
    "groslambert",
    # two channels
    "cross_adev",
    "cross_psd_y",
    "cross_psd_x",
    # drift and power-law models
    "drift_fit",
    "DRIFT_MODELS",
    "powerlaw_fit",
    "powerlaw_sigma",
    "powerlaw_psd",
    # clock models, Kalman filter, timescales
    "ClockModel",
    "clock_model",
    "state_transition",
    "process_noise",
    "NoiseFit",
    "noise_fit",
    "KalmanResult",
    "kalman_filter",
    "kalman_smooth",
    "holdover_error",
    "HoldoverResult",
    "HoldoverPrediction",
    "holdover",
    "SteeringGains",
    "steering_gains",
    "lqg_gains",
    "steering_poles",
    "SteeringResult",
    "SteeringPrediction",
    "steer",
    "innovation_check",
    "prediction_error",
    "tune_q",
    "simulate_clock",
    "ClockComparisons",
    "simulate_ensemble",
    "clock_models",
    "TimescaleResult",
    "timescale",
    "timescale_phase",
    "implicit_mean",
    # statistics
    "calculate_edf",
    "confidence_intervals",
    "bias_correction",
    # noise identification + generation
    "identify_noise",
    "noise_id_nnls",
    "noise_gen",
    # spectral estimators
    "psd_y",
    "psd_x",
    "ssb_phase_noise",
    "integrated_jitter",
    # file IO
    "read_data",
    "start_time",
    "detect_format",
    "read_phase",
    "read_frequency",
    "detrend",
    "fillgaps",
    "find_outliers",
    "remove_outliers",
    "summarize",
    "save",
    "load",
    # runtime introspection
    "julia_version",
    "sigmatau_jl_version",
]

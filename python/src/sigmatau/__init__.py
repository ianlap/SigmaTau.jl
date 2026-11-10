"""sigmatau — clock-stability analysis, Python bindings for SigmaTau.jl.

The public surface mirrors the retired pure-Python port one-to-one (plus the
estimators it never had); all computation runs in an embedded Julia runtime
via juliacall, so results are bit-identical to the Julia package.
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
from .noise import identify_noise, noise_gen
from .powerlaw import powerlaw_fit, powerlaw_psd, powerlaw_sigma
from .spectral import psd_x, psd_y, ssb_phase_noise
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
    # statistics
    "calculate_edf",
    "confidence_intervals",
    "bias_correction",
    # noise identification + generation
    "identify_noise",
    "noise_gen",
    # spectral estimators
    "psd_y",
    "psd_x",
    "ssb_phase_noise",
    # file IO
    "read_data",
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

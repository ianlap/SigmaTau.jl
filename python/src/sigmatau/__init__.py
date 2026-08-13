"""sigmatau — clock-stability analysis, Python bindings for SigmaTau.jl.

The public surface mirrors the retired pure-Python port one-to-one (plus the
estimators it never had); all computation runs in an embedded Julia runtime
via juliacall, so results are bit-identical to the Julia oracle.
"""

__version__ = "0.5.0"

from ._julia import julia_version, sigmatau_jl_version
from .deviations import (
    DEFAULT_CONFIDENCE,
    adev,
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
    detrend,
    fillgaps,
    find_outliers,
    load_result,
    load_suite,
    read_frequency,
    read_phase,
    remove_outliers,
    save_result,
    save_suite,
)
from .noise import identify_noise, noise_gen
from .spectral import L, Sx, Sy
from .suite import DEFAULT_DEVIATIONS, nch, stability
from .types import (
    DynamicStabilityResult,
    FrequencyData,
    PhaseData,
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
    # dynamic deviations
    "dadev",
    "dhdev",
    # compute-all suite + cornered hat
    "stability",
    "DEFAULT_DEVIATIONS",
    "nch",
    # statistics
    "calculate_edf",
    "confidence_intervals",
    "bias_correction",
    # noise identification + generation
    "identify_noise",
    "noise_gen",
    # spectral estimators
    "Sy",
    "Sx",
    "L",
    # file IO
    "read_phase",
    "read_frequency",
    "detrend",
    "fillgaps",
    "find_outliers",
    "remove_outliers",
    "save_result",
    "load_result",
    "save_suite",
    "load_suite",
    # runtime introspection
    "julia_version",
    "sigmatau_jl_version",
]

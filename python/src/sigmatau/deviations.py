"""Deviation wrappers: each calls the identically-named Julia function and
down-converts the returned ``StabilityResult``.

``taus`` accepts a ``TauMode`` grid selector (default ``Octave``) or an
explicit sequence of integer averaging factors; noise ID + χ² confidence
bounds are computed unless ``ci=False``; the total family applies bias
correction by default (``correct_bias=True``); ``mhtotdev`` additionally
removes record-level quadratic drift by default (``remove_drift=True``).
``mtie``/``tierms`` have no CI model — their ``ci``/``confidence`` kwargs are
accepted and ignored, matching the Julia surface.
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import resolve_taus, result_from_jl, to_jl_data
from ._julia import _map_errors, jl
from .grids import Octave, TauMode
from .types import FrequencyData, PhaseData, StabilityResult

DEFAULT_CONFIDENCE = float(jl.DEFAULT_CONFIDENCE)


@_map_errors
def _run(jlfn, data, taus, **kwargs) -> StabilityResult:
    """Shared call path: convert inputs, call Julia, down-convert the result."""
    return result_from_jl(jlfn(to_jl_data(data), resolve_taus(taus), **kwargs))


def adev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Overlapping Allan deviation σ_y(τ)."""
    return _run(jl.adev, data, taus, ci=ci, confidence=float(confidence))


def mdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Allan deviation Mod σ_y(τ)."""
    return _run(jl.mdev, data, taus, ci=ci, confidence=float(confidence))


def tdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Time deviation σ_x(τ) = (τ/√3)·Mod σ_y(τ). Units of seconds."""
    return _run(jl.tdev, data, taus, ci=ci, confidence=float(confidence))


def hdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Overlapping Hadamard deviation."""
    return _run(jl.hdev, data, taus, ci=ci, confidence=float(confidence))


def mhdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Hadamard deviation."""
    return _run(jl.mhdev, data, taus, ci=ci, confidence=float(confidence))


def htdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Hadamard time deviation σ_x,HT(τ) = (τ/√(10/3))·Mod σ_y,MH(τ). Units of seconds."""
    return _run(jl.htdev, data, taus, ci=ci, confidence=float(confidence))


def pdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Parabolic deviation σ_PDEV(τ) (Vernotte). PDEV(τ₀) ≡ ADEV(τ₀)."""
    return _run(jl.pdev, data, taus, ci=ci, confidence=float(confidence))


def theoh(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """ThêoH composite σ_y(τ) (SP1065 §5.2.16): overlapping ADEV at short τ, ThêoBR beyond."""
    return _run(jl.theoh, data, taus, ci=ci, confidence=float(confidence))


def mtie(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Maximum Time Interval Error (ITU-T G.810). Units of seconds.

    MTIE has no published EDF/CI model, so ``ci``/``confidence`` are accepted for
    signature uniformity but are no-ops (the result always has empty CI fields).
    """
    return _run(jl.mtie, data, taus, ci=ci, confidence=float(confidence))


def tierms(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """RMS Time Interval Error (TIE rms, SP1065 §5.2.18). Units of seconds.

    Like ``mtie``, TIE rms has no published EDF/CI model, so ``ci``/``confidence``
    are accepted for signature uniformity but are no-ops (empty CI fields).
    """
    return _run(jl.tierms, data, taus, ci=ci, confidence=float(confidence))


def totdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Total deviation (Howe / SP1065 mean-flip extension)."""
    return _run(
        jl.totdev, data, taus, ci=ci, correct_bias=correct_bias, confidence=float(confidence)
    )


def mtotdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Total deviation (Greenhall extension)."""
    return _run(
        jl.mtotdev, data, taus, ci=ci, correct_bias=correct_bias, confidence=float(confidence)
    )


def ttotdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Time Total deviation σ_x(τ) = (τ/√3)·Mod-Total σ_y(τ). Units of seconds."""
    return _run(
        jl.ttotdev, data, taus, ci=ci, correct_bias=correct_bias, confidence=float(confidence)
    )


def htotdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Hadamard Total deviation (Greenhall extension on y = diff(x))."""
    return _run(
        jl.htotdev, data, taus, ci=ci, correct_bias=correct_bias, confidence=float(confidence)
    )


def theo1(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Thêo1 deviation (SP1065 §5.2.15) — extends the usable τ range to 75 % of the record span."""
    return _run(
        jl.theo1, data, taus, ci=ci, correct_bias=correct_bias, confidence=float(confidence)
    )


def mhtotdev(
    data: PhaseData | FrequencyData,
    taus: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool = True,
    remove_drift: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Hadamard Total deviation (SigmaTau-original; Greenhall methodology)."""
    return _run(
        jl.mhtotdev,
        data,
        taus,
        ci=ci,
        correct_bias=correct_bias,
        remove_drift=remove_drift,
        confidence=float(confidence),
    )


__all__ = [
    "DEFAULT_CONFIDENCE",
    "adev",
    "hdev",
    "htdev",
    "htotdev",
    "mdev",
    "mhdev",
    "mhtotdev",
    "mtie",
    "mtotdev",
    "pdev",
    "tdev",
    "theo1",
    "theoh",
    "tierms",
    "totdev",
    "ttotdev",
]

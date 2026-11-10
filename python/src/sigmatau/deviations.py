"""Deviation wrappers: each calls the identically-named Julia function and
down-converts the returned ``StabilityResult``.

``m`` accepts a ``TauMode`` grid selector (default ``Octave``) or an
explicit sequence of integer averaging factors; noise ID + χ² confidence
bounds are computed unless ``ci=False``; the total family applies bias
correction by default (``correct_bias=True``; ``"published"`` keeps to the
published tables, ``"stable32"`` does what Stable32 does for each estimator,
``False`` returns the raw statistic, as allantools does); the Hadamard family identifies
noise on the record less its frequency drift by default (``remove_drift=True``;
``False`` uses the record as given, as Stable32 does).
``mtie``/``tierms`` have no CI model — their ``ci``/``confidence`` kwargs are
accepted and ignored, matching the Julia surface.
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import resolve_m, result_from_jl, to_jl_data
from ._julia import _map_errors, jl
from .grids import Octave, TauMode
from .types import FrequencyData, PhaseData, StabilityResult

DEFAULT_CONFIDENCE = float(jl.DEFAULT_CONFIDENCE)


def _bias(correct_bias: bool | str):
    """``correct_bias`` for Julia: a bool, or ``"published"`` / ``"stable32"`` as a Symbol."""
    return jl.Symbol(correct_bias) if isinstance(correct_bias, str) else bool(correct_bias)


@_map_errors
def _run(jlfn, data, m, **kwargs) -> StabilityResult:
    """Shared call path: convert inputs, call Julia, down-convert the result."""
    return result_from_jl(jlfn(to_jl_data(data), resolve_m(m), **kwargs))


def adev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    overlapping: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Allan deviation σ_y(τ): overlapping, or the original non-overlapped one with
    ``overlapping=False`` (allantools' ``adev``, Stable32's Normal Allan)."""
    return _run(jl.adev, data, m, ci=ci, overlapping=overlapping, confidence=float(confidence))


def mdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Allan deviation Mod σ_y(τ)."""
    return _run(jl.mdev, data, m, ci=ci, confidence=float(confidence))


def tdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Time deviation σ_x(τ) = (τ/√3)·Mod σ_y(τ). Units of seconds."""
    return _run(jl.tdev, data, m, ci=ci, confidence=float(confidence))


def hdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    remove_drift: bool = True,
    overlapping: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Hadamard deviation: overlapping, or non-overlapped with ``overlapping=False``
    (allantools' ``hdev``, Stable32's Normal Hadamard)."""
    return _run(
        jl.hdev,
        data,
        m,
        ci=ci,
        remove_drift=remove_drift,
        overlapping=overlapping,
        confidence=float(confidence),
    )


def mhdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    remove_drift: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Hadamard deviation."""
    return _run(jl.mhdev, data, m, ci=ci, remove_drift=remove_drift, confidence=float(confidence))


def htdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    remove_drift: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Hadamard time deviation σ_x,HT(τ) = (τ/√(10/3))·Mod σ_y,MH(τ). Units of seconds."""
    return _run(jl.htdev, data, m, ci=ci, remove_drift=remove_drift, confidence=float(confidence))


def pdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Parabolic deviation σ_PDEV(τ) (Vernotte). PDEV(τ₀) ≡ ADEV(τ₀)."""
    return _run(jl.pdev, data, m, ci=ci, confidence=float(confidence))


def theoh(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """ThêoH composite σ_y(τ) (SP1065 §5.2.16): overlapping ADEV at short τ, ThêoBR beyond."""
    return _run(jl.theoh, data, m, ci=ci, confidence=float(confidence))


def mtie(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Maximum Time Interval Error (ITU-T G.810). Units of seconds.

    MTIE has no published EDF/CI model, so ``ci``/``confidence`` are accepted for
    signature uniformity but are no-ops (the result always has empty CI fields).
    """
    return _run(jl.mtie, data, m, ci=ci, confidence=float(confidence))


def tierms(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """RMS Time Interval Error (TIE rms, SP1065 §5.2.18). Units of seconds.

    Like ``mtie``, TIE rms has no published EDF/CI model, so ``ci``/``confidence``
    are accepted for signature uniformity but are no-ops (empty CI fields).
    """
    return _run(jl.tierms, data, m, ci=ci, confidence=float(confidence))


def totdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Total deviation (Howe / SP1065 mean-flip extension)."""
    return _run(jl.totdev, data, m, ci=ci, correct_bias=_bias(correct_bias), confidence=float(confidence))


def mtotdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Total deviation (Greenhall extension)."""
    return _run(jl.mtotdev, data, m, ci=ci, correct_bias=_bias(correct_bias), confidence=float(confidence))


def ttotdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Time Total deviation σ_x(τ) = (τ/√3)·Mod-Total σ_y(τ). Units of seconds."""
    return _run(jl.ttotdev, data, m, ci=ci, correct_bias=_bias(correct_bias), confidence=float(confidence))


def htotdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    remove_drift: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Hadamard Total deviation (Greenhall extension on y = diff(x))."""
    return _run(
        jl.htotdev,
        data,
        m,
        ci=ci,
        correct_bias=_bias(correct_bias),
        remove_drift=remove_drift,
        confidence=float(confidence),
    )


def theo1(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Thêo1 deviation (SP1065 §5.2.15) — extends the usable τ range to 75 % of the record span."""
    return _run(jl.theo1, data, m, ci=ci, correct_bias=_bias(correct_bias), confidence=float(confidence))


def mhtotdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    remove_drift: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Hadamard Total deviation (SigmaTau-original; Greenhall methodology)."""
    return _run(
        jl.mhtotdev,
        data,
        m,
        ci=ci,
        correct_bias=_bias(correct_bias),
        remove_drift=remove_drift,
        confidence=float(confidence),
    )


@_map_errors
def deviation(
    name: str,
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
    correct_bias: bool | str = True,
    remove_drift: bool = True,
) -> StabilityResult:
    """Run the estimator called ``name`` (any name from :func:`estimators`).

    ``correct_bias`` applies to the total family and ``remove_drift`` to the
    Hadamard family; each is ignored elsewhere, so one call shape works for
    every estimator.
    """
    return result_from_jl(
        jl.deviation(
            jl.Symbol(str(name)),
            to_jl_data(data),
            resolve_m(m),
            ci=bool(ci),
            confidence=float(confidence),
            correct_bias=_bias(correct_bias),
            remove_drift=bool(remove_drift),
        )
    )


def estimators() -> list[dict]:
    """The estimators :func:`deviation` accepts, one dict per estimator.

    Keys: ``name``, ``label``, ``family``, ``axis``, ``correct_bias`` and
    ``description``, as in the Julia registry.
    """
    return [
        {
            "name": str(e.name),
            "label": str(e.label),
            "family": str(e.family),
            "axis": str(e.axis),
            "correct_bias": bool(e.correct_bias),
            "description": str(e.description),
        }
        for e in jl.estimators()
    ]


__all__ = [
    "DEFAULT_CONFIDENCE",
    "deviation",
    "estimators",
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

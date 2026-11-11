"""The stability estimators.

Each takes a ``PhaseData`` or ``FrequencyData`` record and ``m``, a ``TauMode``
grid (``Octave`` by default) or a sequence of averaging factors, and returns a
``StabilityResult``. The docstring of each ends with the keywords that it
shares with the others: ``ci`` and ``noise``, ``remove_drift`` in the Hadamard
family, and ``correct_bias`` in the total estimators and ``theo1``.

``mtie`` and ``tierms`` report no interval; they accept ``ci`` and
``confidence`` and ignore them.
"""

from __future__ import annotations

import inspect
import textwrap
from collections.abc import Sequence

from ._convert import powerlawfit_to_jl, resolve_m, result_from_jl, to_jl_data, to_jl_syms
from ._julia import _map_errors, jl
from .grids import Octave, TauMode
from .types import FrequencyData, PhaseData, PowerLawFit, StabilityResult

DEFAULT_CONFIDENCE = float(jl.DEFAULT_CONFIDENCE)


# The keywords that the estimators share. `_estimator` adds to the docstring of
# each estimator those that it takes.
_KEYWORDS = {
    "ci": (
        "Identify the noise type at each τ and report the equivalent degrees of "
        "freedom (EDF) and a χ² confidence interval at ``confidence``. ``False`` "
        "computes the deviation alone."
    ),
    "correct_bias": (
        '``True`` removes the bias of the estimator. ``"published"`` keeps to the '
        'published definitions and tables, ``"stable32"`` does what Stable32 does '
        "for that estimator, and ``False`` returns the uncorrected value."
    ),
    "remove_drift": (
        "Identify the noise on the record less its frequency drift. ``False`` uses "
        "the record as given, as Stable32 does."
    ),
    "noise": (
        'The noise type that sets the EDF and the bias factor. ``"auto"`` identifies '
        'it at each τ. A name such as ``"WHFM"`` states one type for every τ, a '
        "sequence one per averaging factor, and a ``PowerLawFit`` gives its "
        "strongest term at each τ."
    ),
}


def _estimator(fn):
    """Add to the docstring of an estimator the shared keywords that it takes."""
    if fn.__doc__:
        entries = [
            f"``{name}``\n" + textwrap.indent(textwrap.fill(_KEYWORDS[name], 75), "    ")
            for name in inspect.signature(fn).parameters
            if name in _KEYWORDS
        ]
        fn.__doc__ = inspect.cleandoc(fn.__doc__) + "\n\n" + "\n".join(entries)
    return fn


def _noise(noise: str | Sequence[str] | PowerLawFit):
    """``noise`` for Julia: ``"auto"`` or one noise type as a Symbol, one type per
    averaging factor as a ``Vector{Symbol}``, or a ``PowerLawFit``."""
    if isinstance(noise, str):
        return jl.Symbol(noise)
    return powerlawfit_to_jl(noise) if isinstance(noise, PowerLawFit) else to_jl_syms(noise)


def _bias(correct_bias: bool | str):
    """``correct_bias`` for Julia: a bool, or ``"published"`` / ``"stable32"`` as a Symbol."""
    return jl.Symbol(correct_bias) if isinstance(correct_bias, str) else bool(correct_bias)


@_map_errors
def _run(jlfn, data, m, **kwargs) -> StabilityResult:
    """Convert the arguments, call the Julia estimator and convert its result."""
    return result_from_jl(jlfn(to_jl_data(data), resolve_m(m), **kwargs))


@_estimator
def adev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    overlapping: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
    dead_time: float = 0.0,
) -> StabilityResult:
    """Overlapping Allan deviation σ_y(τ) (IEEE Std 1139-2022, Eq. (A.23); SP1065 §5.2.4).

    ``overlapping=False`` gives the non-overlapped estimator of SP1065 §5.2.2
    (Stable32's Normal Allan deviation, allantools' ``adev``), which needs a
    record without gaps.

    ``dead_time`` is for a ``FrequencyData`` record whose samples each average
    over ``tau0`` and start every ``tau0 + dead_time`` seconds. The deviation
    is then divided by √(B2·B3) of :func:`bias_b2` and :func:`bias_b3`, which
    gives σ_y(m·τ₀) (SP1065 §5.15). ``noise`` has to state the noise type, and
    the record must have no gaps. No other estimator takes it.
    """
    return _run(
        jl.adev,
        data,
        m,
        ci=ci,
        overlapping=overlapping,
        noise=_noise(noise),
        confidence=float(confidence),
        dead_time=float(dead_time),
    )


@_estimator
def mdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Allan deviation Mod σ_y(τ) (IEEE Std 1139-2022, Eq. (A.27); SP1065 §5.2.5)."""
    return _run(jl.mdev, data, m, ci=ci, noise=_noise(noise), confidence=float(confidence))


@_estimator
def tdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Time deviation σ_x(τ) = (τ/√3)·Mod σ_y(τ), in seconds (IEEE Std 1139-2022,
    Eq. (A.29); SP1065 §5.2.6)."""
    return _run(jl.tdev, data, m, ci=ci, noise=_noise(noise), confidence=float(confidence))


@_estimator
def hdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    remove_drift: bool = True,
    overlapping: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Overlapping Hadamard deviation (SP1065 §5.2.9), which ignores linear
    frequency drift.

    ``overlapping=False`` gives the non-overlapped estimator of SP1065 §5.2.8
    (Stable32's Normal Hadamard deviation, allantools' ``hdev``), which needs a
    record without gaps.
    """
    return _run(
        jl.hdev,
        data,
        m,
        ci=ci,
        remove_drift=remove_drift,
        overlapping=overlapping,
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
def mhdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    remove_drift: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Hadamard deviation Mod σ_H(τ) (SP1065 §5.2.10; Bregni and Jmoda 2006)."""
    return _run(
        jl.mhdev,
        data,
        m,
        ci=ci,
        remove_drift=remove_drift,
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
def htdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    remove_drift: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Hadamard time deviation (τ/√(10/3))·Mod σ_H(τ), in seconds: the time-stability
    form of ``mhdev``, as ``tdev`` is of ``mdev``.

    The statistic is defined in this package: SP1065 and IEEE Std 1139-2022
    define no Hadamard time deviation.
    """
    return _run(
        jl.htdev,
        data,
        m,
        ci=ci,
        remove_drift=remove_drift,
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
def pdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Parabolic deviation (Vernotte, Lenczner, Bourgeois and Rubiola 2016). At m = 1
    it is the Allan deviation by definition. The record must have no gaps."""
    return _run(jl.pdev, data, m, ci=ci, noise=_noise(noise), confidence=float(confidence))


@_estimator
def theoh(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """ThêoH: the overlapping Allan deviation up to a fifth of the record length and
    ThêoBR, the bias-removed Thêo1, beyond it (Howe 2006; SP1065 §5.2.16).

    The crossover is that of SP1065; Howe (2006) places it at a tenth. The
    record must have no gaps.
    """
    return _run(jl.theoh, data, m, ci=ci, noise=_noise(noise), confidence=float(confidence))


def mtie(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Maximum time interval error, in seconds: the largest peak-to-peak phase
    excursion within any window of m + 1 samples (ITU-T G.810; SP1065 §5.2.17).

    MTIE is a peak statistic and has no χ² interval. ``ci`` and ``confidence``
    are accepted, so that every estimator can be called alike, and ignored.
    """
    return _run(jl.mtie, data, m, ci=ci, confidence=float(confidence))


def tierms(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """RMS time interval error, in seconds: the rms phase change over τ
    (SP1065 §5.2.18; ITU-T G.810).

    It includes any frequency offset in the record. Like ``mtie`` it has no χ²
    interval, and ``ci`` and ``confidence`` are accepted and ignored.
    """
    return _run(jl.tierms, data, m, ci=ci, confidence=float(confidence))


@_estimator
def totdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Total deviation: an estimator of the Allan deviation with more degrees of
    freedom at long τ, from the record extended at each end by reflection about
    the end sample (Greenhall et al. 1999; SP1065 §5.2.11). The record must have
    no gaps."""
    return _run(
        jl.totdev,
        data,
        m,
        ci=ci,
        correct_bias=_bias(correct_bias),
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
def mtotdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified total deviation, the total form of ``mdev`` (Howe and Vernotte 1999;
    Howe and Peppler 2001; SP1065 §5.2.12). The record must have no gaps."""
    return _run(
        jl.mtotdev,
        data,
        m,
        ci=ci,
        correct_bias=_bias(correct_bias),
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
def ttotdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Time total deviation (τ/√3)·Mod σ_total(τ), in seconds: the time-stability
    form of ``mtotdev`` (Howe and Peppler 2001; SP1065 §5.2.13)."""
    return _run(
        jl.ttotdev,
        data,
        m,
        ci=ci,
        correct_bias=_bias(correct_bias),
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
def htotdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    remove_drift: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Hadamard total deviation, the total form of ``hdev`` (Howe et al. 2005;
    SP1065 §5.2.14). The record must have no gaps."""
    return _run(
        jl.htotdev,
        data,
        m,
        ci=ci,
        correct_bias=_bias(correct_bias),
        remove_drift=remove_drift,
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
def theo1(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Thêo1 deviation, with its bias removed by default (ThêoBR): an estimator of
    the Allan deviation out to 75 % of the record length (Howe and Peppler 2003;
    SP1065 §5.2.15 and §5.2.16).

    ``m`` holds even Thêo1 factors, and the result reports τ = 0.75·m·τ₀.
    ``correct_bias=False`` returns Thêo1 uncorrected. The record must have no
    gaps.
    """
    return _run(
        jl.theo1,
        data,
        m,
        ci=ci,
        correct_bias=_bias(correct_bias),
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
def mhtotdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    ci: bool = True,
    correct_bias: bool | str = True,
    remove_drift: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
    confidence: float = DEFAULT_CONFIDENCE,
) -> StabilityResult:
    """Modified Hadamard total deviation, the total form of ``mhdev``.

    The construction, its bias and its EDF are defined in this package.
    Vernotte and Howe (2000) outline a related total form of the same
    variance. The record must have no gaps.
    """
    return _run(
        jl.mhtotdev,
        data,
        m,
        ci=ci,
        correct_bias=_bias(correct_bias),
        remove_drift=remove_drift,
        noise=_noise(noise),
        confidence=float(confidence),
    )


@_estimator
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
    noise: str | Sequence[str] | PowerLawFit = "auto",
    dead_time: float = 0.0,
) -> StabilityResult:
    """Run the estimator called ``name``, any name from :func:`estimators`.

    ``correct_bias`` reaches the total estimators and ``theo1``, and
    ``remove_drift`` the Hadamard family. Each is ignored elsewhere, so one
    call serves every estimator. A ``dead_time`` other than zero is an error
    for any estimator but ``"adev"``. The non-overlapped forms are reached
    through :func:`adev` and :func:`hdev`.
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
            noise=_noise(noise),
            dead_time=float(dead_time),
        )
    )


def estimators() -> list[dict]:
    """The estimators that :func:`deviation` and :func:`stability` run, one dict each.

    Keys: ``name``, ``label``, ``family`` (``"allan"``, ``"total"``, ``"theo"``
    or ``"time"``), ``axis`` (``"sigma_y"``, or ``"sigma_x"`` for seconds),
    ``correct_bias`` (whether that keyword applies) and ``description``. The
    dynamic deviations ``dadev`` and ``dhdev`` are not listed.
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

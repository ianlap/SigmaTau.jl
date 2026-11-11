"""Frequency drift models and estimators."""

from __future__ import annotations

from ._convert import driftfit_from_jl, to_jl_data
from ._julia import _map_errors, jl
from .types import DriftFit, FrequencyData, PhaseData

DRIFT_MODELS = tuple(str(m) for m in jl.seval("SigmaTau.DRIFT_MODELS"))
DRIFT_METHODS = tuple(str(m) for m in jl.seval("SigmaTau.DRIFT_METHODS"))


@_map_errors
def drift_fit(
    data: PhaseData | FrequencyData, *, model: str = "linear", method: str = "frequency_fit"
) -> DriftFit:
    """Fit of a frequency drift model to the record ``data``.

    ``model`` is ``"linear"`` (a + b·t), ``"quadratic"`` (a + b·t + c·t²),
    ``"log"`` (a·ln(1 + b·t) + c, the stabilization of an oscillator) or
    ``"diffusion"`` (a + b·√(t + c)), with t in seconds from the first sample.
    The linear, log and diffusion models are those of SP1065 §5.22 and of
    Stable32's Drift function.

    ``method`` chooses the estimator of the linear model; the other models
    take ``"frequency_fit"``:

    - ``"frequency_fit"``: least-squares line through the frequency (least
      variance on white FM).
    - ``"phase_fit"``: least-squares parabola through the phase (white PM).
    - ``"bisection"``: the mean frequency of the second half of the record
      less that of the first, over half its length; on phase data, the
      three-point fit.
    - ``"second_difference"``: mean of the second differences of the phase
      (random-walk FM).
    - ``"four_point"``: the estimator of Greenhall (1997) from four values of
      the summed phase, close to the least variance on every power-law noise.

    These are the drift methods of SP1065 §5.21 and §5.22 and of Stable32's
    Drift function, and the four-point estimator, which Stable32 does not
    have.

    A phase record is differenced to frequency first, except for
    ``"phase_fit"`` and ``"four_point"``, and ``nan`` samples are skipped.
    ``"four_point"`` needs a record without gaps, and so does ``"phase_fit"``
    of a frequency record. The result gives the parameters, their standard
    errors and the drift rate at the end of the record, per day. The standard
    errors of the least-squares fits assume white residuals and are too small
    otherwise, as on a record with flicker FM or random-walk FM.
    ``"bisection"`` and ``"four_point"`` have none (``nan``): the Julia
    docstring gives their variances in terms of the deviation of the record.
    """
    return driftfit_from_jl(
        jl.drift_fit(to_jl_data(data), model=jl.Symbol(str(model)), method=jl.Symbol(str(method)))
    )


__all__ = ["drift_fit", "DRIFT_MODELS", "DRIFT_METHODS"]

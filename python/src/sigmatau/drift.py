"""Frequency drift models."""

from __future__ import annotations

from ._convert import driftfit_from_jl, to_jl_data
from ._julia import _map_errors, jl
from .types import DriftFit, FrequencyData, PhaseData

DRIFT_MODELS = tuple(str(m) for m in jl.seval("SigmaTau.DRIFT_MODELS"))


@_map_errors
def drift_fit(data: PhaseData | FrequencyData, *, model: str = "linear") -> DriftFit:
    """Least-squares fit of a frequency drift model to the frequency of a record.

    ``model`` is ``"linear"`` (a + b·t), ``"quadratic"`` (a + b·t + c·t²),
    ``"log"`` (a·ln(1 + b·t) + c, the stabilization of an oscillator) or
    ``"diffusion"`` (a + b·√(t + c)), with t in seconds from the first sample.
    The linear, log and diffusion models are those of SP1065 §5.22 and of
    Stable32's Drift function.

    A phase record is differenced to frequency first, and ``nan`` samples are
    skipped. The result gives the parameters, their standard errors and the
    drift rate at the end of the record, per day. The standard errors assume
    white residuals and are too small otherwise, as on a record with flicker
    FM or random-walk FM.
    """
    return driftfit_from_jl(jl.drift_fit(to_jl_data(data), model=jl.Symbol(str(model))))


__all__ = ["drift_fit", "DRIFT_MODELS"]

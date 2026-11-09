"""Frequency drift models — mirrors ``drift.jl``."""

from __future__ import annotations

from ._convert import driftfit_from_jl, to_jl_data
from ._julia import _map_errors, jl
from .types import DriftFit, FrequencyData, PhaseData

DRIFT_MODELS = tuple(str(m) for m in jl.seval("SigmaTau.DRIFT_MODELS"))


@_map_errors
def drift_fit(data: PhaseData | FrequencyData, *, model: str = "linear") -> DriftFit:
    """Fit a frequency drift model to a record, as Stable32's Drift function does.

    ``model`` is ``"linear"``, ``"quadratic"``, ``"log"`` (the MIL-PRF-55310
    aging model a·ln(1 + b·t) + c) or ``"diffusion"`` (a + b·√(t + c)). Phase
    records are differenced to frequency first; ``nan`` samples are skipped.
    The standard errors assume white residuals.
    """
    return driftfit_from_jl(jl.drift_fit(to_jl_data(data), model=jl.Symbol(str(model))))


__all__ = ["drift_fit", "DRIFT_MODELS"]

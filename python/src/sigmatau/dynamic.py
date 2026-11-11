"""Dynamic deviations: the stability of a record as a function of time.

``dadev`` and ``dhdev`` compute the overlapping Allan or Hadamard deviation
inside a window that slides along the record. The map σ_y(t, τ) shows changes
of stability with time, such as a step in noise level or a burst, that a
single σ_y(τ) curve averages away.
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import as_int, dynamic_from_jl, resolve_m, to_jl_data
from ._julia import _map_errors, jl
from .grids import Octave, TauMode
from .types import DynamicStabilityResult, FrequencyData, PhaseData


@_map_errors
def _run(jlfn, data, m, *, window, step) -> DynamicStabilityResult:
    """Convert the arguments, call the Julia estimator and convert its map."""
    kwargs = {}
    if window is not None:
        kwargs["window"] = as_int(window, "window")
    if step is not None:
        kwargs["step"] = as_int(step, "step")
    return dynamic_from_jl(jlfn(to_jl_data(data), resolve_m(m), **kwargs))


def dadev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    window: int | None = None,
    step: int | None = None,
) -> DynamicStabilityResult:
    """Dynamic Allan deviation σ_y(t, τ): the overlapping Allan deviation of a
    window of ``window`` samples that slides along the record in steps of
    ``step`` samples (Galleani and Tavella 2009; IEEE Std 1139-2022, Eq. (A.40)).

    ``window`` defaults to the largest power of two up to an eighth of the
    record, and at least 16, and ``step`` to half the window. A ``TauMode``
    grid is resolved against the window length. The record must have no gaps.

    The result holds one row per window and one column per τ, ``nan`` where the
    window is too short for the averaging factor. No confidence interval is
    reported: Galleani (2011) derives confidence surfaces for the dynamic Allan
    variance, which are not implemented here.
    """
    return _run(jl.dadev, data, m, window=window, step=step)


def dhdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    window: int | None = None,
    step: int | None = None,
) -> DynamicStabilityResult:
    """Dynamic Hadamard deviation: the overlapping Hadamard deviation of a sliding
    window, built as :func:`dadev` is (SP1065 §5.2.20).

    The third difference cancels linear frequency drift within each window, so
    the map of a drifting clock shows its noise and not its drift. Arguments
    and keywords are those of :func:`dadev`. No confidence interval is
    reported.
    """
    return _run(jl.dhdev, data, m, window=window, step=step)


__all__ = ["dadev", "dhdev"]

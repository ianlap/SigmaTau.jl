"""Dynamic (time-resolved) deviations: dadev/dhdev.

Each slides a fixed-length analysis window of ``window`` phase samples across
the record in increments of ``step`` samples and evaluates the corresponding
overlapping kernel (Allan / Hadamard) inside every window, producing the
σ_y(t, τ) map. A ``TauMode`` grid is resolved by Julia against the *window*
length, so pass ``m`` straight through. There is no EDF model for the
time-resolved map, hence no ``ci``/``confidence`` kwargs.
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import dynamic_from_jl, resolve_m, to_jl_data
from ._julia import _map_errors, jl
from .grids import Octave, TauMode
from .types import DynamicStabilityResult, FrequencyData, PhaseData


@_map_errors
def _run(jlfn, data, m, *, window, step) -> DynamicStabilityResult:
    """Shared call path: convert inputs, call Julia, down-convert the dynamic map."""
    kwargs = {"window": int(window)}
    if step is not None:
        kwargs["step"] = int(step)
    return dynamic_from_jl(jlfn(to_jl_data(data), resolve_m(m), **kwargs))


def dadev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    window: int,
    step: int | None = None,
) -> DynamicStabilityResult:
    """Dynamic Allan deviation σ_y(t, τ); rows are windows (centered at t), cols are τ.

    No CI model exists for the time-resolved map, so there are no ci/confidence kwargs.
    """
    return _run(jl.dadev, data, m, window=window, step=step)


def dhdev(
    data: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
    *,
    window: int,
    step: int | None = None,
) -> DynamicStabilityResult:
    """Dynamic Hadamard deviation σ_y(t, τ); rows are windows (centered at t), cols are τ.

    Rejects within-window linear drift. No CI model exists for the time-resolved
    map, so there are no ci/confidence kwargs.
    """
    return _run(jl.dhdev, data, m, window=window, step=step)


__all__ = ["dadev", "dhdev"]

"""Real-time streaming deviations — mirrors ``streaming.jl``."""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import result_from_jl, to_jl_f64, to_jl_ints
from ._julia import _map_errors, jl
from .types import StabilityResult

_push = jl.seval("push!")
_append = jl.seval("append!")


class StreamingStability:
    """Running ADEV, MDEV, HDEV or MHDEV over phase samples fed one at a time.

    ``dev`` is ``"adev"``, ``"mdev"``, ``"hdev"`` or ``"mhdev"``; ``tau0`` the
    sample interval in seconds; ``m_values`` the averaging factors. ``push`` one
    phase sample (seconds) or ``extend`` with a chunk, and read the estimate over
    everything so far with ``snapshot``, which equals the batch call with
    ``ci=False`` on the same samples. The record is not stored: each sample costs
    O(len(m_values)) and memory is a ring buffer of a few ``max(m_values)``.
    """

    @_map_errors
    def __init__(self, dev: str, tau0: float, m_values: Sequence[int]) -> None:
        self._jl = jl.StreamingStability(jl.Symbol(str(dev)), float(tau0), to_jl_ints(m_values))

    def push(self, x: float) -> StreamingStability:
        _push(self._jl, float(x))
        return self

    def extend(self, xs: Sequence[float]) -> StreamingStability:
        _append(self._jl, to_jl_f64(xs))
        return self

    @_map_errors
    def snapshot(self) -> StabilityResult:
        return result_from_jl(jl.snapshot(self._jl))

    def __len__(self) -> int:
        return int(jl.nsamples(self._jl))

    def __repr__(self) -> str:
        return str(jl.repr(self._jl))


__all__ = ["StreamingStability"]

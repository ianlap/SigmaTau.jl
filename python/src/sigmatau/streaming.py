"""Deviations updated as samples arrive."""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import result_from_jl, to_jl_f64, to_jl_ints
from ._julia import _map_errors, jl
from .types import StabilityResult

_push = jl.seval("push!")
_append = jl.seval("append!")


class StreamingStability:
    """ADEV, MDEV, HDEV or MHDEV updated as phase samples arrive, without keeping
    the record.

    ``dev`` is ``"adev"``, ``"mdev"``, ``"hdev"`` or ``"mhdev"``, ``tau0`` the
    sample interval in seconds and ``m_values`` the averaging factors. ``push``
    adds one phase sample, in seconds, and ``extend`` a sequence of them.
    ``snapshot`` gives the deviation of the samples so far, which equals that
    of the estimator on them, as a ``StabilityResult`` without noise type, EDF
    or interval.

    The running sums are those of Dobrogowski and Kasznia (2007), who give
    them for ADEV and TDEV; the same scheme is applied here to the third
    differences of HDEV and MHDEV. A sample costs a fixed number of operations
    per averaging factor, and the last k·max(m) + 1 samples are held, with
    k = 2 for ADEV, 3 for MDEV and HDEV and 4 for MHDEV.

    Every sample must be finite. ``push`` and ``extend`` raise ``ValueError``
    on one that is not; ``extend`` adds the samples before it.
    """

    @_map_errors
    def __init__(self, dev: str, tau0: float, m_values: Sequence[int]) -> None:
        self._jl = jl.StreamingStability(
            jl.Symbol(str(dev)), float(tau0), to_jl_ints(m_values, "an averaging factor")
        )

    @_map_errors
    def push(self, x: float) -> StreamingStability:
        """Add one phase sample ``x``, in seconds, and return the accumulator."""
        _push(self._jl, float(x))
        return self

    @_map_errors
    def extend(self, xs: Sequence[float]) -> StreamingStability:
        """Add the phase samples ``xs``, a sequence in seconds, and return the accumulator."""
        _append(self._jl, to_jl_f64(xs))
        return self

    @_map_errors
    def snapshot(self) -> StabilityResult:
        """The deviation of the samples added so far."""
        return result_from_jl(jl.snapshot(self._jl))

    def __len__(self) -> int:
        return int(jl.nsamples(self._jl))

    def __repr__(self) -> str:
        return str(jl.repr(self._jl))


__all__ = ["StreamingStability"]

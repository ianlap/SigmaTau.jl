"""ITU-T wander masks for MTIE and TDEV.

``TELECOM_MASKS`` names the masks. Twelve limit the wander that a clock
generates in locked mode: ``G811`` and ``G811_1`` (a primary reference clock
and an enhanced one), ``G812_I`` and ``G812_II_III`` (a node clock of type I,
and of types II and III), ``G813_opt1`` and ``G813_opt2`` (an SDH equipment
clock), ``G8262_opt1``, ``G8262_opt2`` and ``G8262_1`` (a synchronous
equipment clock, with the limits of G.813, and an enhanced one),
``G8272_PRTC_A``, ``G8272_PRTC_B`` and ``G8272_1`` (a primary reference time
clock of class A and of class B, and an enhanced one). ``G823_PRC``,
``G823_SSU``, ``G823_SEC`` and ``G823_PDH`` are the network limits for wander
at those synchronization interfaces. The Julia documentation of
``TELECOM_MASKS`` gives the edition and the tables of each recommendation.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np

from ._convert import np_f64, result_to_jl, to_jl_f64
from ._julia import _map_errors, jl
from .types import StabilityResult

TELECOM_MASKS = tuple(str(m) for m in jl.TELECOM_MASKS)


class MaskCheck(NamedTuple):
    """Result of :func:`mask_check`: the name of the mask, the observation intervals
    of the result, the limit at each (``nan`` where the mask does not reach),
    whether the record passed, and the smallest ratio of limit to value."""

    mask: str
    tau: np.ndarray
    limit: np.ndarray
    passed: bool
    margin: float


@_map_errors
def telecom_mask(mask: str, stat: str, tau: float | Sequence[float]) -> float | np.ndarray:
    """Limit of ``mask``, one of ``TELECOM_MASKS``, for ``stat``, ``"mtie"`` or
    ``"tdev"``, at the observation interval ``tau`` in seconds, or at each of a
    sequence of them.

    The limit is in seconds, and ``nan`` outside the range of τ that the
    recommendation covers.
    """
    if np.ndim(tau) == 0:
        return float(jl.telecom_mask(jl.Symbol(mask), jl.Symbol(stat), float(tau)))
    return np_f64(jl.telecom_mask(jl.Symbol(mask), jl.Symbol(stat), to_jl_f64(tau)))


@_map_errors
def mask_check(r: StabilityResult, mask: str) -> MaskCheck:
    """Compare ``r``, an ``mtie`` or ``tdev`` result, with the telecom ``mask``, one of
    ``TELECOM_MASKS``.

    ``passed`` is true when every value is at or under its limit, and ``margin``
    is the smallest ratio of limit to value: a margin above 1 passes. Points
    outside the range of the mask, and points without a finite value, are not
    compared. With none left, ``margin`` is ``nan`` and ``passed`` is false.

    The clock recommendations measure MTIE and TDEV through a 10 Hz low-pass
    filter. No such filter is applied here: the record is compared as it was
    measured.
    """
    c = jl.mask_check(result_to_jl(r), jl.Symbol(mask))
    passed = bool(getattr(c, "pass"))  # a Python keyword, so not c.pass
    return MaskCheck(str(c.mask), np_f64(c.tau), np_f64(c.limit), passed, float(c.margin))


__all__ = ["TELECOM_MASKS", "MaskCheck", "telecom_mask", "mask_check"]

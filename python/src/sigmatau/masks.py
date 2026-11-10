"""ITU-T wander masks for MTIE and TDEV — mirrors ``masks.jl``."""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np

from ._convert import np_f64, result_to_jl, to_jl_f64
from ._julia import _map_errors, jl
from .types import StabilityResult

TELECOM_MASKS = tuple(str(m) for m in jl.TELECOM_MASKS)


class MaskCheck(NamedTuple):
    """``mask_check`` result: the limit at each τ, pass or fail, and the smallest limit/value."""

    mask: str
    tau: np.ndarray
    limit: np.ndarray
    passed: bool
    margin: float


@_map_errors
def telecom_mask(mask: str, stat: str, tau: float | Sequence[float]) -> float | np.ndarray:
    """The limit in seconds of ``mask`` (one of ``TELECOM_MASKS``) for ``stat``
    (``"mtie"`` or ``"tdev"``) at observation interval ``tau`` in seconds, or an
    array for an array of ``tau``; ``nan`` outside the range the recommendation covers.
    """
    if np.ndim(tau) == 0:
        return float(jl.telecom_mask(jl.Symbol(mask), jl.Symbol(stat), float(tau)))
    return np_f64(jl.telecom_mask(jl.Symbol(mask), jl.Symbol(stat), to_jl_f64(tau)))


@_map_errors
def mask_check(r: StabilityResult, mask: str) -> MaskCheck:
    """Compare an ``mtie`` or ``tdev`` result with a telecom mask.

    ``passed`` is true when every value inside the mask's τ range is at or under
    its limit; ``margin`` is the smallest ``limit / value`` (above 1 passes).
    """
    c = jl.mask_check(result_to_jl(r), jl.Symbol(mask))
    passed = bool(getattr(c, "pass"))  # a Python keyword, so not c.pass
    return MaskCheck(str(c.mask), np_f64(c.tau), np_f64(c.limit), passed, float(c.margin))


__all__ = ["TELECOM_MASKS", "MaskCheck", "telecom_mask", "mask_check"]

"""Grids of averaging factors.

A ``TauMode`` names the spacing of a grid, and ``tau_values`` gives its
averaging factors for a record length and an estimator. Julia builds the
grid, so it is the one the estimators use.
"""

from __future__ import annotations

import operator
from enum import Enum

from ._julia import _map_errors, jl


class TauMode(Enum):
    """Spacing of the averaging-factor grid: every integer (``AllTaus``), or
    successive factors in the ratio 2 (``Octave``, the default of every
    estimator), √2 (``HalfOctave``), 2^(1/4) (``QuarterOctave``), 10
    (``Decade``) or √10 (``HalfDecade``).

    The geometric grids are rounded to integers, duplicates are dropped, and
    each grid stops at the largest factor the estimator supports for the
    record length.
    """

    AllTaus = "AllTaus"
    Octave = "Octave"
    HalfOctave = "HalfOctave"
    QuarterOctave = "QuarterOctave"
    Decade = "Decade"
    HalfDecade = "HalfDecade"


# The instances by name, so that a call reads `adev(pd, Octave)` as in Julia.
AllTaus = TauMode.AllTaus
Octave = TauMode.Octave
HalfOctave = TauMode.HalfOctave
QuarterOctave = TauMode.QuarterOctave
Decade = TauMode.Decade
HalfDecade = TauMode.HalfDecade

# The Julia enum instances, resolved once by name.
_JL_TAUMODE = {mode: getattr(jl, mode.name) for mode in TauMode}


@_map_errors
def tau_values(mode: TauMode, n: int, kernel: str) -> list[int]:
    """Averaging factors m (τ = m·τ₀) of the grid ``mode`` for a record of ``n``
    phase samples, up to the largest that the estimator ``kernel`` (``"adev"``,
    ``"mdev"``, …) supports.

    For ``"theo1"``, which takes even factors only, odd entries are rounded up
    to the next even value. A ``"theoh"`` grid is in units of τ = m·τ₀.
    """
    ms = jl.tau_values(_JL_TAUMODE[mode], operator.index(n), jl.Symbol(str(kernel)))
    return [int(m) for m in ms]


__all__ = [
    "TauMode",
    "AllTaus",
    "Octave",
    "HalfOctave",
    "QuarterOctave",
    "Decade",
    "HalfDecade",
    "tau_values",
]

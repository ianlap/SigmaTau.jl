"""Averaging-factor grids — thin delegation to Julia's ``grids.jl``.

``TauMode`` selects the averaging-factor spacing; ``tau_values`` resolves it to
a concrete integer grid bounded by each kernel's algorithmic m-max. Resolution
happens in Julia, so the grids are guaranteed identical to Julia's —
including the theo1 even-rounding and the theoh/tierms m-max rules.
"""

from __future__ import annotations

from enum import Enum

from ._julia import _map_errors, jl


class TauMode(Enum):
    """Averaging-factor grid spacing. Mirrors the Julia ``@enum TauMode``."""

    AllTaus = "AllTaus"
    Octave = "Octave"
    HalfOctave = "HalfOctave"
    QuarterOctave = "QuarterOctave"
    Decade = "Decade"
    HalfDecade = "HalfDecade"


# Module-level singletons so callers write `adev(pd, Octave)`, matching Julia.
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
    """Averaging factors ``m`` (τ = m·τ₀) for ``mode``, bounded by ``kernel``'s m-max."""
    return [int(m) for m in jl.tau_values(_JL_TAUMODE[mode], int(n), jl.Symbol(str(kernel)))]


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

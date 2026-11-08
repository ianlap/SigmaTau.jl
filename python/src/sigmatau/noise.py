"""Power-law noise identification and synthesis — mirrors ``noise.jl``.

``identify_noise`` classifies the dominant power-law noise type per averaging
factor (the lag-1 ACF / B1-ratio classifier the EDF/CI machinery consumes);
``noise_gen`` synthesizes a record with a chosen power-law frequency spectrum,
drawing from Julia's RNG so a seeded draw matches the oracle bit-for-bit.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from ._convert import data_from_jl, np_syms, to_jl_data, to_jl_dict_int_f64, to_jl_f64, to_jl_ints
from ._julia import _map_errors, jl
from .types import FrequencyData, PhaseData

_xoshiro = jl.seval("Random.Xoshiro")


@_map_errors
def identify_noise(
    x: PhaseData | FrequencyData | Sequence[float],
    m_values: Sequence[int],
    *,
    dmin: int = 0,
    dmax: int = 2,
    remove_drift: bool = True,
) -> np.ndarray:
    """Identify the dominant power-law noise type at each averaging factor.

    ``x`` is a ``PhaseData`` or ``FrequencyData`` record (frequency is integrated
    to phase first) or a raw vector of phase samples. Returns an object array of
    symbol strings (``"WHPM"``/``"FLPM"``/``"WHFM"``/``"FLFM"``/``"RWFM"`` or
    ``"unknown"``). ``remove_drift`` (default ``True``, matching allantools)
    quadratically detrends each decimated subseries before classification.
    """
    return np_syms(
        jl.identify_noise(
            to_jl_data(x) if isinstance(x, (PhaseData, FrequencyData)) else to_jl_f64(x),
            to_jl_ints(m_values),
            dmin=int(dmin),
            dmax=int(dmax),
            remove_drift=bool(remove_drift),
        )
    )


@_map_errors
def noise_gen(
    kind: type,
    n: int,
    tau0: float = 1.0,
    *,
    sigma1: Mapping[int, float] | None = None,
    h: Mapping[int, float] | None = None,
    rng: int | None = None,
) -> PhaseData | FrequencyData:
    """Synthesize a length-``n`` clock record with a power-law frequency spectrum.

    ``kind`` is the class ``PhaseData`` (integrated to phase) or ``FrequencyData``
    (raw ``y``). Specify the mixture by ``sigma1[α] = σ_y(τ₀)`` *or* ``h[α] = h_α``
    (not both), with α ∈ {−4,…,2} (white PM down to random-run FM). ``sigma1``
    components are rescaled to hit their target σ exactly for the drawn
    realization; ``h`` components have spectrum h_α (sin(πfτ₀)/(πτ₀))^α.
    Components of different α are independent.

    ``rng`` is an int seed (drawn through Julia's ``Random.Xoshiro``) or ``None``
    (Julia's global RNG). Seed for reproducible output.
    """
    if kind is PhaseData:
        jl_kind = jl.PhaseData
    elif kind is FrequencyData:
        jl_kind = jl.FrequencyData
    else:
        raise TypeError(f"noise_gen: kind must be PhaseData or FrequencyData, got {kind!r}")

    kwargs = {}
    if sigma1 is not None:
        kwargs["sigma1"] = to_jl_dict_int_f64(sigma1)
    if h is not None:
        kwargs["h"] = to_jl_dict_int_f64(h)
    if rng is not None:
        if isinstance(rng, bool) or not isinstance(rng, int):
            raise TypeError("rng: pass an int seed or None; sigmatau draws from Julia's RNG")
        kwargs["rng"] = _xoshiro(int(rng))

    return data_from_jl(jl.noise_gen(jl_kind, int(n), float(tau0), **kwargs))


__all__ = ["identify_noise", "noise_gen"]

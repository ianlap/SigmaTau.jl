"""Identification and simulation of power-law noise.

``identify_noise`` names the dominant power-law noise type at each averaging
factor, ``noise_id_nnls`` fits an amplitude for every type over the whole range
of τ, and ``noise_gen`` simulates a record with a power-law frequency spectrum.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from ._convert import (
    data_from_jl,
    np_f64,
    np_syms,
    resolve_m,
    to_jl_data,
    to_jl_dict_int_f64,
    to_jl_f64,
    to_jl_ints,
)
from ._julia import _map_errors, jl
from .grids import TauMode
from .types import FrequencyData, PhaseData

_xoshiro = jl.seval("Random.Xoshiro")


@_map_errors
def identify_noise(
    x: PhaseData | FrequencyData | Sequence[float],
    m_values: TauMode | Sequence[int],
    *,
    dmin: int = 0,
    dmax: int = 2,
    remove_drift: bool = True,
    b1_length: str = "averaged",
    b1: bool = True,
) -> np.ndarray:
    """Dominant power-law noise type of a record at each averaging factor:
    ``"WHPM"``, ``"FLPM"``, ``"WHFM"``, ``"FLFM"`` or ``"RWFM"``, and with
    ``dmax=3`` also ``"FWFM"`` (flicker-walk FM) and ``"RRFM"`` (random-run FM).

    The test is the lag-1 autocorrelation of the phase decimated by m (Riley
    and Greenhall 2004; SP1065 §5.5.6). Where fewer than 30 decimated points
    remain it gives way to the B₁ ratio, with R(n) to separate the PM types
    (SP1065 §5.5.2), as in Stable32. The limits between the FM labels of a
    decimated record, and the conditions on a PM label at long τ, are those of
    this package. A label is always returned.

    ``x`` is a ``PhaseData`` or ``FrequencyData`` record, which is integrated to
    phase, or a sequence of phase samples. ``m_values`` holds the averaging
    factors; with a record it may be a ``TauMode`` grid.

    ``dmin`` and ``dmax`` are the fewest and the most differences taken: use
    ``dmax=2`` for the Allan family and ``dmax=3`` for the Hadamard family, as
    the estimators do. ``remove_drift`` removes a least-squares quadratic from
    each decimated series before the test, as allantools does. ``b1_length``
    is the record length at which the expected B₁ of each noise type is
    evaluated: ``"averaged"``, the number of m-averaged frequency values, or
    ``"record"``, the number of frequency samples of the whole record, which
    follows Stable32's labels more closely. ``b1=False`` leaves the B₁ ratio
    out, so that an FM result from fewer than 30 points is the type found at
    the longest averaging factor that has 30.
    """
    return np_syms(
        jl.identify_noise(
            to_jl_data(x) if isinstance(x, (PhaseData, FrequencyData)) else to_jl_f64(x),
            resolve_m(m_values),
            dmin=int(dmin),
            dmax=int(dmax),
            remove_drift=bool(remove_drift),
            b1_length=jl.Symbol(str(b1_length)),
            b1=bool(b1),
        )
    )


@_map_errors
def noise_id_nnls(
    data: PhaseData | FrequencyData,
    m_values: Sequence[int] | None = None,
    *,
    deviation: str = "mhdev",
    alphas: Sequence[int] | None = None,
    weights: str | Sequence[float] = "edf",
) -> dict:
    """Amount of each power-law noise type in a record, from one nonnegative
    least-squares fit of the measured variance to σ²(τ) ≈ Σᵢ aᵢ·(τ/τ₀)^μ(αᵢ) over
    the whole range of τ.

    Where :func:`identify_noise` names the dominant type at each τ, this fit
    gives an amplitude for every type in the basis. The basis is the
    asymptotic power law of each type; for coefficients h_α use
    :func:`powerlaw_fit`, which fits the exact response of each component.

    ``m_values`` defaults to the octave grid of the estimator from m = 4.
    ``deviation`` is ``"mdev"``, ``"hdev"``, ``"mhdev"`` or ``"mhtotdev"``.
    ``alphas`` holds the exponents in the basis, by default every type the
    estimator converges for: 2 to −2 for ``"mdev"`` and 2 to −4 for the
    Hadamard family. ``weights`` is ``"edf"``, ``"relative"`` or a sequence of
    weights.

    Returns a dict with ``deviation``, ``tau``, ``var`` (the measured
    variance), ``fit`` (the model at each τ), ``alphas``, ``noise_types``,
    ``amplitudes`` (the variance of each component at τ₀), ``contributions``
    (an array, τ by type), ``dominant`` (the largest contribution at each τ),
    ``residual`` and ``weights``.
    """
    kwargs = {
        "deviation": jl.Symbol(str(deviation)),
        "weights": jl.Symbol(weights) if isinstance(weights, str) else to_jl_f64(weights),
    }
    if alphas is not None:
        kwargs["alphas"] = to_jl_ints(alphas)
    args = (to_jl_data(data),) if m_values is None else (to_jl_data(data), to_jl_ints(m_values))
    r = jl.noise_id_nnls(*args, **kwargs)
    return {
        "deviation": str(r.deviation),
        "tau": np_f64(r.tau),
        "var": np_f64(r.var),
        "fit": np_f64(r.fit),
        "alphas": np.array(r.alphas, dtype=np.int64),
        "noise_types": np_syms(r.noise_types),
        "amplitudes": np_f64(r.amplitudes),
        "contributions": np.array(r.contributions, dtype=np.float64),
        "dominant": np_syms(r.dominant),
        "residual": float(r.residual),
        "weights": np_f64(r.weights),
    }


@_map_errors
def noise_gen(
    kind: type,
    n: int,
    tau0: float = 1.0,
    *,
    sigma1: Mapping[int | str, float] | None = None,
    h: Mapping[int | str, float] | None = None,
    rng: int | None = None,
) -> PhaseData | FrequencyData:
    """Simulated clock record of ``n`` samples at sample interval ``tau0`` whose
    fractional-frequency spectrum is a sum of power laws, S_y(f) = Σ h_α f^α.

    ``kind`` is the class of the record returned, ``PhaseData`` or
    ``FrequencyData``. Each component is white Gaussian noise passed through
    the fractional differencing filter of Kasdin and Walter (1992), and the
    components are independent.

    Give ``sigma1`` or ``h``, not both, keyed by the integer α, from 2 (white
    PM) to −4 (random-run FM), or by the label, as in ``{"WHFM": 1e-12}``.
    ``sigma1`` holds the Allan deviation σ_y(τ₀) of each component, which is
    rescaled so that its overlapping Allan deviation at τ₀ equals the request
    exactly, in every realization. ``h`` holds the coefficient h_α of each:
    the component has S_y(f) = h_α·(sin(πfτ₀)/(πτ₀))^α, which is h_α f^α well
    below the Nyquist frequency.

    ``rng`` is an integer seed for Julia's ``Xoshiro`` generator, for a
    reproducible record: it gives the record that Julia gives with
    ``Xoshiro(seed)``. With ``None`` Julia's global generator is used.
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


__all__ = ["identify_noise", "noise_id_nnls", "noise_gen"]

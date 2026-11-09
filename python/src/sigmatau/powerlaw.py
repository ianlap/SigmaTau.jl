"""Power-law noise model S_y(f) = Σ h_α f^α — mirrors ``powerlaw.jl``.

``powerlaw_fit`` fits the h_α to a record's σ(τ) or to a spectrum;
``powerlaw_sigma`` and ``powerlaw_psd`` turn a model back into σ(τ) or a
spectrum.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from ._convert import (
    np_f64,
    powerlaw_from_jl,
    powerlaw_to_jl,
    resolve_m,
    to_jl_data,
    to_jl_f64,
    to_jl_ints,
)
from ._julia import _map_errors, jl
from .types import FrequencyData, PhaseData, PowerLawFit, SpectralResult

_jl_spectral = jl.seval(
    "(t, f, p, u, np, no, w) -> SpectralResult(Symbol(t), f, p, Symbol(u), np, no, Symbol(w))"
)


def _spectral_to_jl(s: SpectralResult):
    return _jl_spectral(
        s.spectral_type,
        to_jl_f64(s.freq),
        to_jl_f64(s.psd),
        s.units,
        int(s.nperseg),
        int(s.noverlap),
        s.taper,
    )


@_map_errors
def powerlaw_fit(
    data: PhaseData | FrequencyData | SpectralResult | Sequence[float],
    second: Sequence[int] | Sequence[float] | None = None,
    *,
    deviation: str = "mhdev",
    quantity: str = "Sy",
    alphas: Sequence[int] = (2, 1, 0, -1, -2),
    f_carrier: float | None = None,
) -> PowerLawFit:
    """Fit the coefficients h_α of S_y(f) = Σ h_α f^α by nonnegative least squares.

    Three forms, as in Julia:

    - ``powerlaw_fit(record[, m], deviation="mhdev", alphas=...)`` fits a record's
      σ(τ) (``"adev"``, ``"mdev"``, ``"hdev"`` or ``"mhdev"``) on the octave grid
      or at the averaging factors ``m``.
    - ``powerlaw_fit(spectral_result, alphas=..., f_carrier=None)`` fits a
      ``psd_y`` / ``psd_x`` / ``ssb_phase_noise`` result.
    - ``powerlaw_fit(f, S, quantity="Sy", alphas=..., f_carrier=None)`` fits
      spectrum points, for instance a datasheet's ℒ(f) (``quantity="L"``, with the
      carrier in Hz).
    """
    kwargs = {"alphas": to_jl_ints(alphas)}
    if isinstance(data, (PhaseData, FrequencyData)):
        args = (to_jl_data(data),) if second is None else (to_jl_data(data), resolve_m(second))
        return powerlaw_from_jl(
            jl.powerlaw_fit(*args, deviation=jl.Symbol(str(deviation)), **kwargs)
        )
    if f_carrier is not None:
        kwargs["f_carrier"] = float(f_carrier)
    if isinstance(data, SpectralResult):
        return powerlaw_from_jl(jl.powerlaw_fit(_spectral_to_jl(data), **kwargs))
    if second is None:
        raise TypeError("powerlaw_fit: pass a record, a SpectralResult, or frequencies and values")
    return powerlaw_from_jl(
        jl.powerlaw_fit(
            to_jl_f64(data), to_jl_f64(second), quantity=jl.Symbol(str(quantity)), **kwargs
        )
    )


@_map_errors
def powerlaw_sigma(
    h: PowerLawFit | Mapping[int | str, float],
    tau: Sequence[float],
    *,
    deviation: str = "adev",
    tau0: float | None = None,
) -> np.ndarray:
    """σ(τ) the model ``h`` gives in ``deviation`` for a record sampled at ``tau0``.

    ``h`` is a ``PowerLawFit`` (whose ``tau0`` is the default) or an ``{α: h_α}``
    mapping keyed by the integer α or a noise-type name (``"WHFM"``). Every τ
    must be a whole multiple of ``tau0``.
    """
    if tau0 is None:
        if not isinstance(h, PowerLawFit):
            raise ValueError("powerlaw_sigma: pass `tau0`, the sample interval in seconds")
        tau0 = h.tau0
    return np_f64(
        jl.powerlaw_sigma(
            powerlaw_to_jl(h), to_jl_f64(tau), deviation=jl.Symbol(str(deviation)), tau0=float(tau0)
        )
    )


@_map_errors
def powerlaw_psd(
    h: PowerLawFit | Mapping[int | str, float],
    f: Sequence[float],
    *,
    quantity: str = "Sy",
    f_carrier: float | None = None,
) -> np.ndarray:
    """The model ``h`` as a spectrum at ``f`` (Hz).

    ``quantity`` is ``"Sy"`` [1/Hz], ``"Sx"`` [s²/Hz] or ``"L"`` [dBc/Hz, with ``f_carrier``].
    """
    kwargs = {"quantity": jl.Symbol(str(quantity))}
    if f_carrier is not None:
        kwargs["f_carrier"] = float(f_carrier)
    return np_f64(jl.powerlaw_psd(powerlaw_to_jl(h), to_jl_f64(f), **kwargs))


__all__ = ["powerlaw_fit", "powerlaw_sigma", "powerlaw_psd"]

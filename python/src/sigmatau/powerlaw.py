"""The power-law noise model S_y(f) = Σ h_α f^α.

``powerlaw_fit`` fits the coefficients h_α to the deviation of a record or to a
spectrum, and ``powerlaw_sigma`` and ``powerlaw_psd`` turn a model back into a
deviation or a spectrum.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from ._convert import (
    np_f64,
    powerlaw_from_jl,
    resolve_m,
    result_to_jl,
    to_jl_data,
    to_jl_f64,
    to_jl_ints,
)
from ._julia import _map_errors, jl
from .ensemble import ClockModel, NoiseFit, _noise_model_to_jl
from .types import FrequencyData, PhaseData, PowerLawFit, SpectralResult, StabilityResult

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
    data: PhaseData | FrequencyData | StabilityResult | SpectralResult | Sequence[float],
    second: Sequence[int] | Sequence[float] | None = None,
    *,
    deviation: str = "mhdev",
    quantity: str = "Sy",
    alphas: Sequence[int] = (2, 1, 0, -1, -2),
    f_carrier: float | None = None,
) -> PowerLawFit:
    """Coefficients h_α of the power-law model S_y(f) = Σ h_α f^α, fitted by
    nonnegative least squares to the deviation of a record or to a spectrum.

    Four forms, as in Julia, of which ``data`` is the first positional argument
    and ``second`` the averaging factors ``m`` or the values ``S``:

    - ``powerlaw_fit(record[, m], deviation="mhdev", alphas=...)`` computes the
      deviation of the record (``"adev"``, ``"mdev"``, ``"hdev"`` or
      ``"mhdev"``) on the octave grid, or at the averaging factors ``m``, and
      fits it. Each component enters with the variance it has in that
      estimator on a record sampled at τ₀, which is exact from m = 1.
    - ``powerlaw_fit(result, alphas=...)`` fits an existing ADEV, MDEV, HDEV or
      MHDEV result, such as the Allan deviation of one clock from
      :func:`groslambert`.
    - ``powerlaw_fit(spectral_result, alphas=..., f_carrier=None)`` fits a
      result of :func:`psd_y`, :func:`psd_x` or :func:`ssb_phase_noise`.
    - ``powerlaw_fit(f, S, quantity="Sy", alphas=..., f_carrier=None)`` fits
      spectrum points, such as the phase-noise table of a data sheet
      (``quantity="L"``, with the carrier frequency in Hz).

    The Allan pair takes α from 2 to −2. The Hadamard pair takes it to −4 and
    is insensitive to linear frequency drift. A coefficient cannot be
    negative, so a component that the record does not have comes out small
    and positive: compare ``fit`` with ``y`` before trusting a small
    coefficient.
    """
    kwargs = {"alphas": to_jl_ints(alphas)}
    if isinstance(data, (PhaseData, FrequencyData)):
        args = (to_jl_data(data),) if second is None else (to_jl_data(data), resolve_m(second))
        return powerlaw_from_jl(
            jl.powerlaw_fit(*args, deviation=jl.Symbol(str(deviation)), **kwargs)
        )
    if isinstance(data, StabilityResult):
        return powerlaw_from_jl(jl.powerlaw_fit(result_to_jl(data), **kwargs))
    if f_carrier is not None:
        kwargs["f_carrier"] = float(f_carrier)
    if isinstance(data, SpectralResult):
        return powerlaw_from_jl(jl.powerlaw_fit(_spectral_to_jl(data), **kwargs))
    if second is None:
        raise TypeError("powerlaw_fit: pass a record, a result, or frequencies and values")
    return powerlaw_from_jl(
        jl.powerlaw_fit(
            to_jl_f64(data), to_jl_f64(second), quantity=jl.Symbol(str(quantity)), **kwargs
        )
    )


@_map_errors
def powerlaw_sigma(
    h: PowerLawFit | NoiseFit | ClockModel | Mapping[int | str, float],
    tau: Sequence[float],
    *,
    deviation: str = "adev",
    tau0: float | None = None,
) -> np.ndarray:
    """Deviation σ(τ) that the power-law model ``h`` gives in an estimator, for a
    record sampled at ``tau0``.

    The model is taken as the spectrum of a sampled record,
    h_α·(sin(πfτ₀)/(πτ₀))^α up to the Nyquist frequency, which is the spectrum
    that :func:`noise_gen` produces, and it is integrated through the transfer
    function of the estimator. For ADEV at m ≫ 1 this is the usual conversion
    (SP1065 §7.1; IEEE Std 1139-2022, Eq. (B.4)).

    ``h`` is a ``PowerLawFit``, a ``NoiseFit``, a ``ClockModel``, or a mapping
    of h_α keyed by the integer α or by the noise type, as in
    ``{"WHFM": 2e-24}``. ``tau`` holds averaging times in seconds, each a whole
    multiple of ``tau0``. ``deviation`` is ``"adev"``, ``"mdev"``, ``"hdev"`` or
    ``"mhdev"``. A fit or a clock model carries its own ``tau0``; a mapping
    needs the keyword.
    """
    kwargs = {"deviation": jl.Symbol(str(deviation))}
    if tau0 is None and isinstance(h, PowerLawFit):
        tau0 = h.tau0
    if tau0 is not None:
        kwargs["tau0"] = float(tau0)
    elif not isinstance(h, (NoiseFit, ClockModel)):
        raise ValueError("powerlaw_sigma: pass `tau0`, the sample interval in seconds")
    return np_f64(jl.powerlaw_sigma(_noise_model_to_jl(h), to_jl_f64(tau), **kwargs))


@_map_errors
def powerlaw_psd(
    h: PowerLawFit | NoiseFit | ClockModel | Mapping[int | str, float],
    f: Sequence[float],
    *,
    quantity: str = "Sy",
    f_carrier: float | None = None,
) -> np.ndarray:
    """The power-law model ``h`` as a spectrum at the Fourier frequencies ``f``, in Hz.

    ``h`` is a ``PowerLawFit``, a ``NoiseFit``, a ``ClockModel``, or a mapping
    of h_α keyed by α or by noise type. ``quantity`` selects the spectrum
    (IEEE Std 1139-2022, Eqs. (A.13), (A.15) and (A.16)): ``"Sy"``, S_y(f) in
    1/Hz; ``"Sx"``, S_y(f)/(2πf)² in s²/Hz; ``"Sphi"``, S_φ(f) in rad²/Hz; or
    ``"L"``, ℒ(f) in dBc/Hz. The last two need the carrier frequency
    ``f_carrier``, in Hz.
    """
    kwargs = {"quantity": jl.Symbol(str(quantity))}
    if f_carrier is not None:
        kwargs["f_carrier"] = float(f_carrier)
    return np_f64(jl.powerlaw_psd(_noise_model_to_jl(h), to_jl_f64(f), **kwargs))


__all__ = ["powerlaw_fit", "powerlaw_sigma", "powerlaw_psd"]

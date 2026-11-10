"""Spectral densities ``psd_y`` / ``psd_x`` / ``ssb_phase_noise`` — mirrors ``spectral.jl``.

One-sided Welch estimates at sample rate ``fs = 1/τ₀``, returned as a
``SpectralResult``. ``PhaseData`` and ``FrequencyData`` are accepted by every
estimator (converted internally by Julia as needed).
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import spectral_from_jl, to_jl_data, to_jl_f64
from ._julia import _map_errors, jl
from .types import FrequencyData, PhaseData, SpectralResult


def _welch_kwargs(nperseg: int | None, noverlap: int | None, taper: str) -> dict:
    """Welch keyword arguments for Julia: ``None`` segment sizes are omitted."""
    kwargs = {"taper": jl.Symbol(str(taper))}
    if nperseg is not None:
        kwargs["nperseg"] = int(nperseg)
    if noverlap is not None:
        kwargs["noverlap"] = int(noverlap)
    return kwargs


@_map_errors
def psd_y(
    data: PhaseData | FrequencyData,
    *,
    nperseg: int | None = None,
    noverlap: int | None = None,
    taper: str = "hann",
) -> SpectralResult:
    """One-sided PSD of fractional frequency S_y(f) [1/Hz] (Welch, fs=1/τ₀)."""
    return spectral_from_jl(jl.psd_y(to_jl_data(data), **_welch_kwargs(nperseg, noverlap, taper)))


@_map_errors
def psd_x(
    data: PhaseData | FrequencyData,
    *,
    nperseg: int | None = None,
    noverlap: int | None = None,
    taper: str = "hann",
) -> SpectralResult:
    """One-sided PSD of the phase residual S_x(f) [s²/Hz]: the Welch S_y(f) over the
    first-difference response 4 sin²(πfτ₀)/τ₀². NaN at f = 0."""
    return spectral_from_jl(jl.psd_x(to_jl_data(data), **_welch_kwargs(nperseg, noverlap, taper)))


@_map_errors
def ssb_phase_noise(
    data: PhaseData | FrequencyData,
    *,
    f_carrier: float,
    nperseg: int | None = None,
    noverlap: int | None = None,
    taper: str = "hann",
) -> SpectralResult:
    """Single-sideband phase noise ℒ(f) [dBc/Hz] at carrier ``f_carrier`` (IEEE 1139)."""
    kwargs = _welch_kwargs(nperseg, noverlap, taper)
    kwargs["f_carrier"] = float(f_carrier)
    return spectral_from_jl(jl.ssb_phase_noise(to_jl_data(data), **kwargs))


@_map_errors
def integrated_jitter(
    spectrum: SpectralResult | Sequence[float],
    S: Sequence[float] | None = None,
    *,
    f_carrier: float,
    quantity: str = "L",
    band: tuple[float, float] | None = None,
) -> dict:
    """Integrated phase jitter and residual FM over ``band`` (Hz), NIST SP1065 §5.2.19.

    Takes a ``ssb_phase_noise`` / ``psd_x`` / ``psd_y`` result, or spectrum points
    ``(f, S)`` such as a phase-noise analyzer's ℒ(f) in dBc/Hz (``quantity="L"``,
    or ``"Sx"``, ``"Sy"``). Returns ``phase`` (rms rad), ``time`` (rms s) and
    ``residual_fm`` (rms Hz), integrating the spectrum as straight lines on log-log axes.
    """
    kwargs = {"f_carrier": float(f_carrier)}
    if band is not None:
        kwargs["band"] = to_jl_f64(band)
    if isinstance(spectrum, SpectralResult):
        s = jl.SpectralResult(
            jl.Symbol(spectrum.spectral_type),
            to_jl_f64(spectrum.freq),
            to_jl_f64(spectrum.psd),
            jl.Symbol(spectrum.units),
            int(spectrum.nperseg),
            int(spectrum.noverlap),
            jl.Symbol(spectrum.taper),
        )
        j = jl.integrated_jitter(s, **kwargs)
    else:
        kwargs["quantity"] = jl.Symbol(str(quantity))
        j = jl.integrated_jitter(to_jl_f64(spectrum), to_jl_f64(S), **kwargs)
    return {"phase": float(j.phase), "time": float(j.time), "residual_fm": float(j.residual_fm)}


__all__ = ["psd_y", "psd_x", "ssb_phase_noise", "integrated_jitter"]

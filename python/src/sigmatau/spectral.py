"""Spectral densities ``psd_y`` / ``psd_x`` / ``ssb_phase_noise`` — mirrors ``spectral.jl``.

One-sided Welch estimates at sample rate ``fs = 1/τ₀``, returned as a
``SpectralResult``. ``PhaseData`` and ``FrequencyData`` are accepted by every
estimator (converted internally by the oracle as needed).
"""

from __future__ import annotations

from ._convert import spectral_from_jl, to_jl_data
from ._julia import _map_errors, jl
from .types import FrequencyData, PhaseData, SpectralResult


def _welch_kwargs(nperseg: int | None, noverlap: int | None, taper: str) -> dict:
    """Welch keyword arguments for the oracle: ``None`` segment sizes are omitted."""
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
    """One-sided PSD of the phase residual S_x(f) [s²/Hz] (Welch, fs=1/τ₀)."""
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


__all__ = ["psd_y", "psd_x", "ssb_phase_noise"]

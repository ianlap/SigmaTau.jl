"""Spectral densities.

One-sided Welch estimates at the sample rate 1/τ₀, returned as a
``SpectralResult``. Each estimator takes a ``PhaseData`` or ``FrequencyData``
record. A missing sample makes the whole spectrum ``nan``: fill gaps first
with ``fillgaps``.
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import as_int, np_f64, spectral_from_jl, to_jl_data, to_jl_f64
from ._julia import _map_errors, jl
from .types import FrequencyData, PhaseData, SpectralResult


def _welch_kwargs(nperseg: int | None, noverlap: int | None, taper: str) -> dict:
    """Welch keywords for Julia; a segment size of ``None`` is left out."""
    kwargs = {"taper": jl.Symbol(str(taper))}
    if nperseg is not None:
        kwargs["nperseg"] = as_int(nperseg, "nperseg")
    if noverlap is not None:
        kwargs["noverlap"] = as_int(noverlap, "noverlap")
    return kwargs


@_map_errors
def psd_y(
    data: PhaseData | FrequencyData,
    *,
    nperseg: int | None = None,
    noverlap: int | None = None,
    taper: str = "hann",
) -> SpectralResult:
    """One-sided power spectral density of the fractional frequency, S_y(f) in 1/Hz
    (IEEE Std 1139-2022, Eq. (A.8)), estimated by Welch's method (Welch 1967).

    ``data`` is the record, and a phase record is differenced to frequency
    first. ``nperseg`` is the segment length in samples, by default the
    largest power of two up to a quarter of the record and at least 16.
    ``noverlap`` is the number of samples that successive segments share,
    half a segment by default. ``taper`` is ``"hann"``, ``"hamming"`` or
    ``"rectangular"``.

    The mean removed from each segment is the one weighted by the taper, so
    the estimate is zero at f = 0. ``scipy.signal.welch`` with
    ``detrend="constant"`` removes the plain mean; with the Hann and Hamming
    tapers the two estimates differ at f = 0 and in the first bin and are
    equal from the second bin on. The first bin shares the main lobe of the
    taper with the mean that was removed and should not be read. A linear
    frequency drift is not removed and raises the lowest bins: take it out
    first with ``detrend``.
    """
    return spectral_from_jl(jl.psd_y(to_jl_data(data), **_welch_kwargs(nperseg, noverlap, taper)))


@_map_errors
def psd_x(
    data: PhaseData | FrequencyData,
    *,
    nperseg: int | None = None,
    noverlap: int | None = None,
    taper: str = "hann",
) -> SpectralResult:
    """One-sided power spectral density of the time fluctuations, S_x(f) in s²/Hz
    (IEEE Std 1139-2022, Eq. (A.14)): the estimate of :func:`psd_y` divided by
    the response of the first difference, 4·sin²(πfτ₀)/τ₀².

    The value at f = 0 is ``nan``. ``data``, ``noverlap`` and ``taper`` are
    those of :func:`psd_y`, and so is ``nperseg``, except that its default is
    the largest power of two up to a quarter of the number of phase samples,
    and at least 16.
    """
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
    """Phase noise ℒ(f) in dBc/Hz of a signal at the carrier frequency
    ``f_carrier``, in Hz, from its phase or frequency record ``data``.

    The time spectrum of :func:`psd_x` gives the phase spectrum
    S_φ(f) = (2πν₀)²·S_x(f), and the phase noise is half of it,
    ℒ(f) = S_φ(f)/2, reported as 10·log₁₀ ℒ(f) (IEEE Std 1139-2022,
    Eqs. (A.13) and (A.15)). The result runs from the first bin to the Nyquist
    frequency: the bin at f = 0 is dropped. ``nperseg``, ``noverlap`` and
    ``taper`` are those of :func:`psd_x`.
    """
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
    """Phase jitter, time jitter and residual FM of a signal at the carrier
    frequency ``f_carrier``, in Hz: its phase spectrum integrated over ``band``,
    a pair of offset frequencies in Hz, or over the whole spectrum.

    ``spectrum`` is a result of :func:`ssb_phase_noise`, :func:`psd_x` or
    :func:`psd_y`, or the frequencies of spectrum points whose values are in
    ``S``, such as the ℒ(f) table of a phase-noise analyzer. ``quantity`` says
    what ``S`` holds: ``"L"`` (ℒ(f) in dBc/Hz), ``"Sphi"`` (S_φ(f) in rad²/Hz),
    ``"Sx"`` (S_x(f) in s²/Hz) or ``"Sy"`` (S_y(f) in 1/Hz).

    Returns a dict with ``phase`` (radians rms; IEEE Std 1139-2022,
    Eq. (A.17)), ``time`` (seconds rms) and ``residual_fm`` (hertz rms).
    Between neighboring points the spectrum is taken as a power law, a
    straight line on a log-log plot, and each segment is integrated in closed
    form (SP1065 §5.2.19).
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


@_map_errors
def spectrum_sigma(
    f: Sequence[float],
    S: Sequence[float],
    tau: Sequence[float],
    *,
    quantity: str = "L",
    f_carrier: float | None = None,
    deviation: str = "adev",
    band: tuple[float, float] | None = None,
) -> dict:
    """Deviation σ(τ) of the noise whose spectrum is given at the points
    ``(f, S)``, such as the ℒ(f) table of a phase-noise analyzer: S_y(f)
    integrated over the band of the table through the transfer function of
    the deviation (SP1065 §7, Eqs. (63) to (66)).

    ``f`` holds the Fourier frequencies in Hz, rising, and ``tau`` the
    averaging times in seconds. ``quantity`` says what ``S`` holds: ``"L"``
    (ℒ(f) in dBc/Hz), ``"Sphi"`` (S_φ(f) in rad²/Hz), ``"Sx"`` (S_x(f) in
    s²/Hz) or ``"Sy"`` (S_y(f) in 1/Hz); the first two need the carrier
    frequency ``f_carrier``, in Hz. ``deviation`` is ``"adev"``, ``"mdev"``,
    ``"hdev"`` or ``"mhdev"``, the modified ones in the limit of a sample
    interval far below τ. ``band`` is a pair of frequencies inside the table,
    by default all of it. Between neighboring points the spectrum is a power
    law, a straight line on a log-log plot.

    Returns a dict. ``dev`` is the deviation of the noise inside the band at
    each ``tau``. ``below`` and ``above`` are the variance that the spectrum
    would add if it went on below and above the band, on the slope of its
    lowest and of its highest decade, as fractions of ``dev**2``: ``inf``
    where that integral does not converge, as ``above`` does for the Allan
    and Hadamard deviations of white and flicker PM, which rise with the
    bandwidth of the measurement. ``band`` holds the limits of the integral.

    The conversion goes one way: a deviation does not determine the spectrum
    (Greenhall 1998). From σ(τ) to a spectrum, fit a power-law model with
    :func:`powerlaw_fit` and evaluate it with :func:`powerlaw_psd`.
    """
    kwargs = {"quantity": jl.Symbol(str(quantity)), "deviation": jl.Symbol(str(deviation))}
    if f_carrier is not None:
        kwargs["f_carrier"] = float(f_carrier)
    if band is not None:
        kwargs["band"] = to_jl_f64(band)
    r = jl.spectrum_sigma(to_jl_f64(f), to_jl_f64(S), to_jl_f64(tau), **kwargs)
    return {
        "dev": np_f64(r.dev),
        "below": np_f64(r.below),
        "above": np_f64(r.above),
        "band": tuple(float(v) for v in r.band),
    }


__all__ = ["psd_y", "psd_x", "ssb_phase_noise", "integrated_jitter", "spectrum_sigma"]

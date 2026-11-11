"""Estimators for two simultaneous records of one device.

With each record taken against its own reference or through its own channel,
the cross-spectrum and the cross Allan variance keep the noise that the two
share, and the noise in only one averages toward zero.
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import cross_spectral_from_jl, cross_stability_from_jl, resolve_m, to_jl_data
from ._julia import _map_errors, jl
from .grids import Octave, TauMode
from .spectral import _welch_kwargs
from .types import CrossSpectralResult, CrossStabilityResult, FrequencyData, PhaseData


@_map_errors
def cross_psd_y(
    a: PhaseData | FrequencyData,
    b: PhaseData | FrequencyData,
    *,
    nperseg: int | None = None,
    noverlap: int | None = None,
    taper: str = "hann",
) -> CrossSpectralResult:
    """Cross power spectral density of the fractional frequency of two simultaneous
    records of one device, in 1/Hz.

    The estimate is Welch's average of cross-periodograms (Welch 1967): that of
    :func:`psd_y` with the product of the two channels' transforms in place of
    the squared magnitude. Its real part is reported. Noise in only one channel
    averages toward zero and leaves the residual given in ``floor``.

    ``a`` and ``b`` are two ``PhaseData`` or two ``FrequencyData`` records with
    the same sample interval and length and no gaps. ``nperseg``, ``noverlap``
    and ``taper`` are those of :func:`psd_y`.
    """
    return cross_spectral_from_jl(
        jl.cross_psd_y(to_jl_data(a), to_jl_data(b), **_welch_kwargs(nperseg, noverlap, taper))
    )


@_map_errors
def cross_psd_x(
    a: PhaseData | FrequencyData,
    b: PhaseData | FrequencyData,
    *,
    nperseg: int | None = None,
    noverlap: int | None = None,
    taper: str = "hann",
) -> CrossSpectralResult:
    """Cross power spectral density of the time fluctuations of two simultaneous
    records of one device, in s²/Hz.

    It is the estimate of :func:`cross_psd_y`, its floor and the spectrum of
    each channel, each divided by the response of the first difference as in
    :func:`psd_x`. Every field is ``nan`` at f = 0.

    ``a``, ``b``, ``noverlap`` and ``taper`` are those of :func:`cross_psd_y`,
    and ``nperseg`` is counted and its default chosen as in :func:`psd_x`.
    """
    return cross_spectral_from_jl(
        jl.cross_psd_x(to_jl_data(a), to_jl_data(b), **_welch_kwargs(nperseg, noverlap, taper))
    )


@_map_errors
def cross_adev(
    a: PhaseData | FrequencyData,
    b: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
) -> CrossStabilityResult:
    """Cross Allan variance of two simultaneous records of one device: the
    overlapping Allan variance with the square of each second difference
    replaced by the product of the two channels' second differences.

    This is the two-sample covariance that Vernotte, Calosso and Rubiola (2016)
    call the Groslambert covariance. It estimates the Allan variance of the
    noise the channels share; noise in only one averages toward zero, to within
    the ``uncertainty`` given at each τ. The estimate can be negative:
    ``covar`` keeps the sign, and ``dev`` is ``nan`` there.

    ``a`` and ``b`` are two ``PhaseData`` or two ``FrequencyData`` records with
    the same sample interval and length and no gaps, and ``m`` is a
    ``TauMode`` grid or a sequence of averaging factors.
    """
    return cross_stability_from_jl(jl.cross_adev(to_jl_data(a), to_jl_data(b), resolve_m(m)))


__all__ = ["cross_psd_y", "cross_psd_x", "cross_adev"]

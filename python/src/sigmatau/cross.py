"""Two-channel (cross) estimators — mirrors ``cross.jl``.

Two simultaneous records of one device, each against its own reference or
through its own channel: the cross-spectrum and the cross Allan variance keep
the noise the channels share and average away the rest.
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
    """Cross PSD of the fractional frequency of two records of one device [1/Hz]."""
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
    """Cross PSD of the phase of two records of one device [s²/Hz]."""
    return cross_spectral_from_jl(
        jl.cross_psd_x(to_jl_data(a), to_jl_data(b), **_welch_kwargs(nperseg, noverlap, taper))
    )


@_map_errors
def cross_adev(
    a: PhaseData | FrequencyData,
    b: PhaseData | FrequencyData,
    m: TauMode | Sequence[int] = Octave,
) -> CrossStabilityResult:
    """Cross Allan variance of two simultaneous records of one device.

    The overlapping Allan variance with the product of the two channels'
    second differences in place of the square. ``a`` and ``b`` are two
    ``PhaseData`` or two ``FrequencyData`` with the same τ₀ and length and no
    gaps; ``m`` is a ``TauMode`` grid or a sequence of averaging factors.
    """
    return cross_stability_from_jl(jl.cross_adev(to_jl_data(a), to_jl_data(b), resolve_m(m)))


__all__ = ["cross_psd_y", "cross_psd_x", "cross_adev"]

"""Shared data records and result types — mirrors ``types.jl``.

``PhaseData``/``FrequencyData`` are the input records; ``StabilityResult``,
``StabilitySuite``, ``SpectralResult``, and ``DynamicStabilityResult`` are the
outputs. All are frozen dataclasses with NumPy array fields, matching the Julia
package's non-parametric ``Vector{Float64}`` design. This module is pure Python
— importing it does not boot Julia.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import numpy as np


def _as_samples(v: object, name: str) -> np.ndarray:
    """Promote a sample vector to a contiguous 1-D float64 array."""
    arr = np.ascontiguousarray(v, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"{name}: expected a 1-D sample vector, got {arr.ndim}-D")
    return arr


@dataclass(frozen=True)
class PhaseData:
    """Phase residuals ``x(t)`` sampled at uniform interval ``tau0`` (default 1.0 s)."""

    x: np.ndarray
    tau0: float = 1.0
    source: str = "user"

    def __post_init__(self) -> None:
        x = _as_samples(self.x, "PhaseData")
        if not self.tau0 > 0:
            raise ValueError(f"PhaseData: tau0 must be positive, got {self.tau0}")
        if x.size < 2:
            raise ValueError(f"PhaseData: need at least 2 phase samples, got {x.size}")
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "tau0", float(self.tau0))

    def __repr__(self) -> str:
        src = "" if self.source == "user" else f", source={self.source!r}"
        return f"PhaseData(N={self.x.size}, tau0={self.tau0} s{src})"


@dataclass(frozen=True)
class FrequencyData:
    """Fractional-frequency samples ``y(t)`` at uniform interval ``tau0`` (default 1.0 s)."""

    y: np.ndarray
    tau0: float = 1.0
    source: str = "user"

    def __post_init__(self) -> None:
        y = _as_samples(self.y, "FrequencyData")
        if not self.tau0 > 0:
            raise ValueError(f"FrequencyData: tau0 must be positive, got {self.tau0}")
        if y.size < 2:
            raise ValueError(f"FrequencyData: need at least 2 frequency samples, got {y.size}")
        object.__setattr__(self, "y", y)
        object.__setattr__(self, "tau0", float(self.tau0))

    def __repr__(self) -> str:
        src = "" if self.source == "user" else f", source={self.source!r}"
        return f"FrequencyData(N={self.y.size}, tau0={self.tau0} s{src})"


@dataclass(frozen=True)
class StabilityResult:
    """Unified return type for every stability calculation.

    ``noise_type``, ``ci_lower``, ``ci_upper``, and ``edf`` are empty arrays when
    the calculation was invoked with ``ci=False`` (and always for ``mtie``/
    ``tierms``); ``neff`` is always populated. ``tau0`` and ``n`` describe the
    record the result came from (``nan`` and ``0`` when unknown), and
    ``confidence`` is the level of the intervals (``None`` without them).
    """

    deviation_type: str
    tau: np.ndarray
    dev: np.ndarray
    noise_type: np.ndarray
    ci_lower: np.ndarray
    ci_upper: np.ndarray
    edf: np.ndarray
    neff: np.ndarray
    tau0: float = float("nan")
    n: int = 0
    confidence: float | None = None

    def __repr__(self) -> str:
        n = self.tau.size
        ci = "no CI" if self.edf.size == 0 else "with CI"
        if n == 0:
            return f"StabilityResult({self.deviation_type}, 0 pts, {ci})"
        return (
            f"StabilityResult({self.deviation_type}, {n} pts, "
            f"τ∈[{self.tau[0]}, {self.tau[-1]}] s, {ci})"
        )


@dataclass(frozen=True)
class StabilitySuite:
    """Ordered, name-indexable collection of ``StabilityResult``s + session metadata.

    Index by deviation name (``suite["adev"]``) or position (``suite[0]``).
    """

    results: tuple[StabilityResult, ...]
    data_kind: str
    tau0: float
    n: int
    confidence: float | None
    tau_mode: str

    def __len__(self) -> int:
        return len(self.results)

    def __iter__(self) -> Iterator[StabilityResult]:
        return iter(self.results)

    def keys(self) -> tuple[str, ...]:
        return tuple(r.deviation_type for r in self.results)

    def __contains__(self, key: str) -> bool:
        return any(r.deviation_type == key for r in self.results)

    def __getitem__(self, key: int | str) -> StabilityResult:
        if isinstance(key, str):
            for r in self.results:
                if r.deviation_type == key:
                    return r
            raise KeyError(key)
        return self.results[key]

    def __repr__(self) -> str:
        devs = ", ".join(self.keys())
        ci = "no CI" if self.confidence is None else f"CI@{self.confidence}"
        return (
            f"StabilitySuite({len(self.results)} devs [{devs}], "
            f"{self.data_kind}, N={self.n}, τ₀={self.tau0} s, {ci})"
        )


@dataclass(frozen=True)
class SpectralResult:
    """Unified return type for every spectral-density estimate (``Sy``/``Sx``/``L``).

    A flat, non-parametric record: the one-sided frequency grid, the estimated
    spectrum, and the Welch parameters that produced it.
    """

    spectral_type: str
    freq: np.ndarray
    psd: np.ndarray
    units: str
    nperseg: int
    noverlap: int
    taper: str

    def __repr__(self) -> str:
        n = self.freq.size
        rng = "" if n == 0 else f", f∈[{self.freq[0]}, {self.freq[-1]}] Hz"
        return (
            f"SpectralResult({self.spectral_type}, {n} bins{rng}, "
            f"{self.units}, nperseg={self.nperseg})"
        )


@dataclass(frozen=True)
class DynamicStabilityResult:
    """Return type of the dynamic (time-resolved) deviations ``dadev``/``dhdev``.

    ``dev[i, j]`` is the deviation of analysis window ``i`` (centered at
    ``t[i]``) at averaging interval ``tau[j]``; ``NaN`` wherever the window is
    too short to support the averaging factor. No noise-ID/EDF/CI machinery is
    carried — the literature gives no EDF model for the time-resolved map.
    """

    deviation_type: str
    t: np.ndarray
    tau: np.ndarray
    dev: np.ndarray
    window: int
    tau0: float

    def __repr__(self) -> str:
        nt, ntau = self.dev.shape
        t_rng = "" if nt == 0 else f", t∈[{self.t[0]}, {self.t[-1]}] s"
        tau_rng = "" if ntau == 0 else f", τ∈[{self.tau[0]}, {self.tau[-1]}] s"
        return (
            f"DynamicStabilityResult({self.deviation_type}, {nt}×{ntau} map"
            f"{t_rng}{tau_rng}, window={self.window})"
        )


@dataclass(frozen=True)
class DriftFit:
    """Return type of ``drift_fit``: a frequency drift model fitted to y(t).

    ``coef`` holds the model parameters with ``t`` in seconds from the first
    sample (``[a, b]`` for ``"linear"``, ``[a, b, c]`` otherwise; see the Julia
    docstring for each model's form), ``stderr`` their one-σ errors, ``rate``
    the drift dy/dt at the end of the record per day with ``rate_stderr``, and
    ``fit`` the model at each ``t``. The residuals are ``y - fit``.
    """

    model: str
    coef: np.ndarray
    stderr: np.ndarray
    rate: float
    rate_stderr: float
    t: np.ndarray
    y: np.ndarray
    fit: np.ndarray
    rms: float

    def __repr__(self) -> str:
        return f"DriftFit({self.model}, {self.t.size} pts, rate={self.rate:.4g}/day)"


@dataclass(frozen=True)
class PowerLawFit:
    """Return type of ``powerlaw_fit``: the coefficients h_α of S_y(f) = Σ h_α f^α.

    ``source`` names what was fitted (a deviation such as ``"mhdev"``, or a
    spectrum ``"Sy"``/``"Sx"``/``"L"``); ``x`` is τ in seconds or f in Hz
    accordingly, ``y`` the measured values and ``fit`` the model at each ``x``.
    ``tau0`` is the record's sample interval (``nan`` for a spectrum given as
    points).
    """

    alphas: np.ndarray
    h: np.ndarray
    source: str
    tau0: float
    x: np.ndarray
    y: np.ndarray
    fit: np.ndarray

    def __repr__(self) -> str:
        terms = ", ".join(f"h{a}={v:.3g}" for a, v in zip(self.alphas, self.h, strict=True))
        return f"PowerLawFit({self.source}, {self.x.size} pts, {terms})"


@dataclass(frozen=True)
class CrossSpectralResult:
    """Return type of ``cross_psd_y`` / ``cross_psd_x``.

    ``psd`` is the real part of the averaged cross-spectrum (the noise common to
    both channels) and ``floor`` its residual √(S_aa S_bb / K) after ``K``
    independent averages; ``psd_a`` and ``psd_b`` are the channels alone.
    """

    spectral_type: str
    freq: np.ndarray
    psd: np.ndarray
    floor: np.ndarray
    psd_a: np.ndarray
    psd_b: np.ndarray
    nseg: int
    nperseg: int
    noverlap: int
    taper: str

    def __repr__(self) -> str:
        return (
            f"CrossSpectralResult({self.spectral_type}, {self.freq.size} bins, "
            f"{self.nseg} averages, nperseg={self.nperseg})"
        )


@dataclass(frozen=True)
class CrossStabilityResult:
    """Return type of ``cross_adev``: the cross Allan variance of two channels.

    ``covar`` keeps its sign (negative where the shared noise is below the
    floor) and ``dev`` is √covar where positive, else ``nan``; ``uncertainty``
    is the 1σ statistical uncertainty of ``covar``; ``dev_a`` and ``dev_b`` are
    the channels' own Allan deviations.
    """

    tau: np.ndarray
    covar: np.ndarray
    dev: np.ndarray
    uncertainty: np.ndarray
    dev_a: np.ndarray
    dev_b: np.ndarray
    neff: np.ndarray
    tau0: float
    n: int

    def __repr__(self) -> str:
        return f"CrossStabilityResult({self.tau.size} τ, {self.n} samples per channel)"


__all__ = [
    "PhaseData",
    "FrequencyData",
    "StabilityResult",
    "StabilitySuite",
    "SpectralResult",
    "DynamicStabilityResult",
    "DriftFit",
    "PowerLawFit",
    "CrossSpectralResult",
    "CrossStabilityResult",
]

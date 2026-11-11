"""The data records and the result types.

``PhaseData`` and ``FrequencyData`` are the records that the estimators take,
and the other classes are what they return. All are frozen dataclasses with
NumPy array fields. Importing this module does not start Julia.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta

import numpy as np


def _as_samples(v: object, name: str) -> np.ndarray:
    """The samples as a contiguous one-dimensional array of float64."""
    arr = np.ascontiguousarray(v, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"{name}: expected a 1-D sample vector, got {arr.ndim}-D")
    return arr


@dataclass(frozen=True)
class PhaseData:
    """Phase samples x(t) in seconds at the uniform sample interval ``tau0``: the
    time deviation of a clock from its reference (IEEE Std 1139-2022, 3.1).

    ``nan`` marks a missing sample, and a record holds at least two samples.
    ``source`` is ``"user"`` for a record made from an array and the path of
    the file for a record that was read. ``start`` is the time of the first
    sample in seconds since 2000-01-01T00:00Z, or ``None``; read it as a
    ``datetime`` with :func:`start_time`.
    """

    x: np.ndarray
    tau0: float = 1.0
    source: str = "user"
    start: float | None = None

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

    def to_frequency(self) -> FrequencyData:
        """The record as fractional frequency, y[k] = (x[k+1] − x[k])/τ₀: N − 1
        samples (SP1065 §3.4). A missing phase sample makes the two frequency
        samples beside it missing."""
        return FrequencyData(np.diff(self.x) / self.tau0, self.tau0, self.source, self.start)


@dataclass(frozen=True)
class FrequencyData:
    """Fractional-frequency samples y(t) = (ν(t) − ν₀)/ν₀, dimensionless, at the
    uniform sample interval ``tau0`` (IEEE Std 1139-2022, 3.1). Each sample is
    taken as the average over its interval, with no dead time.

    ``nan``, ``source`` and ``start`` are as for :class:`PhaseData`.
    """

    y: np.ndarray
    tau0: float = 1.0
    source: str = "user"
    start: float | None = None

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

    def to_phase(self) -> PhaseData:
        """The record as phase, x[k+1] = x[k] + τ₀·y[k] from x[0] = 0: N + 1 samples
        (SP1065 §3.4). A missing frequency sample leaves every later phase
        unknown and raises ``ValueError``: fill it first with ``fillgaps``."""
        if np.isnan(self.y).any():
            raise ValueError(
                "to_phase: a missing frequency sample leaves the phase after it unknown; "
                "fill the gaps first (fillgaps)"
            )
        x = np.concatenate(([0.0], np.cumsum(self.y))) * self.tau0
        return PhaseData(x, self.tau0, self.source, self.start)


_EPOCH = datetime(2000, 1, 1, tzinfo=UTC)


def start_time(data: PhaseData | FrequencyData) -> datetime | None:
    """UTC time of the first sample of a record, to the millisecond, or ``None``
    when the record does not carry it. The readers take it from ISO 8601 or MJD
    timestamps and from the ``# start:`` line that ``save`` writes. Sample
    ``i``, counted from 0, is at ``start_time(data) + i * tau0``."""
    if data.start is None:
        return None
    return _EPOCH + timedelta(milliseconds=round(data.start * 1000))


@dataclass(frozen=True)
class StabilityResult:
    """Result of every stability estimator: one deviation per averaging time, with
    the noise type, the equivalent degrees of freedom and the confidence
    interval when the estimator was called with ``ci=True``.

    ``dev`` is dimensionless for the σ_y estimators and in seconds for
    ``tdev``, ``htdev``, ``ttotdev``, ``mtie`` and ``tierms``. ``ci_lower``,
    ``ci_upper`` and ``edf`` are empty for a call with ``ci=False``, and always
    for ``mtie`` and ``tierms``; ``noise_type`` is empty too, unless a bias
    correction needed the noise type. ``neff``, the number of terms averaged
    at each τ, is always filled. ``tau0`` and ``n`` are the sample interval
    and the number of phase samples of the record (``nan`` and 0 when
    unknown), and ``confidence`` is the level of the interval (``None``
    without one).
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
    """The results of one ``stability`` call: a ``StabilityResult`` per estimator,
    in the order requested, with the record and grid settings that produced
    them.

    ``suite["adev"]`` returns a result by estimator name and ``suite[0]`` by
    position. The suite iterates in order and supports ``keys``, ``in`` and
    ``len``. ``tau_mode`` is the name of the ``TauMode``, or ``"explicit"``.
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
        """The names of the estimators, in order."""
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
    """Result of the spectral-density estimates ``psd_y``, ``psd_x`` and
    ``ssb_phase_noise``: the one-sided frequency grid ``freq`` in Hz, the
    estimated density ``psd`` and the Welch settings that produced it.

    ``spectral_type`` is ``"Sy"``, ``"Sx"`` or ``"L"``, and ``units`` names the
    unit of ``psd``: ``"per_Hz"``, ``"s2_per_Hz"`` or ``"dBc_per_Hz"``.
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
    """Result of the dynamic deviations ``dadev`` and ``dhdev``: a map σ_y(t, τ) of
    the deviation computed inside a window of ``window`` phase samples that
    slides along the record.

    Row ``i`` of ``dev`` is the window centered at ``t[i]``, in seconds from
    the first sample, and column ``j`` the averaging time ``tau[j]``. An entry
    is ``nan`` where the window is too short for the averaging factor. The
    result carries no noise type, EDF or interval: Galleani (2011) derives
    confidence surfaces for the dynamic Allan variance, which are not
    implemented here.
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
    """Result of ``drift_fit``: a frequency drift model fitted to y(t).

    ``coef`` holds the parameters of the model, with t in seconds from the
    first sample: ``[a, b]`` for ``"linear"`` (a + b·t), and ``[a, b, c]`` for
    ``"quadratic"`` (a + b·t + c·t²), ``"log"`` (a·ln(1 + b·t) + c) and
    ``"diffusion"`` (a + b·√(t + c)). ``method`` names the estimator.
    ``stderr`` holds the standard errors of the parameters, ``rate`` the
    drift rate dy/dt at the end of the record, per day, and ``rate_stderr``
    its standard error; a standard error is ``nan`` where the method gives
    none. ``fit`` is the model at each ``t``, and ``rms`` the rms of the
    residuals ``y - fit``.
    """

    model: str
    method: str
    coef: np.ndarray
    stderr: np.ndarray
    rate: float
    rate_stderr: float
    t: np.ndarray
    y: np.ndarray
    fit: np.ndarray
    rms: float

    def __repr__(self) -> str:
        return f"DriftFit({self.model}, {self.method}, {self.t.size} pts, rate={self.rate:.4g}/day)"


@dataclass(frozen=True)
class PowerLawFit:
    """Result of ``powerlaw_fit``: the coefficients h_α of the power-law model
    S_y(f) = Σ h_α f^α, with S_y in 1/Hz, and the fit they give.

    ``alphas`` holds the exponents, highest first, and ``h`` the coefficient of
    each. ``source`` names what was fitted: a deviation (``"adev"``,
    ``"mdev"``, ``"hdev"`` or ``"mhdev"``) or a spectrum (``"Sy"``, ``"Sx"`` or
    ``"L"``). ``x`` is τ in seconds or f in Hz accordingly, ``y`` the measured
    values and ``fit`` the model at each ``x``. ``tau0`` is the sample
    interval of the record, ``nan`` for a spectrum given as points.
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
    """Result of ``cross_psd_y`` and ``cross_psd_x``.

    ``psd`` is the real part of the averaged cross-spectrum: the spectrum of
    the noise common to both channels. ``floor`` is the rms magnitude of the
    residual that the noise the channels do not share leaves in the average,
    √(S_aa·S_bb/K), with K the number of segments reduced for the correlation
    of overlapping ones (Welch 1967). A bin where ``psd`` is not well above
    the floor is not resolved. ``psd_a`` and ``psd_b`` are the spectra of the
    channels alone, and ``nseg`` is the number of segments averaged.
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
    """Result of ``cross_adev``: the cross Allan variance of two channels.

    ``covar`` keeps its sign, negative where the common noise is below the
    scatter of the estimate, and ``dev`` is √covar where it is positive and
    ``nan`` elsewhere. ``uncertainty`` is the standard uncertainty of
    ``covar``, from the covariance model of Vernotte and Lantz (2019):
    ``covar`` well above it is resolved. ``dev_a`` and ``dev_b`` are the Allan
    deviations of the channels alone, and ``neff`` is the number of second
    differences averaged at each τ.
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

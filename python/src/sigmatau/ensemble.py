"""Clock models, the Kalman filter, clock ensembles and timescales — mirrors
``clockmodel.jl``, ``noisefit.jl``, ``kalman.jl``, ``comparisons.jl`` and
``timescale.jl``.

Results that later calls build on (a filter run, a noise fit) keep their Julia
object, so a smoother or a holdover estimate runs on exactly what the filter
produced.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import juliacall
import numpy as np

from ._convert import (
    data_from_jl,
    np_f64,
    powerlaw_to_jl,
    result_to_jl,
    to_jl_data,
    to_jl_f64,
    to_jl_ints,
)
from ._julia import VecStr, _map_errors, jl
from .noise import _xoshiro
from .types import PhaseData, PowerLawFit, StabilityResult

_MatF64 = jl.seval("Matrix{Float64}")
_VecModel = jl.seval("Vector{ClockModel}")
_VecF64 = jl.seval("Vector{Float64}")
_add = jl.seval(
    "(f, noise, a, b) -> SigmaTau.add_component!(f, noise, a:b)"
)  # a range would come back as a Python range
_component = jl.seval("c -> (c.alpha, c.h, first(c.range), last(c.range))")


def _rng(seed: int | None) -> dict:
    if seed is None:
        return {}
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError("rng: pass an int seed or None; sigmatau draws from Julia's RNG")
    return {"rng": _xoshiro(seed)}


def _mat(a: Any) -> np.ndarray:
    return np.array(a, dtype=np.float64)


# ── Clock model ────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ClockModel:
    """Discrete clock model: phase, frequency and (with ``drift``) drift driven by
    white noise (Zucca & Tavella 2005).

    ``R`` is the white PM phase variance per sample (s²), ``q1`` the white FM
    (s), ``q2`` the random-walk FM (1/s) and ``q3`` the random-run FM (1/s³)
    diffusion coefficient.
    """

    tau0: float
    R: float = 0.0
    q1: float = 0.0
    q2: float = 0.0
    q3: float = 0.0
    drift: bool | None = None

    def __post_init__(self):
        if self.drift is None:
            object.__setattr__(self, "drift", self.q3 > 0)

    def _jl(self) -> Any:
        return jl.ClockModel(
            float(self.tau0),
            float(self.R),
            float(self.q1),
            float(self.q2),
            float(self.q3),
            bool(self.drift),
        )


def _model_from_jl(m: Any) -> ClockModel:
    return ClockModel(
        float(m.tau0), float(m.R), float(m.q1), float(m.q2), float(m.q3), bool(m.drift)
    )


@_map_errors
def clock_model(
    fit: PowerLawFit | NoiseFit | Mapping[int | str, float],
    *,
    tau0: float | None = None,
    drift: bool | None = None,
) -> ClockModel:
    """The ``ClockModel`` with the white PM, white FM, random-walk FM and random-run
    FM levels of a power-law model; flicker components have no finite-state form
    and are left out.

    ``fit`` is a ``PowerLawFit``, a ``NoiseFit`` or an ``{α: h_α}`` mapping (with
    ``tau0``).
    """
    kw = {} if drift is None else {"drift": bool(drift)}
    if isinstance(fit, NoiseFit):
        return _model_from_jl(jl.clock_model(fit._jl, **kw))
    if tau0 is None:
        if not isinstance(fit, PowerLawFit):
            raise ValueError("clock_model: pass `tau0`, the sample interval in seconds")
        tau0 = fit.tau0
    return _model_from_jl(jl.clock_model(powerlaw_to_jl(fit), tau0=float(tau0), **kw))


@_map_errors
def state_transition(m: ClockModel, delta: float) -> np.ndarray:
    """The matrix Φ(δ) that carries the state over a step of ``delta`` seconds."""
    return _mat(jl.state_transition(m._jl(), float(delta)))


@_map_errors
def process_noise(m: ClockModel, delta: float) -> np.ndarray:
    """Covariance Q(δ) of the noise the clock picks up over a step of ``delta`` seconds."""
    return _mat(jl.process_noise(m._jl(), float(delta)))


# ── Interactive noise fit ──────────────────────────────────────────────


class NoiseFit:
    """An interactive power-law fit to a σ(τ) result, built one noise type at a
    time; see ``noise_fit``. ``components`` lists ``(alpha, h, tau_from, tau_to)``
    in the order added.
    """

    def __init__(self, jl_fit: Any):
        self._jl = jl_fit

    @property
    def components(self) -> list[tuple[int, float, float, float]]:
        tau = np_f64(self._jl.result.tau)
        out = []
        for c in self._jl.components:
            a, h, i, j = _component(c)
            out.append((int(a), float(h), float(tau[int(i) - 1]), float(tau[int(j) - 1])))
        return out

    @_map_errors
    def add_component(self, noise: int | str, tau_index: slice | Sequence[int]) -> NoiseFit:
        """Fit one noise type (α or a name such as ``"WHFM"``) on a 0-based range of
        τ indices (``slice(0, 6)`` or ``[0, 5]``, both ends included for a list),
        against what the earlier components leave."""
        if isinstance(tau_index, slice):
            n = len(self._jl.result.tau)
            lo, hi, _ = tau_index.indices(n)
            first, last = lo, hi - 1
        else:
            first, last = int(tau_index[0]), int(tau_index[-1])
        key = noise if isinstance(noise, int) else jl.Symbol(str(noise))
        _add(self._jl, key, first + 1, last + 1)
        return self

    @_map_errors
    def undo(self) -> NoiseFit:
        """Remove the last component."""
        jl.seval("SigmaTau.undo!")(self._jl)
        return self

    def residual(self) -> np.ndarray:
        """The σ² the model leaves at each τ."""
        return np_f64(jl.noise_residual(self._jl))

    def __repr__(self) -> str:
        return f"NoiseFit({len(self.components)} components)"


@_map_errors
def noise_fit(r: StabilityResult) -> NoiseFit:
    """Start an interactive power-law fit on an ADEV, MDEV, HDEV or MHDEV result
    (or a total form); add components with ``add_component``."""
    return NoiseFit(jl.noise_fit(result_to_jl(r)))


# ── Kalman filter ──────────────────────────────────────────────────────


@dataclass(frozen=True)
class KalmanResult:
    """Return type of ``kalman_filter`` / ``kalman_smooth``.

    ``x`` holds one row per sample (phase, frequency, and drift with a drift
    state), ``P`` the covariances with shape (states, states, samples), and
    ``innovation`` / ``innovation_var`` the innovations (``nan`` at gaps).
    """

    model: ClockModel
    t: np.ndarray
    x: np.ndarray
    P: np.ndarray
    innovation: np.ndarray
    innovation_var: np.ndarray
    smoothed: bool
    _jl: Any = field(repr=False, compare=False)


def _kalman_from_jl(r: Any) -> KalmanResult:
    return KalmanResult(
        _model_from_jl(r.model),
        np_f64(r.t),
        _mat(r.x),
        _mat(r.P),
        np_f64(r.innovation),
        np_f64(r.innovation_var),
        bool(r.smoothed),
        r,
    )


@_map_errors
def kalman_filter(
    data: PhaseData, m: ClockModel, *, x0: Sequence[float] | None = None, P0: Any = None
) -> KalmanResult:
    """Run the Kalman filter for clock model ``m`` over a phase record; ``nan``
    samples are gaps the filter predicts across. The start defaults to a
    least-squares polynomial through the first samples."""
    kw = {}
    if x0 is not None:
        kw["x0"] = to_jl_f64(x0)
    if P0 is not None:
        kw["P0"] = juliacall.convert(_MatF64, np.ascontiguousarray(P0, dtype=np.float64))
    return _kalman_from_jl(jl.kalman_filter(to_jl_data(data), m._jl(), **kw))


@_map_errors
def kalman_smooth(kf: KalmanResult) -> KalmanResult:
    """Rauch–Tung–Striebel smoother over a filtered run."""
    return _kalman_from_jl(jl.kalman_smooth(kf._jl))


@_map_errors
def holdover_error(kf: KalmanResult, tau: Sequence[float]) -> np.ndarray:
    """1σ time error (s) after free-running for each τ from the last sample."""
    return np_f64(jl.holdover_error(kf._jl, to_jl_f64(tau)))


@_map_errors
def innovation_check(kf: KalmanResult, *, confidence: float = 0.95) -> dict:
    """Whether the model fits: the NIS mean against its χ² interval and the lag-1
    autocorrelation of the normalized innovations, with ``ok`` when both pass."""
    c = jl.innovation_check(kf._jl, confidence=float(confidence))
    return {
        "nis": float(c.nis),
        "nis_bounds": tuple(float(b) for b in c.nis_bounds),
        "acf": float(c.acf),
        "acf_bound": float(c.acf_bound),
        "ok": bool(c.ok),
    }


@_map_errors
def prediction_error(
    kf: KalmanResult, data: PhaseData, horizons: Sequence[int], *, burnin: int = 0
) -> np.ndarray:
    """RMS error (s) of the filter's phase extrapolation each horizon (samples) ahead."""
    return np_f64(
        jl.prediction_error(kf._jl, to_jl_data(data), to_jl_ints(horizons), burnin=int(burnin))
    )


@_map_errors
def tune_q(
    data: PhaseData,
    m: ClockModel,
    *,
    horizons: Sequence[int] = (10, 100, 1000),
    decades: float = 2,
    points: int = 9,
) -> dict:
    """Retune the white and random-walk FM levels on a log grid for the smallest
    phase prediction error; returns ``model``, the grid axes ``q1``, ``q2`` and
    the ``cost`` surface."""
    r = jl.tune_q(
        to_jl_data(data),
        m._jl(),
        horizons=to_jl_ints(horizons),
        decades=float(decades),
        points=int(points),
    )
    return {
        "model": _model_from_jl(r.model),
        "q1": np_f64(r.q1),
        "q2": np_f64(r.q2),
        "cost": _mat(r.cost),
    }


@_map_errors
def simulate_clock(
    m: ClockModel, n: int, *, rng: int | None = None
) -> tuple[PhaseData, np.ndarray]:
    """A phase record of ``n`` samples from ``m`` and its true states."""
    d, x = jl.simulate_clock(m._jl(), int(n), **_rng(rng))
    return data_from_jl(d), _mat(x)


# ── Ensembles and timescales ───────────────────────────────────────────


@dataclass(frozen=True)
class ClockComparisons:
    """Several clocks read against one reference: column ``i`` of ``x`` is clock i
    minus the reference (s), ``nan`` where a reading is missing."""

    x: np.ndarray
    tau0: float
    names: tuple[str, ...] = ()
    reference: str = ""

    def __post_init__(self):
        x = np.array(self.x, dtype=np.float64)
        if x.ndim != 2:
            raise ValueError("ClockComparisons: x must be a samples × clocks matrix")
        object.__setattr__(self, "x", x)
        object.__setattr__(
            self, "names", tuple(self.names) or tuple(f"clock {i + 1}" for i in range(x.shape[1]))
        )

    def _jl(self) -> Any:
        return jl.ClockComparisons(
            juliacall.convert(_MatF64, np.ascontiguousarray(self.x)),
            float(self.tau0),
            names=juliacall.convert(VecStr, list(self.names)),
            reference=self.reference,
        )

    def __getitem__(self, ij: int | tuple[int, int]) -> PhaseData:
        """``cc[i]`` one clock against the reference, ``cc[i, j]`` clock i − clock j (0-based)."""
        if isinstance(ij, tuple):
            i, j = ij
            return PhaseData(
                self.x[:, i] - self.x[:, j], self.tau0, source=f"{self.names[i]} − {self.names[j]}"
            )
        return PhaseData(self.x[:, ij], self.tau0, source=self.names[ij])


def _comparisons_from_jl(c: Any) -> ClockComparisons:
    return ClockComparisons(
        _mat(c.x), float(c.tau0), tuple(str(s) for s in c.names), str(c.reference)
    )


@_map_errors
def simulate_ensemble(
    models: Sequence[ClockModel], n: int, *, reference: int | None = 0, rng: int | None = None
) -> tuple[ClockComparisons, np.ndarray]:
    """Clocks drawn from ``models`` and read against member ``reference`` (0-based)
    or, with ``reference=None``, a perfect outside reference; returns
    the comparisons and the true states ``truth[k, state, clock]``."""
    ms = juliacall.convert(_VecModel, [m._jl() for m in models])
    cc, truth = jl.simulate_ensemble(
        ms, int(n), reference=0 if reference is None else int(reference) + 1, **_rng(rng)
    )
    return _comparisons_from_jl(cc), _mat(truth)


@dataclass(frozen=True)
class TimescaleResult:
    """Return type of ``timescale``: each clock minus the scale (``x``), its
    estimated frequency (``y``), its weight (rows sum to 1) and its prediction
    error (``nan`` without a reading), one column per clock."""

    algorithm: str
    comparisons: ClockComparisons
    t: np.ndarray
    x: np.ndarray
    y: np.ndarray
    weights: np.ndarray
    error: np.ndarray

    def __repr__(self) -> str:
        n, k = self.x.shape
        return f"TimescaleResult({self.algorithm}, {k} clocks, {n} samples)"


@_map_errors
def timescale(
    cc: ClockComparisons,
    *,
    algorithm: str = "at1",
    tau_y: float | None = None,
    tau_min: float | Sequence[float] | None = None,
    memory: float | None = None,
    bias: str = "paper",
    wmax: float = 1.0,
    outlier: float | None = None,
    models: Sequence[ClockModel] | None = None,
    pivot: int = 0,
) -> TimescaleResult:
    """Form a timescale from the clocks in ``cc``: ``"ewfa"`` weighs every clock the
    same (its frequencies average over ``tau_y`` seconds, default 100 τ₀),
    ``"at1"`` follows Weiss & Weissert 1991, and ``"greenhall"`` runs
    Greenhall's Kalman timescale (reduced Kalman at ``R = 0``, IEMx above) on
    one ``ClockModel`` per clock in ``models``, read at clock ``pivot``
    (0-based, default 0).

    AT1 options: ``tau_min`` (each clock's τ_min in seconds, a number or one per
    clock; default estimated from the record), ``memory`` (the variance time
    constant in seconds, default 20 days, at most a quarter of the record),
    ``bias`` (``"paper"`` or ``"exact"``), ``wmax`` (largest weight) and
    ``outlier`` (reject a clock whose error passes that many ε; default off).
    """
    kw = {"algorithm": jl.Symbol(str(algorithm)), "bias": jl.Symbol(str(bias)), "wmax": float(wmax)}
    if models is not None:
        kw["pivot"] = int(pivot) + 1
        kw["models"] = juliacall.convert(_VecModel, [m._jl() for m in models])
    if tau_min is not None:
        kw["tau_min"] = (
            float(tau_min)
            if np.isscalar(tau_min)
            else juliacall.convert(_VecF64, [float(t) for t in tau_min])
        )
    for name, value in (("tau_y", tau_y), ("memory", memory), ("outlier", outlier)):
        if value is not None:
            kw[name] = float(value)
    r = jl.timescale(cc._jl(), **kw)
    return TimescaleResult(
        str(r.algorithm), cc, np_f64(r.t), _mat(r.x), _mat(r.y), _mat(r.weights), _mat(r.error)
    )


def timescale_phase(r: TimescaleResult) -> PhaseData:
    """The scale against the comparisons' reference, ready for ``adev`` and the
    other estimators; ``nan`` where no clock was read."""
    c = r.comparisons.x
    ok = np.isfinite(c)
    first = ok.argmax(axis=1)
    rows = np.arange(c.shape[0])
    u = np.where(ok.any(axis=1), r.x[rows, first] - c[rows, first], np.nan)
    return PhaseData(u, r.comparisons.tau0, source=r.algorithm.upper())


__all__ = [
    "ClockModel",
    "clock_model",
    "state_transition",
    "process_noise",
    "NoiseFit",
    "noise_fit",
    "KalmanResult",
    "kalman_filter",
    "kalman_smooth",
    "holdover_error",
    "innovation_check",
    "prediction_error",
    "tune_q",
    "simulate_clock",
    "ClockComparisons",
    "simulate_ensemble",
    "TimescaleResult",
    "timescale",
    "timescale_phase",
]

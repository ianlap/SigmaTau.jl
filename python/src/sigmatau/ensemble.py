"""Clock models, the Kalman filter, holdover, steering, clock ensembles and
timescales — mirrors ``clockmodel.jl``, ``noisefit.jl``, ``kalman.jl``,
``holdover.jl``, ``steering.jl``, ``comparisons.jl`` and ``timescale.jl``.

Results that later calls build on (a filter run, a noise fit) keep their Julia
object, so a smoother or a holdover estimate runs on exactly what the filter
produced.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import InitVar, dataclass, field
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
_VecC64 = jl.seval("Vector{ComplexF64}")
_band = jl.seval("b -> b === nothing ? nothing : (Float64(b[1]), Float64(b[2]))")
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
    white noise (Zucca & Tavella 2015), with flicker noise as Markov states.

    ``R`` is the variance of a reading's noise, taken as white (s²), ``q1`` the white FM
    (s), ``q2`` the random-walk FM (1/s) and ``q3`` the random-run FM (1/s³)
    diffusion coefficient. ``flicker_fm`` (h₋₁) and ``flicker_pm`` (h₁) are carried
    by Markov states across ``fm_band`` and ``pm_band``, the first and last τ (s)
    each spans, ``(0, 0)`` for none. Left out, a band covers the τ up to
    ``tau_max`` (default 10⁶ τ₀) where the flicker holds a quarter of the model's
    Allan variance.
    """

    tau0: float
    R: float = 0.0
    q1: float = 0.0
    q2: float = 0.0
    q3: float = 0.0
    drift: bool | None = None
    flicker_fm: float = 0.0
    flicker_pm: float = 0.0
    fm_band: tuple[float, float] | None = None
    pm_band: tuple[float, float] | None = None
    tau_max: InitVar[float | None] = None

    @_map_errors
    def __post_init__(self, tau_max: float | None):
        if self.drift is None:
            object.__setattr__(self, "drift", self.q3 > 0)
        kw = {} if tau_max is None else {"tau_max": float(tau_max)}
        m = jl.ClockModel(
            tau0=float(self.tau0),
            R=float(self.R),
            q1=float(self.q1),
            q2=float(self.q2),
            q3=float(self.q3),
            drift=bool(self.drift),
            flicker_fm=float(self.flicker_fm),
            flicker_pm=float(self.flicker_pm),
            fm_band=_band(self.fm_band),
            pm_band=_band(self.pm_band),
            **kw,
        )
        object.__setattr__(self, "fm_band", tuple(m.fm_band))
        object.__setattr__(self, "pm_band", tuple(m.pm_band))

    def _jl(self) -> Any:
        return jl.ClockModel(
            float(self.tau0),
            float(self.R),
            float(self.q1),
            float(self.q2),
            float(self.q3),
            bool(self.drift),
            float(self.flicker_fm),
            _band(self.fm_band),
            float(self.flicker_pm),
            _band(self.pm_band),
        )


def _model_from_jl(m: Any) -> ClockModel:
    return ClockModel(
        float(m.tau0),
        float(m.R),
        float(m.q1),
        float(m.q2),
        float(m.q3),
        bool(m.drift),
        flicker_fm=float(m.flicker_fm),
        flicker_pm=float(m.flicker_pm),
        fm_band=tuple(m.fm_band),
        pm_band=tuple(m.pm_band),
    )


@_map_errors
def clock_model(
    fit: PowerLawFit | NoiseFit | Mapping[int | str, float],
    *,
    tau0: float | None = None,
    drift: bool | None = None,
    flicker: bool = True,
    tau_max: float | None = None,
) -> ClockModel:
    """The ``ClockModel`` with the white FM, random-walk FM and random-run FM
    levels of a power-law model, and its white PM as the reading noise ``R``.

    Flicker FM and flicker PM are kept as Markov states
    (``ClockModel.flicker_fm``, ``flicker_pm``) across the τ, up to ``tau_max``,
    where each holds a quarter of the model's Allan variance. With
    ``flicker=False`` flicker FM is left out, and flicker PM is folded into ``R``
    as the white PM of the same Allan variance where PM hands over to FM; the
    ``"stein"`` time scale needs such a model.

    ``fit`` is a ``PowerLawFit``, a ``NoiseFit`` or an ``{α: h_α}`` mapping (with
    ``tau0``).
    """
    kw: dict[str, Any] = {"flicker": bool(flicker)}
    if drift is not None:
        kw["drift"] = bool(drift)
    if tau_max is not None:
        kw["tau_max"] = float(tau_max)
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


@dataclass(frozen=True)
class HoldoverResult:
    """Return type of ``holdover``: the statistics of the time error at each
    holdover time ``tau`` (s), over every start the record allows. ``quantile``
    is the |error| that ``confidence`` of the trials stay within; ``claimed`` is
    the standard deviation a Kalman filter gives the error itself (empty for a
    least-squares fit); ``limit`` and ``time_to_limit`` are ``nan`` when no limit
    was asked for. With a stretch
    that grows with ``tau``, ``grow`` is its multiple and ``train`` is ``nan``;
    otherwise ``grow`` is 0."""

    tau: np.ndarray
    rms: np.ndarray
    quantile: np.ndarray
    worst: np.ndarray
    mean: np.ndarray
    n: np.ndarray
    claimed: np.ndarray
    train: float
    grow: int
    drift: bool
    anchor: bool
    confidence: float
    limit: float
    time_to_limit: float


@dataclass(frozen=True)
class HoldoverPrediction:
    """Return type of ``holdover`` on a noise model: the standard deviation of
    the time error at each holdover time ``tau`` (s), its ``quantile`` at
    ``confidence``, and ``budget``, the variance each noise type in ``alphas``
    adds (one row per ``tau``, one column per α)."""

    tau: np.ndarray
    sigma: np.ndarray
    quantile: np.ndarray
    alphas: np.ndarray
    budget: np.ndarray
    train: float
    grow: int
    drift: bool
    anchor: bool
    confidence: float
    limit: float
    time_to_limit: float


@_map_errors
def holdover(
    data: PhaseData | PowerLawFit | ClockModel | Mapping[int | str, float],
    model: ClockModel | None = None,
    *,
    train: float | None = None,
    grow: int | None = None,
    tau: Sequence[float] | None = None,
    tau0: float | None = None,
    drift: bool = False,
    anchor: bool = True,
    limit: float | None = None,
    confidence: float = 0.95,
) -> HoldoverResult | HoldoverPrediction:
    """Holdover: learn the frequency (and drift) on ``train`` seconds, free-run
    from the last reading for each ``tau``, and compare with the clock
    (Delporte et al. 2000).

    On a ``PhaseData`` record the error is measured from every start the record
    allows. On a noise model (a ``PowerLawFit``, a ``ClockModel`` or an
    ``{α: h_α}`` mapping, which needs ``tau0``) it is predicted, with each noise
    type's share in ``budget``; ``train`` is then required.

    ``grow=k`` in place of ``train`` makes the stretch k times each holdover
    time. At 2 (3 with ``drift``) the test is the sibling of ``tdev``
    (``htdev``): √6 TDEV (√20 HTDEV) is the error of the same extrapolation
    made from τ-averages of the phase, which this test equals at τ₀ and
    exceeds by a constant 1.26 to 1.42 on FM noise.

    With a record and a ``ClockModel``, the Kalman filter for that model is the
    predictor: it runs once over the record, and from every sample past
    ``train`` seconds its estimate is carried ``tau`` ahead with no readings.
    ``claimed`` is then the filter's own standard deviation for the error;
    ``drift``, ``anchor`` and ``grow`` do not apply."""
    kw: dict[str, Any] = {"confidence": float(confidence)}
    if model is None:
        kw |= {"drift": bool(drift), "anchor": bool(anchor)}
    if train is not None:
        kw["train"] = float(train)
    if grow is not None:
        kw["grow"] = int(grow)
    if tau is not None:
        kw["tau"] = to_jl_f64(tau)
    if limit is not None:
        kw["limit"] = float(limit)
    if isinstance(data, PhaseData):
        args = (to_jl_data(data),) if model is None else (to_jl_data(data), model._jl())
        r = jl.holdover(*args, **kw)
        return HoldoverResult(
            np_f64(r.tau),
            np_f64(r.rms),
            np_f64(r.quantile),
            np_f64(r.worst),
            np_f64(r.mean),
            np.array(r.n, dtype=np.int64),
            np_f64(r.claimed),
            float(r.train),
            int(r.grow),
            bool(r.drift),
            bool(r.anchor),
            float(r.confidence),
            float(r.limit),
            float(r.time_to_limit),
        )
    if train is None and grow is None:
        raise ValueError("holdover: a noise model needs `train`, the stretch in seconds, or `grow`")
    if tau0 is None:
        tau0 = getattr(data, "tau0", None)
    if tau0 is not None:
        kw["tau0"] = float(tau0)
    r = jl.holdover(data._jl() if isinstance(data, ClockModel) else powerlaw_to_jl(data), **kw)
    return HoldoverPrediction(
        np_f64(r.tau),
        np_f64(r.sigma),
        np_f64(r.quantile),
        np.array(r.alphas, dtype=np.int64),
        _mat(r.budget),
        float(r.train),
        int(r.grow),
        bool(r.drift),
        bool(r.anchor),
        float(r.confidence),
        float(r.limit),
        float(r.time_to_limit),
    )


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
    phase prediction error; returns ``model``, the grid axes ``q1``, ``q2``, the
    ``cost`` surface, the factor ``scale`` on every level of the model (the mean
    square of its normalized innovations) and ``edge``, whether the best pair
    lies on the edge of the grid."""
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
        "scale": float(r.scale),
        "edge": bool(r.edge),
    }


@_map_errors
def simulate_clock(
    m: ClockModel, n: int, *, rng: int | None = None
) -> tuple[PhaseData, np.ndarray]:
    """A phase record of ``n`` samples from ``m`` and its true states."""
    d, x = jl.simulate_clock(m._jl(), int(n), **_rng(rng))
    return data_from_jl(d), _mat(x)


# ── Steering ───────────────────────────────────────────────────────────


@dataclass(frozen=True)
class SteeringGains:
    """Gains of a steering loop (Matsakis 2019). Every ``interval`` seconds the
    frequency of the clock is changed by
    −(``phase``·p̂ + ``frequency``·ŷ + ``integral``·Σp̂), with p̂ and ŷ the Kalman
    estimates of its phase and frequency against the reference and Σp̂ the sum of
    the phase estimates at the steers so far. ``phase`` and ``integral`` are in
    1/s. Build them with ``steering_gains``."""

    phase: float
    frequency: float
    integral: float
    interval: float

    def _jl(self) -> Any:
        return jl.SteeringGains(
            float(self.phase), float(self.frequency), float(self.integral), float(self.interval)
        )


def _gains_from_jl(g: Any) -> SteeringGains:
    return SteeringGains(float(g.phase), float(g.frequency), float(g.integral), float(g.interval))


@_map_errors
def steering_gains(
    T: float | Sequence[complex], interval: float, *, integral: bool = False
) -> SteeringGains:
    """Gains by pole placement. A time constant ``T`` (s) puts the poles together
    at exp(−interval/T): a critically damped loop, with three gains when
    ``integral``. A sequence in place of ``T`` is the poles themselves in the
    z-plane, two or three, complex ones in conjugate pairs."""
    if np.ndim(T) == 0:
        return _gains_from_jl(jl.steering_gains(float(T), float(interval), integral=bool(integral)))
    poles = juliacall.convert(_VecC64, [complex(z) for z in T])
    return _gains_from_jl(jl.steering_gains(poles, float(interval)))


@_map_errors
def lqg_gains(
    interval: float, *, phase: float, frequency: float = 0.0, steer: float
) -> SteeringGains:
    """Gains that minimize the sum over the steers of ``phase``·p² +
    ``frequency``·f² + ``steer``·u², the linear quadratic Gaussian design
    (Koppang & Leland 1999). Only the ratios of the weights matter; a heavier
    ``steer`` weight is a slower loop."""
    return _gains_from_jl(
        jl.lqg_gains(
            float(interval), phase=float(phase), frequency=float(frequency), steer=float(steer)
        )
    )


@_map_errors
def steering_poles(g: SteeringGains) -> dict:
    """Where the gains put the poles of the closed loop: ``poles`` (complex), the
    ``time_constants`` (s) of their decay, the ``frequencies`` (Hz) they ring at,
    and ``stable``, whether all are inside the unit circle."""
    p = jl.steering_poles(g._jl())
    return {
        "poles": np.array(p.poles, dtype=np.complex128),
        "time_constants": np_f64(p.time_constants),
        "frequencies": np_f64(p.frequencies),
        "stable": bool(p.stable),
    }


@dataclass(frozen=True)
class SteeringResult:
    """Return type of ``steer`` on a record: ``steered`` is the record as the
    loop would have left it, ``steers`` the frequency steer made at each sample,
    ``correction`` their running sum, and ``phase`` / ``frequency`` the filter's
    estimates of the steered clock."""

    gains: SteeringGains
    steered: PhaseData
    steers: np.ndarray
    correction: np.ndarray
    phase: np.ndarray
    frequency: np.ndarray


@dataclass(frozen=True)
class SteeringPrediction:
    """Return type of ``steer`` on a clock model: the RMS time error ``phase``
    (s), frequency error ``frequency`` and ``steers`` of the settled loop, and
    the Allan deviation ``adev`` of the steered clock at ``tau`` (s)."""

    gains: SteeringGains
    phase: float
    frequency: float
    steers: float
    tau: np.ndarray
    adev: np.ndarray


@_map_errors
def steer(
    data: PhaseData | ClockModel,
    model: ClockModel | SteeringGains,
    gains: SteeringGains | None = None,
    *,
    tau: Sequence[float] | None = None,
    readings: bool = False,
) -> SteeringResult | SteeringPrediction:
    """Steer a clock to its reference with the Kalman filter for ``model`` and
    the gains of ``steering_gains``.

    ``steer(data, model, gains)`` replays the phase record of the free-running
    clock against the reference through the loop. ``steer(model, gains)``
    predicts the steady state of the loop from the model alone, for a steer at
    every sample; ``tau`` sets the averaging times of its Allan deviation. The
    prediction is for the clock itself; with ``readings`` its ``phase`` and
    ``adev`` are those of the readings, with their noise, as a steered record
    shows them."""
    if isinstance(data, ClockModel):
        if gains is not None or not isinstance(model, SteeringGains):
            raise TypeError("steer: pass (data, model, gains) or (model, gains)")
        kw = {} if tau is None else {"tau": to_jl_f64(tau)}
        r = jl.steer(data._jl(), model._jl(), readings=bool(readings), **kw)
        return SteeringPrediction(
            model,
            float(r.phase),
            float(r.frequency),
            float(r.steers),
            np_f64(r.tau),
            np_f64(r.adev),
        )
    if gains is None or tau is not None or readings:
        raise TypeError("steer: pass (data, model, gains) or (model, gains, tau=...)")
    r = jl.steer(to_jl_data(data), model._jl(), gains._jl())
    return SteeringResult(
        gains,
        data_from_jl(r.steered),
        np_f64(r.steers),
        np_f64(r.correction),
        np_f64(r.phase),
        np_f64(r.frequency),
    )


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


@_map_errors
def clock_models(
    cc: ClockComparisons, *, alphas: Sequence[int] | None = None, flicker: bool = False
) -> list[ClockModel]:
    """One ``ClockModel`` per clock, fitted from the comparisons alone (three or
    more clocks): each clock's Allan deviation from the Groslambert covariance,
    then a power-law fit with the components ``alphas`` (by default Julia's:
    white and flicker PM, white FM and random-walk FM), the flicker PM folded
    into ``R``. With ``flicker=True`` flicker FM joins the fit and both flickers
    are kept as Markov states where they matter."""
    kw: dict[str, Any] = {} if alphas is None else {"alphas": to_jl_ints(alphas)}
    if flicker:
        kw["flicker"] = True
    return [_model_from_jl(m) for m in jl.clock_models(cc._jl(), **kw)]


@dataclass(frozen=True)
class TimescaleResult:
    """Return type of ``timescale``: each clock minus the scale (``x``), its
    estimated frequency (``y``), its weight (rows sum to 1), its prediction
    error (``nan`` without a reading) and its phase as the algorithm estimated
    it (``xhat``: the filter's phase state for the Kalman scales, where ``x`` is
    read at the pivot; equal to ``x`` for the others), one column per clock."""

    algorithm: str
    comparisons: ClockComparisons
    t: np.ndarray
    x: np.ndarray
    y: np.ndarray
    weights: np.ndarray
    error: np.ndarray
    xhat: np.ndarray

    def __repr__(self) -> str:
        n, k = self.x.shape
        return f"TimescaleResult({self.algorithm}, {k} clocks, {n} samples)"


@_map_errors
def timescale(
    cc: ClockComparisons,
    *,
    algorithm: str = "at1",
    tau_min: float | Sequence[float] | None = None,
    memory: float | None = None,
    bias: str = "paper",
    wmax: float = 1.0,
    outlier: float | None = None,
    models: Sequence[ClockModel] | None = None,
    pivot: int = 0,
    steps: bool | str = False,
) -> TimescaleResult:
    """Form a timescale from the clocks in ``cc``: ``"ewfa"`` weighs every clock the
    same (Kettering 2023, section 3.1),
    ``"at1"`` follows Weiss & Weissert 1991, ``"at2"`` adds their Kalman
    frequency update, and ``"greenhall"`` runs
    Greenhall's Kalman timescale (reduced Kalman at ``R = 0``, IEMx above) on
    one ``ClockModel`` per clock in ``models``, read at clock ``pivot``
    (0-based, default 0); ``"jones_tryon"`` is the same filter without the
    covariance reduction, and ``"stein"`` is Stein's Kalman timescale (Stein 2003;
    KAS-2 with ``wmax``) on the same models.

    AT1 and AT2 options: ``tau_min`` (each clock's τ_min in seconds, a number or one per
    clock; default estimated from the record), ``memory`` (the variance time
    constant in seconds, default 20 days, at most a quarter of the record),
    ``bias`` (``"paper"`` or ``"exact"``), ``wmax`` (largest weight) and
    ``outlier`` (reject a clock whose error passes that many ε; default off).
    ``steps`` (AT2 only) detects clock frequency steps: ``"best"`` (or ``True``)
    recomputes the scale from the step, ``"real"`` keeps the values already
    issued and corrects forward.
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
    for name, value in (("memory", memory), ("outlier", outlier)):
        if value is not None:
            kw[name] = float(value)
    if steps:
        kw["steps"] = jl.Symbol("best" if steps is True else str(steps))
    r = jl.timescale(cc._jl(), **kw)
    return TimescaleResult(
        str(r.algorithm), cc, np_f64(r.t), _mat(r.x), _mat(r.y), _mat(r.weights), _mat(r.error),
        _mat(r.xhat),
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


def implicit_mean(r: TimescaleResult, truth: np.ndarray) -> PhaseData:
    """The implicit ensemble mean of a ``"greenhall"`` scale against true time, in
    simulation: the weighted mean ``Σ aᵢ (xᵢ − x̂ᵢ)`` of the corrected clocks, each
    clock's true phase (``truth`` as from ``simulate_ensemble``) less the filter's
    estimate ``r.xhat`` (Greenhall 2006, Eqs. (3.8) and (3.13)). Its stability is
    the scale's own, free of any clock's reading noise."""
    if r.algorithm != "greenhall":
        raise ValueError(f"implicit_mean: the implicit mean belongs to greenhall, not {r.algorithm}")
    truth = np.asarray(truth, dtype=np.float64)
    if truth.ndim != 3 or (truth.shape[0], truth.shape[2]) != r.weights.shape:
        raise ValueError(
            f"implicit_mean: truth is {truth.shape}, the scale has {r.weights.shape[0]} samples "
            f"and {r.weights.shape[1]} clocks"
        )
    mean = (r.weights * (truth[:, 0, :] - r.xhat)).sum(axis=1)
    return PhaseData(mean, r.comparisons.tau0, source="IEM")


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
    "HoldoverResult",
    "HoldoverPrediction",
    "holdover",
    "innovation_check",
    "prediction_error",
    "tune_q",
    "simulate_clock",
    "ClockComparisons",
    "simulate_ensemble",
    "clock_models",
    "TimescaleResult",
    "timescale",
    "timescale_phase",
    "implicit_mean",
]

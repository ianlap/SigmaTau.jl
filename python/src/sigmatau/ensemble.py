"""Clock models, the Kalman filter, holdover, steering and time scales.

A result that a later call builds on, such as a filter run or a noise fit,
keeps its Julia object, so that a smoother or a holdover estimate runs on what
the filter produced.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import InitVar, dataclass, field
from typing import Any

import juliacall
import numpy as np

from ._convert import (
    as_int,
    data_from_jl,
    noise_alpha,
    np_f64,
    powerlaw_to_jl,
    powerlawfit_to_jl,
    result_to_jl,
    to_jl_data,
    to_jl_f64,
    to_jl_ints,
)
from ._julia import VecStr, _map_errors, jl
from .noise import _seeded
from .types import FrequencyData, PhaseData, PowerLawFit, StabilityResult

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
    return {} if seed is None else {"rng": _seeded(seed)}


def _clock(i: object, n: int, name: str) -> int:
    """``i`` as the index, counted from 0, of one of ``n`` clocks."""
    i = as_int(i, name)
    if not 0 <= i < n:
        raise ValueError(f"{name} must be the index of a clock, 0 to {n - 1}, got {i}")
    return i


def _mat(a: Any) -> np.ndarray:
    return np.array(a, dtype=np.float64)


# ── Clock model ────────────────────────────────────────────────────────


@dataclass(frozen=True)
class ClockModel:
    """State-space model of a clock read every ``tau0`` seconds: phase, frequency
    and, with ``drift``, frequency drift driven by white noises, with optional
    Markov states for flicker noise.

    The phase, frequency and drift follow the three-state model of Zucca and
    Tavella (2015), without its deterministic terms. ``R`` is the white PM
    variance of a reading (s²), and ``q1``, ``q2`` and ``q3`` are the white FM
    (s), random-walk FM (1/s) and random-run FM (1/s³) levels, the squares of
    that paper's diffusion coefficients. They stand for the coefficients of
    S_y(f) = Σ h_α f^α as h₂ = 8π²τ₀R, h₀ = 2q₁, h₋₂ = q₂/(2π²) and
    h₋₄ = q₃/(8π⁴). ``drift`` says whether the state carries a frequency
    drift, by default whether ``q3`` is positive.

    ``flicker_fm`` (h₋₁) and ``flicker_pm`` (h₁, in s²) are carried by Markov
    states: the flicker FM model of Davis, Greenhall and Stacey (2005), with
    time constants a decade apart, and the same construction on the readings
    for flicker PM, which is this package's choice. ``fm_band`` and
    ``pm_band`` are the first and last τ, in seconds, that the states span,
    ``(0, 0)`` for no states. Left out, each is the range of τ between τ₀ and
    ``tau_max``, 10⁶ τ₀ by default, at which that flicker holds at least a
    quarter of the Allan variance of the model.
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
    fm_band: tuple[float, float] | None = None,
    pm_band: tuple[float, float] | None = None,
) -> ClockModel:
    """The ``ClockModel`` of a power-law noise model.

    White FM, random-walk FM and random-run FM become ``q1``, ``q2`` and
    ``q3``, and white PM becomes ``R``, by the relations given under
    ``ClockModel``. Flicker FM and flicker PM are kept as its Markov states,
    across the τ where each holds a quarter of the Allan variance of the
    model.

    With ``flicker=False`` flicker FM is left out and flicker PM is folded into
    ``R``: the white PM that has the Allan variance of the white and flicker PM
    together at the last τ where that exceeds the FM terms. The model then has
    two states, or three with drift, as the ``"stein"`` time scale needs, and
    it has less noise than the fit. Flicker-walk FM (α = −3) is left out with
    either setting.

    ``fit`` is a ``PowerLawFit`` that has a sample interval, a ``NoiseFit``, or
    a mapping of h_α keyed by α or by noise type, which needs ``tau0``.
    ``drift`` says whether the model has a drift state, by default whether h₋₄
    is positive. ``tau_max``, ``fm_band`` and ``pm_band`` are those of
    ``ClockModel``; ``tau_max`` defaults to the longest τ of a ``PowerLawFit``
    of a deviation, and to 10⁶ τ₀ otherwise.
    """
    kw: dict[str, Any] = {"flicker": bool(flicker)}
    if drift is not None:
        kw["drift"] = bool(drift)
    if tau_max is not None:
        kw["tau_max"] = float(tau_max)
    if fm_band is not None:
        kw["fm_band"] = _band(fm_band)
    if pm_band is not None:
        kw["pm_band"] = _band(pm_band)
    if isinstance(fit, NoiseFit):
        return _model_from_jl(jl.clock_model(fit._jl, **kw))
    if tau0 is not None:
        kw["tau0"] = float(tau0)
    if isinstance(fit, PowerLawFit) and np.isfinite(fit.tau0):
        # As Julia's own type, so that `tau_max` defaults to the longest τ of the fit.
        return _model_from_jl(jl.clock_model(powerlawfit_to_jl(fit), **kw))
    if tau0 is None:
        raise ValueError("clock_model: pass `tau0`, the sample interval in seconds")
    return _model_from_jl(jl.clock_model(powerlaw_to_jl(fit), **kw))


@_map_errors
def state_transition(m: ClockModel, delta: float) -> np.ndarray:
    """The matrix Φ(δ) that carries the state of ``m`` over a step of ``delta``
    seconds. The Markov states follow the phase, frequency and drift, the
    frequency states first."""
    return _mat(jl.state_transition(m._jl(), float(delta)))


@_map_errors
def process_noise(m: ClockModel, delta: float) -> np.ndarray:
    """Covariance Q(δ) of the noise that the state of ``m`` gathers over a step of
    ``delta`` seconds (Greenhall 2003, Eq. (5); Zucca and Tavella 2015, Eq. (3))."""
    return _mat(jl.process_noise(m._jl(), float(delta)))


# ── Interactive noise fit ──────────────────────────────────────────────


class NoiseFit:
    """A power-law model built one component at a time on a measured σ(τ), by
    :func:`noise_fit`, ``add_component`` and ``undo``. ``powerlaw_sigma``,
    ``powerlaw_psd``, ``clock_model`` and ``holdover`` accept it as a model.

    ``components`` lists ``(alpha, h, tau_from, tau_to)`` in the order added:
    the exponent, the level h_α and the first and last τ, in seconds, that the
    component was fitted on.
    """

    def __init__(self, jl_fit: Any):
        self._jl = jl_fit

    @property
    def components(self) -> list[tuple[int, float, float, float]]:
        """The components so far, as ``(alpha, h, tau_from, tau_to)`` in the order added."""
        tau = np_f64(self._jl.result.tau)
        out = []
        for c in self._jl.components:
            a, h, i, j = _component(c)
            out.append((int(a), float(h), float(tau[int(i) - 1]), float(tau[int(j) - 1])))
        return out

    @_map_errors
    def add_component(self, noise: int | str, tau_index: slice | Sequence[int]) -> NoiseFit:
        """Add the power law ``noise``, an α or a name such as ``"WHFM"``, fitted at
        the τ indices ``tau_index``, and refit the levels already there.

        ``tau_index`` counts from 0: a slice such as ``slice(0, 6)``, or a
        sequence whose first and last entries are the first and last index,
        both included. Points where the residual is not positive are skipped.
        """
        if isinstance(tau_index, slice):
            n = len(self._jl.result.tau)
            lo, hi, _ = tau_index.indices(n)
            first, last = lo, hi - 1
        else:
            first, last = as_int(tau_index[0], "tau_index"), as_int(tau_index[-1], "tau_index")
        _add(self._jl, noise_alpha(noise), first + 1, last + 1)
        return self

    @_map_errors
    def undo(self) -> NoiseFit:
        """Remove the component added last and refit the rest."""
        jl.seval("SigmaTau.undo!")(self._jl)
        return self

    def residual(self) -> np.ndarray:
        """The measured σ²(τ) less the variance of every component so far. It can go
        negative where the model overshoots."""
        return np_f64(jl.noise_residual(self._jl))

    def __repr__(self) -> str:
        return f"NoiseFit({len(self.components)} components)"


@_map_errors
def noise_fit(r: StabilityResult) -> NoiseFit:
    """Start an interactive power-law fit on ``r``, an ADEV, MDEV, HDEV or MHDEV
    result or one of their total forms, which are fitted with the response of
    the parent. Add components with ``add_component`` and take the last one
    back with ``undo``.

    The level of a component is fitted over its own range of τ against what
    the other components leave. After each change every level is refitted
    until they settle, so the model does not depend on the order in which the
    components were added.
    """
    return NoiseFit(jl.noise_fit(result_to_jl(r)))


def _noise_model_to_jl(h: Any) -> Any:
    """A noise model for Julia: a ``ClockModel`` or a ``NoiseFit`` as Julia's own
    type, a ``PowerLawFit`` or a mapping as a ``Dict`` of h_α."""
    if isinstance(h, ClockModel):
        return h._jl()
    if isinstance(h, NoiseFit):
        return h._jl
    return powerlaw_to_jl(h)


# ── Kalman filter ──────────────────────────────────────────────────────


@dataclass(frozen=True)
class KalmanResult:
    """Result of :func:`kalman_filter` or :func:`kalman_smooth`: the state estimate
    at every sample, its covariance and the innovations.

    ``x`` holds one row per sample: phase (s), frequency, drift (1/s) for a
    model with a drift state, then the Markov frequency states, which add to
    the frequency, and the Markov phase states, which add to the readings.
    ``P`` holds the covariances, ``P[:, :, k]`` at sample k. ``innovation`` is
    each reading less its prediction (s) and ``innovation_var`` its variance
    (s²), both ``nan`` at gaps. ``t`` is in seconds from the first sample, and
    ``smoothed`` says whether ``x`` and ``P`` are smoothed.
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
    """Kalman filter of clock model ``m`` over a phase record.

    A ``nan`` sample is a gap: the state is predicted across it without an
    update. The covariance is updated in the Joseph form. The ``tau0`` of the
    record must be that of the model, and a model with no noise at all is
    refused.

    ``x0`` and ``P0`` are the state and its covariance at the first sample,
    before its reading is used. By default the filter starts from a
    least-squares line through the first 16 readings, a parabola for a model
    with drift, with the covariance of that fit multiplied by ten, and with
    the Markov states at zero with their stationary variances.
    """
    kw = {}
    if x0 is not None:
        kw["x0"] = to_jl_f64(x0)
    if P0 is not None:
        kw["P0"] = juliacall.convert(_MatF64, np.ascontiguousarray(P0, dtype=np.float64))
    return _kalman_from_jl(jl.kalman_filter(to_jl_data(data), m._jl(), **kw))


@_map_errors
def kalman_smooth(kf: KalmanResult) -> KalmanResult:
    """Rauch–Tung–Striebel smoother over a filtered run: every estimate then uses
    the whole record. The innovations of the result are those of the filter."""
    return _kalman_from_jl(jl.kalman_smooth(kf._jl))


@_map_errors
def holdover_error(kf: KalmanResult, tau: Sequence[float]) -> np.ndarray:
    """Standard deviation, in seconds, of the time error of the clock after running
    free for each τ in ``tau`` (seconds) from the last sample, as the model
    gives it: the error of the last estimate carried forward and the noise that
    the clock gathers on the way. The noise of a reading is not in it.
    :func:`holdover` measures the same error on the record."""
    return np_f64(jl.holdover_error(kf._jl, to_jl_f64(tau)))


@dataclass(frozen=True)
class HoldoverResult:
    """Result of :func:`holdover` on a record: the statistics of the time error at
    each holdover time ``tau`` (s), over every start the record allows.

    ``rms``, ``worst`` and ``mean`` are the RMS, the largest magnitude and the
    mean of the error (s), ``quantile`` the magnitude that ``confidence`` of
    the trials stay within, and ``n`` the number of trials. ``claimed`` is the
    standard deviation that a Kalman filter itself gives the error, empty for
    a least-squares fit. ``train`` is the length of the training stretch in
    seconds, ``nan`` when it grows with ``tau``; ``grow`` is then its multiple
    of each holdover time, and 0 otherwise. ``limit`` and ``time_to_limit``
    are ``nan`` when no limit was asked for.
    """

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
    """Result of :func:`holdover` on a noise model: the standard deviation ``sigma``
    of the time error at each holdover time ``tau`` (s), its ``quantile`` at
    ``confidence``, and ``budget``, the variance (s²) that each noise type in
    ``alphas`` adds, with one row per ``tau`` and one column per α. The rows
    of ``budget`` sum to ``sigma**2``. The other fields are those of
    ``HoldoverResult``."""

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
    data: PhaseData | PowerLawFit | NoiseFit | ClockModel | Mapping[int | str, float],
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
    """The time error of a clock that runs on without its reference.

    ``holdover(record, ...)`` measures it on the record, at every start the
    record allows. A line is fitted to ``train`` seconds of phase, or a
    parabola with ``drift=True``. From the last reading of the stretch the
    clock is left to run on the fitted frequency and drift for each ``tau``,
    and the error is what the record then shows. The stretch moves on by one
    sample, as in the test of Delporte et al. (2000). With ``anchor=False``
    the value of the fit takes the place of the last reading, which is the
    extrapolated fit of that paper. ``nan`` samples are gaps.

    ``grow=k`` in place of ``train`` makes the stretch k times each holdover
    time. At 2, or 3 with ``drift``, the test is the counterpart of the time
    deviations: √6 ``tdev`` is the RMS error of extrapolating the next
    τ-average of the phase from the two before it, and √20 ``htdev`` that of
    the next from the three before it. The test equals them at τ₀, and on FM
    noise it reads 1.22 to 1.42 times as much at long τ.

    ``holdover(record, model)`` is the same test with the Kalman filter for a
    ``ClockModel`` as the predictor. The filter runs once over the record, and
    at every sample past the first ``train`` seconds its estimate is carried
    ``tau`` ahead without further readings. ``claimed`` is the standard
    deviation that the filter itself gives the error, and ``rms`` agrees with
    it when the model fits the clock. ``drift``, ``anchor`` and ``grow`` do
    not apply.

    ``holdover(model, train=...)`` predicts the test from a noise model: a
    ``PowerLawFit``, a ``NoiseFit``, a ``ClockModel``, or a mapping of h_α keyed
    by α or by noise type, which needs ``tau0``. ``budget`` holds the share of
    each noise type. ``train`` has no default here, and ``grow`` needs
    ``tau``.

    ``train`` is a whole multiple of τ₀, by default a quarter of the record.
    ``tau`` holds the holdover times in seconds, rising whole multiples of τ₀,
    by default octave multiples of τ₀. ``limit`` is a time error in seconds:
    ``time_to_limit`` is then the holdover time at which the quantile at
    ``confidence`` reaches it.
    """
    if isinstance(data, FrequencyData):
        raise TypeError("holdover: the record must be PhaseData, as `data.to_phase()` gives")
    kw: dict[str, Any] = {"confidence": float(confidence)}
    if model is None:
        kw |= {"drift": bool(drift), "anchor": bool(anchor)}
    elif grow is not None or drift or not anchor:
        raise ValueError("holdover: `grow`, `drift` and `anchor` do not apply with a ClockModel")
    if train is not None:
        kw["train"] = float(train)
    if grow is not None:
        kw["grow"] = as_int(grow, "grow")
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
    r = jl.holdover(_noise_model_to_jl(data), **kw)
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
    """Test of the model against the record, from the innovations of the filter.

    Under the model the normalized innovations are independent, with zero mean
    and unit variance. Returns a dict with ``nis``, the mean of their squares,
    and ``nis_bounds``, its χ² interval at ``confidence``; ``acf``, their
    lag-1 autocorrelation, and ``acf_bound``, the bound on its magnitude; and
    ``ok``, whether both are inside. A mean above the interval says that the
    noise of the model is too low. Both intervals narrow as 1/√n, so on a long
    record the test fails a model that departs a little from the noise and
    still predicts well.
    """
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
    """RMS error, in seconds, of the predictions of the readings ``horizons`` samples
    ahead. From the estimate at each sample past ``burnin`` the state is carried
    forward, and the reading it predicts is compared with the reading there, so
    the error holds the noise of that reading as well as the clock's. ``kf`` is
    the run of the filter over ``data``, not a smoothed one."""
    return np_f64(
        jl.prediction_error(
            kf._jl,
            to_jl_data(data),
            to_jl_ints(horizons, "a horizon"),
            burnin=as_int(burnin, "burnin"),
        )
    )


@_map_errors
def tune_q(
    data: PhaseData,
    m: ClockModel,
    *,
    burnin: int | None = None,
    maxiter: int = 300,
) -> dict:
    """Retune the white FM and random-walk FM levels of ``m`` by maximum likelihood.

    ``q1`` and ``q2`` are searched (Nelder-Mead on their logarithms, started
    from those of ``m``) for the smallest -2 ln L of the innovations of the
    filter (Tryon and Jones 1983, Eq. (5.1)), with the other levels held
    fixed. ``burnin`` is the number of samples at the start left out of the
    sum, a tenth of the record by default, and ``maxiter`` limits the
    iterations of the search.

    Returns a dict with ``model``, the tuned model; ``cost``, -2 ln L at its
    levels; ``converged``, whether the search ended within tolerance before
    ``maxiter``; and ``evaluations``, the number of filter runs.
    """
    kw = {} if burnin is None else {"burnin": as_int(burnin, "burnin")}
    r = jl.tune_q(to_jl_data(data), m._jl(), maxiter=as_int(maxiter, "maxiter"), **kw)
    return {
        "model": _model_from_jl(r.model),
        "cost": float(r.cost),
        "converged": bool(r.converged),
        "evaluations": int(r.evaluations),
    }


@_map_errors
def simulate_clock(
    m: ClockModel, n: int, *, rng: int | None = None
) -> tuple[PhaseData, np.ndarray]:
    """A phase record of ``n`` samples drawn from clock model ``m``, with the true
    states behind it: one row per sample and the columns of ``KalmanResult.x``.
    ``rng`` is an integer seed for a reproducible record."""
    d, x = jl.simulate_clock(m._jl(), as_int(n, "n"), **_rng(rng))
    return data_from_jl(d), _mat(x)


# ── Steering ───────────────────────────────────────────────────────────


@dataclass(frozen=True)
class SteeringGains:
    """Gains of a steering loop. Every ``interval`` seconds the frequency of the
    clock is changed by −(``phase``·p̂ + ``frequency``·ŷ + ``integral``·Σp̂), with
    p̂ and ŷ the Kalman estimates of its phase and frequency against the
    reference and Σp̂ the sum of the phase estimates at the steers so far, this
    one included.

    ``phase`` and ``frequency`` are the gains g₁ and g₂ of Matsakis (2019), and
    ``integral`` is the integral gain of Matsakis (2020). ``phase`` and
    ``integral`` are in 1/s. Build the gains with :func:`steering_gains` or
    :func:`lqg_gains`.
    """

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
    """Gains of a steering loop by pole placement.

    A time constant ``T``, in seconds, gives a critically damped loop: its
    poles sit together at exp(−interval/T), and a disturbance dies away as
    exp(−t/T) without ringing. The gains are those of Matsakis (2019,
    Eqs. (3.10)–(3.11)) and, with ``integral=True`` and three poles, of
    Matsakis (2020, Eqs. (3.15)–(3.17)). A sequence in place of ``T`` holds
    the poles themselves in the z-plane: two for phase and frequency gains,
    three with an integral gain, complex ones in conjugate pairs.
    """
    if np.ndim(T) == 0:
        return _gains_from_jl(jl.steering_gains(float(T), float(interval), integral=bool(integral)))
    poles = juliacall.convert(_VecC64, [complex(z) for z in T])
    return _gains_from_jl(jl.steering_gains(poles, float(interval)))


@_map_errors
def lqg_gains(
    interval: float, *, phase: float, frequency: float = 0.0, steer: float
) -> SteeringGains:
    """Gains of the steering loop that minimizes the sum over the steers of
    ``phase``·p² + ``frequency``·f² + ``steer``·u²: the linear quadratic Gaussian
    design of Koppang and Leland (1999), with p and f the phase and frequency
    of the clock against its reference and u the steer.

    The gains do not depend on the noise of the clock or of the readings, and
    only the ratios of the weights matter. A small ``steer`` weight takes out
    the whole estimated error in one steer, and a large one gives a slow loop.
    The design has no integral gain.
    """
    return _gains_from_jl(
        jl.lqg_gains(
            float(interval), phase=float(phase), frequency=float(frequency), steer=float(steer)
        )
    )


@_map_errors
def steering_poles(g: SteeringGains) -> dict:
    """Where the gains put the poles of the closed loop (Matsakis 2019; Matsakis
    2020). Returns a dict with ``poles``, complex; ``time_constants``, for each
    pole the T (s) of its decay exp(−t/T), negative for a pole that grows;
    ``frequencies``, the frequency (Hz) that each pole rings at; and
    ``stable``, whether every pole is inside the unit circle."""
    p = jl.steering_poles(g._jl())
    return {
        "poles": np.array(p.poles, dtype=np.complex128),
        "time_constants": np_f64(p.time_constants),
        "frequencies": np_f64(p.frequencies),
        "stable": bool(p.stable),
    }


@dataclass(frozen=True)
class SteeringResult:
    """Result of :func:`steer` on a record. ``steered`` is the record as the loop
    would have left it: each reading plus the phase that the steers had added
    by then. ``steers`` is the frequency steer made at each sample, zero
    between steers, and ``correction`` their sum so far. ``phase`` and
    ``frequency`` are the filter's estimates of the steered clock."""

    gains: SteeringGains
    steered: PhaseData
    steers: np.ndarray
    correction: np.ndarray
    phase: np.ndarray
    frequency: np.ndarray


@dataclass(frozen=True)
class SteeringPrediction:
    """Result of :func:`steer` on a clock model, the steady state of the loop: the
    RMS time error ``phase`` (s) of the steered clock against its reference,
    its RMS frequency error ``frequency``, the RMS frequency steer ``steers``,
    and the Allan deviation ``adev`` of the steered clock at ``tau`` (s)."""

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
    """Steer a clock to its reference with the Kalman filter for ``model`` and the
    gains of :func:`steering_gains` or :func:`lqg_gains`.

    ``steer(data, model, gains)`` replays a record through the loop. ``data``
    is the phase of the free-running clock against the reference. The filter
    follows the clock, and every ``gains.interval`` seconds, a whole multiple
    of the sample interval, the frequency of the clock is changed as
    ``SteeringGains`` states (Koppang and Leland 1999; Matsakis 2019).

    ``steer(model, gains)`` gives the steady state that the loop settles to on
    a clock with the noise of ``model``, for a steer at every sample
    (Matsakis 2019; Matsakis 2020). ``tau`` sets the averaging times of its
    Allan deviation. The figures are for the clock itself; with
    ``readings=True`` its ``phase`` and ``adev`` are those of the readings,
    with their noise, as a steered record shows them. Random-run FM has no
    steady state.
    """
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
    """Phase readings of several clocks against one reference, taken at the same
    instants: the input of :func:`timescale`.

    Column ``i`` of ``x`` is clock i minus the reference, in seconds, with one
    row per sample and ``nan`` where a reading is missing. The reference may
    be one of the clocks, whose column is then zero, or an outside one.
    ``names`` holds one name per clock, all different, ``"clock 1"``,
    ``"clock 2"``, … by default, and ``reference`` says what the clocks were
    read against. ``cc[i]`` is the record of clock i as ``PhaseData``, and
    ``cc[i, j]`` the difference clock i − clock j, in which the reference
    cancels; both count from 0. :meth:`from_records` takes one phase record
    per clock in place of the matrix.
    """

    x: np.ndarray
    tau0: float
    names: tuple[str, ...] = ()
    reference: str = ""

    def __post_init__(self):
        x = np.array(self.x, dtype=np.float64)
        if x.ndim != 2:
            raise ValueError("ClockComparisons: x must be a samples × clocks matrix")
        if x.shape[1] < 2:
            raise ValueError("ClockComparisons: needs two or more clocks")
        if not self.tau0 > 0:
            raise ValueError(f"ClockComparisons: tau0 must be positive, got {self.tau0}")
        names = tuple(self.names) or tuple(f"clock {i + 1}" for i in range(x.shape[1]))
        if len(names) != x.shape[1]:
            raise ValueError("ClockComparisons: one name per clock")
        if len(set(names)) != len(names):
            raise ValueError("ClockComparisons: clock names must differ")
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "names", names)

    @classmethod
    def from_records(
        cls, records: Sequence[PhaseData], *, names: Sequence[str] = (), reference: str = ""
    ) -> ClockComparisons:
        """Comparisons from one phase record per clock, each read against the same
        reference at the same instants: the records share one length and one τ₀."""
        records = list(records)
        if not records:
            raise ValueError("ClockComparisons: no records")
        tau0 = records[0].tau0
        if not all(np.isclose(r.tau0, tau0, rtol=1e-9, atol=0.0) for r in records):
            raise ValueError("ClockComparisons: the records must share one τ₀")
        if len({len(r.x) for r in records}) != 1:
            raise ValueError("ClockComparisons: the records must have the same length")
        return cls(np.column_stack([r.x for r in records]), tau0, tuple(names), reference)

    def _jl(self) -> Any:
        return jl.ClockComparisons(
            juliacall.convert(_MatF64, np.ascontiguousarray(self.x)),
            float(self.tau0),
            names=juliacall.convert(VecStr, list(self.names)),
            reference=self.reference,
        )

    def __getitem__(self, ij: int | tuple[int, int]) -> PhaseData:
        """``cc[i]``, one clock against the reference, or ``cc[i, j]``, clock i − clock j."""
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
    """``n`` samples of an ensemble in which clock i follows ``models[i]``, with the
    true states behind the readings.

    Each clock is read against clock ``reference``, counted from 0, whose own
    column is zero, or with ``reference=None`` against a perfect outside
    reference. Returns the comparisons and the true states ``truth[k, s, i]``
    of clock i at sample k against true time: phase (s = 0), frequency (1)
    and drift (2), which is zero for a model without drift. ``rng`` is an
    integer seed for a reproducible ensemble.
    """
    ms = juliacall.convert(_VecModel, [m._jl() for m in models])
    ref = 0 if reference is None else 1 + _clock(reference, len(ms), "simulate_ensemble: reference")
    cc, truth = jl.simulate_ensemble(ms, as_int(n, "n"), reference=ref, **_rng(rng))
    return _comparisons_from_jl(cc), _mat(truth)


@_map_errors
def clock_models(
    cc: ClockComparisons, *, alphas: Sequence[int] | None = None, flicker: bool = False
) -> list[ClockModel]:
    """One ``ClockModel`` per clock, fitted from the comparisons alone, for the
    Kalman time scales of :func:`timescale`. This takes three or more clocks.

    The Allan deviation of each clock is separated from the others' by the
    Groslambert covariance of the clock differences (``groslambert``), in
    which the reference cancels. ``powerlaw_fit`` fits the components
    ``alphas`` to that curve, and :func:`clock_model` turns the fit into a
    model. The default components are white PM, flicker PM, white FM and
    random-walk FM, and the flicker PM is folded into ``R``. With
    ``flicker=True``, flicker FM joins the fit and each model keeps its
    flicker FM and flicker PM as Markov states.
    """
    kw: dict[str, Any] = {} if alphas is None else {"alphas": to_jl_ints(alphas, "an exponent α")}
    if flicker:
        kw["flicker"] = True
    return [_model_from_jl(m) for m in jl.clock_models(cc._jl(), **kw)]


@dataclass(frozen=True)
class TimescaleResult:
    """Result of :func:`timescale`, with one row per sample and one column per
    clock.

    ``x`` is each clock minus the scale (s), ``y`` its frequency against the
    scale as the algorithm estimated it, ``weights`` its weight and ``error``
    its prediction error (s), ``nan`` where it had no reading. ``xhat`` is its
    phase against the scale as the algorithm estimated it: the phase state of
    the filter for the Kalman scales, where ``x`` is read at the pivot, and
    ``x`` itself for the others.

    The weights of a sample sum to 1, or are zero where no clock of the scale
    is read. ``"jones_tryon"`` defines none, and its weights are ``nan``.
    Before a clock's first reading its ``x``, ``y`` and ``xhat`` are ``nan``.
    :func:`timescale_phase` gives the scale itself.
    """

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
    step_threshold: float = 5.0,
    lmax: int | None = None,
    step_inflate: float = 4.0,
) -> TimescaleResult:
    """A time scale from the clocks of ``cc``: the phase and frequency of each
    clock against the scale, and its weight, at every sample.

    ``algorithm`` is one of

    - ``"ewfa"``: equal-weights frequency averaging (Kettering 2023, §3.1);
    - ``"at1"``: the NIST AT1 scale as Weiss and Weissert (1991) give it, with
      weights from the prediction errors of the clocks;
    - ``"at2"``: AT1 with the frequency of each clock from a Kalman filter,
      from the same paper;
    - ``"greenhall"``: one Kalman filter over all the clocks with the
      covariance reduction of Greenhall (2003), and with noise on the readings
      the IEMx algorithm of Greenhall (2006);
    - ``"jones_tryon"``: the same filter without the reduction (Jones and
      Tryon 1983), for which no weights are reported;
    - ``"stein"``: the scale of Stein (2003), which Levine (2012) describes as
      KAS-2, with a Kalman filter for the difference of each clock from the
      pivot.

    The first three learn what they need from the comparisons. The last three
    take one ``ClockModel`` per clock in ``models``, which :func:`clock_models`
    fits from ``cc``, and are read at clock ``pivot``, counted from 0.

    For ``"at1"`` and ``"at2"``, ``tau_min`` is the averaging time in seconds
    at which the Allan deviation of a clock is smallest, one number or one per
    clock, estimated from the record by default. ``memory`` is the time
    constant in seconds of the average behind the weights, 20 days by default
    and a quarter of the record where that is shorter. ``bias`` is ``"paper"``
    or ``"exact"``, the correction for a clock being part of the scale that it
    is compared with. ``outlier`` leaves a clock out of a sample when its
    prediction error is more than that many times its ε from the median error
    of the clocks read, and is off by default. ``wmax`` is the largest weight
    that a clock may take in ``"at1"``, ``"at2"`` and ``"stein"``.

    ``steps``, for ``"at2"``, looks for steps in the frequency of each clock,
    after Weiss and Weissert (1991, §4): ``"best"``, or ``True``, forms the
    scale again from the step on, and ``"real"`` keeps the values already
    issued and corrects those that follow. ``step_threshold`` is the number of
    standard deviations that a frequency difference must pass to count as a
    step. ``lmax`` is one more than the longest window of the test, in
    intervals, by default τ_min over the interval for each clock and at least
    4. ``step_inflate`` is the factor on the error variance of a clock when it
    returns after a step.

    The documentation of the Julia function gives the equations of each scale,
    what each does with a clock that drops out and returns, and where each
    departs from its source.
    """
    kw = {
        "algorithm": jl.Symbol(str(algorithm)),
        "bias": jl.Symbol(str(bias)),
        "wmax": float(wmax),
        "step_threshold": float(step_threshold),
        "step_inflate": float(step_inflate),
    }
    if lmax is not None:
        kw["lmax"] = as_int(lmax, "lmax")
    if models is not None:
        kw["pivot"] = 1 + _clock(pivot, cc.x.shape[1], "timescale: pivot")
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
    """The reference of the comparisons against the scale, as a phase record for
    the usual stability analysis; ``nan`` where no clock was read.

    Its deviations are those of the scale as far as the reference is the more
    stable of the two. Against the ideal reference of ``simulate_ensemble``
    with ``reference=None`` it is the error of the scale, with the sign
    reversed.
    """
    c = r.comparisons.x
    ok = np.isfinite(c)
    first = ok.argmax(axis=1)
    rows = np.arange(c.shape[0])
    u = np.where(ok.any(axis=1), r.x[rows, first] - c[rows, first], np.nan)
    return PhaseData(u, r.comparisons.tau0, source=r.algorithm.upper())


def implicit_mean(r: TimescaleResult, truth: np.ndarray) -> PhaseData:
    """The implicit ensemble mean of a ``"greenhall"`` scale against true time, for
    a simulated ensemble whose true phases are known.

    It is the weighted mean Σ aᵢ·(xᵢ − x̂ᵢ) of the corrected clocks: the true
    phase of each clock, from the ``truth`` of :func:`simulate_ensemble`, less
    the filter's estimate of it, ``r.xhat``, with the weights of the scale
    (Greenhall 2006, Eqs. (3.8) and (3.13)). The mean carries the reading
    noise of no clock, so its stability is that of the scale itself. A
    measurement cannot give it: :func:`timescale_phase` reads the scale at one
    clock and adds the noise of that clock's reading.
    """
    if r.algorithm != "greenhall":
        raise ValueError(
            f"implicit_mean: the implicit mean belongs to greenhall, not {r.algorithm}"
        )
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
    "SteeringGains",
    "steering_gains",
    "lqg_gains",
    "steering_poles",
    "SteeringResult",
    "SteeringPrediction",
    "steer",
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

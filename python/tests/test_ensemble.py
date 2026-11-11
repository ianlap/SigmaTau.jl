"""Clock models, the noise fit, the Kalman filter and timescales: each wrapper
agrees with a direct Julia call."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import _convert


@pytest.fixture(scope="module")
def clock():
    m = st.ClockModel(1.0, R=1e-22, q1=1e-24, q2=1e-30)
    data, truth = st.simulate_clock(m, 2000, rng=1)
    return m, data, truth


def test_clock_model_matches_julia(jl):
    m = st.clock_model({"WHPM": 8e-22, 0: 2e-24, -2: 1e-30}, tau0=1.0)
    direct = jl.clock_model(_convert.to_jl_dict_int_f64({2: 8e-22, 0: 2e-24, -2: 1e-30}), tau0=1.0)
    assert (m.R, m.q1, m.q2, m.drift) == (
        float(direct.R),
        float(direct.q1),
        float(direct.q2),
        False,
    )
    np.testing.assert_array_equal(
        st.state_transition(st.ClockModel(1.0, q3=1e-40), 2.0), [[1, 2, 2], [0, 1, 2], [0, 0, 1]]
    )
    np.testing.assert_allclose(
        st.process_noise(st.ClockModel(1.0, q1=1.0, q2=1.0), 1.0), [[4 / 3, 0.5], [0.5, 1]]
    )
    with pytest.raises(ValueError):
        st.clock_model({0: 1e-24})  # no tau0


def test_flicker_states_match_julia(jl):
    h = {1: 2.2e-22, 0: 2e-24, -1: 1.6e-26, -2: 1.7e-31}
    m = st.clock_model(h, tau0=1.0, flicker=True)
    direct = jl.clock_model(_convert.to_jl_dict_int_f64(h), tau0=1.0, flicker=True)
    assert (m.flicker_fm, m.fm_band, m.flicker_pm, m.pm_band) == (
        float(direct.flicker_fm),
        tuple(direct.fm_band),
        float(direct.flicker_pm),
        tuple(direct.pm_band),
    )
    np.testing.assert_array_equal(st.process_noise(m, 1.0), np.array(jl.process_noise(direct, 1.0)))
    assert st.clock_model(h, tau0=1.0).flicker_fm == 0  # folded into R unless asked for
    assert st.ClockModel(1.0, flicker_fm=1e-26, tau_max=1e3).fm_band == (1.0, 1000.0)
    with pytest.raises(ValueError):
        st.ClockModel(1.0, fm_band=(10.0, 100.0))  # a band without its flicker
    data, truth = st.simulate_clock(m, 500, rng=1)
    n = st.process_noise(m, 1.0).shape[0]
    assert n > 2 and st.kalman_filter(data, m).x.shape == truth.shape == (500, n)
    # An ensemble with flicker PM on two channels: the fit, and the scales that take it.
    models = [st.ClockModel(1.0, q1=1e-24, q2=1e-30, flicker_pm=p) for p in (0.0, 2e-22, 2e-21)]
    cc, _ = st.simulate_ensemble(models, 8192, reference=None, rng=3)
    fitted = st.clock_models(cc, flicker=True)
    direct = jl.seval("c -> clock_models(c; flicker=true)")(cc._jl())
    assert [(f.flicker_pm, f.pm_band) for f in fitted] == [
        (float(d.flicker_pm), tuple(d.pm_band)) for d in direct
    ]
    assert fitted[2].pm_band[0] == 1.0 and st.clock_models(cc)[2].flicker_pm == 0
    r = st.timescale(cc, algorithm="greenhall", models=fitted)
    jr = jl.seval("(c, ms) -> timescale(c; algorithm=:greenhall, models=ms)")(cc._jl(), direct)
    np.testing.assert_array_equal(r.weights, np.array(jr.weights))
    with pytest.raises(ValueError):
        st.timescale(cc, algorithm="stein", models=fitted)


def test_noise_fit_builds_and_undoes():
    r = st.adev(st.noise_gen(st.PhaseData, 4096, 1.0, h={0: 1e-24}, rng=3))
    f = st.noise_fit(r).add_component("WHFM", slice(0, 6))
    ((alpha, h, lo, hi),) = f.components
    assert (alpha, lo, hi) == (0, r.tau[0], r.tau[5]) and h == pytest.approx(1e-24, rel=0.2)
    assert st.clock_model(f).q1 == pytest.approx(h / 2)
    assert f.residual().shape == r.tau.shape
    assert f.undo().components == []
    with pytest.raises(ValueError):
        f.add_component("WHFM", slice(3, 3))


def test_kalman_filter_matches_julia(clock, jl):
    m, data, truth = clock
    kf = st.kalman_filter(data, m)
    direct = jl.kalman_filter(_convert.to_jl_data(data), m._jl())
    np.testing.assert_array_equal(kf.x, np.array(direct.x))
    assert kf.P.shape == (2, 2, 2000) and not kf.smoothed
    assert np.sqrt(np.mean((kf.x[500:, 1] - truth[500:, 1]) ** 2)) < 3 * np.sqrt(
        np.mean(kf.P[1, 1, 500:])
    )
    assert st.kalman_smooth(kf).smoothed
    assert st.innovation_check(kf)["ok"]
    assert np.all(np.diff(st.holdover_error(kf, [10.0, 100.0, 1000.0])) > 0)
    assert st.prediction_error(kf, data, [1, 10]).shape == (2,)
    t = st.tune_q(data, m, points=3)
    assert t["cost"].shape == (3, 3) and t["model"].R == m.R


def test_holdover_matches_julia(clock, jl):
    m, data, _ = clock
    h = st.holdover(data, train=200, tau=[10, 100], limit=1e-10)
    direct = jl.holdover(
        _convert.to_jl_data(data), train=200, tau=_convert.to_jl_f64([10, 100]), limit=1e-10
    )
    np.testing.assert_array_equal(h.rms, np.array(direct.rms))
    assert h.n.tolist() == [1791, 1701] and h.time_to_limit == float(direct.time_to_limit)
    assert (h.train, h.anchor, h.drift, h.confidence) == (200.0, True, False, 0.95)
    assert np.isnan(st.holdover(data).limit)
    with pytest.raises(ValueError):
        st.holdover(data, train=0.5)
    # A stretch of twice the holdover time: at τ₀ the second difference behind TDEV.
    g = st.holdover(data, grow=2, tau=[1, 10])
    assert g.grow == 2 and np.isnan(g.train) and h.grow == 0
    np.testing.assert_allclose(g.rms[0], np.sqrt(6) * st.tdev(data, [1], ci=False).dev[0])
    assert st.holdover(m, grow=2, tau=[1, 10]).grow == 2


def test_holdover_with_a_kalman_filter(clock, jl):
    m, data, _ = clock
    k = st.holdover(data, m, train=200, tau=[10, 100])
    direct = jl.holdover(
        _convert.to_jl_data(data), m._jl(), train=200, tau=_convert.to_jl_f64([10, 100])
    )
    np.testing.assert_array_equal(k.rms, np.array(direct.rms))
    np.testing.assert_array_equal(k.claimed, np.array(direct.claimed))
    assert k.train == 200.0 and st.holdover(data, train=200).claimed.size == 0


def test_holdover_prediction_matches_julia(clock, jl):
    m, _, _ = clock
    p = st.holdover(m, train=200, tau=[10, 100])
    direct = jl.holdover(m._jl(), train=200, tau=_convert.to_jl_f64([10, 100]))
    np.testing.assert_array_equal(p.sigma, np.array(direct.sigma))
    assert p.alphas.tolist() == [2, 0, -2] and p.budget.shape == (2, 3)
    np.testing.assert_allclose(p.budget.sum(axis=1), p.sigma**2)
    np.testing.assert_allclose(
        st.holdover(
            {"WHPM": 8 * np.pi**2 * m.R, 0: 2 * m.q1, -2: m.q2 / (2 * np.pi**2)},
            tau0=1.0,
            train=200,
            tau=[10, 100],
        ).sigma,
        p.sigma,
    )
    with pytest.raises(ValueError):
        st.holdover({0: 1e-24}, train=200)  # no tau0


def test_steering_matches_julia(clock, jl):
    m, data, _ = clock
    g = st.steering_gains(30.0, 1.0)
    r = np.exp(-1 / 30)
    np.testing.assert_allclose([g.phase, g.frequency, g.integral], [(1 - r) ** 2, 1 - r**2, 0])
    p = st.steering_poles(g)
    assert p["stable"] and p["poles"].dtype == np.complex128
    np.testing.assert_allclose(p["time_constants"], [30.0, 30.0], rtol=1e-6)
    k = st.steering_gains([0.9, 0.8 + 0.1j, 0.8 - 0.1j], 2.0)
    np.testing.assert_allclose(
        np.sort_complex(st.steering_poles(k)["poles"]), [0.8 - 0.1j, 0.8 + 0.1j, 0.9]
    )
    s = st.steer(data, m, g)
    direct = jl.steer(_convert.to_jl_data(data), m._jl(), g._jl())
    np.testing.assert_array_equal(s.steered.x, np.array(direct.steered.x))
    np.testing.assert_array_equal(s.steers, np.array(direct.steers))
    assert s.steered.tau0 == 1.0 and s.gains == g
    q = st.steer(m, g, tau=[1, 16])
    assert (q.phase, q.steers) == tuple(
        float(v) for v in (jl.steer(m._jl(), g._jl()).phase, jl.steer(m._jl(), g._jl()).steers)
    )
    assert q.tau.tolist() == [1.0, 16.0] and q.adev.shape == (2,)
    z = st.steer(m, g, readings=True)  # the readings carry their own noise
    assert z.phase == float(jl.steer(m._jl(), g._jl(), readings=True).phase) > q.phase
    with pytest.raises(ValueError):
        st.steer(m, st.steering_gains(30.0, 10.0))
    with pytest.raises(TypeError):
        st.steer(data, m)
    # Koppang & Leland (1999), §IV, first design, as printed.
    k = st.lqg_gains(3600.0, phase=1e-3, frequency=1e-3, steer=1e5)
    np.testing.assert_allclose([k.phase / 6.503e-5, k.frequency / 0.5771], [1, 1], rtol=5e-4)
    assert k.integral == 0


def test_timescale_matches_julia(jl):
    models = [st.ClockModel(1.0, q1=q) for q in (1e-24, 4e-24, 9e-24)]
    cc, truth = st.simulate_ensemble(models, 3000, reference=None, rng=2)
    assert cc.names == ("clock 1", "clock 2", "clock 3") and truth.shape[0] == 3000
    np.testing.assert_array_equal(cc[0, 1].x, cc.x[:, 0] - cc.x[:, 1])
    r = st.timescale(cc)
    direct = jl.timescale(cc._jl())
    np.testing.assert_array_equal(r.weights, np.array(direct.weights))
    assert r.weights[-1, 0] > r.weights[-1, 2]
    np.testing.assert_array_equal(st.timescale_phase(r).x, np.array(jl.timescale_phase(direct).x))
    assert np.all(st.timescale(cc, algorithm="ewfa").weights == pytest.approx(1 / 3))
    for mode in ("best", "real"):
        assert st.timescale(cc, algorithm="at2", steps=mode).weights.shape == (3000, 3)
    with pytest.raises(Exception):
        st.timescale(cc, steps="best")
    ex = st.timescale(cc, bias="exact", tau_min=[50, 100, 100])
    direct = jl.seval("c -> timescale(c; bias=:exact, tau_min=[50.0, 100.0, 100.0])")(cc._jl())
    np.testing.assert_array_equal(ex.weights, np.array(direct.weights))
    rk = st.timescale(cc, algorithm="greenhall", models=models)
    w = 1 / np.array([1.0, 4.0, 9.0])
    np.testing.assert_allclose(rk.weights[-1], w / w.sum(), rtol=0.02)
    with pytest.raises(ValueError):
        st.timescale(cc, algorithm="at3")


def test_clock_models_and_implicit_mean_match_julia(jl):
    models = [st.ClockModel(1.0, R=1e-22, q1=q) for q in (1e-24, 4e-24, 9e-24)]
    cc, truth = st.simulate_ensemble(models, 4000, reference=None, rng=5)
    fitted = st.clock_models(cc)
    direct = jl.clock_models(cc._jl())
    assert [m.q1 for m in fitted] == [float(m.q1) for m in direct]
    r = st.timescale(cc, algorithm="greenhall", models=fitted)
    jr = jl.seval("(c, ms) -> timescale(c; algorithm=:greenhall, models=ms)")(cc._jl(), direct)
    iem = st.implicit_mean(r, truth)
    np.testing.assert_allclose(iem.x, np.array(jl.implicit_mean(jr, truth).x), rtol=1e-12)
    with pytest.raises(ValueError):
        st.implicit_mean(st.timescale(cc), truth)

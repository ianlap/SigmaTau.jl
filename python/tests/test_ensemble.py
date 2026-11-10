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
    _, data, _ = clock
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

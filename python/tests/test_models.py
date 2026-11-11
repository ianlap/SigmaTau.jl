"""Drift and power-law models, the two-channel estimators and the Groslambert
covariance: each wrapper agrees bitwise with a direct Julia call."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import FrequencyData, PhaseData, _convert


@pytest.mark.parametrize("model", st.DRIFT_MODELS)
def test_drift_fit_parity(small_freq, model, jl):
    f = st.drift_fit(small_freq, model=model)
    direct = _convert.driftfit_from_jl(
        jl.drift_fit(_convert.to_jl_data(small_freq), model=jl.Symbol(model))
    )
    assert f.model == model and f.coef.size == (2 if model == "linear" else 3)
    np.testing.assert_array_equal(f.coef, direct.coef)
    np.testing.assert_array_equal(f.fit, direct.fit)
    assert f.rate == direct.rate and f.rms == direct.rms


@pytest.mark.parametrize("method", st.DRIFT_METHODS)
def test_drift_fit_methods(method, jl):
    y = 1e-9 + 2e-15 * np.arange(1000.0) + 1e-14 * np.sin(np.arange(1000.0) / 9)
    fd = FrequencyData(y, 1.0)
    f = st.drift_fit(fd, method=method)
    direct = jl.drift_fit(_convert.to_jl_data(fd), method=jl.Symbol(method))
    assert f.method == method and f.rate == float(direct.rate)
    assert f.rate == pytest.approx(2e-15 * 86400, rel=0.05)
    # No standard error from the record for these two.
    assert np.isnan(f.rate_stderr) == (method in ("bisection", "four_point"))
    # detrend by the fit removes the model: the residuals of the fit remain.
    np.testing.assert_allclose(st.detrend(fd, f).y, f.y - f.fit, rtol=0, atol=1e-25)
    with pytest.raises(ValueError):
        st.drift_fit(fd, model="log", method="bisection")


def test_drift_fit_recovers_a_ramp():
    y = 1e-9 + 2e-15 * np.arange(1000.0)
    f = st.drift_fit(FrequencyData(y, 1.0))
    assert f.rate == pytest.approx(2e-15 * 86400, rel=1e-9)
    assert st.drift_fit(
        PhaseData(np.concatenate([[0.0], np.cumsum(y)]), 1.0)
    ).rate == pytest.approx(f.rate)
    with pytest.raises(ValueError):
        st.drift_fit(FrequencyData(y, 1.0), model="cubic")


def test_powerlaw_fit_record_and_roundtrip(jl):
    pd = st.noise_gen(PhaseData, 4096, 1.0, h={0: 1e-24}, rng=11)
    f = st.powerlaw_fit(pd, deviation="adev", alphas=[2, 0, -2])
    direct = _convert.powerlaw_from_jl(
        jl.powerlaw_fit(
            _convert.to_jl_data(pd),
            deviation=jl.Symbol("adev"),
            alphas=_convert.to_jl_ints([2, 0, -2]),
        )
    )
    assert f.source == "adev" and list(f.alphas) == [2, 0, -2]
    np.testing.assert_array_equal(f.h, direct.h)
    np.testing.assert_array_equal(f.fit, direct.fit)
    assert f.h[1] == pytest.approx(1e-24, rel=0.3)
    # The model's σ(τ) at the fitted τ is the fit itself.
    np.testing.assert_allclose(st.powerlaw_sigma(f, f.x), f.fit, rtol=1e-12)
    np.testing.assert_allclose(
        st.powerlaw_sigma({"WHFM": 2e-24}, [1.0, 4.0], tau0=1.0),
        st.powerlaw_sigma({0: 2e-24}, [1.0, 4.0], tau0=1.0),
    )
    with pytest.raises(ValueError):
        st.powerlaw_sigma({0: 1e-24}, [1.0])  # no tau0
    psd = st.powerlaw_psd({0: 1e-24, -2: 1e-26}, [0.01, 0.1, 1.0])
    np.testing.assert_allclose(psd, 1e-24 + 1e-26 / np.array([0.01, 0.1, 1.0]) ** 2)
    assert st.powerlaw_psd({0: 1e-24}, [1.0], quantity="L", f_carrier=1e7).shape == (1,)


def test_powerlaw_fit_spectrum_forms(small_phase, jl):
    s = st.psd_y(small_phase)
    f = st.powerlaw_fit(s, alphas=[2, 0])
    jl_s = jl.psd_y(_convert.to_jl_data(small_phase), taper=jl.Symbol("hann"))
    direct = _convert.powerlaw_from_jl(jl.powerlaw_fit(jl_s, alphas=_convert.to_jl_ints([2, 0])))
    assert f.source == "Sy" and f.tau0 == direct.tau0
    np.testing.assert_array_equal(f.h, direct.h)
    pts = st.powerlaw_fit(
        [1, 10, 100, 1e3], [-95, -125, -140, -150], quantity="L", f_carrier=1e7, alphas=[2, 0, -2]
    )
    assert pts.source == "L" and np.isnan(pts.tau0) and pts.h.size == 3
    with pytest.raises(TypeError):
        st.powerlaw_fit([1.0, 2.0])


def test_powerlaw_fit_of_a_result(jl):
    pd = st.noise_gen(PhaseData, 4096, 1.0, h={0: 1e-24, -2: 1e-30}, rng=11)
    f = st.powerlaw_fit(st.mhdev(pd))
    direct = _convert.powerlaw_from_jl(
        jl.seval("d -> powerlaw_fit(mhdev(d))")(_convert.to_jl_data(pd))
    )
    assert f.source == "mhdev" and f.tau0 == 1.0
    np.testing.assert_array_equal(f.h, direct.h)
    np.testing.assert_array_equal(f.fit, direct.fit)
    with pytest.raises(ValueError):
        st.powerlaw_fit(st.totdev(pd, [1, 2, 4]))  # no model of the total deviation


def test_powerlaw_curves_of_a_model_and_a_noise_fit(jl):
    tau, f = [1.0, 16.0, 256.0], [1e-3, 1e-2, 1e-1]
    jtau, jf = _convert.to_jl_f64(tau), _convert.to_jl_f64(f)
    m = st.ClockModel(1.0, R=1e-22, q1=1e-24, q2=1e-30, flicker_fm=1.6e-26)
    np.testing.assert_array_equal(
        st.powerlaw_sigma(m, tau), np.array(jl.powerlaw_sigma(m._jl(), jtau))
    )
    np.testing.assert_array_equal(st.powerlaw_psd(m, f), np.array(jl.powerlaw_psd(m._jl(), jf)))
    # White FM alone: σ_y(τ) = √(q₁/τ) and S_y = 2q₁.
    white = st.ClockModel(1.0, q1=1e-24)
    np.testing.assert_allclose(st.powerlaw_sigma(white, [1.0, 4.0]), [1e-12, 5e-13], rtol=1e-12)
    np.testing.assert_allclose(st.powerlaw_psd(white, f), 2e-24, rtol=1e-12)
    # A noise fit is read as its components, at the sample interval of its result.
    nf = st.noise_fit(st.adev(st.noise_gen(PhaseData, 4096, 1.0, h={0: 1e-24}, rng=3)))
    h = nf.add_component("WHFM", slice(0, 6)).components[0][1]
    np.testing.assert_allclose(
        st.powerlaw_sigma(nf, tau), st.powerlaw_sigma({0: h}, tau, tau0=1.0), rtol=1e-12
    )
    np.testing.assert_array_equal(st.powerlaw_psd(nf, f), np.array(jl.powerlaw_psd(nf._jl, jf)))


@pytest.fixture(scope="module")
def two_channels():
    dut = st.noise_gen(FrequencyData, 2048, 1.0, h={0: 1e-24}, rng=2).y
    r1 = st.noise_gen(FrequencyData, 2048, 1.0, h={0: 1e-23}, rng=3).y
    r2 = st.noise_gen(FrequencyData, 2048, 1.0, h={0: 1e-23}, rng=4).y
    return FrequencyData(dut - r1, 1.0), FrequencyData(dut - r2, 1.0)


def test_cross_adev_parity(two_channels, jl):
    a, b = two_channels
    ms = [1, 2, 4, 8]
    x = st.cross_adev(a, b, ms)
    direct = _convert.cross_stability_from_jl(
        jl.cross_adev(_convert.to_jl_data(a), _convert.to_jl_data(b), _convert.to_jl_ints(ms))
    )
    np.testing.assert_array_equal(x.covar, direct.covar)
    np.testing.assert_array_equal(x.uncertainty, direct.uncertainty)
    assert list(x.tau) == ms and x.n == 2049
    # The device's own ADEV at 1 s sits well below either channel's.
    assert x.dev[0] < 0.6 * x.dev_a[0]
    assert st.cross_adev(a, b).tau[0] == 1.0
    with pytest.raises(ValueError):
        st.cross_adev(a, FrequencyData(b.y[:-1], 1.0), ms)


@pytest.mark.parametrize("name", ["cross_psd_y", "cross_psd_x"])
def test_cross_psd_parity(two_channels, name, jl):
    a, b = two_channels
    s = getattr(st, name)(a, b, nperseg=256)
    direct = _convert.cross_spectral_from_jl(
        getattr(jl, name)(
            _convert.to_jl_data(a), _convert.to_jl_data(b), nperseg=256, taper=jl.Symbol("hann")
        )
    )
    np.testing.assert_array_equal(s.psd, direct.psd)
    np.testing.assert_array_equal(s.floor, direct.floor)
    assert s.nperseg == 256 and s.nseg == direct.nseg
    assert s.spectral_type == ("Sy" if name.endswith("y") else "Sx")


def test_groslambert_forms(jl):
    rng = np.random.default_rng(7)
    x = [s * rng.standard_normal(1024) for s in (1.0, 2.0, 0.5)]
    d = [PhaseData(x[0] - x[1], 1.0), PhaseData(x[1] - x[2], 1.0), PhaseData(x[2] - x[0], 1.0)]
    ms = [1, 2, 4]
    g = st.groslambert(*d, m=ms)
    direct = [
        _convert.result_from_jl(r)
        for r in jl.groslambert(*(_convert.to_jl_data(r) for r in d), _convert.to_jl_ints(ms))
    ]
    assert len(g) == 3 and g[0].deviation_type == "adev" and g[0].edf.size == 3
    for a, b in zip(g, direct, strict=True):
        np.testing.assert_array_equal(a.dev, b.dev)
        np.testing.assert_array_equal(a.ci_lower, b.ci_lower)
    matrix = [[None, d[0], None], [None, None, d[1]], [d[2], None, None]]
    for a, b in zip(g, st.groslambert(matrix, m=ms), strict=True):
        np.testing.assert_array_equal(a.dev, b.dev)
    assert st.groslambert(*d, m=ms, confidence=0.95)[0].confidence == 0.95
    k = st.groslambert(*d, m=ms, interval="klts")
    np.testing.assert_array_equal(k[1].dev, g[1].dev)
    assert np.all(k[1].ci_lower <= k[1].dev) and np.all(k[1].dev <= k[1].ci_upper)
    with pytest.raises(TypeError):
        st.groslambert(d[0], d[1])

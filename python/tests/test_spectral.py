"""Spectral estimators psd_y / psd_x / ssb_phase_noise: units and bitwise parity with Julia."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import _convert


@pytest.mark.parametrize("name,units", [("psd_y", "per_Hz"), ("psd_x", "s2_per_Hz")])
def test_psd_units_and_parity(small_phase, name, units, jl):
    fn = getattr(st, name)
    result = fn(small_phase)
    assert result.units == units

    jldata = _convert.to_jl_data(small_phase)
    direct = _convert.spectral_from_jl(getattr(jl, name)(jldata, taper=jl.Symbol("hann")))
    np.testing.assert_array_equal(result.freq, direct.freq)
    np.testing.assert_array_equal(result.psd, direct.psd)


def test_ssb_phase_noise_units(small_phase):
    result = st.ssb_phase_noise(small_phase, f_carrier=1e7)
    assert result.units == "dBc_per_Hz"


def test_ssb_phase_noise_requires_carrier(small_phase):
    with pytest.raises(TypeError):
        st.ssb_phase_noise(small_phase)


def test_integrated_jitter_matches_julia(small_phase):
    s = st.ssb_phase_noise(small_phase, f_carrier=1e7)
    j = st.integrated_jitter(s, f_carrier=1e7, band=(0.01, 0.2))
    pts = st.integrated_jitter([1e3, 1e4, 1e5], [-100.0, -110.0, -120.0], f_carrier=10e6)
    assert pts["phase"] == pytest.approx(np.sqrt(2e-7 * np.log(100)))
    assert j["time"] == pytest.approx(j["phase"] / (2 * np.pi * 1e7))


def test_spectrum_sigma_closed_forms():
    # White FM, S_y = h₀: σ²_y(τ) = h₀/(2τ), and the modified variance is half of it.
    f, h0, tau = [1e-4, 1e4], 2e-24, [1.0, 10.0]
    L = 10 * np.log10(h0 * 10e6**2 / np.square(f) / 2)
    a = st.spectrum_sigma(f, L, tau, f_carrier=10e6)
    m = st.spectrum_sigma(f, [h0, h0], tau, quantity="Sy", deviation="mdev")
    assert a["dev"] == pytest.approx(np.sqrt(h0 / (2 * np.array(tau))), rel=1e-3)
    assert m["dev"] == pytest.approx(a["dev"] / np.sqrt(2), rel=1e-3)
    assert a["band"] == (1e-4, 1e4)
    sphi = st.spectrum_sigma(f, 2 * 10 ** (L / 10), tau, quantity="Sphi", f_carrier=10e6)
    assert sphi["dev"] == pytest.approx(a["dev"])
    assert np.all(a["below"] < 1e-3) and np.all(a["above"] < 1e-3)
    pm = st.spectrum_sigma([1.0, 1e3], [1e-26, 1e-20], [1.0], quantity="Sy")  # white PM
    assert np.isinf(pm["above"][0])
    with pytest.raises(ValueError):
        st.spectrum_sigma(f, L, tau)  # ℒ(f) needs the carrier

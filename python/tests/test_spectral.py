"""Spectral estimators Sy / Sx / L: units and bitwise parity with Julia."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import _convert


@pytest.mark.parametrize("name,units", [("Sy", "per_Hz"), ("Sx", "s2_per_Hz")])
def test_psd_units_and_parity(small_phase, name, units, jl):
    fn = getattr(st, name)
    result = fn(small_phase)
    assert result.units == units

    jldata = _convert.to_jl_data(small_phase)
    direct = _convert.spectral_from_jl(getattr(jl, name)(jldata, window=jl.Symbol("hann")))
    np.testing.assert_array_equal(result.freq, direct.freq)
    np.testing.assert_array_equal(result.psd, direct.psd)


def test_L_units(small_phase):
    result = st.L(small_phase, f_carrier=1e7)
    assert result.units == "dBc_per_Hz"


def test_L_requires_carrier(small_phase):
    with pytest.raises(TypeError):
        st.L(small_phase)

"""Telecom masks and the record check match the Julia package."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st


def test_telecom_mask_scalar_and_vector():
    assert st.telecom_mask("G8262_opt1", "tdev", 10.0) == pytest.approx(3.2e-9)
    lim = st.telecom_mask("G8262_opt1", "tdev", [10.0, 1e9])
    assert lim[0] == pytest.approx(3.2e-9) and np.isnan(lim[1])
    assert "G811" in st.TELECOM_MASKS


def test_mask_check(small_phase):
    c = st.mask_check(st.mtie(small_phase), "G811")
    assert c.mask == "G811" and c.tau.shape == c.limit.shape
    assert c.passed == (c.margin >= 1)
    with pytest.raises(ValueError):
        st.mask_check(st.adev(small_phase), "G811")


def test_summarize(small_phase):
    s = st.summarize(small_phase)
    assert s["n"] == small_phase.x.size and s["span"] == s["n"] * s["tau0"]
    np.testing.assert_array_equal(s["outliers"], st.find_outliers(small_phase))

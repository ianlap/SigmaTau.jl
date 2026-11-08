"""StabilitySuite behaviour: keys, indexing, iteration, ci=False, and parity."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st

MS = [1, 2, 4]


@pytest.fixture(scope="module")
def suite(small_phase):
    return st.stability(small_phase, m=MS)


def test_default_keys(suite):
    assert suite.keys() == st.DEFAULT_DEVIATIONS
    assert st.DEFAULT_DEVIATIONS == ("adev", "mdev", "hdev", "mhdev")


def test_index_by_name_is_index_by_position(suite):
    assert suite["adev"] is suite[0]


def test_iteration_order(suite):
    assert [r.deviation_type for r in suite] == list(st.DEFAULT_DEVIATIONS)


def test_contains(suite):
    assert "adev" in suite
    assert "nope" not in suite


def test_missing_key_raises(suite):
    with pytest.raises(KeyError):
        suite["nope"]


def test_ci_false_drops_confidence_and_edf(small_phase):
    suite = st.stability(small_phase, m=MS, ci=False)
    assert suite.confidence is None
    for r in suite:
        assert r.edf.size == 0


def test_devs_argument_honored_in_order(small_phase):
    suite = st.stability(small_phase, devs=("adev", "tdev"), m=MS)
    assert suite.keys() == ("adev", "tdev")


def test_suite_matches_standalone(suite, small_phase):
    standalone = st.adev(small_phase, MS)
    np.testing.assert_array_equal(suite["adev"].dev, standalone.dev)
    np.testing.assert_array_equal(suite["adev"].edf, standalone.edf)
    np.testing.assert_array_equal(suite["adev"].ci_lower, standalone.ci_lower)
    np.testing.assert_array_equal(suite["adev"].ci_upper, standalone.ci_upper)

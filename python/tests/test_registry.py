"""The estimator registry: deviation(name, ...) matches the named functions."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st


def test_estimators_lists_every_static_deviation():
    names = {e["name"] for e in st.estimators()}
    for name in ("adev", "mdev", "hdev", "mhdev", "totdev", "mhtotdev", "mtie", "theo1"):
        assert name in names
    first = st.estimators()[0]
    assert set(first) == {"name", "label", "family", "axis", "correct_bias", "description"}


@pytest.mark.parametrize("name", ["adev", "mdev", "hdev", "totdev"])
def test_deviation_matches_named_function(small_phase, name):
    by_name = st.deviation(name, small_phase, [1, 2, 4])
    direct = getattr(st, name)(small_phase, [1, 2, 4])
    np.testing.assert_array_equal(by_name.dev, direct.dev)
    np.testing.assert_array_equal(by_name.edf, direct.edf)


def test_deviation_unknown_name(small_phase):
    with pytest.raises(ValueError):
        st.deviation("nosuchdev", small_phase)

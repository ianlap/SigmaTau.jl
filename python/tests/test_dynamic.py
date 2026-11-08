"""Dynamic deviations dadev / dhdev: shapes, NaN columns, and parity.

window > N raises ValueError (surfaced from Julia, carrying .julia_error).
"""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import _convert

WINDOW = 128


def test_dadev_default_grid_shape(small_phase):
    result = st.dadev(small_phase, window=WINDOW)
    assert result.dev.shape == (result.t.size, result.tau.size)


def test_dadev_explicit_taus_column_count(small_phase):
    result = st.dadev(small_phase, [1, 2, 4], window=WINDOW)
    assert result.dev.shape[1] == 3


def test_dadev_smaller_step_gives_more_rows(small_phase):
    default = st.dadev(small_phase, [1, 2, 4], window=WINDOW)
    stepped = st.dadev(small_phase, [1, 2, 4], window=WINDOW, step=32)
    assert stepped.t.size > default.t.size


def test_dadev_unsupported_m_is_nan_column(small_phase):
    # 128 - 2*64 = 0 < 2, so the m=64 column is entirely NaN.
    result = st.dadev(small_phase, [1, 2, 4, 64], window=WINDOW)
    assert np.isnan(result.dev[:, 3]).all()
    assert not np.isnan(result.dev[:, 0]).any()


def test_dadev_parity_with_julia(small_phase, jl):
    ms = [1, 2, 4]
    result = st.dadev(small_phase, ms, window=WINDOW)
    jldata = _convert.to_jl_data(small_phase)
    direct = _convert.dynamic_from_jl(jl.dadev(jldata, _convert.to_jl_ints(ms), window=WINDOW))
    np.testing.assert_array_equal(result.dev, direct.dev)
    np.testing.assert_array_equal(result.t, direct.t)
    np.testing.assert_array_equal(result.tau, direct.tau)


def test_dhdev_shape(small_phase):
    result = st.dhdev(small_phase, [1, 2, 4], window=WINDOW)
    assert result.dev.shape == (result.t.size, 3)


def test_window_larger_than_record_raises(small_phase):
    with pytest.raises(ValueError):
        st.dadev(small_phase, [1, 2, 4], window=1000)

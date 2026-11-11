"""Error mapping: Julia exceptions surface as native Python errors.

read_phase on a missing path raises ValueError (an ArgumentError from Julia) and
keeps the original JuliaError reachable as ``.julia_error``.
"""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import PhaseData


def test_invalid_averaging_factor(small_phase):
    with pytest.raises(ValueError):
        st.adev(small_phase, [0])


def test_single_sample_rejected():
    with pytest.raises(ValueError):
        PhaseData(np.array([1.0]), 1.0)


def test_read_phase_missing_file():
    with pytest.raises((OSError, ValueError)) as exc:
        st.read_phase("/nonexistent/x.csv")
    # Surfaced from Julia, so the original JuliaError is preserved.
    assert hasattr(exc.value, "julia_error")


def test_suite_missing_key(small_phase):
    suite = st.stability(small_phase, m=[1, 2, 4])
    with pytest.raises(KeyError):
        suite["does_not_exist"]


def test_a_float_is_refused_where_julia_counts(small_phase):
    with pytest.raises(TypeError, match="averaging factor"):
        st.adev(small_phase, [2.7])
    with pytest.raises(TypeError, match="n must be an integer, got 16.5"):
        st.noise_gen(PhaseData, 16.5, 1.0, sigma1={0: 1e-12})
    with pytest.raises(TypeError, match="lags must be an integer, got 4.5"):
        st.acf(small_phase, lags=4.5)
    with pytest.raises(TypeError, match="value_col must be an integer, got 2.0"):
        st.read_phase("clock.csv", value_col=2.0)
    assert st.adev(small_phase, np.array([1, 2]), ci=False).tau.tolist() == [1.0, 2.0]


def test_numpy_integer_seed():
    by_numpy = st.noise_gen(PhaseData, 16, 1.0, sigma1={0: 1e-12}, rng=np.int64(3))
    np.testing.assert_array_equal(
        by_numpy.x, st.noise_gen(PhaseData, 16, 1.0, sigma1={0: 1e-12}, rng=3).x
    )
    with pytest.raises(TypeError, match="integer seed"):
        st.noise_gen(PhaseData, 16, 1.0, sigma1={0: 1e-12}, rng=3.0)


def test_unknown_noise_type_names_the_key(small_phase):
    with pytest.raises(ValueError, match="unknown noise type 'XX'"):
        st.clock_model({"XX": 1.0}, tau0=1.0)
    with pytest.raises(ValueError, match="unknown noise type 'XX'"):
        st.noise_fit(st.adev(small_phase)).add_component("XX", slice(0, 3))


def test_delimiter_that_is_no_string(tmp_path):
    with pytest.raises(ValueError, match="one character"):
        st.read_phase(tmp_path / "x.csv", delim=1)

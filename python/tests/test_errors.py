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
    suite = st.stability(small_phase, taus=[1, 2, 4])
    with pytest.raises(KeyError):
        suite["does_not_exist"]

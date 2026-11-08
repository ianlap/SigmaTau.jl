"""Averaging-factor grids and averaging-factor validation."""

from __future__ import annotations

import pytest

import sigmatau as st
from sigmatau import AllTaus, Octave


def test_octave_adev_grid_is_powers_of_two():
    assert st.tau_values(Octave, 8192, "adev") == [2**k for k in range(12)]


def test_alltaus_small_record():
    assert st.tau_values(AllTaus, 16, "adev") == [1, 2, 3, 4, 5, 6, 7]


def test_theo1_grid_entries_all_even():
    grid = st.tau_values(Octave, 8192, "theo1")
    assert grid
    assert all(m % 2 == 0 for m in grid)


def test_unknown_kernel_name_raises():
    with pytest.raises(ValueError):
        st.tau_values(Octave, 8192, "not_a_kernel")


@pytest.mark.parametrize("bad", [[0], [-1]], ids=["zero", "negative"])
def test_invalid_averaging_factors_raise(small_phase, bad):
    with pytest.raises(ValueError):
        st.adev(small_phase, bad)

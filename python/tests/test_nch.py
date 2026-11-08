"""N-cornered-hat separation: three-argument form and equivalent matrix form."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import PhaseData

MS = [1, 2, 4]


@pytest.fixture(scope="module")
def pairwise():
    """Three independent clocks -> the three cyclic pairwise adev results."""
    clocks = [st.noise_gen(PhaseData, 512, 1.0, sigma1={0: 1e-12}, rng=s) for s in (1, 2, 3)]
    a, b, c = (clk.x for clk in clocks)
    r_ab = st.adev(PhaseData(a - b, 1.0), MS)
    r_bc = st.adev(PhaseData(b - c, 1.0), MS)
    r_ca = st.adev(PhaseData(c - a, 1.0), MS)
    return r_ab, r_bc, r_ca


def test_three_argument_form(pairwise):
    separated = st.ncornered_hat(*pairwise)
    assert len(separated) == 3
    for clock in separated:
        assert clock.dev.size == len(MS)
        assert np.all(np.isfinite(clock.dev))


def test_matrix_form_matches_three_argument(pairwise):
    r_ab, r_bc, r_ca = pairwise
    three = st.ncornered_hat(r_ab, r_bc, r_ca)
    matrix = [[None, r_ab, None], [None, None, r_bc], [r_ca, None, None]]
    from_matrix = st.ncornered_hat(matrix)
    assert len(from_matrix) == 3
    for a, b in zip(three, from_matrix, strict=True):
        np.testing.assert_array_equal(a.dev, b.dev)

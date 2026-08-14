"""Noise identification, EDF, confidence intervals, and bias correction.

The EDF/CI helpers must reproduce a full adev(ci=True) result *bitwise* when fed
the raw pieces — the record duration is (N-1)*tau0, matching the Julia oracle.
"""

from __future__ import annotations

import numpy as np

import sigmatau as st
from sigmatau import PhaseData

MS = [1, 2, 4]


def test_identify_noise_white_fm():
    rec = st.noise_gen(PhaseData, 2048, 1.0, sigma1={0: 1e-12}, rng=7)
    noises = st.identify_noise(rec.x, MS)
    assert noises[0] == "WHFM"
    assert all(str(n) for n in noises)


def test_calculate_edf_and_ci_reproduce_adev(small_phase):
    ms = [1, 2, 4, 8]
    full = st.adev(small_phase, ms, ci=True, confidence=0.683)
    bare = st.adev(small_phase, ms, ci=False)
    noises = st.identify_noise(small_phase.x, ms)
    n = small_phase.x.size
    t = (n - 1) * small_phase.tau0

    edf = st.calculate_edf("adev", bare.dev, noises, ms, bare.tau, n, t)
    np.testing.assert_array_equal(edf, full.edf)

    lower, upper = st.confidence_intervals(bare.dev, edf, noises, n, 0.683)
    np.testing.assert_array_equal(lower, full.ci_lower)
    np.testing.assert_array_equal(upper, full.ci_upper)


def test_bias_correction_unknown_var_type_is_identity():
    b = st.bias_correction(["WHFM"], "unknown", [1.0], 512.0)
    assert np.all(b == 1.0)

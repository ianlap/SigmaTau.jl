"""Noise identification, EDF, confidence intervals and bias correction.

The EDF and interval functions must reproduce a full adev(ci=True) result
bitwise from its pieces. The record duration is (N-1)*tau0, as in Julia.
"""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import FrequencyData, PhaseData, _convert

MS = [1, 2, 4]


def test_identify_noise_white_fm():
    rec = st.noise_gen(PhaseData, 2048, 1.0, sigma1={0: 1e-12}, rng=7)
    noises = st.identify_noise(rec.x, MS)
    assert noises[0] == "WHFM"
    assert all(str(n) for n in noises)


def test_identify_noise_accepts_records():
    rec = st.noise_gen(PhaseData, 2048, 1.0, sigma1={0: 1e-12}, rng=7)
    assert list(st.identify_noise(rec, MS)) == list(st.identify_noise(rec.x, MS))
    # Frequency records are integrated to phase, not read as phase.
    y = np.diff(rec.x) / rec.tau0
    from_freq = st.identify_noise(FrequencyData(y, rec.tau0), MS)
    assert from_freq[0] == "WHFM"


def test_identify_noise_keywords_match_julia(jl):
    rec = st.noise_gen(PhaseData, 4096, 1.0, h={0: 1e-24, -1: 1.6e-26, -2: 1e-30}, rng=1)
    ms = [1, 8, 64, 256, 1024]
    direct = jl.seval("(d, m; kw...) -> string.(identify_noise(d, m; kw...))")
    jd, jm = _convert.to_jl_data(rec), _convert.to_jl_ints(ms)
    default = list(st.identify_noise(rec, ms))
    no_b1 = list(st.identify_noise(rec, ms, b1=False))
    whole = list(st.identify_noise(rec, ms, b1_length="record"))
    assert no_b1 == list(direct(jd, jm, b1=False)) != default
    assert whole == list(direct(jd, jm, b1_length=jl.Symbol("record"))) != default
    assert list(st.identify_noise(rec, st.Octave)) == list(direct(jd, jl.Octave))
    with pytest.raises(ValueError):
        st.identify_noise(rec, [0])
    with pytest.raises(TypeError):
        st.identify_noise(rec.x, st.Octave)  # bare samples carry no grid


def test_calculate_edf_and_ci_reproduce_adev(small_phase):
    ms = [1, 2, 4, 8]
    full = st.adev(small_phase, ms, ci=True, confidence=0.683)
    bare = st.adev(small_phase, ms, ci=False)
    noises = st.identify_noise(small_phase.x, ms)
    n = small_phase.x.size
    t = (n - 1) * small_phase.tau0

    edf = st.calculate_edf("adev", bare.dev, noises, ms, bare.tau, n, t)
    np.testing.assert_array_equal(edf, full.edf)

    lower, upper = st.confidence_intervals(bare.dev, edf, 0.683)
    np.testing.assert_array_equal(lower, full.ci_lower)
    np.testing.assert_array_equal(upper, full.ci_upper)


def test_bias_correction_unknown_var_type_is_identity():
    b = st.bias_correction(["WHFM"], "unknown", [1], 512)
    assert np.all(b == 1.0)


def test_bias_correction_published_leaves_pm_alone():
    assert st.bias_correction(["WHPM"], "htot", [16], 512)[0] > 1.3
    assert st.bias_correction(["WHPM"], "htot", [16], 512, published=True)[0] == 1.0

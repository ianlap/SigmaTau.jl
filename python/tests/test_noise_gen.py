"""noise_gen: seeded determinism, target sigma, and argument validation."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import PhaseData


def test_same_seed_is_deterministic():
    a = st.noise_gen(PhaseData, 512, 1.0, sigma1={0: 1e-12}, rng=42)
    b = st.noise_gen(PhaseData, 512, 1.0, sigma1={0: 1e-12}, rng=42)
    np.testing.assert_array_equal(a.x, b.x)


def test_different_seeds_differ():
    a = st.noise_gen(PhaseData, 512, 1.0, sigma1={0: 1e-12}, rng=42)
    b = st.noise_gen(PhaseData, 512, 1.0, sigma1={0: 1e-12}, rng=43)
    assert not np.array_equal(a.x, b.x)


def test_target_sigma_is_hit():
    rec = st.noise_gen(PhaseData, 512, 1.0, sigma1={0: 1e-12}, rng=7)
    dev = st.adev(rec, [1], ci=False).dev[0]
    assert dev == pytest.approx(1e-12, rel=0.3)


def test_sigma1_and_h_both_rejected():
    with pytest.raises(ValueError):
        st.noise_gen(PhaseData, 512, 1.0, sigma1={0: 1e-12}, h={0: 1e-24})


def test_numpy_generator_rejected():
    with pytest.raises(TypeError):
        st.noise_gen(PhaseData, 512, 1.0, sigma1={0: 1e-12}, rng=np.random.default_rng(1))


def test_noise_id_nnls_finds_the_mixture():
    rec = st.noise_gen(PhaseData, 8192, 1.0, sigma1={0: 1e-12, -2: 1e-14}, rng=3)
    r = st.noise_id_nnls(rec)
    assert r["deviation"] == "mhdev"
    assert r["contributions"].shape == (r["tau"].size, r["alphas"].size)
    np.testing.assert_allclose(r["fit"], r["contributions"].sum(axis=1))
    assert r["dominant"][0] == "WHFM" and r["dominant"][-1] == "RWFM"

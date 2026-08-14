"""File IO, detrend, gap filling, outliers, and result/suite round-trips."""

from __future__ import annotations

import numpy as np

import sigmatau as st
from sigmatau import PhaseData


def _write_two_column_csv(path, values):
    with open(path, "w") as f:
        for i, v in enumerate(values, start=1):
            f.write(f"{i},{v}\n")


def test_read_phase_carries_source(tmp_path):
    csv_path = tmp_path / "data.csv"
    _write_two_column_csv(csv_path, np.sin(np.arange(1, 21) / 3.0))
    pd = st.read_phase(str(csv_path), tau0=1.0, value_col=2)
    assert pd.source == str(csv_path)
    assert pd.x.size == 20


def test_detrend_changes_values_preserves_n_and_source(small_phase):
    det = st.detrend(small_phase, method="linear")
    assert not np.array_equal(det.x, small_phase.x)
    assert det.x.size == small_phase.x.size
    assert det.source == small_phase.source


def test_fillgaps_interpolates_nan(small_phase):
    gapped = small_phase.x.copy()
    gapped[100] = np.nan
    filled = st.fillgaps(PhaseData(gapped, 1.0))
    assert np.all(np.isfinite(filled.x))


def test_find_outliers_flags_planted_glitch(small_phase):
    glitched = small_phase.x.copy()
    glitched[150] += 1e-6
    idx = st.find_outliers(PhaseData(glitched, 1.0))
    assert 150 in idx


def test_remove_outliers_returns_finite(small_phase):
    glitched = small_phase.x.copy()
    glitched[150] += 1e-6
    cleaned = st.remove_outliers(PhaseData(glitched, 1.0))
    assert cleaned.x.size <= small_phase.x.size
    assert np.all(np.isfinite(cleaned.x))


def test_result_round_trip(tmp_path, small_phase):
    result = st.adev(small_phase, [1, 2, 4], ci=True)
    path = tmp_path / "result.txt"
    st.save_result(str(path), result)
    loaded = st.load_result(str(path))
    for field in ("tau", "dev", "ci_lower", "ci_upper", "edf", "neff"):
        np.testing.assert_array_equal(getattr(loaded, field), getattr(result, field))
    np.testing.assert_array_equal(loaded.noise_type, result.noise_type)


def test_suite_round_trip(tmp_path, small_phase):
    suite = st.stability(small_phase, taus=[1, 2, 4])
    path = tmp_path / "suite.txt"
    st.save_suite(str(path), suite)
    loaded = st.load_suite(str(path))
    assert loaded.keys() == suite.keys()
    assert loaded.confidence == suite.confidence
    for key in suite.keys():
        np.testing.assert_array_equal(loaded[key].dev, suite[key].dev)

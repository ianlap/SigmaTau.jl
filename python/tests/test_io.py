"""File IO, detrend, gap filling, outliers, and result/suite round-trips."""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pytest

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


FORMATS = Path(__file__).resolve().parents[2] / "test" / "fixtures" / "formats"


def test_read_data_instrument_files():
    tim = str(FORMATS / "timelab.tim")
    assert st.detect_format(tim) == "timelab"
    pd = st.read_data(tim)
    assert isinstance(pd, PhaseData) and pd.tau0 == 0.1 and pd.x.size == 64
    gps = st.read_data(str(FORMATS / "gpsdo.csv"), column="OSC offset (ppb)")
    assert isinstance(gps, st.FrequencyData) and gps.tau0 == 1.0


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
    st.save(str(path), result)
    loaded = st.load(str(path))
    for field in ("tau", "dev", "ci_lower", "ci_upper", "edf", "neff"):
        np.testing.assert_array_equal(getattr(loaded, field), getattr(result, field))
    np.testing.assert_array_equal(loaded.noise_type, result.noise_type)


def test_suite_round_trip(tmp_path, small_phase):
    suite = st.stability(small_phase, m=[1, 2, 4])
    path = tmp_path / "suite.txt"
    st.save(str(path), suite)
    loaded = st.load(str(path))
    assert loaded.keys() == suite.keys()
    assert loaded.confidence == suite.confidence
    for key in suite.keys():
        np.testing.assert_array_equal(loaded[key].dev, suite[key].dev)


def test_result_round_trip_keeps_record_metadata(tmp_path, small_phase):
    result = st.adev(small_phase, [1, 2, 4], ci=True, confidence=0.95)
    assert result.tau0 == small_phase.tau0
    assert result.n == small_phase.x.size
    assert result.confidence == 0.95
    path = tmp_path / "result.txt"
    st.save(str(path), result)
    loaded = st.load(str(path))
    assert (loaded.tau0, loaded.n, loaded.confidence) == (result.tau0, result.n, 0.95)


def test_no_ci_result_has_no_confidence(small_phase):
    assert st.adev(small_phase, [1, 2], ci=False).confidence is None


def test_data_round_trip(tmp_path, small_phase):
    path = tmp_path / "phase.txt"
    st.save(str(path), small_phase)
    loaded = st.load(str(path))
    assert isinstance(loaded, PhaseData)
    assert loaded.tau0 == small_phase.tau0
    np.testing.assert_allclose(loaded.x, small_phase.x, rtol=1e-15, atol=0)


def test_save_rejects_other_types(tmp_path):
    with pytest.raises(TypeError):
        st.save(str(tmp_path / "x.txt"), [1, 2, 3])


def test_start_time_from_timestamps(tmp_path):
    # ISO timestamps set the start, which survives conversion and a save/load.
    p = tmp_path / "ts.csv"
    rows = ("2026-11-08T12:00:00Z,0", "2026-11-08T12:00:01Z,1e-9", "2026-11-08T12:00:02Z,3e-9")
    p.write_text("time,phase\n" + "\n".join(rows) + "\n")
    rec = st.read_data(str(p))
    t0 = datetime(2026, 11, 8, 12, tzinfo=UTC)
    assert st.start_time(rec) == t0
    assert st.start_time(rec.to_frequency()) == t0
    q = tmp_path / "rec.txt"
    st.save(str(q), rec)
    assert st.start_time(st.load(str(q))) == t0
    assert st.start_time(PhaseData(np.zeros(4))) is None

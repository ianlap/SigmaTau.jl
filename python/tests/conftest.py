"""Shared session-scoped fixtures for the sigmatau test suite.

The Julia runtime boots once per session (``import sigmatau`` is slow), so every
fixture here is ``scope="session"``. Records are deterministic — no RNG — unless
a test seeds ``noise_gen`` explicitly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sigmatau import FrequencyData, PhaseData
from sigmatau._julia import jl as _jl

REPO_ROOT = Path(__file__).resolve().parents[2]
_DAT = REPO_ROOT / "test" / "fixtures" / "validation" / "stable32gen.DAT"
_CSV = REPO_ROOT / "test" / "fixtures" / "validation" / "stable32out" / "stable32_data_full.csv"


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def stable32_csv() -> Path:
    """Path to the Stable32 reference CSV (skips if the fixture is absent)."""
    if not _CSV.exists():
        pytest.skip(f"Stable32 reference CSV not found: {_CSV}")
    return _CSV


@pytest.fixture(scope="session")
def phase8192() -> PhaseData:
    """The 8192-sample Stable32 phase record (10 header lines, tau0=1.0)."""
    if not _DAT.exists():
        pytest.skip(f"Stable32 phase fixture not found: {_DAT}")
    lines = _DAT.read_text().splitlines()
    samples = np.array([float(v) for v in lines[10:]], dtype=np.float64)
    return PhaseData(samples, 1.0)


@pytest.fixture(scope="session")
def small_phase() -> PhaseData:
    """Deterministic 512-sample phase record (no RNG)."""
    x = np.cumsum(np.sin(np.pi * np.arange(1, 513) / 7.3)) / 1e9
    return PhaseData(x, 1.0)


@pytest.fixture(scope="session")
def small_freq() -> FrequencyData:
    """Deterministic 512-sample fractional-frequency record (no RNG)."""
    y = np.sin(np.arange(1, 513) / 5.0) / 1e10
    return FrequencyData(y, 1.0)


@pytest.fixture(scope="session")
def jl():
    """The embedded Julia module with SigmaTau loaded (for direct-call parity)."""
    return _jl

"""StreamingStability: the running estimate equals the batch call."""

from __future__ import annotations

import numpy as np
import pytest

import sigmatau as st
from sigmatau import PhaseData


@pytest.mark.parametrize("dev", ["adev", "mdev", "hdev", "mhdev"])
def test_snapshot_matches_batch(small_phase, dev):
    x = small_phase.x
    acc = st.StreamingStability(dev, small_phase.tau0, [1, 2, 4])
    acc.push(x[0]).extend(x[1:])
    assert len(acc) == x.size
    batch = getattr(st, dev)(PhaseData(x, small_phase.tau0), [1, 2, 4], ci=False)
    np.testing.assert_allclose(acc.snapshot().dev, batch.dev, rtol=1e-12)


def test_rejects_total_family():
    with pytest.raises(ValueError):
        st.StreamingStability("totdev", 1.0, [1])


def test_refuses_a_sample_that_is_not_finite(small_phase):
    acc = st.StreamingStability("adev", 1.0, [1, 2]).extend(small_phase.x[:10])
    with pytest.raises(ValueError, match="finite"):
        acc.push(np.nan)
    with pytest.raises(ValueError, match="finite"):
        acc.extend([0.0, np.inf, 1.0])
    assert len(acc) == 11  # the samples before the refused one stay

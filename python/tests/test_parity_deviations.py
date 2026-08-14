"""The core contract: every static deviation, computed through the Python wrapper,
is *bitwise* identical to a direct Julia call on the same converted input.

Also a single golden test against the Stable32 reference CSV, verifying the whole
conversion layer against an independent oracle.
"""

from __future__ import annotations

import csv

import numpy as np
import pytest

import sigmatau as st
from sigmatau import Octave, _convert

# All 16 static deviations exposed by the package.
STATIC_DEVIATIONS = [
    "adev",
    "mdev",
    "tdev",
    "hdev",
    "mhdev",
    "htdev",
    "totdev",
    "mtotdev",
    "ttotdev",
    "htotdev",
    "mhtotdev",
    "mtie",
    "tierms",
    "pdev",
    "theo1",
    "theoh",
]

_CI_FIELDS = ("tau", "dev", "edf", "ci_lower", "ci_upper")


def _explicit_ms(name: str) -> list[int]:
    # theo1 requires even averaging factors; tau = 0.75 * m * tau0.
    return [2, 4, 8] if name == "theo1" else [1, 2, 4, 8]


@pytest.mark.parametrize("name", STATIC_DEVIATIONS, ids=STATIC_DEVIATIONS)
@pytest.mark.parametrize("fixture", ["small_phase", "small_freq"])
@pytest.mark.parametrize("tau_kind", ["explicit", "octave"])
def test_wrapper_matches_direct_julia(name, fixture, tau_kind, request, jl):
    data = request.getfixturevalue(fixture)
    wrapper_fn = getattr(st, name)
    jl_fn = getattr(jl, name)
    jldata = _convert.to_jl_data(data)

    if tau_kind == "explicit":
        ms = _explicit_ms(name)
        wrapped = wrapper_fn(data, ms, ci=True, confidence=0.683)
        jl_taus = _convert.to_jl_ints(ms)
    else:
        wrapped = wrapper_fn(data, Octave, ci=True, confidence=0.683)
        jl_taus = _convert.resolve_taus(Octave)

    direct = _convert.result_from_jl(jl_fn(jldata, jl_taus, ci=True, confidence=0.683))

    assert wrapped.deviation_type == direct.deviation_type
    for field in _CI_FIELDS:
        # assert_array_equal is NaN-safe (NaNs in matching positions compare equal).
        np.testing.assert_array_equal(getattr(wrapped, field), getattr(direct, field))
    np.testing.assert_array_equal(wrapped.noise_type, direct.noise_type)
    np.testing.assert_array_equal(wrapped.neff, direct.neff)


def test_adev_golden_against_stable32(phase8192, stable32_csv):
    """adev on the 8192-sample record reproduces Stable32's Overlapping Allan sigmas."""
    ms = [1, 2, 4, 8]
    result = st.adev(phase8192, ms, ci=True)

    reference = {}
    with open(stable32_csv, newline="") as f:
        for row in csv.DictReader(f):
            if row["Type"] == "Overlapping Allan":
                reference[int(row["AF"])] = float(row["Sigma"])

    expected = np.array([reference[m] for m in ms])
    np.testing.assert_allclose(result.dev, expected, rtol=2e-4)

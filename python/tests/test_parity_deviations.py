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
        jl_taus = _convert.resolve_m(Octave)

    direct = _convert.result_from_jl(jl_fn(jldata, jl_taus, ci=True, confidence=0.683))

    assert wrapped.deviation_type == direct.deviation_type
    for field in _CI_FIELDS:
        # assert_array_equal is NaN-safe (NaNs in matching positions compare equal).
        np.testing.assert_array_equal(getattr(wrapped, field), getattr(direct, field))
    np.testing.assert_array_equal(wrapped.noise_type, direct.noise_type)
    np.testing.assert_array_equal(wrapped.neff, direct.neff)


def test_htotdev_remove_drift_sets_the_noise_type_and_the_bias(phase8192):
    """A frequency drift HTOTDEV ignores does not reach its noise type by
    default; with remove_drift=False the type at τ = 1024 s is the drift's,
    and the bias factor follows."""
    drifting = st.PhaseData(phase8192.x + 1e-5 * np.arange(len(phase8192.x)) ** 2, 1.0)
    default = st.htotdev(drifting, [1024])
    assert default.noise_type[0] == "FLFM"
    assert default.dev[0] == pytest.approx(st.htotdev(phase8192, [1024]).dev[0], rel=1e-8)
    as_given = st.htotdev(drifting, [1024], remove_drift=False)
    assert as_given.noise_type[0] == "RWFM"
    assert as_given.dev[0] == pytest.approx(2.74011e-03, rel=1e-5)
    assert st.deviation("htotdev", drifting, [1024], remove_drift=False).dev[0] == as_given.dev[0]


def test_correct_bias_takes_a_reference_name(phase8192):
    """At τ = 16 s the record reads white PM, where Stable32 and the published
    tables apply no factor to HTOTDEV and the default lowers it."""
    raw = st.htotdev(phase8192, [16], correct_bias=False).dev[0]
    assert raw == pytest.approx(7.7624e-02, rel=1e-4)
    assert st.htotdev(phase8192, [16], correct_bias="stable32").dev[0] == raw
    assert st.deviation("htotdev", phase8192, [16], correct_bias="published").dev[0] == raw
    assert st.htotdev(phase8192, [16]).dev[0] < 0.9 * raw


def test_noise_states_the_type_behind_the_edf(small_phase):
    """A stated noise type stands in for the identified one: one name, one per
    averaging factor, or a fitted model (flicker PM is this record's strongest
    term at every τ)."""
    ms = [1, 4, 16]
    auto = st.adev(small_phase, ms)
    stated = st.adev(small_phase, ms, noise=list(auto.noise_type))
    np.testing.assert_array_equal(stated.edf, auto.edf)
    flfm = st.adev(small_phase, ms, noise="FLFM")
    assert list(flfm.noise_type) == ["FLFM"] * 3
    assert flfm.edf[0] != auto.edf[0]
    suite = st.stability(small_phase, m=ms, noise="FLFM")
    np.testing.assert_array_equal(suite["adev"].edf, flfm.edf)
    fitted = st.deviation("hdev", small_phase, ms, noise=st.powerlaw_fit(small_phase))
    assert list(fitted.noise_type) == ["FLPM"] * 3
    with pytest.raises(ValueError, match="noise"):
        st.adev(small_phase, ms, noise="PINK")


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


@pytest.mark.parametrize("name", ["adev", "hdev"])
def test_non_overlapped(small_phase, name):
    fn, m = getattr(st, name), [1, 4, 16]
    r, o = fn(small_phase, m, overlapping=False), fn(small_phase, m)
    assert r.dev[0] == o.dev[0] and (r.neff[1:] < o.neff[1:]).all()

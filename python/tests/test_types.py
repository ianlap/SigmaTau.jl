"""PhaseData / FrequencyData validation, immutability, defaults, and reprs."""

from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from sigmatau import FrequencyData, PhaseData


@pytest.mark.parametrize("cls", [PhaseData, FrequencyData], ids=["phase", "freq"])
@pytest.mark.parametrize("bad_tau0", [0.0, -1.0], ids=["zero", "negative"])
def test_tau0_must_be_positive(cls, bad_tau0):
    with pytest.raises(ValueError):
        cls(np.arange(10.0), bad_tau0)


@pytest.mark.parametrize("cls", [PhaseData, FrequencyData], ids=["phase", "freq"])
def test_too_few_samples(cls):
    with pytest.raises(ValueError):
        cls(np.array([1.0]), 1.0)


@pytest.mark.parametrize("cls", [PhaseData, FrequencyData], ids=["phase", "freq"])
def test_two_dimensional_input_rejected(cls):
    with pytest.raises(ValueError):
        cls(np.ones((4, 2)), 1.0)


@pytest.mark.parametrize("cls", [PhaseData, FrequencyData], ids=["phase", "freq"])
def test_frozen(cls):
    d = cls(np.arange(10.0), 1.0)
    with pytest.raises(FrozenInstanceError):
        d.tau0 = 2.0


@pytest.mark.parametrize("cls", [PhaseData, FrequencyData], ids=["phase", "freq"])
def test_source_defaults_to_user(cls):
    assert cls(np.arange(10.0), 1.0).source == "user"


@pytest.mark.parametrize("cls", [PhaseData, FrequencyData], ids=["phase", "freq"])
def test_repr_contains_n_and_tau0(cls):
    r = repr(cls(np.arange(10.0), 2.5))
    assert "N=10" in r
    assert "2.5" in r

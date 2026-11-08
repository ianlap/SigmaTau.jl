"""Public surface resolves and the version/introspection helpers work."""

from __future__ import annotations

import pytest

import sigmatau as st
from sigmatau import julia_version, sigmatau_jl_version


@pytest.mark.parametrize("name", st.__all__, ids=st.__all__)
def test_public_name_resolves(name):
    assert getattr(st, name, None) is not None


def test_all_is_complete():
    # 56 public names, per the packaged contract.
    assert len(st.__all__) == 56


def test_version():
    assert st.__version__ == "0.5.0"


def test_julia_version_nonempty():
    assert isinstance(julia_version(), str)
    assert julia_version()


def test_sigmatau_jl_version_nonempty():
    assert isinstance(sigmatau_jl_version(), str)
    assert sigmatau_jl_version()

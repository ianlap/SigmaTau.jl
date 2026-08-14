"""The Python↔Julia boundary: eager, copying conversions in both directions.

Down-conversion (Julia → Python) copies every array into a fresh NumPy array so
results never hold views into Julia-owned memory; up-conversion (Python →
Julia) goes through ``juliacall.convert`` against cached Julia type objects
because SigmaTau's kernels dispatch on concrete ``Vector`` types (a raw
``PyArray`` is an ``AbstractVector``, not a ``Vector``).
"""

from __future__ import annotations

from typing import Any

import juliacall
import numpy as np

from ._julia import (
    DictIntF64,
    VecF64,
    VecInt,
    VecResult,
    VecStr,
    _make_result,
    _to_syms,
    jl,
)
from .grids import _JL_TAUMODE, TauMode
from .types import (
    DynamicStabilityResult,
    FrequencyData,
    PhaseData,
    SpectralResult,
    StabilityResult,
    StabilitySuite,
)

# ── Python → Julia ─────────────────────────────────────────────────────


def to_jl_f64(a: object) -> Any:
    """Sample/float vector → ``Vector{Float64}`` (copy)."""
    return juliacall.convert(VecF64, np.ascontiguousarray(a, dtype=np.float64))


def to_jl_ints(ms: object) -> Any:
    """Integer sequence → ``Vector{Int}``."""
    return juliacall.convert(VecInt, [int(m) for m in ms])


def to_jl_syms(ss: object) -> Any:
    """String sequence → ``Vector{Symbol}``."""
    return _to_syms(juliacall.convert(VecStr, [str(s) for s in ss]))


def to_jl_dict_int_f64(d: dict) -> Any:
    """``{int: float}`` mapping → ``Dict{Int,Float64}`` (noise_gen exponents)."""
    return juliacall.convert(DictIntF64, {int(k): float(v) for k, v in d.items()})


def resolve_taus(taus: TauMode | object) -> Any:
    """The ``taus`` argument as Julia sees it: a TauMode enum instance or ``Vector{Int}``.

    TauMode passes through as the enum so Julia's own dispatch resolves the
    grid (dadev/dhdev clamp it to the analysis window, not the record).
    Averaging factors are validated here — the kernels assume m ≥ 1.
    """
    if isinstance(taus, TauMode):
        return _JL_TAUMODE[taus]
    ms = [int(m) for m in taus]
    if any(m < 1 for m in ms):
        raise ValueError(f"averaging factors must be ≥ 1, got {min(ms)}")
    return to_jl_ints(ms)


def to_jl_data(data: PhaseData | FrequencyData) -> Any:
    """Python data record → Julia ``PhaseData``/``FrequencyData``."""
    if isinstance(data, PhaseData):
        return jl.PhaseData(to_jl_f64(data.x), data.tau0, source=data.source)
    if isinstance(data, FrequencyData):
        return jl.FrequencyData(to_jl_f64(data.y), data.tau0, source=data.source)
    raise TypeError(f"expected PhaseData or FrequencyData, got {type(data).__name__}")


def result_to_jl(r: StabilityResult) -> Any:
    """Python ``StabilityResult`` → Julia (for ``save_result``/``save_suite``/``nch``)."""
    return _make_result(
        r.deviation_type,
        to_jl_f64(r.tau),
        to_jl_f64(r.dev),
        to_jl_syms(r.noise_type),
        to_jl_f64(r.ci_lower),
        to_jl_f64(r.ci_upper),
        to_jl_f64(r.edf),
        to_jl_ints(r.neff),
    )


def suite_to_jl(s: StabilitySuite) -> Any:
    """Python ``StabilitySuite`` → Julia (for ``save_suite``)."""
    results = juliacall.convert(VecResult, [result_to_jl(r) for r in s.results])
    confidence = float("nan") if s.confidence is None else float(s.confidence)
    return jl.StabilitySuite(
        results,
        jl.Symbol(s.data_kind),
        float(s.tau0),
        int(s.n),
        confidence,
        jl.Symbol(s.tau_mode),
    )


# ── Julia → Python ─────────────────────────────────────────────────────


def np_f64(v: Any) -> np.ndarray:
    return np.array(v, dtype=np.float64)


def np_int(v: Any) -> np.ndarray:
    return np.array(v, dtype=np.int64)


def np_syms(v: Any) -> np.ndarray:
    """``Vector{Symbol}`` → object array of plain strings (e.g. ``"WHFM"``)."""
    return np.array([str(s) for s in v], dtype=object)


def data_from_jl(d: Any) -> PhaseData | FrequencyData:
    """Julia data record → Python (read_phase/detrend/fillgaps round-trips)."""
    if jl.isa(d, jl.PhaseData):
        return PhaseData(np_f64(d.x), float(d.tau0), source=str(d.source))
    return FrequencyData(np_f64(d.y), float(d.tau0), source=str(d.source))


def result_from_jl(r: Any) -> StabilityResult:
    return StabilityResult(
        deviation_type=str(r.deviation_type),
        tau=np_f64(r.tau),
        dev=np_f64(r.dev),
        noise_type=np_syms(r.noise_type),
        ci_lower=np_f64(jl.ci_lower(r)),
        ci_upper=np_f64(jl.ci_upper(r)),
        edf=np_f64(r.edf),
        neff=np_int(r.neff),
    )


def suite_from_jl(s: Any) -> StabilitySuite:
    confidence = float(s.confidence)
    return StabilitySuite(
        results=tuple(result_from_jl(r) for r in s.results),
        data_kind=str(s.data_kind),
        tau0=float(s.tau0),
        n=int(s.n),
        confidence=None if np.isnan(confidence) else confidence,
        tau_mode=str(s.tau_mode),
    )


def spectral_from_jl(r: Any) -> SpectralResult:
    return SpectralResult(
        spectral_type=str(r.spectral_type),
        freq=np_f64(r.freq),
        psd=np_f64(r.psd),
        units=str(r.units),
        nperseg=int(r.nperseg),
        noverlap=int(r.noverlap),
        window=str(r.window),
    )


def dynamic_from_jl(r: Any) -> DynamicStabilityResult:
    return DynamicStabilityResult(
        deviation_type=str(r.deviation_type),
        t=np_f64(r.t),
        tau=np_f64(r.tau),
        dev=np.array(r.dev, dtype=np.float64),
        window=int(r.window),
        tau0=float(r.tau0),
    )

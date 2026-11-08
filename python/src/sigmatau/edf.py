"""Equivalent degrees of freedom, confidence intervals, and bias correction.

Thin wrappers over ``edf.jl``. These operate on the raw pieces of a
``StabilityResult`` (deviations, noise types, averaging factors, τ grid), so
they can be composed independently of the deviation call that produced them.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ._convert import np_f64, to_jl_f64, to_jl_ints, to_jl_syms
from ._julia import _map_errors, jl


@_map_errors
def calculate_edf(
    method: str,
    devs: Sequence[float],
    noises: Sequence[str],
    m_values: Sequence[int],
    taus: Sequence[float],
    n: int,
    t: float,
) -> np.ndarray:
    """Equivalent degrees of freedom per τ for the given deviation ``method``.

    ``method`` is a deviation name (e.g. ``"adev"``, ``"mhtotdev"``); ``noises``
    are the per-τ noise-type symbols (from :func:`identify_noise`); ``n`` is the
    record length and ``t`` the record duration ``(N-1)·τ₀`` in seconds.
    """
    return np_f64(
        jl.calculate_edf(
            jl.Symbol(str(method)),
            to_jl_f64(devs),
            to_jl_syms(noises),
            to_jl_ints(m_values),
            to_jl_f64(taus),
            int(n),
            float(t),
        )
    )


@_map_errors
def confidence_intervals(
    devs: Sequence[float],
    edfs: Sequence[float],
    noises: Sequence[str],
    n: int,
    confidence: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Lower/upper confidence limits — χ² when edf ≥ 1, Gaussian fallback otherwise.

    Returns a ``(lower, upper)`` tuple of NumPy arrays.
    """
    result = jl.confidence_intervals(
        to_jl_f64(devs),
        to_jl_f64(edfs),
        to_jl_syms(noises),
        int(n),
        float(confidence),
    )
    return np_f64(jl.getindex(result, 1)), np_f64(jl.getindex(result, 2))


@_map_errors
def bias_correction(
    noises: Sequence[str],
    var_type: str,
    taus: Sequence[float],
    t: float,
) -> np.ndarray:
    """Variance-scale bias factor B(α) = E[estimator²]/true_variance (apply as σ/√B).

    ``var_type`` is one of ``"totvar"``, ``"mtot"``, ``"htot"``, ``"mhtot"``;
    anything else returns all ones.
    """
    return np_f64(
        jl.bias_correction(
            to_jl_syms(noises),
            jl.Symbol(str(var_type)),
            to_jl_f64(taus),
            float(t),
        )
    )


__all__ = ["calculate_edf", "confidence_intervals", "bias_correction"]

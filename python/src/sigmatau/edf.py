"""Equivalent degrees of freedom, confidence intervals and bias factors.

The estimators use these for the interval they report. The functions take the
pieces of a ``StabilityResult`` (deviations, noise types, averaging factors and
averaging times), so that each can be used on its own.
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
    """Equivalent degrees of freedom (EDF) of a variance estimate at each averaging
    factor: the ν of the χ² distribution behind :func:`confidence_intervals`.

    ``method`` names the estimator: ``"adev"``, ``"mdev"``, ``"hdev"``,
    ``"mhdev"``, ``"totdev"``, ``"mtotdev"``, ``"htotdev"``, ``"mhtotdev"``,
    ``"theo1"``, ``"theobr"`` (Thêo1 with its bias removed, as :func:`theo1`
    returns it by default) or ``"pdev"``. Any other gives ``nan``. ``noises``
    holds the noise type at each point, as :func:`identify_noise` returns it,
    ``m_values`` and ``taus`` the averaging factor and averaging time of each,
    ``n`` the number of phase samples and ``t`` the record duration (N − 1)·τ₀
    in seconds. Of ``devs`` only the number is used.

    The models are those of Greenhall and Riley (2003) for the Allan and
    Hadamard families, of Greenhall et al. (1999), Howe and Peppler (2001) and
    Howe et al. (2005) for the total estimators, of SP1065 §5.2.15 for Thêo1
    and of Vernotte et al. (2021) for the parabolic deviation. The EDF of
    ``"mhtotdev"``, and that of ``"totdev"`` and ``"htotdev"`` on PM noise,
    are models of this package.
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
    confidence: float,
) -> tuple[np.ndarray, np.ndarray]:
    """χ² confidence limits of each deviation from its equivalent degrees of
    freedom (SP1065 §5.3.2).

    ``confidence`` is the probability that the interval covers, such as 0.683.
    Returns the pair ``(lower, upper)``; a limit is ``nan`` where the deviation
    is not finite or the EDF is not positive.
    """
    result = jl.confidence_intervals(to_jl_f64(devs), to_jl_f64(edfs), float(confidence))
    return np_f64(jl.getindex(result, 1)), np_f64(jl.getindex(result, 2))


@_map_errors
def bias_correction(
    noises: Sequence[str],
    var_type: str,
    m_values: Sequence[int],
    n: int,
    *,
    published: bool = False,
) -> np.ndarray:
    """Bias factor B = E[total variance]/E[parent variance] of a total estimator at
    each averaging factor, for a record of ``n`` phase samples.

    The estimators divide a total deviation by √B. ``var_type`` is ``"totvar"``
    (total against Allan variance), ``"mtot"`` (modified total against modified
    Allan), ``"htot"`` (Hadamard total against Hadamard) or ``"mhtot"``; any
    other gives B = 1. ``noises`` gives the noise type at each entry of
    ``m_values``.

    The published values are those of Greenhall et al. (1999) for ``"totvar"``
    on FM noise, and of Howe and Peppler (2001) for ``"mtot"`` and Howe et al.
    (2005) for ``"htot"`` on FM noise, both from τ = 16τ₀. The values of
    ``"totvar"`` and ``"htot"`` on PM noise, those of ``"mtot"`` and ``"htot"``
    below m = 16 and all of ``"mhtot"`` are the expectation of each statistic,
    computed in this package. ``published=True`` leaves them out: it gives
    B = 1 on PM for ``"totvar"`` and ``"htot"``, and the published value of
    ``"mtot"`` and ``"htot"`` at every m > 1. ``"mhtot"`` is unchanged.
    """
    return np_f64(
        jl.bias_correction(
            to_jl_syms(noises),
            jl.Symbol(str(var_type)),
            to_jl_ints(m_values),
            int(n),
            published=bool(published),
        )
    )


__all__ = ["calculate_edf", "confidence_intervals", "bias_correction"]

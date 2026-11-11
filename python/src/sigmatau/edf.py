"""Equivalent degrees of freedom, confidence intervals, bias factors and the
bias functions for samples with dead time.

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


def _mu(mu: float | str):
    """``mu`` for Julia: the exponent as a float, or a noise type as a Symbol."""
    return jl.Symbol(mu) if isinstance(mu, str) else float(mu)


def _count(n: float) -> int | float:
    """A sample count for Julia: an int, or ``inf`` kept as a float."""
    return float(n) if n == float("inf") else int(n)


@_map_errors
def bias_b1(n: float, r: float, mu: float | str) -> float:
    """Barnes's bias function B1(N, r, μ): the expected ``n``-sample variance of
    frequency averages over their expected two-sample variance (NBS TN 375,
    Eq. (8); NIST TN 1318, Eq. (11)).

    The averages are over τ and start every T = ``r``·τ, so that ``r`` > 1 is
    dead time, on noise whose Allan variance goes as τ^μ. ``n`` is an integer
    ≥ 2, or ``math.inf`` for the limit of a long record (infinite for μ ≥ 0).
    ``mu`` is the exponent, from −2 to 2, or a noise type: ``"WHPM"`` and
    ``"FLPM"`` (μ = −2), ``"WHFM"`` (−1), ``"FLFM"`` (0) or ``"RWFM"`` (1).
    At μ = −2 the function is that of white PM.
    """
    return float(jl.bias_b1(_count(n), float(r), _mu(mu)))


@_map_errors
def bias_b2(r: float, mu: float | str) -> float:
    """Barnes's bias function B2(r, μ): the expected two-sample variance of
    frequency averages with dead time over the Allan variance at the same
    averaging time (NBS TN 375, Eq. (10); NIST TN 1318, Eq. (12)).

    ``r`` ≥ 0 is the ratio T/τ of the measurement period to the averaging
    time, and ``mu`` the exponent or noise type, as in :func:`bias_b1`.
    """
    return float(jl.bias_b2(float(r), _mu(mu)))


@_map_errors
def bias_b3(n: int, m: int, r: float, mu: float | str) -> float:
    """Barnes and Allan's bias function B3(N, M, r, μ): the expected ``n``-sample
    variance of averages of ``m`` measurements, with the dead time spread among
    them, over that of single measurements ``m`` times as long with the dead
    time at the end (NIST TN 1318, Eqs. (7) and (13)).

    B2(r, μ)·B3(2, M, r, μ) is the ratio of the two-sample variance of
    M-averages of dead-time data to the Allan variance at τ = M·τ₀, which
    :func:`adev` divides by under ``dead_time``. The reference gives the
    expression for ``n`` = 2; for more samples the function evaluates the
    definition by the same steps. ``r`` and ``mu`` are as in :func:`bias_b1`.
    """
    return float(jl.bias_b3(int(n), int(m), float(r), _mu(mu)))


__all__ = [
    "calculate_edf",
    "confidence_intervals",
    "bias_correction",
    "bias_b1",
    "bias_b2",
    "bias_b3",
]

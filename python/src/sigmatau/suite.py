"""Several deviations in one call, and the separation of clocks measured against
each other.

``stability`` runs a set of estimators on one record and returns a
``StabilitySuite``. ``ncornered_hat`` gives the deviation of each clock from the
deviations measured between every pair of them, and ``groslambert`` gives its
Allan deviation from the difference records themselves.
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import (
    resolve_m,
    result_from_jl,
    result_to_jl,
    suite_from_jl,
    to_jl_data,
    to_jl_syms,
)
from ._julia import _map_errors, jl
from .deviations import DEFAULT_CONFIDENCE, _bias, _noise
from .grids import Octave, TauMode
from .types import FrequencyData, PhaseData, PowerLawFit, StabilityResult, StabilitySuite

DEFAULT_DEVIATIONS = tuple(str(s) for s in jl.DEFAULT_DEVIATIONS)

# A Julia matrix for the matrix forms, with the cells that are not given left
# unassigned: Julia reads either triangle, so only the measured pairs are set.
_new_result_matrix = jl.seval("n -> Matrix{StabilityResult}(undef, n, n)")
_new_data_matrix = jl.seval("n -> Matrix{AbstractTimingData}(undef, n, n)")
_assign_result = jl.seval("(M, r, i, j) -> (M[i, j] = r; nothing)")


def _fill_matrix(new, convert, matrix, name):
    """A part-filled square Julia matrix from a list of lists with ``None`` holes."""
    n = len(matrix)
    jl_matrix = new(n)
    for i, row in enumerate(matrix):
        if len(row) != n:
            raise ValueError(
                f"{name}: pairwise matrix must be square, got a {n}-row matrix "
                f"with a length-{len(row)} row"
            )
        for j, entry in enumerate(row):
            if entry is not None:
                _assign_result(jl_matrix, convert(entry), i + 1, j + 1)
    return jl_matrix


@_map_errors
def stability(
    data: PhaseData | FrequencyData,
    *,
    devs: str | Sequence[str] = DEFAULT_DEVIATIONS,
    m: TauMode | Sequence[int] = Octave,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
    correct_bias: bool | str = True,
    remove_drift: bool = True,
    noise: str | Sequence[str] | PowerLawFit = "auto",
) -> StabilitySuite:
    """Several deviations of the record ``data`` in a single call, as a
    ``StabilitySuite`` indexed by name.

    ``devs`` names the estimators to run, from :func:`estimators`; one name
    runs one estimator. The default
    is ``DEFAULT_DEVIATIONS``, ``("adev", "mdev", "hdev", "mhdev")``. ``m`` is a
    ``TauMode`` grid, cut for each estimator at the largest averaging factor it
    supports, or a sequence of averaging factors that every estimator takes as
    it stands. ``theo1`` reads such a sequence as even Thêo1 factors and
    ``theoh`` as factors of τ = m·τ₀, so prefer a ``TauMode`` when the suite
    mixes them with other estimators.

    ``ci``, ``confidence`` and ``noise`` are those of the estimators.
    ``correct_bias`` reaches the total estimators and ``theo1``, and
    ``remove_drift`` the Hadamard family.
    """
    return suite_from_jl(
        jl.stability(
            to_jl_data(data),
            devs=to_jl_syms([devs] if isinstance(devs, str) else devs),
            m=resolve_m(m),
            ci=bool(ci),
            confidence=float(confidence),
            correct_bias=_bias(correct_bias),
            remove_drift=bool(remove_drift),
            noise=_noise(noise),
        )
    )


@_map_errors
def ncornered_hat(
    *results: StabilityResult | Sequence[Sequence[StabilityResult | None]],
) -> list[StabilityResult]:
    """Deviation of each of N ≥ 3 clocks from the deviations measured between every
    pair of them: the three-cornered hat of Gray and Allan (1974) and its
    extension to N clocks (SP1065 §10.14).

    Two call forms for ``*results``:

    - ``ncornered_hat(r_ab, r_bc, r_ca)``: for three clocks, the results for
      A−B, B−C and C−A. The clocks are returned in the order A, B, C.
    - ``ncornered_hat(matrix)``: a square list of lists whose entry ``[i][j]``
      is the result of an estimator on the difference record of clocks i and j,
      counted from 0, with ``None`` elsewhere. One entry per pair is enough, in
      either triangle.

    The sign of each difference does not matter. All results must come from
    one estimator on one τ grid, and the estimator must be one whose variance
    adds for independent clocks, such as ``adev``, ``mdev`` or ``hdev``.

    The estimate can be negative for a clock much more stable than the others,
    for clocks that are correlated and where few terms are averaged; ``dev``
    is ``nan`` there. The EDF and the interval are filled when every input has
    an EDF (Vernotte and Lantz 2019, Appendix B). The interval is too narrow
    for a clock that is no longer resolved under the others.
    """
    if len(results) == 3:
        separated = jl.ncornered_hat(
            result_to_jl(results[0]),
            result_to_jl(results[1]),
            result_to_jl(results[2]),
        )
        return [result_from_jl(r) for r in separated]

    if len(results) == 1:
        jl_matrix = _fill_matrix(_new_result_matrix, result_to_jl, results[0], "ncornered_hat")
        return [result_from_jl(r) for r in jl.ncornered_hat(jl_matrix)]

    raise TypeError(
        "ncornered_hat: pass three pairwise results (r_ab, r_bc, r_ca) "
        "or one square matrix of results"
    )


@_map_errors
def groslambert(
    *records: PhaseData | FrequencyData | Sequence[Sequence[PhaseData | FrequencyData | None]],
    m: TauMode | Sequence[int] = Octave,
    confidence: float = DEFAULT_CONFIDENCE,
    interval: str = "chi2",
    weighted: bool = False,
) -> list[StabilityResult]:
    """Allan deviation of each of N ≥ 3 clocks from the difference records between
    them, by the Groslambert covariance: the cross Allan variance of two records
    that share the clock (Vernotte, Calosso and Rubiola 2016, after Groslambert
    et al. 1981).

    Two call forms for ``*records``, which are ``PhaseData`` or
    ``FrequencyData`` records with one sample interval and length and no gaps:

    - ``groslambert(x_ab, x_bc, x_ca)``: for three clocks, the records A−B, B−C
      and C−A. The clocks are returned in the order A, B, C.
    - ``groslambert(matrix)``: a square list of lists whose entry ``[i][j]`` is
      the record x_i − x_j, counted from 0, with ``None`` elsewhere. One entry
      per pair is enough, in either triangle: the entry across the diagonal
      stands in, with its sign reversed.

    Noise that each record carries on its own, such as that of the counter
    that measured it, biases the hat of :func:`ncornered_hat` and averages out
    of the covariance. Like the hat, the estimate can be negative for a clock
    much more stable than the others; ``dev`` is ``nan`` there.

    ``m`` is a ``TauMode`` grid or a sequence of averaging factors. Each result
    is an ``"adev"`` ``StabilityResult`` without noise types.

    ``interval="chi2"`` gives the χ² interval at ``confidence`` on the EDF of
    the estimate, from the variance model of Vernotte and Lantz (2019,
    Appendix B). It holds for a clock that dominates its records and is too
    narrow for a clock that is no longer resolved under the others.
    ``interval="klts"``, for three clocks, gives the interval of Lantz et al.
    (2019), which is given also where the estimate is negative.

    With ``weighted=True`` and four clocks or more, each pair of the other
    clocks counts by the inverse of the product of its two records' Allan
    variances, so that a quiet clock is read against the other quiet clocks.
    """
    kwargs = {
        "confidence": float(confidence),
        "interval": jl.Symbol(interval),
        "weighted": bool(weighted),
    }
    if len(records) == 3:
        separated = jl.groslambert(
            to_jl_data(records[0]),
            to_jl_data(records[1]),
            to_jl_data(records[2]),
            resolve_m(m),
            **kwargs,
        )
        return [result_from_jl(r) for r in separated]

    if len(records) == 1:
        jl_matrix = _fill_matrix(_new_data_matrix, to_jl_data, records[0], "groslambert")
        return [result_from_jl(r) for r in jl.groslambert(jl_matrix, resolve_m(m), **kwargs)]

    raise TypeError(
        "groslambert: pass three difference records (x_ab, x_bc, x_ca) "
        "or one square matrix of records"
    )


__all__ = ["stability", "DEFAULT_DEVIATIONS", "ncornered_hat", "groslambert"]

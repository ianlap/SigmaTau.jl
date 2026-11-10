"""Compute-all entry point and the clock-ensemble separations — mirrors ``suite.jl``
and ``ensemble.jl``.

``stability`` runs a set of deviations on one record in a single Julia call and
collects them into a ``StabilitySuite``; ``ncornered_hat`` separates per-clock
deviations from pairwise results (the classic three-cornered hat at N=3) and
``groslambert`` does the same from the pairwise difference records without
subtracting variances.
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
from .deviations import DEFAULT_CONFIDENCE
from .grids import Octave, TauMode
from .types import FrequencyData, PhaseData, StabilityResult, StabilitySuite

DEFAULT_DEVIATIONS = tuple(str(s) for s in jl.DEFAULT_DEVIATIONS)

# Build a possibly part-filled Julia StabilityResult matrix for the nch matrix
# form: unassigned cells stay ``undef`` (Julia reads the upper or lower
# triangle, so only the measured pairs need filling).
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
    devs: Sequence[str] = DEFAULT_DEVIATIONS,
    m: TauMode | Sequence[int] = Octave,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
    correct_bias: bool = True,
    remove_drift: bool = True,
) -> StabilitySuite:
    """Compute several deviations on one record and collect them into a suite.

    ``devs`` is a sequence of deviation names (default ``DEFAULT_DEVIATIONS``,
    the four σ_y quantities ``("adev", "mdev", "hdev", "mhdev")``). ``m`` is a
    ``TauMode`` grid selector or an explicit sequence of averaging factors;
    ``ci``/``confidence``, ``correct_bias`` (total family) and ``remove_drift``
    (Hadamard family) flow through to each deviation.
    """
    return suite_from_jl(
        jl.stability(
            to_jl_data(data),
            devs=to_jl_syms(devs),
            m=resolve_m(m),
            ci=bool(ci),
            confidence=float(confidence),
            correct_bias=bool(correct_bias),
            remove_drift=bool(remove_drift),
        )
    )


@_map_errors
def ncornered_hat(
    *results: StabilityResult | Sequence[Sequence[StabilityResult | None]],
) -> list[StabilityResult]:
    """N-cornered-hat noise separation (Gray–Allan three-cornered hat at N=3).

    Two call forms, mirroring the Julia package:

    - ``ncornered_hat(r_ab, r_bc, r_ca)`` — the three-cornered convenience form: the three
      pairwise deviations in cyclic order A−B, B−C, C−A. Returns ``[A, B, C]``.
    - ``ncornered_hat(matrix)`` — a square list-of-lists ``matrix[i][j]`` of the pairwise
      deviation between clocks ``i`` and ``j`` (0-based), with ``None`` where a
      pair is not filled. Fill either triangle (or both); the orientation of
      each difference record does not matter.

    Every pairwise result must share one ``deviation_type`` and one τ grid.
    Returns one ``StabilityResult`` per clock, in clock order.
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
) -> list[StabilityResult]:
    """Groslambert covariance: each clock's Allan deviation from the pairwise records.

    Two call forms, as for ``ncornered_hat`` but on the difference records
    themselves (``PhaseData`` or ``FrequencyData``, same τ₀ and length):

    - ``groslambert(x_ab, x_bc, x_ca)`` — the records A−B, B−C, C−A. Returns ``[A, B, C]``.
    - ``groslambert(matrix)`` — ``matrix[i][j]`` is the record x_i − x_j (0-based),
      ``None`` where a pair is missing; the transposed entry, negated, stands in.

    Each result is an ``"adev"`` ``StabilityResult`` with χ² intervals at
    ``confidence`` from the covariance model of Vernotte & Lantz (2019).
    """
    kwargs = {"confidence": float(confidence)}
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

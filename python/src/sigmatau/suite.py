"""Compute-all entry point and the N-cornered-hat separation — mirrors ``suite.jl``.

``stability`` runs a set of deviations on one record in a single Julia call and
collects them into a ``StabilitySuite``; ``nch`` separates per-clock deviations
from a set of pairwise measurements (the classic three-cornered hat at N=3).
"""

from __future__ import annotations

from collections.abc import Sequence

from ._convert import (
    resolve_taus,
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
# form: unassigned cells stay ``undef`` (the oracle reads the upper or lower
# triangle, so only the measured pairs need filling).
_new_result_matrix = jl.seval("n -> Matrix{StabilityResult}(undef, n, n)")
_assign_result = jl.seval("(M, r, i, j) -> (M[i, j] = r; nothing)")


@_map_errors
def stability(
    data: PhaseData | FrequencyData,
    *,
    devs: Sequence[str] = DEFAULT_DEVIATIONS,
    taus: TauMode | Sequence[int] = Octave,
    ci: bool = True,
    confidence: float = DEFAULT_CONFIDENCE,
    correct_bias: bool = True,
) -> StabilitySuite:
    """Compute several deviations on one record and collect them into a suite.

    ``devs`` is a sequence of deviation names (default ``DEFAULT_DEVIATIONS``,
    the four σ_y quantities ``("adev", "mdev", "hdev", "mhdev")``). ``taus`` is a
    ``TauMode`` grid selector or an explicit sequence of averaging factors;
    ``ci``/``confidence`` and (for the total family) ``correct_bias`` flow
    through to each deviation.
    """
    return suite_from_jl(
        jl.stability(
            to_jl_data(data),
            devs=to_jl_syms(devs),
            taus=resolve_taus(taus),
            ci=bool(ci),
            confidence=float(confidence),
            correct_bias=bool(correct_bias),
        )
    )


@_map_errors
def nch(
    *results: StabilityResult | Sequence[Sequence[StabilityResult | None]],
) -> list[StabilityResult]:
    """N-cornered-hat noise separation (Gray–Allan three-cornered hat at N=3).

    Two call forms, mirroring the Julia oracle:

    - ``nch(r_ab, r_bc, r_ca)`` — the three-cornered convenience form: the three
      pairwise deviations in cyclic order A−B, B−C, C−A. Returns ``[A, B, C]``.
    - ``nch(matrix)`` — a square list-of-lists ``matrix[i][j]`` of the pairwise
      deviation between clocks ``i`` and ``j`` (0-based), with ``None`` where a
      pair is not filled. Fill either triangle (or both); the orientation of
      each difference record does not matter.

    Every pairwise result must share one ``deviation_type`` and one τ grid.
    Returns one ``StabilityResult`` per clock, in clock order.
    """
    if len(results) == 3:
        separated = jl.nch(
            result_to_jl(results[0]),
            result_to_jl(results[1]),
            result_to_jl(results[2]),
        )
        return [result_from_jl(r) for r in separated]

    if len(results) == 1:
        matrix = results[0]
        n = len(matrix)
        jl_matrix = _new_result_matrix(n)
        for i, row in enumerate(matrix):
            if len(row) != n:
                raise ValueError(
                    f"nch: pairwise matrix must be square, got a {n}-row matrix "
                    f"with a length-{len(row)} row"
                )
            for j, entry in enumerate(row):
                if entry is not None:
                    _assign_result(jl_matrix, result_to_jl(entry), i + 1, j + 1)
        return [result_from_jl(r) for r in jl.nch(jl_matrix)]

    raise TypeError(
        "nch: pass three pairwise results (r_ab, r_bc, r_ca) or one square matrix of results"
    )


__all__ = ["stability", "DEFAULT_DEVIATIONS", "nch"]

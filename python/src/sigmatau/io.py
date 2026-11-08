"""File IO, detrending, gap filling, outliers, and result round-trip — mirrors ``io/``.

Readers ingest raw phase/frequency tables into ``PhaseData``/``FrequencyData``
(carrying the file path in ``source``); ``detrend``/``fillgaps``/``remove_outliers``
return cleaned records; ``save``/``load`` round-trip results, suites and records
through the oracle's self-describing text format.
"""

from __future__ import annotations

import numpy as np

from ._convert import (
    data_from_jl,
    np_int,
    result_from_jl,
    result_to_jl,
    suite_from_jl,
    suite_to_jl,
    to_jl_data,
)
from ._julia import _char, _map_errors, jl
from .types import FrequencyData, PhaseData, StabilityResult, StabilitySuite


def _read_kwargs(tau0, time_col, value_col, header, delim, scaling, detrend, fillgaps) -> dict:
    """Reader keyword arguments for the oracle: ``None`` ``tau0``/``delim`` are omitted."""
    kwargs = {
        "time_col": int(time_col),
        "value_col": int(value_col),
        "header": int(header),
        "scaling": float(scaling),
        "detrend": jl.Symbol(str(detrend)),
        "fillgaps": bool(fillgaps),
    }
    if tau0 is not None:
        kwargs["tau0"] = float(tau0)
    if delim is not None:
        kwargs["delim"] = _char(str(delim))
    return kwargs


@_map_errors
def read_phase(
    path: str,
    *,
    tau0: float | None = None,
    time_col: int = 1,
    value_col: int = 2,
    header: int = 0,
    delim: str | None = None,
    scaling: float = 1.0,
    detrend: str = "none",
    fillgaps: bool = False,
) -> PhaseData:
    """Read a phase-data file into a :class:`PhaseData` (carrying ``path`` in ``source``).

    ``tau0`` is auto-resolved from a ``# tau0:`` comment header or the time column
    (median Δt) when omitted, and required when ``time_col=0`` with neither
    fallback. ``time_col``/``value_col`` are 1-based; ``delim`` is auto-detected
    (``,`` for .csv, tab for .tsv, whitespace otherwise). ``scaling`` multiplies
    the samples; ``detrend`` ∈ {none, mean, endpoint, linear, quadratic};
    ``fillgaps`` equispaces on the min spacing and Howe-imputes NaNs.
    """
    return data_from_jl(
        jl.read_phase(
            str(path),
            **_read_kwargs(tau0, time_col, value_col, header, delim, scaling, detrend, fillgaps),
        )
    )


@_map_errors
def read_frequency(
    path: str,
    *,
    tau0: float | None = None,
    time_col: int = 1,
    value_col: int = 2,
    header: int = 0,
    delim: str | None = None,
    scaling: float = 1.0,
    detrend: str = "none",
    fillgaps: bool = False,
) -> FrequencyData:
    """Read a fractional-frequency file into a :class:`FrequencyData`.

    Same keyword arguments and comment-header handling as :func:`read_phase`;
    only the return type differs. The record carries ``path`` in ``source``.
    """
    return data_from_jl(
        jl.read_frequency(
            str(path),
            **_read_kwargs(tau0, time_col, value_col, header, delim, scaling, detrend, fillgaps),
        )
    )


@_map_errors
def detrend(
    data: PhaseData | FrequencyData, *, method: str = "linear"
) -> PhaseData | FrequencyData:
    """Return a new record with a detrended sample vector (the original is untouched).

    ``method`` ∈ {``linear`` (default), ``endpoint``, ``mean``, ``quadratic``, ``none``}.
    """
    return data_from_jl(jl.detrend(to_jl_data(data), method=jl.Symbol(str(method))))


@_map_errors
def fillgaps(data: PhaseData | FrequencyData) -> PhaseData | FrequencyData:
    """Impute ``NaN`` samples with Howe & Schlossberger's reflect-and-filter algorithm.

    Returns a new record with the same ``tau0``; the input is treated as already
    equispaced and returned unchanged when it has no ``NaN`` samples.
    """
    return data_from_jl(jl.fillgaps(to_jl_data(data)))


@_map_errors
def find_outliers(data: PhaseData | FrequencyData, *, n_sigma: float = 5.0) -> np.ndarray:
    """Indices of outlying samples via the Stable32-style scaled-MAD test.

    Returns a **0-based** integer array of sample indices (the oracle's 1-based
    indices minus one). For ``PhaseData`` the test runs on first differences, so a
    single phase glitch flags the two samples the step spans.
    """
    return np_int(jl.find_outliers(to_jl_data(data), nsigma=float(n_sigma))) - 1


@_map_errors
def remove_outliers(
    data: PhaseData | FrequencyData, *, n_sigma: float = 5.0
) -> PhaseData | FrequencyData:
    """Return a new record with scaled-MAD outliers gapped and Howe-imputed.

    The Stable32 deglitching idiom: flag samples more than ``n_sigma`` scaled MADs
    from the median, gap them, and fill. The input is returned as-is when nothing
    is flagged.
    """
    return data_from_jl(jl.remove_outliers(to_jl_data(data), nsigma=float(n_sigma)))


@_map_errors
def save(
    path: str,
    x: StabilityResult | StabilitySuite | PhaseData | FrequencyData,
    *,
    source_file: str = "",
) -> str:
    """Write a result, a whole suite, or a data record to a text file.

    Results and suites go to a tab-delimited file with a metadata header; a
    ``PhaseData``/``FrequencyData`` record goes to a SigmaTau data file that
    :func:`load` reads back.

    ``source_file`` is recorded in a suite's header (ignored for a single
    result). Returns ``path``. Read it back with :func:`load`.
    """
    if isinstance(x, StabilitySuite):
        jl.save(str(path), suite_to_jl(x), source_file=str(source_file))
    elif isinstance(x, StabilityResult):
        jl.save(str(path), result_to_jl(x))
    elif isinstance(x, (PhaseData, FrequencyData)):
        jl.save(str(path), to_jl_data(x))
    else:
        raise TypeError(
            f"expected StabilityResult, StabilitySuite, PhaseData or FrequencyData, "
            f"got {type(x).__name__}"
        )
    return str(path)


@_map_errors
def load(path: str) -> StabilityResult | StabilitySuite | PhaseData | FrequencyData:
    """Read a file written by :func:`save`; the file's header says which type comes back."""
    r = jl.load(str(path))
    if jl.isa(r, jl.StabilitySuite):
        return suite_from_jl(r)
    if jl.isa(r, jl.StabilityResult):
        return result_from_jl(r)
    return data_from_jl(r)


__all__ = [
    "read_phase",
    "read_frequency",
    "detrend",
    "fillgaps",
    "find_outliers",
    "remove_outliers",
    "save",
    "load",
]

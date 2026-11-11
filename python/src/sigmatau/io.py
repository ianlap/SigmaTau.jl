"""Reading and writing files, and cleaning records.

The readers return a ``PhaseData`` or ``FrequencyData`` record that carries the
path of its file as ``source``. ``detrend``, ``fillgaps`` and
``remove_outliers`` return a cleaned copy of a record, and ``save`` and
``load`` write and read records, results and suites as text.
"""

from __future__ import annotations

import numpy as np

from ._convert import (
    data_from_jl,
    driftfit_to_jl,
    np_f64,
    np_int,
    result_from_jl,
    result_to_jl,
    suite_from_jl,
    suite_to_jl,
    to_jl_char,
    to_jl_data,
)
from ._julia import _map_errors, jl
from .types import DriftFit, FrequencyData, PhaseData, StabilityResult, StabilitySuite


def _read_kwargs(tau0, time_col, value_col, header, delim, scaling, detrend, fillgaps) -> dict:
    """Reader keywords for Julia; a ``tau0`` or ``delim`` of ``None`` is left out."""
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
        kwargs["delim"] = to_jl_char(delim)
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
    """Read a text file of numeric columns as a ``PhaseData`` record: one sample on
    each line, in seconds, in column ``value_col``, with its time in column
    ``time_col``. A file with a single column is read as the samples alone.

    Lines that begin with ``#`` at the head of the file are skipped, and the
    ``# tau0:`` and ``# start:`` lines that :func:`save` writes are read. The
    sample interval is ``tau0``, else that of the header, else the step of the
    time column; an error is raised when the file gives none. The time column
    also places the samples: its steps are fitted to an even grid, and a
    reading that is missing from the file leaves a ``nan`` sample in its place.

    ``time_col`` and ``value_col`` count from 1, and ``time_col=0`` reads a
    file without a time column. ``header`` is the number of lines to skip
    after the leading ``#`` lines, such as a line of column names. ``delim``
    is the character that separates the fields: by default a tab in a ``.tsv``
    file whose first line of data holds one, else a comma when that line holds
    one, else runs of blanks. ``scaling`` multiplies the samples, ``1e-9`` for
    a file in nanoseconds. ``detrend`` names a trend to remove, one of the
    methods of :func:`detrend`. ``fillgaps=True`` fills the missing samples
    with :func:`fillgaps`.
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
    """Read a text file of numeric columns as a ``FrequencyData`` record of
    fractional-frequency samples. The header, the time column and the keywords
    are those of :func:`read_phase`."""
    return data_from_jl(
        jl.read_frequency(
            str(path),
            **_read_kwargs(tau0, time_col, value_col, header, delim, scaling, detrend, fillgaps),
        )
    )


@_map_errors
def read_data(
    path: str,
    *,
    format: str = "auto",
    tau0: float | None = None,
    column: int | str | None = None,
    nominal: float | None = None,
    kind: str = "auto",
    scaling: float = 1.0,
    delim: str | None = None,
    fillgaps: bool = False,
    keep_zeros: bool = False,
) -> PhaseData | FrequencyData:
    """Read a measurement file and return its record, phase or frequency as the
    file has it.

    ``format`` is one of the formats below, or ``"auto"`` to find it with
    :func:`detect_format`:

    - ``"columns"``: numeric columns, read as :func:`read_phase` reads them by
      default. The record is phase, or frequency when the median size of the
      samples is under 10⁻¹¹.
    - ``"stable32"``: a Stable32 ``.phd`` or ``.frd`` file. A zero is a missing
      sample, as in Stable32, except the first and the last sample of phase
      data.
    - ``"timelab"``: a TimeLab ``.tim`` file with a ``TIC`` block.
    - ``"counter"``: readings of a counter, one on each line, alone or after a
      timestamp. Readings of median size under 1 are taken as time intervals
      in seconds, and larger ones as frequencies in hertz about ``nominal``.
    - ``"table"``: columns under a line of names, as in a CSV export.
    - ``"bipm_utc"``: a UTC − UTC(k) file of the BIPM.
    - ``"rinex_clk"``: a RINEX clock file, of which ``column`` chooses the clock.

    ``tau0`` is the sample interval in seconds, by default that of the header
    of the file, else that of its timestamps. ``column`` is the clock to read:
    in a table, the number of its column, counted from 1, or its name; in a
    RINEX clock file, its position in order of first appearance or its name,
    such as ``"G08"``; in a plain column file, the number of the column of
    samples, counted from 1 (2 by default), and a name is an error. ``nominal`` is the nominal frequency in hertz for
    readings in hertz, by default the mean reading rounded to four significant
    figures. ``kind``, ``"phase"`` or ``"frequency"``, sets the type of the
    record whatever the file suggests. ``scaling`` multiplies the samples
    after any conversion of units. ``delim`` is the character that separates
    the fields of a columns, counter or table file. ``fillgaps=True`` fills
    the missing samples with :func:`fillgaps`, and ``keep_zeros=True`` reads
    the zeros of a Stable32 file as samples.

    The documentation of the Julia function says how each format is read.
    """
    kwargs = {
        "format": jl.Symbol(str(format)),
        "kind": jl.Symbol(str(kind)),
        "scaling": float(scaling),
        "fillgaps": bool(fillgaps),
        "keep_zeros": bool(keep_zeros),
    }
    if delim is not None:
        kwargs["delim"] = to_jl_char(delim)
    if tau0 is not None:
        kwargs["tau0"] = float(tau0)
    if column is not None:
        kwargs["column"] = column if isinstance(column, str) else int(column)
    if nominal is not None:
        kwargs["nominal"] = float(nominal)
    return data_from_jl(jl.read_data(str(path), **kwargs))


@_map_errors
def detect_format(path: str, *, delim: str | None = None) -> str:
    """Name the format of a measurement file from its extension and its first 200
    lines: ``"timelab"``, ``"stable32"``, ``"bipm_utc"``, ``"rinex_clk"``,
    ``"counter"``, ``"table"`` or ``"columns"``.

    The tests are made in that order, and the first that holds decides.
    ``delim`` is the character at which the counter and table tests split a
    line: by default a comma, or a tab, when the first line of data holds one,
    and runs of blanks otherwise.
    """
    kwargs = {} if delim is None else {"delim": to_jl_char(delim)}
    return str(jl.detect_format(str(path), **kwargs))


@_map_errors
def detrend(
    data: PhaseData | FrequencyData, fit: DriftFit | None = None, *, method: str = "linear"
) -> PhaseData | FrequencyData:
    """Copy of a record with a trend removed from its samples.

    ``method`` is ``"mean"``, ``"endpoint"`` (the line through the first and
    last samples), ``"linear"`` (a least-squares line), ``"quadratic"`` (a
    least-squares parabola) or ``"none"``. In phase data a line is a frequency
    offset and a parabola a linear frequency drift; in frequency data a line
    is a linear frequency drift (SP1065 §5.18, §5.21 and §5.22). Missing
    samples stay in place and do not enter the fit.

    With ``fit``, the result of :func:`drift_fit` on the same record, the
    frequency model of that fit is removed instead: ``fit.fit`` from the
    frequency samples, or its running sum times tau0 from the phase samples.
    """
    if fit is not None:
        return data_from_jl(jl.detrend(to_jl_data(data), driftfit_to_jl(fit)))
    return data_from_jl(jl.detrend(to_jl_data(data), method=jl.Symbol(str(method))))


@_map_errors
def fillgaps(data: PhaseData | FrequencyData, *, form: str = "auto") -> PhaseData | FrequencyData:
    """Copy of a record with its missing samples (``nan``) filled by the reflection
    method of Howe and Schlossberger (2022).

    ``form`` is the series that is filled: ``"phase"``, ``"frequency"``, or
    ``"auto"`` to choose it from the noise of the record. The measured samples
    are returned as they came. The paper prescribes neither the series nor the
    length of its filter, and this package sets both from the record.

    Fill only for a function that needs a complete record, such as the
    spectra, the total, Thêo, parabolic and dynamic deviations and the cross
    estimators: the Allan and Hadamard deviations leave out the terms that a
    missing sample enters. A filled record is not the record, and the
    documentation of the Julia function states what a fill does not recover.
    """
    return data_from_jl(jl.fillgaps(to_jl_data(data), form=jl.Symbol(str(form))))


@_map_errors
def find_outliers(data: PhaseData | FrequencyData, *, n_sigma: float = 5.0) -> np.ndarray:
    """Indices of the outlying samples of a record, counted from 0, by the median
    absolute deviation test of SP1065 §10.11 and of Stable32's Check function.

    A frequency sample is an outlier when it lies more than ``n_sigma`` times
    MAD/0.6745 from the median. For a ``PhaseData`` record the test is made on
    the first differences, so a single displaced phase sample x_j marks
    samples j and j + 1. Missing samples are never marked, and the test is
    applied once.
    """
    return np_int(jl.find_outliers(to_jl_data(data), nsigma=float(n_sigma))) - 1


@_map_errors
def summarize(data: PhaseData | FrequencyData) -> dict:
    """A quick look at a record before analysis, as the explorer's Data page shows
    it.

    Keys: ``n`` (number of samples), ``tau0`` and ``span`` (n·τ₀), in seconds,
    ``gaps`` (number of ``nan`` samples), ``outliers`` (the indices that
    :func:`find_outliers` returns at ``n_sigma=5``) and ``drift_per_day``
    (least-squares slope of the fractional frequency against time, per day).
    """
    s = jl.summarize(to_jl_data(data))
    return {
        "n": int(s.n),
        "tau0": float(s.tau0),
        "span": float(s.span),
        "gaps": int(s.gaps),
        "outliers": np_int(s.outliers) - 1,
        "drift_per_day": float(s.drift_per_day),
    }


@_map_errors
def sample_statistics(data: PhaseData | FrequencyData) -> dict:
    """Basic statistics of a record's samples, those of Stable32's Statistics
    function.

    Keys: ``n`` (samples measured), ``gaps`` (``nan`` samples), ``min``,
    ``max``, ``mean``, ``median`` and ``std`` (with n − 1) of the samples
    measured, phase in seconds or fractional frequency, ``adev`` (the Allan
    deviation at τ₀) and ``noise`` (the noise type at τ₀ that
    :func:`identify_noise` gives, on the record with each gap bridged by a
    straight line). The record needs three samples measured.

    ``std_y`` is the standard deviation of the fractional frequency, the
    Std Dev that Stable32 shows for either kind of record; ``std`` is that of
    the samples as given.
    """
    s = jl.sample_statistics(to_jl_data(data))
    out = {k: float(getattr(s, k)) for k in ("min", "max", "mean", "median", "std", "std_y", "adev")}
    return {"n": int(s.n), "gaps": int(s.gaps), **out, "noise": str(s.noise)}


@_map_errors
def sample_histogram(data: PhaseData | FrequencyData, *, bins: int | None = None) -> dict:
    """Histogram of a record's samples, as Stable32's Histogram function plots
    it.

    The bins have equal widths and span the samples measured; ``bins=None``
    takes their number from Scott's rule, h = 3.49·s·n^(−1/3) (Scott 1979),
    with at most 1000. Keys: ``edges`` (the bins + 1 edges, rising),
    ``counts`` (each bin holds its lower edge, the last also its upper one),
    and ``n``, ``mean`` and ``std`` of the samples measured, for the normal
    curve n·h·φ((v − mean)/std)/std over a bin of width h.
    """
    kw = {} if bins is None else {"bins": int(bins)}
    h = jl.sample_histogram(to_jl_data(data), **kw)
    return {
        "edges": np_f64(h.edges),
        "counts": np_int(h.counts),
        "n": int(h.n),
        "mean": float(h.mean),
        "std": float(h.std),
    }


@_map_errors
def acf(data: PhaseData | FrequencyData, *, m: int = 1, lags: int | None = None) -> dict:
    """Autocorrelation of a record at averaging factor ``m``, as Stable32's
    Autocorrelation function plots it (SP1065 §5.5.3, Eq. 47).

    The series is the phase taken every m samples, or the frequency averaged
    over m samples without overlap. ``lags`` is the largest lag, by default
    half the series. Keys: ``lag`` and ``r`` (lags 0 to K and their
    autocorrelation), ``m``, ``tau`` (m·τ₀) and ``n`` (points in the series),
    ``alpha`` (the lag-1 estimate of the power-law exponent, ``nan`` with
    fewer than 30 points) and ``noise`` (the type :func:`identify_noise`
    gives at m). A record with ``nan`` gaps raises an error: fill them first.
    """
    kw = {"m": int(m)} if lags is None else {"m": int(m), "lags": int(lags)}
    a = jl.acf(to_jl_data(data), **kw)
    return {
        "lag": np_int(a.lag),
        "r": np_f64(a.r),
        "m": int(a.m),
        "tau": float(a.tau),
        "n": int(a.n),
        "alpha": float(a.alpha),
        "noise": str(a.noise),
    }


@_map_errors
def remove_outliers(
    data: PhaseData | FrequencyData, *, n_sigma: float = 5.0
) -> PhaseData | FrequencyData:
    """Copy of a record with the outliers that :func:`find_outliers` finds at
    ``n_sigma`` replaced by filled values.

    An isolated outlier takes the mean of its two neighbors, and a run of them
    is filled by reflection, as in :func:`fillgaps`. A ``PhaseData`` record is
    cleaned as frequency, which takes a phase jump out. Samples that were
    missing stay missing, and a record without outliers is returned as it is.

    SP1065 §10.11 and Stable32 replace an outlier by a gap, which the
    estimators then skip. For that, set the samples that :func:`find_outliers`
    marks in a frequency record to ``nan``.
    """
    return data_from_jl(jl.remove_outliers(to_jl_data(data), nsigma=float(n_sigma)))


@_map_errors
def save(
    path: str,
    x: StabilityResult | StabilitySuite | PhaseData | FrequencyData,
    *,
    source_file: str = "",
) -> str:
    """Write a record, a result or a suite to a text file and return ``path``.
    :func:`load` reads the file back.

    A ``PhaseData`` or ``FrequencyData`` record is written as two columns, the
    sample time and the sample, after a comment header with its ``tau0`` and
    start time. The columns are separated by a comma when ``path`` ends in
    ``.csv`` and by a tab otherwise, and :func:`read_phase` and
    :func:`read_frequency` read the file as well. A ``StabilityResult`` is
    written as a tab-separated table with one row for each τ, and a
    ``StabilitySuite`` as one such table for each result, with ``source_file``
    in its header.
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
    """Read a file written by :func:`save`. Its first line says what the file holds:
    a result, a suite, or a phase or frequency record."""
    r = jl.load(str(path))
    if jl.isa(r, jl.StabilitySuite):
        return suite_from_jl(r)
    if jl.isa(r, jl.StabilityResult):
        return result_from_jl(r)
    return data_from_jl(r)


__all__ = [
    "read_data",
    "detect_format",
    "read_phase",
    "read_frequency",
    "detrend",
    "fillgaps",
    "find_outliers",
    "remove_outliers",
    "summarize",
    "save",
    "load",
]

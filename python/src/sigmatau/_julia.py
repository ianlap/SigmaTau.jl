"""Julia bootstrap: one embedded Julia runtime with SigmaTau.jl loaded.

Importing this module starts (or attaches to) the in-process Julia runtime via
juliacall, loads SigmaTau, and caches the conversion targets used at the
Python↔Julia boundary. Every other module imports ``jl`` from here.
"""

from __future__ import annotations

import functools
import os
import sys
import warnings

# juliacall reads its configuration from the environment at import time: signal
# handling must be Julia's (or Ctrl-C during a long kernel can crash the
# process), and the thread count must be >1 or the threaded total-family
# kernels (mtotdev/htotdev/mhtotdev) run serial, ~3-4x slower. Explicit user
# settings are respected.
if "juliacall" in sys.modules:
    warnings.warn(
        "juliacall was imported before sigmatau, so sigmatau could not "
        "configure it (PYTHON_JULIACALL_HANDLE_SIGNALS=yes, "
        "PYTHON_JULIACALL_THREADS=auto); Ctrl-C during a Julia call may crash "
        "the process, and the threaded kernels may run serial.",
        stacklevel=2,
    )
else:
    if os.environ.get("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes") != "yes":
        warnings.warn(
            "PYTHON_JULIACALL_HANDLE_SIGNALS="
            f"{os.environ['PYTHON_JULIACALL_HANDLE_SIGNALS']!r}; sigmatau "
            "expects 'yes' — Ctrl-C during a Julia call may crash the process.",
            stacklevel=2,
        )
    os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")
    os.environ.setdefault("PYTHON_JULIACALL_THREADS", "auto")

import juliacall
from juliacall import JuliaError

# A private module keeps SigmaTau's exports out of the process-shared Main
# namespace (another embedded-Julia package may also be running here).
jl = juliacall.newmodule("SigmaTauPy")
jl.seval("using SigmaTau")
jl.seval("import Random")

# Conversion targets and helpers, seval'd once. `juliacall.convert` needs the
# Julia type objects; the lambdas cover constructs attribute access can't reach.
VecF64 = jl.seval("Vector{Float64}")
VecInt = jl.seval("Vector{Int}")
VecStr = jl.seval("Vector{String}")
VecResult = jl.seval("Vector{StabilityResult}")
DictIntF64 = jl.seval("Dict{Int,Float64}")
_to_syms = jl.seval("v -> Symbol.(v)")  # Vector{String} -> Vector{Symbol}
_char = jl.seval("only")  # 1-char str -> Char at the io boundary

# Rebuild a Julia StabilityResult from unpacked Python fields (save_*/nch paths).
_make_result = jl.seval("""
    (dt, tau, dev, noise, lo, hi, edf, neff) -> StabilityResult(
        Symbol(dt), tau, dev, noise,
        [(lo = l, hi = h) for (l, h) in zip(lo, hi)], edf, neff)
""")

# Classify a Julia exception for _map_errors: (python-error-kind, message).
_err_info = jl.seval("""
    e -> e isa ArgumentError ? ("value", e.msg) :
         e isa DomainError   ? ("value", sprint(showerror, e)) :
         e isa KeyError      ? ("key",   sprint(showerror, e)) :
         e isa MethodError   ? ("type",  sprint(showerror, e)) :
         e isa SystemError   ? ("os",    sprint(showerror, e)) :
         e isa Base.IOError  ? ("os",    sprint(showerror, e)) :
                               ("other", "")
""")

_PY_ERR = {"value": ValueError, "key": KeyError, "type": TypeError, "os": OSError}


def _map_errors(fn):
    """Re-raise Julia exceptions as native Python ones (ValueError, KeyError, …).

    The original JuliaError stays reachable as ``exc.julia_error``; the Julia
    backtrace is dropped from the traceback (``raise … from None``).
    """

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except JuliaError as e:
            kind, msg = _err_info(e.exception)
            if kind == "other":
                raise
            err = _PY_ERR[kind](str(msg))
            err.julia_error = e
            raise err from None

    return wrapper


def julia_version() -> str:
    """Version of the embedded Julia runtime (for bug reports)."""
    return str(jl.seval("string(VERSION)"))


def sigmatau_jl_version() -> str:
    """Version of the underlying SigmaTau.jl package."""
    return str(jl.seval("string(pkgversion(SigmaTau))"))

# SigmaTau Python Port Outline

Last updated: 2026-06-07.

Do not start the Python package until the Julia `0.4.0` release is tagged and
can serve as the frozen reference.

## Goal

Build a Python-native package that keeps SigmaTau's validated math while being
easier to use and faster on long records than allantools.

Positioning:

- Faster than allantools on long records.
- More complete result objects than allantools' tuple returns.
- Integrated per-tau noise ID, EDF, confidence intervals, and bias correction.
- Calibrated composite power-law noise generation.
- Julia-compatible TSV result and suite round trips.

## Proposed Package Layout

Mirrors the flattened Julia `src/` — one Python module per Julia source file,
same names — so the two libraries stay structurally aligned.

```text
pyproject.toml
src/
  sigmatau/
    __init__.py
    types.py        # <- types.jl       (PhaseData, FrequencyData, StabilityResult, StabilitySuite, SpectralResult)
    grids.py        # <- grids.jl        (TauMode, tau_values, _freq_to_phase, _f64, _default_m_values)
    kernels.py      # <- kernels.jl      (internal _*_core array kernels)
    deviations.py   # <- deviations.jl   (public adev … pdev)
    suite.py        # <- suite.jl        (stability compute-all entry point)
    edf.py          # <- edf.jl          (EDF / CI / bias)
    noise.py        # <- noise.jl        (identify_noise, synthesis, noise_gen)
    spectral.py     # <- spectral.jl     (Welch core + Sy/Sx/L)
    io/             # <- io/             (read, detrend, fillgaps, results)
    plotting.py     # (no Julia analogue; mirrors the RecipesBase extension)
tests/
  fixtures/julia/   # golden outputs + shared input records exported from the oracle
docs/
```

The MVP (milestone 2) ships only `types.py`, `grids.py`, `kernels.py`, and
`deviations.py`; the remaining modules land as their Julia counterparts get
ported.

## API Sketch

```python
from sigmatau import PhaseData, Octave, stability, pdev

data = PhaseData(x, tau0=1.0)
result = pdev(data, taus=Octave, ci=True, confidence=0.683)
suite = stability(data, devs=("adev", "mdev", "hdev", "pdev"))
suite["pdev"].dev
```

The port mirrors the Julia oracle's names: `stability`, `devs`, `taus`, `ci`,
`tau0` (default `1.0`), and `confidence` are identical across both libraries —
only the surface casing differs (Python keyword args vs Julia's). Result fields
(`tau`, `dev`, `noise_type`, `ci_lower`, `ci_upper`, `edf`) match one-to-one.

Core result types:

```python
@dataclass(frozen=True)
class StabilityResult:
    deviation_type: str
    tau: np.ndarray
    dev: np.ndarray
    noise_type: np.ndarray
    ci_lower: np.ndarray
    ci_upper: np.ndarray
    edf: np.ndarray

@dataclass(frozen=True)
class StabilitySuite:
    results: tuple[StabilityResult, ...]
    data_kind: str
    tau0: float
    n: int
    confidence: float | None
    tau_mode: str
```

## Implementation Strategy

Use NumPy plus Numba for kernels.

Kernel rules:

- Convert inputs to contiguous `float64` arrays at the API boundary.
- Keep kernels array based and free of Python objects.
- Use `@njit(cache=True)` for scalar-loop kernels.
- Use `parallel=True` only where work is large enough and reduction order is
  controlled.
- Do not use `fastmath=True` until parity has been established and measured.
- For modified-total kernels, parallelize over fixed subsequence chunks, not
  over `m`.
- Implement PDEV with SigmaTau.jl's rolling O(N) recurrence.
- Implement MTIE with the monotonic deque algorithm.

Recommended dependencies:

- Required: `numpy`, `numba`, `scipy`.
- Optional: `matplotlib`, `pandas`.
- Developer: `pytest`, `ruff`, `mypy` or `pyright`, `pytest-benchmark`.

## Validation Plan

Phase 1 fixtures:

- Generate golden SigmaTau.jl TSV fixtures for each deviation on deterministic
  phase and frequency records.
- Include Stable32 and allantools fixtures already used by Julia.
- Test Python results against Julia before comparing speed.

Phase 2 properties:

- Constant and linear phase records produce expected zeros where applicable.
- `FrequencyData` dispatch matches phase integration.
- `tdev = mdev * tau / sqrt(3)`.
- `htdev = mhdev * tau / sqrt(10/3)`.
- `pdev(tau0) = adev(tau0)`.
- `ci=False` leaves CI/EDF/noise arrays empty.

Phase 3 Monte Carlo:

- Keep Monte Carlo tests small in CI.
- Store long Monte Carlo artifacts outside normal unit tests.
- Document seed, generator, N grid, tau grid, and fit method.

## Milestones

1. Freeze Julia `0.4.0` and export golden fixtures.
2. Build Python MVP: data/result types, tau grids, ADEV/MDEV/TDEV/HDEV/MHDEV.
3. Add long-record kernels: total family, MTIE, PDEV.
4. Add stats/noise: EDF, confidence intervals, bias correction, `noise_gen`.
5. Add IO, spectral estimators, plotting helpers, and benchmark docs.

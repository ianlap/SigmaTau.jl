# sigmatau

Python bindings for [SigmaTau.jl](https://github.com/ianlap/SigmaTau.jl). `sigmatau`
is a [juliacall](https://github.com/JuliaPy/PythonCall.jl)-based wrapper — the same
architecture as [PySR](https://github.com/MilesCranmer/PySR) — where every
computation runs in an embedded Julia runtime with SigmaTau.jl loaded, so results are
bit-identical to the Julia package, not a reimplementation. It supersedes the earlier
pure-Python/NumPy port ([sigmatau.py](https://github.com/ianlap/sigmatau.py)), which
is being deprecated; this package mirrors that port's public API and extends it.

## Install

Not yet on PyPI. For now, install from a checkout of this repo:

```bash
pip install -e "python/[dev]"
```

The first `import sigmatau` performs a one-time
[juliapkg](https://github.com/JuliaPy/pyjuliapkg) resolve: it installs Julia 1.11
(via juliaup if needed), creates a private Julia environment under the venv,
dev-installs this repo into it, and precompiles — about two minutes. Warm imports
after that take a few seconds (Julia boot time).

Import `sigmatau` before `juliacall`. On import, `sigmatau` sets
`PYTHON_JULIACALL_HANDLE_SIGNALS=yes` (so Ctrl-C during a long-running call
interrupts Julia cleanly instead of risking a crash) and
`PYTHON_JULIACALL_THREADS=auto` (so the threaded total-family kernels get all
cores — the juliacall default is a single thread). Both respect values you set
yourself; if `juliacall` is imported first, neither can take effect and
`sigmatau` warns.

## 30-second example

```python
import sigmatau as st

pd = st.noise_gen(st.PhaseData, 4096, 1.0, sigma1={0: 1e-12}, rng=42)
r = st.adev(pd)                      # Octave grid, noise ID + 68.3% CIs
print(r.tau, r.dev, r.ci_lower, r.ci_upper, r.noise_type)

suite = st.stability(pd)             # adev+mdev+hdev+mhdev in one call
suite["adev"].dev
```

Explicit τ (averaging factors instead of a named grid):

```python
st.mhtotdev(pd, [1, 2, 4, 8])
```

Reading a data file:

```python
st.read_phase("clock.dat", tau0=1.0)
```

Dynamic (time-resolved) deviation:

```python
st.dadev(pd, window=512)
```

## Development configuration vs release

`python/src/sigmatau/juliapkg.json` currently points at this checkout:

```json
{
    "julia": "1.11",
    "packages": {
        "SigmaTau": {
            "uuid": "c34468cf-613a-4a89-b990-738073237867",
            "path": "../../..",
            "dev": true
        }
    }
}
```

The `"dev": true` is load-bearing: without it, juliapkg `Pkg.add`s a git snapshot of
the path instead of dev-installing it, and local Julia edits become invisible to the
Python package until the next resolve.

For a PyPI release, this entry must instead pin a released SigmaTau.jl. Until
SigmaTau.jl is registered in Julia's General registry:

```json
{ "url": "https://github.com/ianlap/SigmaTau.jl", "rev": "v0.5.0" }
```

and once it is registered:

```json
{ "version": "=0.5.0" }
```

The wrapper's version stays in lockstep with SigmaTau.jl's; wrapper-only fixes append
a `.postN` suffix rather than bumping the SigmaTau.jl version they track.

## Differences from sigmatau.py (the retired pure-Python port)

1. `DEFAULT_DEVIATIONS` is `("adev", "mdev", "hdev", "mhdev")` — the Julia truth. The
   port's fourth slot was `tdev`.
2. `StabilityResult` gains a `neff` field.
3. `PhaseData`/`FrequencyData` gain a `source` field; the file readers record the
   path they were read from.
4. `noise_gen(rng=int)` seeds Julia's Xoshiro RNG. The same seed produces different
   samples than the port's NumPy stream — draws are reproducible within this wrapper,
   not bit-compatible with the port.
5. New API the port never had: `tierms`, `theo1`, `theoh`, `dadev`, `dhdev`, `nch`,
   `find_outliers` (returns 0-based indices), `remove_outliers`, plus `tau_values`
   grids for the new kernels.
6. `theo1` reports τ = 0.75·m·τ₀ (the Julia/SP1065 convention), not τ = m·τ₀.

## Testing

From the repo root, with the package installed into the active venv:

```bash
python -m pytest python/tests
```

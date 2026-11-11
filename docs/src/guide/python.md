# Python

`sigmatau` is the Python package of this project. It is a wrapper, not a
second implementation: each function calls the Julia function of the same
name in a Julia runtime inside the Python process, through
[juliacall](https://github.com/JuliaPy/PythonCall.jl). The numbers are
those of the Julia package, and this documentation applies to both.

## Install

```bash
pip install sigmatau        # the latest release on PyPI
pip install -e python       # from a checkout of the repository, as it stands
```

The package needs Python 3.11 or later. From a checkout the install has to
be an editable one (`-e`), because the package finds the Julia code by its
path in the checkout. The first `import sigmatau` finds a Julia of version
1.11 or later, or installs one, creates a Julia environment inside the
Python environment and precompiles the package, which takes a few minutes.
Later imports take a few seconds.

Import `sigmatau` before `juliacall`. On import it asks juliacall to let
Julia handle signals, so that Ctrl-C interrupts a long calculation, and to
start Julia with all threads, for the total and Thêo estimators. It cannot
do either once `juliacall` is loaded, and warns. Leave Python's
`faulthandler` off: it takes the signals that Julia's threads use among
themselves for crashes.

## A first analysis

```python
import sigmatau as st

p = st.noise_gen(st.PhaseData, 4096, 1.0, sigma1={0: 1e-12}, rng=42)

r = st.adev(p)                    # octave spacing, noise type and intervals
r.tau, r.dev, r.noise_type, r.ci_lower, r.ci_upper, r.edf

r = st.adev(p, [1, 10, 100], confidence=0.95, noise="WHFM")

suite = st.stability(p, devs=("adev", "mdev", "hdev"))
suite["mdev"].dev

st.save("adev.tsv", r)            # st.load reads it back
```

A record of your own is `st.PhaseData(x, tau0)` or
`st.FrequencyData(y, tau0)` from a NumPy array, or the result of
`st.read_phase`, `st.read_frequency` or `st.read_data` on a file.

## From Julia to Python

The names and keywords are those of the Julia functions, so the
[Reference](../reference/estimators.md) serves for both. The differences
follow from the languages:

| Julia | Python |
|:---|:---|
| a symbol, `:WHFM`, `:stable32` | a string, `"WHFM"`, `"stable32"` |
| a named grid, `adev(p, Decade)` | `st.adev(p, st.Decade)` |
| `r.ci`, a vector of `(lo, hi)` | `r.ci_lower`, `r.ci_upper`, two arrays |
| `suite[:mdev]` | `suite["mdev"]` |
| a named tuple, `summarize(p).gaps` | a dict, `st.summarize(p)["gaps"]` |
| indices that count from 1 | indices that count from 0: `find_outliers`, the `outliers` of `summarize`, the `reference` and `pivot` of the ensemble functions |
| `find_outliers(p; nsigma = 5)` | `st.find_outliers(p, n_sigma=5)` |
| `mask_check(r, :G811).pass` | `st.mask_check(r, "G811").passed` |
| `FrequencyData(p)`, `PhaseData(f)` | `p.to_frequency()`, `f.to_phase()` |
| `Dict(0 => 1e-12)` | `{0: 1e-12}` |
| `rng = Xoshiro(42)` | `rng=42`, which seeds the same generator |
| `push!(s, x)`, `append!(s, xs)`, `snapshot(s)` on a [`StreamingStability`](@ref) | `s.push(x)`, `s.extend(xs)`, `s.snapshot()` |
| `add_component!(fit, :WHFM, 1:6)`, `undo!(fit)`, `noise_residual(fit)` on a [`NoiseFit`](@ref) | `fit.add_component("WHFM", slice(0, 6))`, `fit.undo()`, `fit.residual()` |
| `ClockComparisons(records)`, from one phase record per clock | `st.ClockComparisons.from_records(records)` |
| an `ArgumentError` | a `ValueError` |
| a `MethodError`, as for a float where an integer counts samples | a `TypeError` |

A float is refused even when it is a whole number, so a grid of averaging
factors built with NumPy is converted first:
`st.adev(p, np.unique(np.round(np.logspace(0, 3, 30)).astype(int)))`.

Arrays come back as NumPy arrays and results as Python objects with the
fields of the Julia types.

Not wrapped: live measurement ([`LiveStability`](@ref) and its sources)
and the plot recipes. Plot a result with any Python library from `r.tau`
and `r.dev`.

## More of the package

```python
st.deviation("mhtotdev", p)                 # any estimator by name
st.psd_y(p)                                 # spectrum of fractional frequency
st.drift_fit(p).rate                        # frequency drift, per day
st.mask_check(st.mtie(p), "G8262_opt1").passed
st.bias_b2(36, "FLFM")                      # dead-time bias function B2(r, μ)

fit = st.powerlaw_fit(p)                    # power-law model of the noise
model = st.clock_model(fit)                 # clock model with those levels
kf = st.kalman_filter(p, model)             # phase and frequency estimates
gains = st.steering_gains(3600.0, p.tau0)   # a steering loop, time constant 1 h
st.steer(p, model, gains).steered           # the record as that loop would leave it

cc, truth = st.simulate_ensemble([model] * 4, 10_000, rng=1)
st.timescale(cc).weights                    # AT1 weights against time
```

The package's own
[README](https://github.com/ianlap/SigmaTau.jl/blob/main/python/README.md)
covers development installs and the release procedure.

# sigmatau

The Python package of [SigmaTau.jl](https://github.com/ianlap/SigmaTau.jl),
for the analysis of clock and oscillator stability. It is a wrapper, not a
second implementation: each function calls the Julia function of the same
name in a Julia runtime inside the Python process, through
[juliacall](https://github.com/JuliaPy/PythonCall.jl). The numbers are
those of the Julia package, and the
[documentation](https://ianlap.github.io/SigmaTau.jl/dev/) serves for
both. Its [Python page](https://ianlap.github.io/SigmaTau.jl/dev/guide/python/)
lists what differs between the two languages.

## Install

```bash
pip install sigmatau        # the latest release on PyPI
pip install -e python       # from a checkout of the repository, as it stands
```

The package needs Python 3.11 or later. From a checkout the install has to
be an editable one (`-e`), because the package finds the Julia code by its
path in the checkout.

The first `import sigmatau` sets up Julia through
[juliapkg](https://github.com/JuliaPy/pyjuliapkg): it finds a Julia of
version 1.11 or later, or installs one, creates a Julia environment inside
the Python environment, installs SigmaTau.jl there and precompiles it. That
takes a few minutes. Later imports take a few seconds.

Import `sigmatau` before `juliacall`. On import it sets
`PYTHON_JULIACALL_HANDLE_SIGNALS=yes`, so that Ctrl-C interrupts a long
calculation in Julia, and `PYTHON_JULIACALL_THREADS=auto`, so that the
estimators that use threads run on all cores. It leaves a value that you
have set yourself. Once `juliacall` is loaded neither can take effect, and
`sigmatau` warns.

Leave Python's `faulthandler` off: it takes the signals that Julia's
threads use among themselves for crashes. Under pytest that is
`-p no:faulthandler`, as in the configuration of this package.

## Example

```python
import sigmatau as st

p = st.noise_gen(st.PhaseData, 8192, 1.0, sigma1={0: 1e-12, -2: 3e-14}, rng=3)

r = st.adev(p)                     # octave spacing, with intervals
r.tau, r.dev, r.ci_lower, r.ci_upper, r.noise_type

suite = st.stability(p)            # ADEV, MDEV, HDEV and MHDEV in one call
suite["mdev"].dev

st.mhtotdev(p, [1, 2, 4, 8])       # averaging factors of your own
st.deviation("theo1", p)           # any estimator by name: st.estimators()
st.save("adev.tsv", r)             # st.load("adev.tsv") reads it back
```

A record of your own comes from a file or from a NumPy array:

```python
p = st.read_data("clock.csv")      # also st.read_phase, st.read_frequency
p = st.PhaseData(x, 1.0)           # phase in seconds, sample interval in seconds
f = st.FrequencyData(y, 1.0)       # fractional frequency
```

The rest of the Julia package is there under the same names:

```python
st.summarize(p)                                  # samples, gaps, outliers, drift
st.psd_y(p)                                      # spectrum of the frequency
st.dadev(p, window=512)                          # deviation in a sliding window
st.drift_fit(p, model="log").rate                # a drift model, rate per day
st.mask_check(st.tdev(p), "G8262_opt1").passed   # an ITU-T wander mask

fit = st.powerlaw_fit(p)                         # power-law model of the noise
model = st.clock_model(fit)                      # clock model with those levels
kf = st.kalman_smooth(st.kalman_filter(p, model))
st.holdover(p, train=1024, limit=1e-9).time_to_limit

gains = st.steering_gains(3600.0, p.tau0)        # a steering loop with a time constant of 1 h
st.adev(st.steer(p, model, gains).steered)       # the record as that loop would leave it

cc, truth = st.simulate_ensemble([model] * 3, 10_000, rng=1)
st.timescale(cc, algorithm="at1").weights        # a time scale from three clocks
```

Of the Julia package, the Python package leaves out live measurement
(`LiveStability`, the instrument sources and `read_reading`), the plot
recipes, and the helpers `nsamples`, `tau_mode` and `INPUT_FORMATS`.

## Development

`python/src/sigmatau/juliapkg.json` names the checkout:

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

With `"dev": true` juliapkg installs the path itself, as `Pkg.develop`
does, so the Python package runs the Julia source of the checkout as it
stands.

```bash
pip install -e "python/[dev]"
pytest python/tests -q
```

`sigmatau.__version__` equals the version of SigmaTau.jl, and a test
checks it.

## Releases

A wheel cannot name a path. `python/tools/pin_juliapkg.py` replaces the
entry above with the release of SigmaTau.jl of the same version in the
General registry, `"version": "=X.Y.Z"`, and stops when that version is
not registered.

`.github/workflows/PythonPublish.yml` runs on a tag `python-vX.Y.Z`. It
checks that the tag, `__version__` and `Project.toml` agree, runs the pin
script, builds the source distribution and the wheel, installs the wheel
on a clean runner, runs the tests against it and publishes through PyPI
trusted publishing. Push the tag only after that version of SigmaTau.jl is
in the General registry. The header of the workflow lists the one-time
setup on PyPI and on GitHub.

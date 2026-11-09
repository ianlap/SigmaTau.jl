# Python bindings

SigmaTau ships an official Python package, `sigmatau`, in `python/`. It is a
[juliacall](https://github.com/JuliaPy/PythonCall.jl)-based wrapper: every
computation runs in an embedded Julia runtime with SigmaTau.jl loaded, so results
match the Julia API on this site bit-for-bit — there is one numerics source, not a
parallel Python implementation. It supersedes the earlier pure-Python/NumPy port
([sigmatau.py](https://github.com/ianlap/sigmatau.py)), which is being deprecated;
the wrapper mirrors that port's public API and extends it.

## Install

```bash
pip install sigmatau                 # or, from a checkout: pip install -e "python/[dev]"
```

The first `import sigmatau` performs a one-time juliapkg resolve — installing Julia
1.11, creating a private Julia environment under the venv, and precompiling — which
takes about two minutes. Warm imports after that take a few seconds.

## Example

```python
import sigmatau as st

pd = st.noise_gen(st.PhaseData, 4096, 1.0, sigma1={0: 1e-12}, rng=42)
r = st.adev(pd)                      # Octave grid, noise ID + 68.3% CIs
print(r.tau, r.dev, r.ci_lower, r.ci_upper, r.noise_type)

suite = st.stability(pd)             # adev+mdev+hdev+mhdev in one call
suite["adev"].dev

st.deviation("mhtotdev", pd)         # any estimator by name; see st.estimators()
st.psd_y(pd)                         # Welch S_y(f)
st.drift_fit(pd).rate                # frequency drift per day
st.powerlaw_fit(pd)                  # h_α of the power-law noise model
st.save("adev.tsv", r)               # st.load("adev.tsv") reads it back

m = st.clock_model(st.powerlaw_fit(pd), tau0=1.0)
kf = st.kalman_filter(pd, m)         # phase and frequency estimates
cc, _ = st.simulate_ensemble([m] * 4, 10_000)
st.timescale(cc).weights             # AT1 weights over time
```

The Python names and keywords match the Julia ones on this site (`psd_y`,
`ncornered_hat`, `groslambert`, `cross_adev`, `save`/`load`, `m=` for the τ
grid). Indices into and out of Python, such as `find_outliers` and
`simulate_ensemble(reference=)`, are 0-based. Live
measurement ([`LiveStability`](@ref) and `sigmatau live`), `noise_id_nnls` and
[`StreamingStability`](@ref) are Julia only.

## More

`python/README.md` is the canonical reference: full install notes, more examples
(explicit τ grids, file IO, dynamic deviations), the development-vs-release
`juliapkg.json` configuration, the divergence ledger against the retired
`sigmatau.py` port, and how to run the test suite.

# Coming from allantools

allantools and SigmaTau implement the same SP1065 estimators; on the
shipped fixture their raw deviations agree to 1e-13 or better
([Validation](../explanation/validation.md)). SigmaTau adds noise
identification and confidence intervals to every call. This page maps
allantools calls onto SigmaTau from Python or Julia, and lists where the
numbers differ.

## The call

allantools takes an array, a sample rate and τ in seconds, and returns four
arrays. SigmaTau takes a typed record with its τ₀ and averaging factors
`m` (τ = m·τ₀), and returns one result:

```python
import allantools
taus, dev, err, n = allantools.oadev(x, rate=10.0, data_type="phase", taus="octave")

import sigmatau as st
r = st.adev(st.PhaseData(x, tau0=0.1))          # Octave grid by default
r.tau, r.dev, r.ci_lower, r.ci_upper, r.noise_type, r.edf
```

```julia
using SigmaTau
r = adev(PhaseData(x, 0.1))
```

| allantools | SigmaTau |
|---|---|
| `data_type="phase"` / `"freq"` | `PhaseData(x, tau0)` / `FrequencyData(y, tau0)` |
| `rate` (Hz) | `tau0 = 1 / rate` (s) |
| `taus="octave"`, `"decade"`, `"all"` | `Octave` (default), `Decade`, `AllTaus`; also `HalfOctave`, `QuarterOctave`, `HalfDecade` |
| `taus=[1.0, 10.0]` (s) | `m=[10, 100]` in Python, `adev(p, [10, 100])` in Julia (factors of τ₀) |
| `err = dev / √n` | `ci_lower`, `ci_upper`: χ² intervals from the Greenhall–Riley EDF, 68.3 % by default (`confidence=`) |
| `n` | `neff` |

## Function names

allantools' plain `adev` and `hdev` are the non-overlapping estimators;
SigmaTau's are overlapping (IEEE 1139), so the names shift:

| allantools | SigmaTau |
|---|---|
| `adev`, `hdev` (non-overlapping) | not provided; the overlapping forms have more degrees of freedom |
| `oadev` | [`adev`](@ref) |
| `mdev`, `tdev` | [`mdev`](@ref), [`tdev`](@ref) |
| `ohdev` | [`hdev`](@ref) |
| `totdev`, `mtotdev`, `ttotdev`, `htotdev` | the same names; see the bias note below |
| `theo1` | [`theo1`](@ref); see the τ note below |
| `mtie`, `tierms` | [`mtie`](@ref), [`tierms`](@ref) |
| `mask.py` (ITU-T TDEV/MTIE masks) | [`telecom_mask`](@ref), [`mask_check`](@ref) |
| `gradev` | [`adev`](@ref): `NaN` gaps are skipped, as `gradev` does |
| `gcodev` | [`groslambert`](@ref), [`cross_adev`](@ref) |
| `three_cornered_hat_phase` | [`ncornered_hat`](@ref) (any number of clocks ≥ 3) |
| `noise.white`, `noise.brown`, `noise_kasdin.Noise` | [`noise_gen`](@ref), any power law and mixture |
| `autocorr_noise_id` | [`identify_noise`](@ref), run inside every estimator |
| `edf_greenhall`, `confidence_interval` | automatic with `ci=true`; [`calculate_edf`](@ref) and [`confidence_intervals`](@ref) on their own |
| `realtime.*` | [`StreamingStability`](@ref) (Julia only) |
| `Dataset`, `Plot` | [`stability`](@ref) for several estimators at once; [Plotting and tables](plotting.md) |

SigmaTau also has [`mhdev`](@ref), [`htdev`](@ref), [`mhtotdev`](@ref),
[`theoh`](@ref), [`pdev`](@ref), the dynamic deviations, Welch spectra and
[`powerlaw_fit`](@ref), none of which allantools provides.

## Where the numbers differ

| | allantools | SigmaTau default | To match allantools |
|---|---|---|---|
| TOTDEV, MTOTDEV, TTOTDEV, HTOTDEV | raw estimators | divided by the SP1065 bias factor, so unbiased | `correct_bias=False` |
| Thêo1 | raw Thêo1, reported at τ = m·τ₀ | ThêoBR, reported at τ = 0.75·m·τ₀ (Thêo1's effective τ) | `correct_bias=False`, and compare at 0.75 × allantools' τ |
| Error bars | dev/√n | χ² intervals that follow the noise type | none; the intervals are wider and asymmetric where the EDF is small |

[Validation](../explanation/validation.md) compares the two at every τ of
the fixture with these settings.

## Reading files

allantools reads arrays you load yourself. SigmaTau reads phase and
frequency text files, Stable32 files, TimeLab exports and counter logs
directly, with their τ₀ and time columns:

```python
p = st.read_data("clock.dat")                 # or read_phase / read_frequency
```

[Reading and saving data](../tutorials/02_reading_your_data.md) covers the
formats, scaling to seconds, detrending and gap filling.

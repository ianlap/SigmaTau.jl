# Coming from allantools

[allantools](https://github.com/aewallin/allantools) and this package
compute the same estimators from the same definitions, and on the same
record their uncorrected deviations agree to eight digits or more
([Related software](../explanation/related.md)). They differ in what a call
takes and returns, in a few names, and in what is applied to a result by
default. This page maps allantools 2024.6 onto this package, in Python and
in Julia.

## The call

An allantools function takes an array, a sample rate and averaging times in
seconds, and returns four arrays. Here a function takes a record, which
carries its sample interval, and averaging factors ``m``, and returns one
result:

```python
import allantools
taus, dev, err, n = allantools.oadev(x, rate=10.0, data_type="phase", taus=[1.0, 10.0])

import sigmatau as st
r = st.adev(st.PhaseData(x, 0.1), [10, 100])
r.tau, r.dev, r.noise_type, r.ci_lower, r.ci_upper, r.edf, r.neff
```

```julia
using SigmaTau
r = adev(PhaseData(x, 0.1), [10, 100])
```

| allantools | SigmaTau |
|:---|:---|
| `data_type="phase"`, `"freq"` | `PhaseData(x, tau0)`, `FrequencyData(y, tau0)` |
| `rate`, in Hz | `tau0 = 1/rate`, in seconds |
| `taus=[1.0, 10.0]`, in seconds | the averaging factors, `[10, 100]` |
| `taus="octave"`, the default | `Octave`, the default |
| `taus="all"` | `AllTaus` |
| `taus="decade"`, which is 1, 2, 4, 10, 20, 40, … | the factors themselves. `Decade` is 1, 10, 100, … |
| `taus="log10"`, about ten values per decade | `QuarterOctave`, about thirteen, or the factors themselves |
| the `taus` returned | `r.tau` |
| the `ns` returned | `r.neff` |
| the `errs` returned, ``σ/\sqrt{n}`` | `r.ci`, a χ² interval ([below](@ref allantools-differences)) |

allantools rounds each requested time to a multiple of ``τ_0`` and leaves
out the ones the record cannot support. This package keeps a row for every
factor asked for, with `NaN` where the record is too short.

## Names

The plain `adev` and `hdev` of allantools are the non-overlapping
estimators. Here the plain names are the overlapping ones, which the
standard recommends over them for their better confidence
([IEEE Std 1139-2022, A.3.2](@cite ieee1139-2022-definitions)).

| allantools | SigmaTau |
|:---|:---|
| `oadev`, `ohdev` | [`adev`](@ref), [`hdev`](@ref) |
| `adev`, `hdev` | `adev(p; overlapping = false)`, `hdev(p; overlapping = false)` |
| `mdev`, `tdev`, `pdev`, `mtie`, `tierms` | the same names |
| `totdev`, `mtotdev`, `ttotdev`, `htotdev` | the same names, bias-corrected by default |
| `theo1` | [`theo1`](@ref), corrected by default and at another ``τ`` |
| `mtie_phase_fast` | [`mtie`](@ref) |
| `gradev` | [`adev`](@ref), which skips the terms that a `NaN` enters |
| `gcodev(a, b)` | [`cross_adev`](@ref) |
| `three_cornered_hat_phase` | [`ncornered_hat`](@ref) on three results, or [`groslambert`](@ref) on three records |
| `frequency2phase`, `phase2frequency` | `PhaseData(f)`, `FrequencyData(p)`; `f.to_phase()`, `p.to_frequency()` in Python |
| `noise.white`, `brown`, `violet`, `pink`; `noise_kasdin.Noise` | [`noise_gen`](@ref) |
| `autocorr_noise_id` | [`identify_noise`](@ref) |
| `edf_greenhall`, `edf_totdev`, `edf_mtotdev` | [`calculate_edf`](@ref) |
| `confidence_interval` | [`confidence_intervals`](@ref) |
| `confidence_interval_noiseID` | the estimator, which identifies the noise and computes the EDF and the interval in one call |
| `mask.prc_tdev`, `mask.prc_mtie` | `telecom_mask(:G811, :tdev, τ)`, `telecom_mask(:G811, :mtie, τ)`; [`mask_check`](@ref) |
| `realtime.oadev_realtime`, `ohdev_realtime`, `tdev_realtime` | [`StreamingStability`](@ref) |
| `Dataset` | [`stability`](@ref) for several estimators, [`save`](@ref) for the file |
| `Plot` | the plot recipes ([Plotting and tables](plotting.md)) |

Three of these return more here. `cross_adev` keeps the sign of the
covariance, of which `gcodev` returns the square root of the magnitude.
`ncornered_hat` returns all the clocks, with an interval for each and `NaN`
where an estimate is negative. `three_cornered_hat_phase` returns the first
clock and sets a negative variance to zero.

## [Where the numbers differ](@id allantools-differences)

| | allantools | SigmaTau by default | To reproduce allantools |
|:---|:---|:---|:---|
| `totdev`, `mtotdev`, `ttotdev`, `htotdev` | uncorrected | divided by the bias factor of the noise type identified ([Bias correction](../explanation/bias.md)) | `correct_bias = false` |
| `theo1` | uncorrected, at ``τ = mτ_0`` | bias removed by ThêoBR, at ``τ = 0.75\,mτ_0`` | `correct_bias = false`, and 0.75 times the ``τ`` of allantools |
| `mtie`, `tierms` of a frequency record | the mean frequency is subtracted before the record is integrated | the record is integrated as given, so a frequency offset counts | subtract the mean of the samples first |
| Error estimate | ``σ/\sqrt{n}`` | χ² interval from the EDF of the noise type identified ([Confidence intervals](../explanation/confidence.md)) | `r.dev ./ sqrt.(r.neff)` |
| EDF of the total deviations | `edf_totdev`, `edf_mtotdev` | [`calculate_edf`](@ref), which takes the modified total deviation from another table | not offered; the [cross-check report](../explanation/crosscheck.md) lists the cases |

With `correct_bias = false` the deviations are those of allantools.
[Validation](../explanation/validation.md) compares the two at every
averaging time of a reference record, and a job of the test suite repeats
the comparison with allantools 2024.6 on every change.

## Not in this package

- `psd2allan`, which integrates a tabulated spectrum to an Allan or modified
  Allan deviation. [`powerlaw_sigma`](@ref) gives the deviations of a
  power-law model.
- The masks of G.8272 and G.8272.1 for primary reference time clocks
  (`prtcA_*`, `prtcB_*`, `eprtc_*`). [`TELECOM_MASKS`](@ref) lists the
  masks here.
- `phase2radians` and `frequency2fractional`. [`read_data`](@ref) converts
  counter readings in hertz to fractional frequency as it reads them.

## Not in allantools

The modified Hadamard and Hadamard time deviations ([`mhdev`](@ref),
[`htdev`](@ref)), the modified Hadamard total deviation
([`mhtotdev`](@ref)), ThêoH ([`theoh`](@ref)), the dynamic deviations
([`dadev`](@ref), [`dhdev`](@ref)), bias correction of the total
deviations, drift fits, readers for instrument files, and the clock model
with its Kalman filter and time scales. The
[home page](../index.md) has the whole list.

## Reading a file

allantools takes arrays that you load yourself. This package reads text
columns, Stable32 and TimeLab files and counter logs, with their sample
interval and their gaps:

```python
p = st.read_data("clock.dat")
```

[Reading and saving data](../tutorials/02_reading_your_data.md) covers the
formats, and [Python](python.md) the rest of the Python interface.

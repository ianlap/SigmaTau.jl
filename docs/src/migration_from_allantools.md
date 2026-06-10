# Migrating from allantools

If you come from [allantools](https://github.com/aewallin/allantools), the
statistics are the same — phase or frequency input, the overlapping
estimators, the power-law noise model, an octave-spaced τ grid by default.
SigmaTau is numerically cross-validated against allantools (see
[Validation](validation/methodology.md)). What changes is the calling
pattern and the result container: instead of passing `data=` + `rate=` +
`data_type=` on every call, you wrap the samples once in a typed record;
instead of a `(taus, devs, errs, ns)` tuple, you get a single
`StabilityResult` that also carries per-τ noise identification, χ²
confidence bounds, and equivalent degrees of freedom.

## Deviation name mapping

| allantools function | SigmaTau function | Notes |
|---------------------|-------------------|-------|
| `adev`              | `adev`            | SigmaTau's `adev` is always overlapping — it is allantools' `oadev`. There is no non-overlapping call. |
| `oadev`             | `adev`            | one-to-one |
| `mdev`              | `mdev`            | |
| `tdev`              | `tdev`            | |
| `hdev`              | `hdev`            | always overlapping — it is allantools' `ohdev` |
| `ohdev`             | `hdev`            | one-to-one |
| `totdev`            | `totdev`          | SigmaTau bias-corrects by default; see [Differences that change numbers](@ref) |
| `mtotdev`           | `mtotdev`         | same |
| `ttotdev`           | `ttotdev`         | same |
| `htotdev`           | `htotdev`         | same |
| `mtie`              | `mtie`            | |
| `pdev`              | `pdev`            | |
| `theo1`             | `theo1`           | same eq. 30 values at equal `m` (cross-validated to machine precision), but allantools reports them at `τ = m·τ0` while SigmaTau reports the SP1065/Stable32 effective `τ = 0.75·m·τ0`; SigmaTau also applies ThêoBR bias removal by default — pass `correct_bias=false` to match allantools' raw values. See [Theory: Allan family](theory/allan_family.md) |
| `tierms`            | `tierms`          | same SP1065 §5.2.18 eq. 37 statistic; cross-validated to machine precision |
| `gradev`            | — none            | no gap-robust ADEV; the SigmaTau route is `fillgaps` (Howe–Schlossberger imputation), then `adev` |
| `three_cornered_hat_phase` | `nch` | operates on pairwise `StabilityResult`s rather than raw phase, and generalizes to N ≥ 3 clocks; see the [three-cornered-hat tutorial](tutorials/06_three_cornered_hat.md) |
| `gcodev`            | — none            | no Groslambert codeviation |
| `noise.*`, `noise_kasdin` | `noise_gen` | power-law noise synthesis from `sigma1` or `h` coefficient targets |
| `autocorr_noise_id` | `identify_noise`  | runs automatically inside every `ci=true` deviation call; see below |

Four SigmaTau deviations have **no allantools counterpart**:

- `theoh` — the SP1065 ThêoH composite (overlapping ADEV below the 20 %-of-T
  crossover, bias-removed ThêoBR above); allantools has no ThêoBR/ThêoH
  helper. See [Theory: Allan family](theory/allan_family.md).
- `mhdev` — modified Hadamard deviation, a phase-averaged third difference:
  drift-insensitive like HDEV while still separating white-PM from
  flicker-PM like MDEV. See [Theory: Allan family](theory/allan_family.md).
- `htdev` — Hadamard time deviation, the σ_x form of `mhdev` (scaled by
  τ/√(10/3)); it is to MHDEV what TDEV is to MDEV. Original to this
  package. Same theory page.
- `mhtotdev` — modified Hadamard total deviation, the boundary-extended
  (total) form of `mhdev` for long-τ estimates. Original to this
  package; its bias and EDF have no published model and are measured by
  Monte Carlo, documented in
  [MHTOTDEV bias and EDF](theory/mhtotdev_bias_edf.md).

## Calling-convention mapping

allantools passes the sampling metadata into every call; SigmaTau attaches
it to the data once:

| allantools argument            | SigmaTau equivalent |
|--------------------------------|---------------------|
| `data=x, data_type="phase"`    | `PhaseData(x, tau0)` |
| `data=y, data_type="freq"`     | `FrequencyData(y, tau0)` |
| `rate=r` (sampling rate, Hz)   | `tau0 = 1/r` (sample interval, seconds) |
| `taus="octave"` (the default)  | `Octave` (also the default) |
| `taus="all"`                   | `AllTaus` |
| `taus="decade"`                | `Decade` |
| `taus=array` (seconds)         | `m_values::Vector{Int}` (averaging factors, τ = m·τ₀) |

Two unit conversions to watch:

- `rate` is a frequency in Hz; `tau0` is the sample interval in seconds, so
  `tau0 = 1/rate`. Data taken at 10 Hz is `PhaseData(x, 0.1)`.
- An explicit allantools `taus` array is in **seconds**; an explicit
  SigmaTau grid is integer **averaging factors** `m`, with τ = m·τ₀.
  `tau_values(Octave, length(p.x), :adev)` materializes a `TauMode` grid as
  that integer vector if you want to inspect or edit it.

`Octave`, `AllTaus`, and `Decade` are instances of the `TauMode` selector,
passed in place of the explicit vector — `adev(p, Decade)`. SigmaTau adds
`HalfOctave`, `QuarterOctave`, and `HalfDecade` spacings, which have no
allantools string equivalent.

### Return values

Every allantools deviation returns a tuple `(taus, ad, ade, ns)`. Every
SigmaTau deviation returns a [`StabilityResult`](reference/types.md):

| allantools tuple element | `StabilityResult` field |
|--------------------------|-------------------------|
| `taus`                   | `r.tau` (seconds)       |
| `ad`                     | `r.dev`                 |
| `ade` (= `ad/√ns`)       | `r.ci` — χ²-based absolute `(lo, hi)` bounds per τ, not a ± half-width; `ci_lower(r)` / `ci_upper(r)` give plain vectors |
| `ns`                     | `r.neff` — the number of analysis windows; `r.edf` separately carries the equivalent degrees of freedom |

The result also carries `r.noise_type` (the per-τ power-law identification,
as SP1065 symbols `:WHPM`, `:FLPM`, `:WHFM`, `:FLFM`, `:RWFM`) and
`r.deviation_type` (which estimator produced it). The `noise_type`, `ci`,
and `edf` vectors are populated when `ci=true` (the default) and empty when
`ci=false`; `r.neff` is populated either way. `mtie` and `tierms` are the
exceptions: they have no published CI model (`mtie` is a deterministic
envelope), so they return `noise_type`, `ci`, and `edf` empty even when
`ci=true`.

## The same analysis in both libraries

allantools:

```python
import numpy as np
import allantools

x = np.loadtxt("clock.dat")          # phase residuals in seconds, 1 Hz
taus, devs, errs, ns = allantools.oadev(
    x, rate=1.0, data_type="phase", taus="octave")
mtaus, mdevs, merrs, mns = allantools.mdev(
    x, rate=1.0, data_type="phase", taus="octave")
```

SigmaTau:

```julia
using SigmaTau

p = read_phase("clock.dat"; time_col=0, value_col=1, tau0=1.0)

r  = adev(p)        # octave grid by default, like taus="octave"
rm = mdev(p)

r.tau               # τ in seconds            (allantools taus)
r.dev               # σ_y(τ)                  (allantools ad)
ci_lower(r)         # χ² lower bounds, 68.3 % (no allantools analog in the tuple)
ci_upper(r)         # χ² upper bounds — or per-τ tuples via r.ci[i].lo / r.ci[i].hi
r.neff              # number of analysis windows (allantools ns)
r.noise_type        # identified power-law type at each τ
r.edf               # equivalent degrees of freedom
```

An explicit grid is `adev(p, [1, 2, 4, 8, 16])`. Where in allantools you
would loop over deviation functions (or `Dataset.compute` one statistic at a
time), `stability(p; devs=(:adev, :mdev, :hdev, :tdev))` computes several in
one call and returns a `StabilitySuite` indexed by symbol
(`suite[:adev]`).

## Differences that change numbers

For the same input and the same τ grid, the σ values agree to the
documented validation tolerance. Three things legitimately differ:

**Error bars.** allantools' deviation functions return the simple estimate
`ade = ad/√ns`, where `ns` is the number of analysis windows (SigmaTau's
`r.neff`). SigmaTau instead identifies the noise type at each τ, computes
the Greenhall–Riley equivalent degrees of freedom
[greenhall-2003-edf-stability](@cite), and maps the deviation through the
χ² distribution into asymmetric `(lo, hi)` bounds in `r.ci` at the 68.3 %
level (`DEFAULT_CONFIDENCE = 0.683`;
override per call with `confidence=0.95`). The two are different statistics
and will not match numerically — allantools can produce the χ² interval
too, but only through separate manual calls to `edf_greenhall` and
`confidence_interval`. See [Theory: Confidence intervals](theory/confidence.md).

**Bias correction on the total family.** `totdev`, `mtotdev`, `ttotdev`,
`htotdev`, and `mhtotdev` apply noise-type-dependent unbias corrections by
default (`correct_bias=true`), because the raw total estimators are biased
relative to their parent statistics — MTOT reads high by roughly 3 – 14 %
in σ depending on noise type, HTOT reads low. That percentage is the
estimator's own bias against the parent statistic, not a kernel
disagreement between libraries: with `correct_bias=false` SigmaTau
reproduces allantools' totals to machine precision on the validation
fixture. allantools always returns the raw estimators (its source marks
the corrections as to-do); Stable32's policy is per-estimator — it
corrects TOTDEV and HTOTDEV but reports MTOTDEV raw, so
`correct_bias=false` matches Stable32 only for MTOTDEV. See
[Theory: Total family](theory/total_family.md) and
[Validation: Stable32](validation/stable32.md).

**Noise identification is automatic and feeds everything.** Both packages
use lag-1-autocorrelation noise identification
[riley-2004-lag1-acf](@cite); allantools exposes it as the separate
`autocorr_noise_id` function, while SigmaTau runs it per τ inside every
`ci=true` call except `mtie`/`tierms` (with a B1-ratio fallback) and reports the
result in `r.noise_type`. Because the EDF — and, for the total family, the bias
factor — depend on the identified α, a different identification at some τ
moves the confidence interval (and the corrected σ) at that τ. See
[Theory: Noise identification](theory/noise_id.md).

## Performance

On identical inputs the measured speedup over allantools ranges from a few
times (`totdev`) through roughly 20× (`adev`) to several thousand times
(`mtotdev` / `htotdev`); the timings and method are on the
[Performance](performance.md) page.

## See also

- [Getting Started](getting_started.md) — installation and a first run.
- [Migrating from Stable32](migration_from_stable32.md) — the same guide
  for Stable32 users.
- [Validation Methodology](validation/methodology.md) — how the
  cross-validation against allantools and Stable32 is done.

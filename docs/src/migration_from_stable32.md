# Migrating from Stable32

If you already think in Stable32 terms, this page maps that workflow onto
SigmaTau.jl. The vocabulary is the same — phase/frequency data, averaging
factors, the Allan/Hadamard/Total families, χ² confidence intervals, lag-1
noise identification — only the interface changes from a desktop GUI to a Julia
API you can script and embed.

SigmaTau is numerically cross-validated against Stable32 (see
[Validation](validation/methodology.md)); for the same input and the same τ
grid you should get the same σ to within the documented tolerance.

## Deviation name mapping

Stable32's run-types map to SigmaTau functions one-to-one:

| Stable32 run / column        | SigmaTau function |
|------------------------------|-------------------|
| Allan (ADEV)                 | `adev`            |
| Overlapping Allan (OADEV)    | `adev` (always overlapping) |
| Modified Allan (MDEV)        | `mdev`            |
| Time (TDEV)                  | `tdev`            |
| Hadamard (HDEV)              | `hdev`            |
| Overlapping Hadamard (OHDEV) | `hdev` (always overlapping) |
| Total (TOTDEV)               | `totdev`          |
| Modified Total (MTOT)        | `mtotdev`         |
| Time Total (TTOT)            | `ttotdev`         |
| Hadamard Total (HTOT)        | `htotdev`         |
| MTIE                         | `mtie`            |
| Thêo1 / ThêoBR               | `theo1` (`correct_bias=false` is the raw Thêo1 run; the default `correct_bias=true` is ThêoBR) |
| ThêoH                        | `theoh`           |

SigmaTau's Allan and Hadamard deviations are always the overlapping estimators
(the Stable32 "O…" variants), which is the standard choice for stability
analysis. There is no separate non-overlapping call.

Four SigmaTau deviations have **no Stable32 run type**:

- `mhdev` — the modified Hadamard deviation (phase-averaged third
  difference); see [Theory: Allan family](theory/allan_family.md).
- `htdev` — the Hadamard time deviation, the σ_x form of `mhdev`;
  original to this package. Same theory page.
- `mhtotdev` — the modified Hadamard total deviation, the
  boundary-extended form of `mhdev`; original to this package, with
  Monte-Carlo-calibrated bias and EDF
  ([Theory: MHTOTDEV bias and EDF](theory/mhtotdev_bias_edf.md)).
- `pdev` — the parabolic deviation (PVAR); see
  [Theory: Allan family](theory/allan_family.md).

## Result-field mapping

A single call returns a [`StabilityResult`](reference/types.md). Its fields line
up with the columns Stable32 prints in a deviation run:

| Stable32 column     | `StabilityResult` field |
|---------------------|-------------------------|
| Tau                 | `r.tau` (seconds)       |
| Sigma / Dev         | `r.dev`                 |
| Min Sigma / Max Sigma (CI) | `r.ci` — one `(lo, hi)` tuple per τ (`r.ci[1].lo`); `ci_lower(r)` / `ci_upper(r)` give plain vectors |
| # (analysis points) | `r.neff`                |
| EDF                 | `r.edf`                 |
| Noise type (α)      | `r.noise_type`          |

Noise-type symbols follow the SP1065 power-law labels: `:WHPM`, `:FLPM`,
`:WHFM`, `:FLFM`, `:RWFM` (α = +2, +1, 0, −1, −2). The confidence and EDF fields
are populated when `ci=true` (the default) and left empty when
`ci=false`; `r.neff` is populated either way.

## A complete session

### 1. Load a data file

Stable32 reads a column of phase (or frequency) values. So does SigmaTau, via
`read_phase` / `read_frequency`. A typical `.DAT` file is a single column of
phase values; give the sample interval τ₀ explicitly:

```julia
using SigmaTau

# Single column of phase residuals (seconds), τ₀ = 1 s.
# time_col=0 says "no time column"; value_col=1 reads the lone column.
p = read_phase("clock.DAT"; time_col=0, value_col=1, tau0=1.0)
```

If the file has a time column and a value column (the common two-column
layout), point at the value column — this is the default, so `tau0` is inferred
from the time column if you omit it:

```julia
p = read_phase("clock.DAT"; time_col=1, value_col=2)
```

`read_frequency` is identical but returns a `FrequencyData`. Both accept
`scaling` (multiply the values, e.g. to convert units), `detrend`
(`:none`, `:mean`, `:endpoint`, `:linear`), and `fillgaps=true` for the
Howe–Schlossberger gap imputation — the preprocessing you'd otherwise reach for
in Stable32's data menu, done in the same call.

You can also build data in memory without a file:

```julia
p = PhaseData(phase_vector, 1.0)        # τ₀ = 1 s
f = FrequencyData(freq_vector, 1.0)
```

### 2. Run a deviation

A Stable32 "Run" on a single deviation is a single function call. Pass the
averaging factors `m` (τ = m·τ₀):

```julia
r = adev(p, [1, 2, 4, 8, 16, 32, 64])
```

Or let SigmaTau pick the grid for you. The default is an octave-spaced grid (the
Stable32 default); other spacings are available through the `TauMode` selector:

```julia
r = adev(p)                 # octave-spaced, like Stable32's default
r = adev(p, AllTaus)        # every averaging factor
r = adev(p, Decade)         # decade-spaced
```

Confidence intervals are on by default (`ci=true`) and use the
Greenhall–Riley EDF with a χ² mapping — the same machinery behind Stable32's
error bars. Turn them off for a faster bare-σ run with `ci=false`.

### 3. Run a whole suite at once

Stable32 runs one deviation per pass. SigmaTau can compute several in one call
with `stability`, returning a `StabilitySuite` you index by deviation symbol:

```julia
suite = stability(p; devs=(:adev, :mdev, :hdev, :tdev))

suite[:adev].dev      # the ADEV curve
suite[:hdev].tau      # its τ grid
keys(suite)           # which deviations are present, in order
```

`devs` defaults to `(:adev, :mdev, :hdev, :tdev)`. Detrend and bias options
passed to `stability` are forwarded only to the total family.

### 4. Save and plot

Where Stable32 writes a results table and a plot, SigmaTau round-trips results
to a self-describing text file and plots through any `Plots`-compatible backend:

```julia
using Plots

save_result("adev.tsv", r)          # single result
suite_path = save_suite("run.tsv", suite)   # a whole suite, with session metadata

plot(r)                             # log-log τ–σ with CI error bars
plot(suite)                         # overlay every deviation in the suite
```

`load_result` / `load_suite` read them back. The suite file also records the
package version, a timestamp, τ₀, N, the confidence level, and the τ-mode, so a
saved run is self-documenting.

## Quick API reference

| Task                          | Stable32             | SigmaTau                          |
|-------------------------------|----------------------|-----------------------------------|
| Load phase file               | File → Open          | `read_phase(path; tau0=…)`        |
| Load frequency file           | File → Open          | `read_frequency(path; tau0=…)`    |
| Detrend / remove slope        | Edit menu            | `detrend(data; method=:linear)` or `detrend=` kwarg on read |
| Fill gaps                     | Edit menu            | `fillgaps(data)` or `fillgaps=true` |
| Overlapping Allan             | OADEV run            | `adev(data, m)`                   |
| Confidence intervals          | Sigma error bars     | `ci=true` (default)          |
| Noise identification          | automatic            | `r.noise_type` (per τ)            |
| Run several deviations        | one pass each        | `stability(data; devs=…)`         |
| Export results table          | Save results         | `save_result` / `save_suite`      |
| Plot                          | Plot window          | `plot(r)` / `plot(suite)`         |

## See also

- [Getting Started](getting_started.md) — installation and a first run.
- [Tutorial 0: Julia for metrologists](tutorials/00_julia_for_metrologists.md) —
  a longer walkthrough aimed at Stable32 users new to Julia.
- [Performance](performance.md) — how the runtime compares to allantools.

# SigmaTau.jl

[![CI](https://github.com/ianlap/SigmaTau.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ianlap/SigmaTau.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://ianlap.github.io/SigmaTau.jl/stable/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Time-and-frequency stability analysis in Julia. One call gives the deviation,
the identified noise type, the equivalent degrees of freedom and the χ²
confidence interval at every τ. The estimators agree with Stable32 to its
printed precision and with allantools to machine precision, and the heavy
total-family kernels run thousands of times faster than allantools.

```julia
using SigmaTau

p = read_phase("clock.DAT"; time_col=0, value_col=1, tau0=1.0)

r = adev(p)                 # octave τ grid, overlapping estimator
r.tau, r.dev                # σ_y(τ)
r.noise_type                # :WHPM / :FLPM / :WHFM / :FLFM / :RWFM per τ
r.ci[1].lo, r.ci[1].hi      # 68.3 % χ² bounds (Greenhall–Riley EDF)
r.neff, r.edf

suite = stability(p; devs=(:adev, :mdev, :hdev, :mhtotdev))
suite[:hdev].dev

using Plots
plot(suite)                 # log-log σ(τ) with CI error bars
```

## Install

```julia
pkg> add SigmaTau
```

Julia 1.11 or later. `Plots` is optional and enables the plot recipes.

Each [release](https://github.com/ianlap/SigmaTau.jl/releases) also carries
prebuilt downloads that need no Julia install: the desktop explorer (macOS
`.dmg` for Apple silicon and Intel, Linux `.AppImage` and `.deb`) and the
`sigmatau` command line as a `.tar.gz` per platform (unpack it and run
`bin/sigmatau`).

## Estimators

| Family | Functions |
|---|---|
| Allan | `adev`, `mdev`, `tdev` |
| Hadamard | `hdev`, `mhdev`, `htdev` |
| Total | `totdev`, `mtotdev`, `ttotdev`, `htotdev`, `mhtotdev` |
| Thêo | `theo1` (ThêoBR by default), `theoh` |
| Time error | `mtie`, `tierms`, `pdev` |
| Time-resolved | `dadev`, `dhdev` |
| Real time | `StreamingStability`, `LiveStability` (`adev`, `mdev`, `hdev`, `mhdev` from an instrument) |
| Multi-clock | `ncornered_hat`, `groslambert` |
| Two channels | `cross_adev`, `cross_psd_y`, `cross_psd_x` |
| Noise and drift models | `identify_noise`, `noise_id_nnls`, `powerlaw_fit` (`powerlaw_sigma`, `powerlaw_psd`), `drift_fit` |
| Spectra and synthesis | `psd_y`, `psd_x`, `ssb_phase_noise`, `noise_gen` |

Every function accepts `PhaseData` or `FrequencyData`, an explicit vector of
averaging factors or a `TauMode` grid (`Octave`, `Decade`, `AllTaus`, ...).
`read_phase` / `read_frequency` read plain columns and `read_data` the
instrument formats (TimeLab, Stable32, counter logs); `detrend`, `fillgaps`
and `remove_outliers` cover the usual preprocessing; `save` / `load` round-trip records, results
and suites as plain text. `estimators()` lists the estimators and
`deviation(:adev, data, m)` runs one chosen by name. The
[documentation](https://ianlap.github.io/SigmaTau.jl/stable/) has a
[getting-started tutorial](https://ianlap.github.io/SigmaTau.jl/stable/getting_started/),
a guide to
[choosing an estimator](https://ianlap.github.io/SigmaTau.jl/stable/guide/choosing/)
and the full API reference.

## Python

`python/` holds the `sigmatau` Python package, a juliacall wrapper with the
same names (`st.adev(pd)`, `st.psd_y(pd)`, `st.deviation("mhtotdev", pd)`).
Install it with `pip install sigmatau`; see
[`python/README.md`](python/README.md).

## Compared with Stable32 and allantools

| | Stable32 | allantools | SigmaTau |
|---|---|---|---|
| Platform | Windows GUI, closed source | Python library | Julia library + desktop GUI |
| χ² confidence intervals with per-τ noise ID | yes | manual (`edf_greenhall` + `confidence_interval`) | every call |
| Total-family bias correction | TOTDEV, HTOTDEV | none | all five, switchable |
| Modified Hadamard Class (`mhdev`, `htdev`, `mhtotdev`) | no | no | yes |
| Thêo1 / ThêoBR / ThêoH | yes | Thêo1 only | yes |
| Dynamic (time-resolved) deviation | yes | no | yes |
| Real-time accumulators | no | ADEV, HDEV, TDEV | ADEV, MDEV, HDEV, MHDEV |
| N-cornered hat | three clocks | three clocks | N ≥ 3 |
| Multithreaded kernels | n/a | no | yes |

Validation: the core estimators match Stable32 to ≤ 0.005 % (its
five-significant-figure output) and allantools to ≤ 1e-13 relative on the
shipped fixture; raw total kernels match allantools to ≤ 4.4e-15. The
cross-check suite in `test/validation/` repeats this for every estimator on
records from white PM to random-run FM, plus EDF, intervals, noise ID and
spectra (about 7,000 comparisons;
[`tools/validation/REPORT.md`](tools/validation/REPORT.md)). Where a
default differs from a reference (MTOTDEV bias correction, noise-ID
detrending), the per-τ tables in the
[validation page](https://ianlap.github.io/SigmaTau.jl/stable/explanation/validation/)
name the switch that reproduces the reference.

Performance: on identical inputs, `adev`/`mdev`/`hdev`/`tdev` run roughly 10
to 40 times faster than allantools, `totdev` about 70 times, and
`mtotdev`/`htotdev` 1,600 times faster on 4 threads to 4,200 times faster on 12
(0.46 s vs 30 min on a 407,000-sample record). The harness, method and both measurements are in
[`bench/`](bench/README.md).

## Command line

`sigmatau` is a compiled command line over the same code (about 30 ms to
start, no Julia needed to run it):

```bash
sigmatau analyze clock.dat --devs adev,mdev,tdev            # CSV on stdout
sigmatau analyze clock.dat --format json --out result.json  # versioned JSON
sigmatau check clock.dat                                    # samples, gaps, outliers, drift
sigmatau spectrum clock.dat                                 # S_y(f), S_x(f)
sigmatau dynamic clock.dat --dev dhdev                      # sliding-window HDEV
sigmatau drift clock.dat --model log                        # aging model and rate per day
sigmatau cross a.dat --with b.dat                           # two channels, one device
sigmatau live tcp://192.168.1.20:5025 --hz                  # running σ(τ) from a counter
sigmatau report clock.dat --out clock.html                  # self-contained HTML report
```

Build it from `cli/` with JuliaC; see [the docs page](docs/src/guide/cli.md).

## Explorer

`app/web/` is a point-and-click explorer for the whole workflow. You can load a
file (Stable32 `.DAT`, CSV, TSV, whitespace columns or an instrument format,
several clocks per file), synthesize a power-law record, or stream readings
from an instrument on the Live page. It checks the record for gaps, outliers
and drift, lets you pick the range to analyse, then plots σ_y and σ_x with χ²
intervals and the fitted power-law model, the S_y and S_x spectra and the
dynamic σ(t, τ) map; pinned runs stay on the plot for comparison. The File
menu writes a self-contained HTML report and saves or restores a session. The
pages are plain HTML and JavaScript over the compiled `sigmatau` command
line, which answers in tens of milliseconds.

```bash
python3 app/web/serve.py --open          # http://127.0.0.1:9385
```

[`app/web/README.md`](app/web/README.md) shows how to build the core, and
[`app/desktop/README.md`](app/desktop/README.md) wraps the same pages in a
Tauri desktop window.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Citing

```bibtex
@software{lapinski-sigmatau,
  author  = {Lapinski, Ian},
  title   = {{SigmaTau.jl}: time-and-frequency stability analysis in {Julia}},
  url     = {https://github.com/ianlap/SigmaTau.jl},
  version = {0.6.0},
  year    = {2026},
}
```

Definitions follow NIST SP1065 (Riley & Howe), IEEE Std 1139-2022 and
Greenhall & Riley 2003. MIT license.

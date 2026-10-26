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
pkg> add https://github.com/ianlap/SigmaTau.jl
```

Julia 1.11 or later. Optional: `Plots` for the plot recipes, `Pluto` for the
interactive explorer in `app/`.

## Estimators

| Family | Functions |
|---|---|
| Allan | `adev`, `mdev`, `tdev` |
| Hadamard | `hdev`, `mhdev`, `htdev` |
| Total | `totdev`, `mtotdev`, `ttotdev`, `htotdev`, `mhtotdev` |
| Thêo | `theo1` (ThêoBR by default), `theoh` |
| Time error | `mtie`, `tierms`, `pdev` |
| Time-resolved | `dadev`, `dhdev` |
| Real time | `StreamingStability` (`adev`, `mdev`, `hdev`, `mhdev`) |
| Multi-clock | `nch` (N-cornered hat) |
| Spectra and synthesis | `Sy`, `Sx`, `L`, `noise_gen` |

Every function accepts `PhaseData` or `FrequencyData`, an explicit vector of
averaging factors or a `TauMode` grid (`Octave`, `Decade`, `AllTaus`, ...).
`read_phase` / `read_frequency`, `detrend`, `fillgaps` and `remove_outliers`
cover the usual preprocessing; `save` / `load_result` / `load_suite`
round-trip results as plain text.

## Compared with Stable32 and allantools

| | Stable32 | allantools | SigmaTau |
|---|---|---|---|
| Platform | Windows GUI, closed source | Python library | Julia library + Pluto GUI |
| χ² confidence intervals with per-τ noise ID | yes | manual (`edf_greenhall` + `confidence_interval`) | every call |
| Total-family bias correction | TOTDEV, HTOTDEV | none | all five, switchable |
| MDEV-class Hadamard (`mhdev`, `htdev`, `mhtotdev`) | no | no | yes |
| Thêo1 / ThêoBR / ThêoH | yes | Thêo1 only | yes |
| Dynamic (time-resolved) deviation | yes | no | yes |
| Real-time accumulators | no | ADEV, HDEV, TDEV | ADEV, MDEV, HDEV, MHDEV |
| N-cornered hat | three clocks | three clocks | N ≥ 3 |
| Multithreaded kernels | n/a | no | yes |

Validation: the core estimators match Stable32 to ≤ 0.005 % (its
five-significant-figure output) and allantools to ≤ 1e-13 relative on the
shipped fixture; raw total kernels match allantools to ≤ 4.4e-15. Where a
default differs from a reference (MTOTDEV bias correction, noise-ID
detrending), the per-τ tables in the
[validation page](https://ianlap.github.io/SigmaTau.jl/stable/validation/)
name the switch that reproduces the reference.

Performance: on identical inputs, `adev`/`mdev`/`hdev`/`tdev` run roughly 10
to 40 times faster than allantools, and `mtotdev`/`htotdev` 1,600 times faster
on 4 cores to 4,200 times faster on 12 (0.47 s vs 30 min on a 407,000-sample
record). The harness, method and both measurements are in
[`bench/`](bench/README.md).

## Interactive explorer

`app/StabilityExplorer.jl` is a Pluto notebook: load a file or synthesize a
record, pick deviations, τ grid, confidence level and preprocessing, and get
the plot, the table and a TSV export.

```bash
julia --project=app -e 'using Pkg; Pkg.instantiate(); using Pluto; Pluto.run(notebook="app/StabilityExplorer.jl")'
```

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
  version = {0.5.0},
  year    = {2026},
}
```

Definitions follow NIST SP1065 (Riley & Howe), IEEE Std 1139-2022 and
Greenhall & Riley 2003. MIT license.

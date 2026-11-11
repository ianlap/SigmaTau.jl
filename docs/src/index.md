# SigmaTau.jl

SigmaTau.jl measures the stability of clocks and oscillators from records of
phase or fractional frequency. It computes the Allan and Hadamard deviations
with their modified, time and total forms, Thêo1, ThêoH and the parabolic
deviation, and at each averaging time it identifies the noise type and gives
a χ² confidence interval. Around them it has spectra and phase noise, MTIE
with telecom masks, tools to prepare a record, a clock model with a Kalman
filter, steering loops, and time scales formed from several clocks.

```@example home
using SigmaTau, Random

# 8192 phase samples at 1 s: white FM of 1e-12 at 1 s, and random-walk FM
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 3e-14),
              rng = Xoshiro(5))
adev(p)
```

## What it computes

| | Functions |
|:---|:---|
| Allan family | [`adev`](@ref), [`mdev`](@ref), [`tdev`](@ref) |
| Hadamard family | [`hdev`](@ref), [`mhdev`](@ref), [`htdev`](@ref) |
| Total estimators | [`totdev`](@ref), [`mtotdev`](@ref), [`ttotdev`](@ref), [`htotdev`](@ref), [`mhtotdev`](@ref) |
| Thêo and parabolic | [`theo1`](@ref), [`theoh`](@ref), [`pdev`](@ref) |
| Time interval error | [`mtie`](@ref), [`tierms`](@ref), [`mask_check`](@ref) |
| Dynamic and live stability | [`dadev`](@ref), [`dhdev`](@ref), [`StreamingStability`](@ref), [`LiveStability`](@ref) |
| Several clocks | [`ncornered_hat`](@ref), [`groslambert`](@ref), [`cross_adev`](@ref) |
| Spectra | [`psd_y`](@ref), [`psd_x`](@ref), [`ssb_phase_noise`](@ref), [`integrated_jitter`](@ref), [`cross_psd_y`](@ref) |
| Noise | [`identify_noise`](@ref), [`powerlaw_fit`](@ref), [`noise_gen`](@ref) |
| Preparing a record | [`read_data`](@ref), [`fillgaps`](@ref), [`remove_outliers`](@ref), [`detrend`](@ref), [`drift_fit`](@ref) |
| Clock model and Kalman filter | [`noise_fit`](@ref), [`clock_model`](@ref), [`kalman_filter`](@ref), [`holdover`](@ref) |
| Steering | [`steering_gains`](@ref), [`lqg_gains`](@ref), [`steer`](@ref) |
| Time scales | [`timescale`](@ref), [`simulate_ensemble`](@ref) |

The estimators take a [`PhaseData`](@ref) or [`FrequencyData`](@ref)
record. The [Reference](reference/estimators.md) documents every function.

## Ways to use it

| | For | Start with |
|:---|:---|:---|
| Julia | scripts, notebooks and batch work | [Getting started](getting_started.md) |
| Explorer | reading a file and its curves without writing code | [The explorer](explorer/index.md) |
| Command line | shell scripts and other languages, with CSV or JSON output | [Command line](guide/cli.md) |
| Python | the same functions from Python | [Python](guide/python.md) |

All four run the same code. Python calls the Julia functions themselves.
The compiled program, which the explorer also runs, is tested to give their
results within a relative difference of ``10^{-12}``.

## Install

```julia-repl
pkg> add SigmaTau
```

installs the latest registered release, and
`pkg> add https://github.com/ianlap/SigmaTau.jl` the repository as it
stands. The package needs Julia 1.11 or later. `Plots` is optional: loading
it turns on the plot recipes.

A release carries the [explorer](explorer/install.md) and the
[command-line program](guide/cli.md#Install) as programs that need no
Julia. Their pages say how to install them, and how to run or build them
from the repository.

## This documentation

- **Tutorials** work through an analysis from start to finish: a first
  result, a drifting clock, a file of your own, three clocks compared, an
  ensemble.
- **How-to guides** each answer one question: which estimator to use, how
  to prepare a record, how to plot, how to run the command line or Python.
- **The explorer** is the manual of the graphical program.
- **Background** explains what the numbers rest on: the conventions, the
  noise types and how they are identified, confidence intervals and bias
  factors, and the evidence that the results are right.
- **Reference** documents every exported function and type.

## Sources

The definitions follow the NIST handbook of frequency stability analysis
([Riley, 2008](@cite riley-2008-sp1065)), cited as SP1065 throughout, and
[IEEE Std 1139-2022](@cite ieee1139-2022-definitions). Confidence intervals
use the equivalent degrees of freedom of
[Greenhall and Riley (2003)](@cite greenhall-2003-edf-stability). Each page
cites what it states, with the section or equation, and the
[Bibliography](bibliography.md) lists every source.

Where the package departs from a source, or goes beyond one, the page that
covers the subject says so and says what the choice rests on.
[Validation](explanation/validation.md) describes the evidence for the
results, and [Related software](explanation/related.md) the programs they
were compared with.

## Citing and license

SigmaTau.jl is released under the MIT license. To cite it, use the entry in
[`CITATION.cff`](https://github.com/ianlap/SigmaTau.jl/blob/main/CITATION.cff)
of the repository.

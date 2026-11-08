# SigmaTau.jl

Time-and-frequency stability analysis in Julia. SigmaTau computes the Allan,
modified Allan, Hadamard, modified Hadamard and total deviation families,
Thêo1 and ThêoH, MTIE, TIE rms and PDEV. Every deviation call (except MTIE
and TIE rms, which have no χ² model) also identifies the noise type at each τ
and returns the equivalent degrees of freedom and a χ² confidence interval.
The estimators that Stable32 and allantools also implement agree with Stable32
to its printed precision and with allantools to machine precision; the others
are checked against reference kernels and identities. See
[Validation](explanation/validation.md).

```@example home
using SigmaTau, Random
Random.seed!(42)

# 8192 samples at τ₀ = 1 s: white FM at σ_y(1 s) = 1e-12 plus random-walk FM
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 3e-14))

adev(p)
```

## Ways to use it

| | Use it for | Start here |
|:---|:---|:---|
| **Julia REPL** | scripts, notebooks, custom analysis, batch jobs | [Getting started](getting_started.md) |
| **Explorer** | loading a file and reading the curves without writing code | [Using the explorer](guide/explorer.md) |
| **Command line** | one-off runs, shell scripts and other languages, with CSV or JSON output | [Command line](guide/cli.md) |
| **Python** | the same API from Python, through an embedded Julia | [Python](guide/python.md) |

All of them run the same library code, so a deviation computed in the
explorer, on the command line or from Python equals the one `adev(p)` returns
in the REPL.

## Install

```julia-repl
pkg> add https://github.com/ianlap/SigmaTau.jl
```

Julia 1.11 or later. `Plots` is optional: loading it turns on the plot
recipes. The explorer and the command-line tool are built separately; their
pages say how.

## How these docs are organised

- **Tutorials** teach by example, from a first analysis to reading your own
  files and separating clocks with the N-cornered hat.
- **How-to guides** answer one question each: which estimator to use, how to
  plot and export, how to run the explorer, the
  command-line tool or the Python package.
- **Background** explains the ideas the numbers rest on: noise types, EDF,
  confidence intervals and bias correction, plus the validation against
  Stable32 and allantools and the derivation of the MHTOTDEV coefficients.
- **Reference** lists every exported function and type, grouped by task.

## Coming from Stable32

| Stable32 | SigmaTau |
|:---|:---|
| OADEV, MDEV, TDEV | `adev`, `mdev`, `tdev` (the Allan estimators are always overlapping) |
| OHDEV | `hdev` |
| TOTDEV, MTOT, TTOT, HTOT | `totdev`, `mtotdev`, `ttotdev`, `htotdev` |
| Thêo1 / ThêoBR, ThêoH | `theo1` (`correct_bias=false` gives raw Thêo1), `theoh` |
| MTIE, TIE rms | `mtie`, `tierms` |
| Dynamic stability | `dadev`, `dhdev` |
| Three-cornered hat | `ncornered_hat` |
| Noise generation | `noise_gen` |
| Drift, outlier and gap tools | `detrend`, `remove_outliers`, `fillgaps` |
| PSD and phase noise | `psd_y`, `psd_x`, `ssb_phase_noise` |

Not in Stable32: `mhdev`, `htdev`, `mhtotdev`, `pdev` and the real-time
`StreamingStability` accumulators. Two defaults differ from Stable32
(MTOTDEV bias correction and noise-ID detrending); the
[validation page](explanation/validation.md) names the switch that reproduces
Stable32 in each case.

## References

Definitions follow NIST SP1065 [riley-2008-sp1065](@cite) and
IEEE Std 1139-2022 [ieee1139-2022-definitions](@cite); confidence intervals
use the Greenhall–Riley EDF algorithm [greenhall-2003-edf-stability](@cite).

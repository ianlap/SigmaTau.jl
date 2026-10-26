# SigmaTau.jl

Time-and-frequency stability analysis in Julia: the Allan, modified Allan,
Hadamard, modified Hadamard and total deviation families, Thêo1/ThêoH, MTIE,
TIE rms and PDEV, with per-τ noise identification and χ² confidence
intervals computed in the same call. The core estimators are validated
against Stable32 and allantools; see [Validation](validation.md).

## Install

```julia-repl
pkg> add https://github.com/ianlap/SigmaTau.jl
```

Julia 1.11 or later.

## First run

```@example basic
using SigmaTau, Random
Random.seed!(42)

# White FM at σ_y(1 s) = 1e-12 plus random-walk FM, 8192 samples at τ₀ = 1 s
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 3e-14))

r = adev(p)                      # octave τ grid, χ² CI, per-τ noise ID
```

Every deviation returns a [`StabilityResult`](@ref) with `tau`, `dev`,
`neff`, and, when `ci=true` (the default), `noise_type`, `ci` (one `(lo, hi)`
pair per τ) and `edf`.

```@example basic
suite = stability(p; devs=(:adev, :mdev, :hdev, :mhdev))   # several at once
suite[:mdev]
```

Reading a file, choosing a τ grid and plotting:

```julia
p = read_phase("clock.DAT"; time_col=0, value_col=1, tau0=1.0)
adev(p, [1, 2, 4, 8, 16, 32])    # explicit averaging factors
adev(p, Decade)                  # decade-spaced grid

using Plots
plot(adev(p))                    # log-log σ(τ) with CI error bars
plot(stability(p))               # overlay a suite
```

## Contents

- [Plotting and tables](plotting.md): the plot recipes, confidence-band styles,
  dynamic deviation waterfalls, and table export.
- [Validation](validation.md): per-τ agreement with Stable32 and allantools.
- [MHTOTDEV bias and EDF](mhtotdev.md): how the coefficients of the one
  estimator with no published reference were measured.
- Tutorials: a [drifting clock and the Hadamard family](tutorials/01_drifting_clock_hadamard.md),
  [reading and saving data](tutorials/02_reading_your_data.md), and the
  [three-cornered hat](tutorials/03_three_cornered_hat.md).
- API reference: [types](reference/types.md) and [stability functions](reference/stab.md).

## Stable32 run types

| Stable32 | SigmaTau |
|---|---|
| OADEV, MDEV, TDEV | `adev`, `mdev`, `tdev` (the Allan estimators are always overlapping) |
| OHDEV | `hdev` |
| TOTDEV, MTOT, TTOT, HTOT | `totdev`, `mtotdev`, `ttotdev`, `htotdev` |
| Thêo1 / ThêoBR, ThêoH | `theo1` (`correct_bias=false` gives raw Thêo1), `theoh` |
| MTIE, TIE rms | `mtie`, `tierms` |
| Dynamic stability | `dadev`, `dhdev` |
| Three-cornered hat | `nch` |
| Noise generation | `noise_gen` |
| Drift / outlier / gap tools | `detrend`, `remove_outliers`, `fillgaps` |

Not in Stable32: `mhdev`, `htdev`, `mhtotdev`, `pdev`, `StreamingStability`.

## References

Definitions follow NIST SP1065 [riley-2008-sp1065](@cite) and
IEEE Std 1139-2022 [ieee1139-2022-definitions](@cite); confidence intervals
use the Greenhall–Riley EDF algorithm [greenhall-2003-edf-stability](@cite).

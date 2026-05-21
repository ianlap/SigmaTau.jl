# Getting Started

## Installation

SigmaTau.jl is a single registerable package. From the Pkg REPL:

```julia-repl
pkg> add https://github.com/ianlap/SigmaTau.jl
```

`using SigmaTau` brings the shared types (`PhaseData`, `FrequencyData`,
`StabilityResult`) and every public deviation, noise-ID, and IO function
into scope under one flat namespace.

For clock state-space estimation and Kalman steering, see the sister
package [ClockEnsemble.jl](https://github.com/ianlap/ClockEnsemble.jl).

## A minimal example

```@example basic
using SigmaTau
using Random
Random.seed!(42)

# 10 minutes of phase residuals at 1 Hz, white-PM-like
phase = randn(600)
p = PhaseData(phase, 1.0)

# Compute overlapping Allan deviation at three τ values
r = adev(p, [1, 4, 16]; calc_ci=true)

# Round for stable display
round.(r.dev; sigdigits=4)
```

The result is a [`StabilityResult`](reference/types.md#SigmaTau.StabilityResult)
populated with `tau`, `dev`, per-τ `noise_type`, χ²-based confidence bounds
`ci_lower` / `ci_upper`, and `edf` per Greenhall–Riley.

## Where to next

- [Tutorial 0: Julia for metrologists](tutorials/00_julia_for_metrologists.md) — coming from Stable32? Start here. Covers installation, loading `.DAT` files, your first `adev`, overlaying plots, and saving results to disk.
- [Tutorial 1: Phase data](tutorials/01_phase_data.md) — `PhaseData` / `FrequencyData` basics.
- [Tutorial 2: Computing Allan deviation](tutorials/02_compute_adev.md) — `adev` and the `StabilityResult` it returns (incl. EDF, χ² CI, noise type).
- [Tutorial 6: Three-cornered hat](tutorials/06_three_cornered_hat.md) — separating clock noise from a three-clock comparison.
- [Theory: Overview](theory/overview.md) — what σ_y(τ) means.

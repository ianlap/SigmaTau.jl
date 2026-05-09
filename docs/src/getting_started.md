# Getting Started

## Installation

SigmaTau.jl is a single registerable package. From the Pkg REPL:

```julia-repl
pkg> add https://github.com/ianlap/SigmaTau.jl
```

`using SigmaTau` brings the shared types (`PhaseData`, `FrequencyData`,
`StabilityResult`) and every export of the `SigmaTau.Stab` and
`SigmaTau.Est` submodules into scope.

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

- [Tutorial 1: Phase Data](tutorials/01_phase_data.md) — load and inspect a record.
- [Tutorial 2: Computing Allan Deviation](tutorials/02_compute_adev.md) — `adev`, `mdev`, `hdev`.
- [Tutorial 3: Identifying Noise Type](tutorials/03_identify_noise.md) — lag-1 ACF.
- [Theory: Overview](theory/overview.md) — what σ_y(τ) means.

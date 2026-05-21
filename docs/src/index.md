# SigmaTau.jl

A Julia package for time-and-frequency stability analysis: Allan / Modified
Allan / Hadamard / Total deviations, χ²-based confidence intervals, lag-1
ACF noise identification, MTIE, parabolic deviation, and calibrated
power-law noise generation.

For clock state-space estimation and Kalman steering, see the sister
package [ClockEnsemble.jl](https://github.com/ianlap/ClockEnsemble.jl).

## Package layout

SigmaTau.jl is a single, flat Julia 1.11 package. Shared types,
file IO, deviation kernels, noise identification, and EDF/CI all live
under one namespace:

| Surface | Responsibility |
|---|---|
| [Shared types](reference/types.md) | `PhaseData`, `FrequencyData`, `StabilityResult`. |
| [Stability API](reference/stab.md) | Deviation kernels, EDF / CI, noise identification, MTIE, PDEV, `noise_gen`. |

## Where to next

- **First time?** See [Getting Started](getting_started.md).
- **Need theory?** Start with [Theory: Overview](theory/overview.md).
- **Want a worked example?** Pick a [Tutorial](tutorials/01_phase_data.md).
- **Looking up an API?** Browse the [API Reference](reference/types.md).
- **Verifying numerics?** See [Validation](validation/methodology.md).

## Reference math

Numerical conventions and confidence-interval formulas follow:

- NIST SP1065 [Riley & Howe 2008](@cite riley-2008-sp1065)
- Greenhall–Riley 2003 EDF approximations [Greenhall 2003](@cite greenhall-2003-edf-stability)
- [IEEE Std 1139-2022](@cite ieee1139-2022-definitions)

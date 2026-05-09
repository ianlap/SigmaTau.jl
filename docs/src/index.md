# SigmaTau.jl

A Julia package for time-and-frequency stability analysis: Allan / Modified
Allan / Hadamard / Total deviations, χ²-based confidence intervals, lag-1
ACF noise identification, and Kalman-based clock state estimation with
PID steering.

## Package layout

SigmaTau.jl is a Julia 1.11 workspace monorepo. Three subpackages are
re-exported through an umbrella `SigmaTau`:

| Subpackage | Responsibility |
|---|---|
| [`SigmaTauBase`](reference/base.md)         | Core data types: `PhaseData`, `FrequencyData`, `StabilityResult`. |
| [`SigmaTauStability`](reference/stability.md) | Deviation kernels, EDF / CI, noise identification. |
| [`SigmaTauEnsemble`](reference/ensemble.md)   | Clock state-space models, Kalman filter, PID steering. |

## Where to next

- **First time?** See [Getting Started](getting_started.md).
- **Need theory?** Start with [Theory: Overview](theory/overview.md).
- **Want a worked example?** Pick a [Tutorial](tutorials/01_phase_data.md).
- **Looking up an API?** Browse the [API Reference](reference/base.md).
- **Verifying numerics?** See [Validation](validation/methodology.md).

## Reference math

Numerical conventions and confidence-interval formulas follow:

- NIST SP1065 [Riley & Howe 2008](@cite riley-2008-sp1065)
- Greenhall–Riley 2003 EDF approximations [Greenhall 2003](@cite greenhall-2003-edf-stability)
- [IEEE Std 1139-2022](@cite ieee1139-2022-definitions)

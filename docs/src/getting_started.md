# Getting Started

## Installation

SigmaTau.jl is a single registerable package. From the Pkg REPL:

```julia-repl
pkg> add https://github.com/ianlap/SigmaTau.jl
```

`using SigmaTau` brings the shared types (`PhaseData`, `FrequencyData`,
`StabilityResult`) and every public deviation, noise-ID, and IO function
into scope under one flat namespace.

## A minimal example

```@example basic
using SigmaTau
using Random
Random.seed!(42)

# 10 minutes of phase residuals at 1 Hz, white-PM-like
phase = randn(600)
p = PhaseData(phase, 1.0)

# Compute overlapping Allan deviation at three τ values
r = adev(p, [1, 4, 16]; ci=true)

# Round for stable display
round.(r.dev; sigdigits=4)
```

The result is a [`StabilityResult`](reference/types.md#SigmaTau.StabilityResult)
populated with `tau`, `dev`, per-τ `noise_type`, χ²-based confidence bounds
`ci` (one `(lo, hi)` tuple per τ — `r.ci[1].lo`, or plain vectors via the
`ci_lower` / `ci_upper` accessors), `edf` per Greenhall–Riley, and `neff`,
the number of analysis windows averaged at each τ.

## Where to next

- [Migrating from Stable32](migration_from_stable32.md) — coming from Stable32? Start here for the run-type and result-column mapping, then the tutorial below for a longer walkthrough.
- [Migrating from allantools](migration_from_allantools.md) — the function-name and calling-convention mapping for Python allantools users.
- [Tutorial 0: Julia for metrologists](tutorials/00_julia_for_metrologists.md) — covers installation, loading `.DAT` files, your first `adev`, overlaying plots, and saving results to disk.
- [Tutorial 1: Phase data](tutorials/01_phase_data.md) — `PhaseData` / `FrequencyData` basics.
- [Tutorial 2: Computing Allan deviation](tutorials/02_compute_adev.md) — `adev` and the `StabilityResult` it returns (incl. EDF, χ² CI, noise type).
- [Tutorial 3: A drifting clock and the Hadamard family](tutorials/03_drifting_clock_hadamard.md) — `hdev`, `mhdev`, `htdev`, and `mhtotdev` on a record with linear frequency drift.
- [Tutorial 4: Reading and saving your data](tutorials/04_reading_your_data.md) — `read_phase` / `read_frequency` for `.DAT` files and counter CSVs, `fillgaps` / `detrend` data prep, and `save_result` / `save_suite` for sharing results.
- [Tutorial 6: Three-cornered hat](tutorials/06_three_cornered_hat.md) — separating clock noise from a three-clock comparison.
- [Theory: Overview](theory/overview.md) — what σ_y(τ) means.
- [FAQ](faq.md) — common stumbling blocks, including the `Package Plots not found` error.

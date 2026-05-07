# Changelog

All notable changes to **SigmaTau.jl** are tracked here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- `tdev(::PhaseData, m_values; …)` API wrapper. Wraps `mdev` and scales by
  `τ/√3`; CI bounds inherit MDEV's χ²/Gaussian limits scaled by the same
  factor. Now exported from `SigmaTauStability`.
- `FrequencyData` entry points for every stability deviation (`adev`, `mdev`,
  `tdev`, `hdev`, `mhdev`, `ldev`, `totdev`, `mtotdev`, `htotdev`,
  `mhtotdev`). Each converts via `_freq_to_phase` (`cumsum(y) · τ₀`) and
  delegates to the existing `PhaseData` method. Lives in
  `lib/SigmaTauStability/src/utils.jl`.
- `edf::Vector{Float64}` field on `StabilityResult`. Populated when
  `calc_ci=true`; empty when `calc_ci=false`. Lets users compute custom
  confidence intervals without re-running noise identification.
- `SigmaTauRecipesBaseExt` package extension at
  `ext/SigmaTauRecipesBaseExt.jl`. Provides a `RecipesBase.@recipe` for
  `StabilityResult` that draws a log-log τ–σ plot with optional
  `ci_lower`/`ci_upper` error bars; activates automatically when
  `RecipesBase` (or `Plots`) is loaded alongside `SigmaTau`.
- Boundary test for `identify_noise` at `N_eff ∈ {29, 31}` to lock in the
  `NEFF_RELIABLE` threshold.

### Changed

- `NEFF_RELIABLE` lowered from 50 → 30 in
  `lib/SigmaTauStability/src/noise/lag1.jl` per the legacy GEMINI.md §2
  mandate. `identify_noise` now switches from the B1-ratio fallback to the
  lag-1 ACF path at `N_eff = 30` instead of `N_eff = 50`.
- Plot recipes are no longer compiled into the umbrella module. The previous
  `src/PlotRecipes.jl` stub has been deleted; recipe code lives entirely in
  the `RecipesBase` package extension. The umbrella package no longer pulls
  any plotting dependency unless one is explicitly loaded.

### Fixed / Verified

- Kalman filter parity scaffolding (legacy `safe_sqrt` clamping, opt-in
  `legacy_compat=true` flag, `clamp_covariance_diag` helper) is in place.
  Stale `test_output.log` from the prior session predates the fix and shows
  failures that no longer reflect the source. **Pending**: re-run
  `Pkg.test()` against the current tree to confirm 4/4 parity.

## [0.1.0] — 2026-05-07

Initial commit. Modularised legacy SigmaTau into three subpackages
(`SigmaTauBase`, `SigmaTauStability`, `SigmaTauEnsemble`) plus an umbrella
`SigmaTau` package. Workspace + path-source wiring per Julia 1.11 monorepo
semantics. MIT licensed.

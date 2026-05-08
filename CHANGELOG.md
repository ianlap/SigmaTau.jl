# Changelog

All notable changes to **SigmaTau.jl** are tracked here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- `detrend::Symbol` kwarg on the four total-family kernels and API
  wrappers (`totdev`, `mtotdev`, `htotdev`, `mhtotdev`). Each kernel
  exposes its canonical recipe and a `:linear` alternative that swaps
  the kernel's natural detrending for a per-window full LS slope+intercept
  while preserving the kernel's natural extension shape. `:legacy`
  aliases the pre-1.0 default of each kernel for backward compat.
- `_totdev_howe` helper implementing NIST SP1065 eqn 25 verbatim
  (Greenhall–Howe–Percival 1998 eq 3): no detrend, mean-flip endpoint
  reflection, sum over centers `n=2..N-1`. Cross-checked against
  Stable32's TOTDEV at `rtol=1e-4` (12 of 13 m values) and against
  allantools' `totdev` to ~7 significant figures on the same fixture.
  The single (m=512) Stable32 row that misses the rtol=1e-4 floor is
  identified by Stable32 as FLFM and carries an alpha-aware correction
  that diverges from the raw SP1065 value; documented in the test
  comment and `TODO.md`.

### Changed

- **Breaking:** `totdev` default detrend recipe is now `:howe` (SP1065
  eqn 25: no detrend, mean-flip endpoint reflection). Previous behavior
  (global LS detrend on top of the same reflection) is available via
  `totdev(...; detrend=:legacy)` or the `:linear` alias. Output values
  change for all τ; the new default matches allantools' raw `totdev` to
  ~7 significant figures and Stable32's `Total` column at `rtol=1e-4`
  (m=512 outlier excepted, see Added section). Bias correction
  `bias_correction(:totvar, ...)` is now correctly calibrated against
  the new `:howe` default — the previous `:legacy` path silently
  removed drift before the SP1065 bias factor was applied.
- **Breaking:** `mhtotdev` default detrend recipe is now `:greenhall`
  (per-window half-mean slope removal), aligning with the convention of
  MTOT and HTOT in the Hadamard-modified family. Previous behavior
  (per-window full LS) is available via `mhtotdev(...; detrend=:linear)`
  or `:legacy`. MHTOTDEV is novel to SigmaTau; no external numerical
  reference exists, so the recipe choice is a methodology decision
  rather than a parity contract.
- MTOTDEV and HTOTDEV default outputs are unchanged. Both keep
  `:greenhall` (per-window half-mean) as the default; the new `:linear`
  recipe is opt-in.
- MHTOTDEV bias-correction policy made explicit: `bias_correction`
  now short-circuits to `B = 1` for `var_type = :mhtot`, with the
  rationale promoted into the function docstring (FCS 2001 and
  NIST SP1065 publish no bias model for MHTOTDEV; Stable32 and
  AllanLab also treat it as unbiased). `mhtotdev` docstring
  rewritten to state this contract instead of the inline
  comment it carried before.
- `calculate_edf` docstring expanded with the WPM/FLPM (α ∈ {1, 2})
  fallback rationale for TOTDEV and HTOTDEV — the substitute
  ADEV/HDEV-style Greenhall–Riley formula is a documented policy
  choice, not a derivation, since SP1065 Table 9 and FCS 2001
  cover α ∈ {0, -1, -2} only.
- `Random` is no longer a hard dep of `SigmaTauStability`. The
  `_gen_powerlaw_phase` helper now uses the global RNG (callers seed via
  `Random.seed!`) so the convention "Random in [extras] for both
  subpackages" is preserved. `Random` stays in `[extras]` and the test
  target.
- Test runner under `lib/SigmaTauEnsemble/test/runtests.jl` now wraps the
  legacy `clock_model.jl` / `filter.jl` includes in a `LegacyKF` module
  and selectively imports the symbols it needs. This stops the legacy
  `step!`/`PIDController` from shadowing the new ones at test scope.
- `.gitignore` now excludes the `CLAUDE.md` and `AGENTS.md` agent-context
  briefs. They are local-only working documents and should not appear on
  the public repo.
- **Project workflow**: every shipped change must remove its TODO entry and
  add a CHANGELOG line under `## [Unreleased]` in the same commit. TODO.md
  rewritten to drop the now-stale "Author README", "Maintain CHANGELOG"
  bullets and re-prioritise remaining work into Critical / High / Medium /
  Low / Documentation buckets.

### Added

- PID steering controller in `SigmaTauEnsemble`. New `PIDController`
  struct (g_p / g_i / g_d gains + integral state) plus `step!(pid, x)` and
  `steer_to_correction(steer, ns, dt)`. `predict!(est, model, dt;
  steering=…)` now accepts an optional steering vector that's added to the
  propagated state mean — matches the legacy `filter_step!` semantics where
  a PID's last steer is folded in after Φ propagation.
- `examples/clock_steering.jl` — PID + Kalman walkthrough on a drifting
  clock; demonstrates ~600× residual-phase reduction vs the unsteered case.
- Power-law phase-noise synthesizer in `SigmaTauStability` at
  `src/noise/synth.jl`. `_gen_powerlaw_phase(α, N; tau0, rng)` generates an
  N-sample phase residual whose fractional-frequency PSD ∝ `f^α` via
  `f^(α/2)` shaping of white Gaussian noise (DC zeroed, integrated to phase).
  Non-exported helper. Requires an `AbstractFFTs` backend (e.g. `using
  FFTW`) loaded by the caller; `AbstractFFTs` is now a hard dep, `FFTW`
  is in test extras.
- MTOTDEV multi-noise validation extended from 3 → 5 power-law types
  (added FLPM and FLFM via the new synthesizer). New
  "ADEV/HDEV across all 5 power-law noise types" testset extends
  `_adev_core` / `_mdev_core` / `_hdev_core` / `_mhdev_core` kernel parity
  to the same alpha range. 120 new assertions; full Stability test count
  is now `339/339` passing.
- Stable32 cross-validation testset against the
  `reference/validation/stable32gen.DAT` fixture (8192 phase samples) and
  Stable32's published outputs (`stable32_data_full.csv`). Verifies the new
  kernels match Stable32's reported sigmas at:
  - `rtol=1e-4` for OADEV / Modified Allan / Overlapping Hadamard / Time
  - `rtol=0.05` for raw MTOTDEV (Stable32 reports unbiased; our API applies
    SP1065 B≈1.27 — kernel match is ~3% per the legacy comparison report)
  - `rtol=0.10` for raw HTOTDEV (~0.5% bias offset + boundary effects)
  - `rtol=0.15` for TOTDEV (close at short τ; documented boundary-reflection
    discrepancy at long τ)

  86 new assertions (85 row-driven + 1 sanity-count).
- Strict numerical parity testset for the eight stability kernels
  (`_adev_core`, `_mdev_core`, `_hdev_core`, `_mhdev_core`, `_totdev_core`,
  `_mtotdev_core`, `_htotdev_core`, `_mhtotdev_core`) against the legacy
  SigmaTau Julia reference. The legacy kernels are inlined verbatim under
  `lib/SigmaTauStability/test/legacy_kernels.jl` so the tests run on CI
  without depending on the gitignored `legacy/` tree. 52 new assertions
  at `rtol=1e-12`.
- MTOTDEV multi-noise validation: kernel + end-to-end pipeline checks
  against WPM, WHFM, and RWFM synthetic fixtures.
- TOTDEV/HTOTDEV WPM/FLPM EDF fallback: when `_coeff_totvar` /
  `_coeff_htot` return `(NaN, NaN)` for `α∈{1,2}` (since SP1065 Table 9 /
  FCS 2001 only cover `α∈{0,-1,-2}`), `calculate_edf` now falls back to
  the ADEV-style (`d=2`) or HDEV-style (`d=3`) Greenhall–Riley formula,
  yielding finite EDFs for every noise type instead of NaN.
- GitHub Actions CI: matrix on Julia 1.11 × {ubuntu, macOS} × all three
  subpackages, plus an umbrella `using SigmaTau` smoke job.
- `examples/quickstart.jl`: end-to-end walkthrough exercising
  `adev`/`mdev`/`tdev`/`hdev`/`totdev`, `FrequencyData` input, and a
  `ThreeStateClock` Kalman filter. Verified to run.

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

- Kalman filter parity now confirmed end-to-end: `Pkg.test()` against
  `SigmaTauEnsemble` reports **15/15 pass**, including the 4 `legacy_compat`
  parity assertions (phase, frequency, drift, P-history). The stale
  `test_output.log` from the prior session predated the
  `clamp_covariance_diag` wiring and is no longer representative.
- `Pkg.test()` against `SigmaTauStability` reports **46/46 pass**, covering
  the new `tdev` wrapper, `FrequencyData ↔ PhaseData` equivalence, the
  `edf` field, and the `NEFF_RELIABLE = 30` boundary at `N_eff ∈ {29, 31}`.
- `Random` was missing from both subpackages' `[extras]`; added so
  `Pkg.test()` resolves the standard-library dependency. Without this,
  `Random.seed!` calls inside test files raised `Package Random not found`
  before any assertion could run.

## [0.1.0] — 2026-05-07

Initial commit. Modularised legacy SigmaTau into three subpackages
(`SigmaTauBase`, `SigmaTauStability`, `SigmaTauEnsemble`) plus an umbrella
`SigmaTau` package. Workspace + path-source wiring per Julia 1.11 monorepo
semantics. MIT licensed.

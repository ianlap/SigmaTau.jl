# SigmaTau.jl — Next-Step Agent Prompts

> Prioritized by impact: correctness bugs > missing core math > missing tests > polish > new features

---

## PROMPT 1 — 🔴 Fix Kalman Filter Parity Failure

**Priority**: Critical — blocking  
**What**: Debug and fix the `StandardKalmanFilter` in `SigmaTauEnsemble` to achieve bit-identical parity with the legacy `kalman_filter` function.

**Context**: All 4 parity tests fail (`test_output.log`). The new filter diverges from legacy starting at step 2. Key symptoms:
- Legacy `P_history[:,:,1]` is all zeros; new filter shows `[1e-22, 0, 0; 0, 1e-12, 0; 0, 0, 1e-12]`
- Phase/freq/drift estimates diverge by orders of magnitude by step 5
- The legacy `safe_sqrt` clamping and predict/update ordering are likely mismatched

**Instructions**:
1. Read `legacy/julia/src/filter.jl` line-by-line to understand the exact predict→update→store ordering
2. Read `legacy/julia/src/clock_model.jl` to verify Φ, Q, H, R construction matches `SigmaTauEnsemble/src/models/clocks.jl`
3. Compare the `predict!` skip-on-k==0 logic against the legacy initialization convention
4. Check whether legacy applies `safe_sqrt` or similar P-clamping that changes the covariance trajectory
5. Fix `lib/SigmaTauEnsemble/src/estimators/filters.jl` to match legacy behavior exactly
6. If `safe_sqrt` is the cause, implement it as an *optional* kwarg (default off) to preserve the AD-clean path while enabling legacy parity

**Reference files**:
- [legacy/julia/src/filter.jl](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/julia/src/filter.jl) (8.2 KB)
- [legacy/julia/src/clock_model.jl](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/julia/src/clock_model.jl) (8.0 KB)
- [legacy/docs/equations/kalman.md](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/docs/equations/kalman.md)
- [lib/SigmaTauEnsemble/test_output.log](file:///Users/ianlapinski/Downloads/SigmaTau-dev/lib/SigmaTauEnsemble/test_output.log)

**Target files**:
- `lib/SigmaTauEnsemble/src/estimators/filters.jl`
- `lib/SigmaTauEnsemble/test/runtests.jl`

**Acceptance criteria**:
- All 4 parity tests pass: `phase_est`, `freq_est`, `drift_est`, `P_history` match legacy within `rtol=1e-12`
- AD-friendliness preserved: no in-place array mutation in hot path
- `StaticArrays` usage retained

---

## PROMPT 2 — 🔴 Fix Package Resolution (Project.toml Wiring)

**Priority**: Critical — blocking  
**What**: Fix all `Project.toml` files so the monorepo resolves correctly with `Pkg.resolve()` from the root.

**Context**: `err.log` shows `expected package SigmaTauBase [e7b0a8c4] to be registered`. Local packages need `[sources]` sections pointing to relative paths (Julia 1.11 workspace feature).

**Instructions**:
1. Add `StaticArrays` to `lib/SigmaTauEnsemble/Project.toml` `[deps]` section
2. Add `SigmaTauBase` to root `Project.toml` `[deps]` section
3. Add `[sources]` sections to each subpackage's `Project.toml` pointing to `SigmaTauBase`:
   ```toml
   [sources.SigmaTauBase]
   path = "../SigmaTauBase"
   ```
4. Add `[sources]` to root `Project.toml` for all three subpackages
5. Run `julia --project=. -e 'using Pkg; Pkg.resolve()'` from root to verify
6. Run `julia --project=lib/SigmaTauStability -e 'using Pkg; Pkg.resolve()'` to verify subpackage

**Target files**:
- `Project.toml` (root)
- `lib/SigmaTauStability/Project.toml`
- `lib/SigmaTauEnsemble/Project.toml`

**Acceptance criteria**:
- `Pkg.resolve()` succeeds from root
- `using SigmaTau` loads without errors in a clean Julia session
- `Pkg.test()` can be invoked on each subpackage independently

---

## PROMPT 3 — 🟡 Add NIST SP1065 Parity Tests for All Deviation Kernels

**Priority**: High — correctness  
**What**: Replace the weak "is finite" tests with rigorous numerical parity tests against known SP1065 reference values or legacy Julia outputs.

**Context**: Current tests only check `isfinite` and `length`. The walkthrough claims "30/30 pass" but these tests would pass even with wrong numerical values. Legacy `legacy/julia/src/engine.jl` and `legacy/julia/src/deviations/` contain the reference implementations.

**Instructions**:
1. Read `legacy/julia/src/engine.jl` to understand the legacy calling convention
2. Generate reference data: create a fixture of ~1000 points of known noise type (e.g., white FM, random walk FM)
3. Run the legacy `_adev_kernel`, `_mdev_kernel`, `_hdev_kernel`, etc. on this fixture
4. Save expected outputs as hardcoded vectors in the test file
5. Compare `_adev_core`, `_mdev_core`, etc. against these reference values at `rtol=1e-12`
6. Test at multiple m-values including m=1, m=10, m=100
7. Test edge cases: N exactly at boundary, m that makes L=1

**Reference files**:
- [legacy/julia/src/engine.jl](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/julia/src/engine.jl) (5.7 KB)
- `legacy/julia/src/deviations/` directory
- [legacy/docs/equations/allan.md](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/docs/equations/allan.md)
- [legacy/docs/equations/hadamard.md](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/docs/equations/hadamard.md)
- [legacy/docs/equations/total.md](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/docs/equations/total.md)

**Target files**:
- `lib/SigmaTauStability/test/runtests.jl`

**Acceptance criteria**:
- At least 5 m-values tested per kernel (9 kernels × 5 = 45 minimum assertions)
- All deviations match legacy within `rtol=1e-12` for at least 2 noise types
- Noise ID categorical matches tested for N_eff above and below `NEFF_RELIABLE`
- Tests reproducible (seeded RNG)

---

## PROMPT 4 — 🟡 Add Missing `tdev` API Wrapper

**Priority**: Medium — completeness  
**What**: Add a user-facing `tdev(::PhaseData)` function that returns a `StabilityResult`, matching the pattern of `adev`, `mdev`, etc.

**Context**: `_tdev_core` exists in `core/allan.jl` and is exported, but there is no `tdev(data::PhaseData, m_values)` wrapper. TDEV = τ·MDEV/√3, so it should follow the same pattern as `ldev` wrapping `mhdev`.

**Instructions**:
1. Add `tdev` to `api/allan.jl` following the exact pattern of `ldev` in `api/hadamard.jl`
2. It should call `mdev(data, m_values; ...)` and scale: `TDEV = τ · MDEV / √3`
3. CI bounds should scale proportionally (same as ldev approach)
4. Export `tdev` from `SigmaTauStability.jl`
5. Add tests

**Reference files**:
- [legacy/docs/equations/allan.md](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/docs/equations/allan.md) (TDEV formula)
- [api/hadamard.jl](file:///Users/ianlapinski/Downloads/SigmaTau-dev/lib/SigmaTauStability/src/api/hadamard.jl) (ldev pattern to follow)

**Target files**:
- `lib/SigmaTauStability/src/api/allan.jl`
- `lib/SigmaTauStability/src/SigmaTauStability.jl` (add export)
- `lib/SigmaTauStability/test/runtests.jl`

**Acceptance criteria**:
- `tdev(PhaseData(...), [1,2,4])` returns a `StabilityResult` with `deviation_type == :tdev`
- Values match `_tdev_core` output
- CI bounds present when `calc_ci=true`
- Test passes

---

## PROMPT 5 — 🟡 Validate MTOTDEV Across All 5 Noise Types

**Priority**: Medium — correctness  
**What**: Implement the "Mtot multi-noise validation" task from `legacy/TODO.md`.

**Context**: MTOTDEV and its bias correction have not been validated across WPM, FLPM, WHFM, FLFM, and RWFM noise types. The bias correction table in `_coeff_mtot` uses approximate coefficients that need verification.

**Instructions**:
1. Generate synthetic phase data for each of the 5 standard noise types using `legacy/julia/src/noise_gen.jl`
2. Run `_mtotdev_core` on each dataset and compare against legacy `_mtotdev_kernel`
3. Verify `bias_correction(:mtot, ...)` produces values consistent with Stable32 output
4. Document any discrepancies found

**Reference files**:
- [legacy/julia/src/noise_gen.jl](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/julia/src/noise_gen.jl) (3.9 KB)
- [legacy/docs/equations/total.md](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/docs/equations/total.md)

**Target files**:
- `lib/SigmaTauStability/test/runtests.jl` (add noise-type validation testset)

**Acceptance criteria**:
- 5 noise types × 3+ m-values tested for mtotdev
- Bias correction direction verified (divided, not multiplied)
- All pass within `rtol=1e-10`

---

## PROMPT 6 — 🟡 Update NEFF_RELIABLE Threshold

**Priority**: Medium — standards compliance  
**What**: Update `NEFF_RELIABLE` from 50 to 30 in `noise/lag1.jl` per GEMINI.md §2 mandate.

**Context**: Legacy TODO.md says: "Update `NEFF_RELIABLE` to 30 (from 50) in both MATLAB and Julia to match GEMINI.md mandate §2." The new monorepo copied the legacy value of 50.

**Instructions**:
1. Change `const NEFF_RELIABLE = 50` to `const NEFF_RELIABLE = 30` in `lib/SigmaTauStability/src/noise/lag1.jl`
2. Update any tests that depend on the threshold boundary (current test uses N=100, m=4 → N_eff=25 which is below both 30 and 50, so the test should still work)
3. Add a test specifically at the boundary: N_eff=29 (should use B1-ratio) and N_eff=31 (should use lag-1 ACF)

**Target files**:
- `lib/SigmaTauStability/src/noise/lag1.jl` (line 4)
- `lib/SigmaTauStability/test/runtests.jl`

**Acceptance criteria**:
- `NEFF_RELIABLE == 30`
- Boundary test at N_eff=29 and N_eff=31 passes

---

## PROMPT 7 — 🟡 Add `_coeff_totvar` Entries for α=2 and α=1

**Priority**: Medium — correctness  
**What**: The `_coeff_totvar` function in `edf.jl` only handles α ∈ {0, -1, -2} and returns `(NaN, NaN)` for WPM (α=2) and FLPM (α=1) noise. This means `calculate_edf(:totdev, ...)` produces NaN for these noise types.

**Instructions**:
1. Research the correct Greenhall/Riley totvar EDF coefficients for α=2 and α=1 from SP1065 or GHP99
2. Add entries to `_coeff_totvar`
3. Similarly verify `_coeff_htot` (currently only handles α ∈ {0, -1, -2})
4. Add test coverage for totdev EDF at each noise type

**Reference files**:
- [legacy/julia/src/stats.jl](file:///Users/ianlapinski/Downloads/SigmaTau-dev/legacy/julia/src/stats.jl) (11.5 KB)
- SP1065 tables (in `legacy/docs/papers/`)

**Target files**:
- `lib/SigmaTauStability/src/stats/edf.jl`
- `lib/SigmaTauStability/test/runtests.jl`

**Acceptance criteria**:
- `_coeff_totvar(2)` and `_coeff_totvar(1)` return valid coefficients
- `calculate_edf(:totdev, ...)` produces finite values for all 5 noise types
- Tests pass

---

## PROMPT 8 — 🟡 Add FrequencyData Support to All API Functions

**Priority**: Medium — completeness  
**What**: Add method dispatches so all user-facing API functions (`adev`, `mdev`, `hdev`, ...) accept `FrequencyData` in addition to `PhaseData`.

**Context**: `FrequencyData` is defined in `SigmaTauBase` but no API function accepts it. The conversion is `x[k] = τ₀ · Σⱼ₌₁ᵏ y[j]` (phase = cumulative sum of fractional frequency × τ₀).

**Instructions**:
1. Add a helper function `_freq_to_phase(data::FrequencyData) → PhaseData` in a shared utility file
2. Add `adev(data::FrequencyData, ...)` dispatches that convert to PhaseData and delegate
3. Repeat for all 10 API functions
4. Add tests for at least `adev` and `hdev` with FrequencyData input

**Target files**:
- New file: `lib/SigmaTauStability/src/utils.jl`
- All files in `lib/SigmaTauStability/src/api/`
- `lib/SigmaTauStability/src/SigmaTauStability.jl` (include + export)
- `lib/SigmaTauStability/test/runtests.jl`

**Acceptance criteria**:
- `adev(FrequencyData(y, tau0), m_values)` returns identical results to `adev(PhaseData(cumsum(y)*tau0, tau0), m_values)`
- At least 2 API functions tested with FrequencyData

---

## PROMPT 9 — 🟢 Implement PlotRecipes for StabilityResult

**Priority**: Low — polish  
**What**: Uncomment and finalize the `PlotRecipes.jl` stub to provide automatic sigma-tau plotting.

**Instructions**:
1. Add `RecipesBase` to root `Project.toml` deps (or make it a package extension for `Plots.jl`)
2. Implement `@recipe function f(res::StabilityResult)` with:
   - Log-log scaling
   - Error bars from `ci_lower`/`ci_upper`
   - Proper axis labels
   - Legend showing deviation type
3. Add a basic visual smoke test

**Reference files**:
- [src/PlotRecipes.jl](file:///Users/ianlapinski/Downloads/SigmaTau-dev/src/PlotRecipes.jl) (existing stub)

**Target files**:
- `src/PlotRecipes.jl`
- `Project.toml`

**Acceptance criteria**:
- `plot(adev_result)` produces a log-log sigma-tau plot
- Error bars displayed correctly
- No errors when plotting a result with empty CI vectors (`calc_ci=false`)

---

## PROMPT 10 — 🟢 Add StabilityResult `edf` Field

**Priority**: Low — usability  
**What**: Add an `edf::Vector{Float64}` field to `StabilityResult` so users can access the computed equivalent degrees of freedom.

**Context**: EDF is computed in every API function but discarded. Power users need EDF for custom confidence interval calculations.

**Instructions**:
1. Add `edf::Vector{Float64}` field to `StabilityResult` in `SigmaTauBase`
2. Update all API functions in `SigmaTauStability` to pass EDF into the constructor
3. When `calc_ci=false`, pass `Float64[]` for edf
4. Update tests

**Target files**:
- `lib/SigmaTauBase/src/SigmaTauBase.jl`
- All files in `lib/SigmaTauStability/src/api/`
- `lib/SigmaTauStability/test/runtests.jl`

**Acceptance criteria**:
- `result.edf` is accessible and contains finite values when `calc_ci=true`
- `result.edf` is empty when `calc_ci=false`
- No breaking changes to existing API

# SigmaTau.jl Restructure Plan — Workspace → Single Package

> **For agentic workers:** Use `superpowers:executing-plans` (inline) to
> execute. Steps use checkbox (`- [ ]`) syntax. Plan author: 2026-05-09.

**Goal:** Collapse the three-subpackage workspace
(`SigmaTauBase` / `SigmaTauStability` / `SigmaTauEnsemble`) into a single
registerable `SigmaTau` package whose internals are organized as two
submodules — `SigmaTau.Stab` and `SigmaTau.Est` — with shared types
(`PhaseData`, `FrequencyData`, `StabilityResult`) at the top level. All
existing exports are re-exported from the parent so casual user code
(`using SigmaTau; adev(data, [1,2,4])`) keeps working unchanged.

**Architecture:** Hybrid pattern (Option B3). Top-level package owns the
shared types and the umbrella module; two child modules (`Stab`, `Est`)
are defined inline in `src/SigmaTau.jl` via `module Stab … end` blocks
that `include` from `src/stab/…` and `src/est/…`. Both children import
shared types from the parent via `..SigmaTau`. The umbrella does
`using .Stab; using .Est` so every public symbol stays unqualified at
the top level. Plot recipes remain in `ext/SigmaTauRecipesBaseExt.jl`.

**Tech Stack:** Julia 1.11, single Project.toml, weakdep extension for
RecipesBase, vanilla Documenter for docs.

---

## File Move Map (source → destination)

### Shared types (Base subpackage → top-level `src/types/`)

`lib/SigmaTauBase/src/SigmaTauBase.jl` is one file holding four
declarations. Split it into four files:

| Old (single file, declaration) | New |
|---|---|
| `lib/SigmaTauBase/src/SigmaTauBase.jl::AbstractTimingData` | `src/types/abstract.jl` |
| `lib/SigmaTauBase/src/SigmaTauBase.jl::PhaseData` | `src/types/phase_data.jl` |
| `lib/SigmaTauBase/src/SigmaTauBase.jl::FrequencyData` | `src/types/frequency_data.jl` |
| `lib/SigmaTauBase/src/SigmaTauBase.jl::StabilityResult` | `src/types/stability_result.jl` |

The `module SigmaTauBase` / `end` wrapping and `using DocStringExtensions`
header are dropped from these split files; the parent `SigmaTau.jl` will
declare `using DocStringExtensions` once at the top so the `$(TYPEDFIELDS)`
docstring macros keep resolving.

### Stab submodule (was `SigmaTauStability`)

| Old | New |
|---|---|
| `lib/SigmaTauStability/src/utils.jl` | `src/stab/utils.jl` |
| `lib/SigmaTauStability/src/core/allan.jl` | `src/stab/core/allan.jl` |
| `lib/SigmaTauStability/src/core/hadamard.jl` | `src/stab/core/hadamard.jl` |
| `lib/SigmaTauStability/src/core/total.jl` | `src/stab/core/total.jl` |
| `lib/SigmaTauStability/src/noise/lag1.jl` | `src/stab/noise/lag1.jl` |
| `lib/SigmaTauStability/src/noise/synth.jl` | `src/stab/noise/synth.jl` |
| `lib/SigmaTauStability/src/stats/edf.jl` | `src/stab/stats/edf.jl` |
| `lib/SigmaTauStability/src/api/allan.jl` | `src/stab/api/allan.jl` |
| `lib/SigmaTauStability/src/api/hadamard.jl` | `src/stab/api/hadamard.jl` |
| `lib/SigmaTauStability/src/api/total.jl` | `src/stab/api/total.jl` |

The wrapper file `lib/SigmaTauStability/src/SigmaTauStability.jl` is
**not** moved — its `module SigmaTauStability … end`, `using` lines,
`include` ordering, and exports collapse into the `module Stab` block
inside the new `src/SigmaTau.jl`.

### Est submodule (was `SigmaTauEnsemble`)

| Old | New |
|---|---|
| `lib/SigmaTauEnsemble/src/models/clocks.jl` | `src/est/models/clocks.jl` |
| `lib/SigmaTauEnsemble/src/estimators/filters.jl` | `src/est/estimators/filters.jl` |

The wrapper `lib/SigmaTauEnsemble/src/SigmaTauEnsemble.jl` collapses
into the `module Est` block.

### Tests

| Old | New |
|---|---|
| `lib/SigmaTauStability/test/runtests.jl` | `test/stab/runtests.jl` |
| `lib/SigmaTauStability/test/legacy_kernels.jl` | `test/stab/legacy_kernels.jl` |
| `lib/SigmaTauStability/test/allantools_cross_validation.jl` | `test/stab/allantools_cross_validation.jl` |
| `lib/SigmaTauEnsemble/test/runtests.jl` | `test/est/runtests.jl` |

`SigmaTauBase` ships no tests of its own today (its `Project.toml` has no
`[extras]` test target and no `test/` directory). The plan therefore
creates `test/types/runtests.jl` from scratch with three smoke checks
(constructor of each type, basic field access). A new top-level
`test/runtests.jl` `include`s all three under `@testset`.

### Top-level files

| Old | New / change |
|---|---|
| `src/SigmaTau.jl` (8-line umbrella) | rewritten in place — see Template below |
| `Project.toml` (workspace root) | rewritten — single package, no `[sources]`, no `[workspace]` |
| `ext/SigmaTauRecipesBaseExt.jl` | edited — `using SigmaTauBase: StabilityResult` → `using SigmaTau: StabilityResult` |

### Vault Code-link rename

| Old folder | New folder |
|---|---|
| `legdocs/obsidian/SigmaTauVault/Code links/base/` | `legdocs/obsidian/SigmaTauVault/Code links/types/` |
| `legdocs/obsidian/SigmaTauVault/Code links/stability/` | `legdocs/obsidian/SigmaTauVault/Code links/stab/` |
| `legdocs/obsidian/SigmaTauVault/Code links/ensemble/` | `legdocs/obsidian/SigmaTauVault/Code links/est/` |

47 Code-link `.md` files have YAML metadata that needs updating. The
prose body of every note is left alone — only `module:`, `file:`, and
`subpackage:` fields change.

### Recovery point

| Old | New |
|---|---|
| `lib/` (entire workspace tree) | `lib.bak/` (gitignored, kept on disk) |

---

## Phase A — Plan only (this file)

- [x] **Step A1: Author `docs/_restructure_plan.md`** (this document).
- [ ] **Step A2: Pause for user approval before any file move.**

User: review this plan in full and reply *approved* (or with edits).
Phase B does not start until approval.

---

## Phase B — Execute restructure

### Task B1 — Stage new directory tree (no deletions yet)

**Files:** create `src/types/`, `src/stab/`, `src/est/`, `test/types/`,
`test/stab/`, `test/est/`. Copy (not move) the source files from `lib/`
into the new locations.

- [ ] **Step B1.1: Make new src directories.**

```bash
mkdir -p src/types src/stab/{api,core,noise,stats} src/est/{models,estimators}
mkdir -p test/types test/stab test/est
```

- [ ] **Step B1.2: Copy stab sources into `src/stab/`.**

```bash
cp lib/SigmaTauStability/src/utils.jl                src/stab/utils.jl
cp lib/SigmaTauStability/src/core/allan.jl           src/stab/core/allan.jl
cp lib/SigmaTauStability/src/core/hadamard.jl        src/stab/core/hadamard.jl
cp lib/SigmaTauStability/src/core/total.jl           src/stab/core/total.jl
cp lib/SigmaTauStability/src/noise/lag1.jl           src/stab/noise/lag1.jl
cp lib/SigmaTauStability/src/noise/synth.jl          src/stab/noise/synth.jl
cp lib/SigmaTauStability/src/stats/edf.jl            src/stab/stats/edf.jl
cp lib/SigmaTauStability/src/api/allan.jl            src/stab/api/allan.jl
cp lib/SigmaTauStability/src/api/hadamard.jl         src/stab/api/hadamard.jl
cp lib/SigmaTauStability/src/api/total.jl            src/stab/api/total.jl
```

- [ ] **Step B1.3: Copy est sources into `src/est/`.**

```bash
cp lib/SigmaTauEnsemble/src/models/clocks.jl         src/est/models/clocks.jl
cp lib/SigmaTauEnsemble/src/estimators/filters.jl    src/est/estimators/filters.jl
```

- [ ] **Step B1.4: Copy stab tests into `test/stab/`.**

```bash
cp lib/SigmaTauStability/test/runtests.jl                  test/stab/runtests.jl
cp lib/SigmaTauStability/test/legacy_kernels.jl            test/stab/legacy_kernels.jl
cp lib/SigmaTauStability/test/allantools_cross_validation.jl test/stab/allantools_cross_validation.jl
```

- [ ] **Step B1.5: Copy est test into `test/est/`.**

```bash
cp lib/SigmaTauEnsemble/test/runtests.jl test/est/runtests.jl
```

### Task B2 — Split `SigmaTauBase` into four `src/types/` files

**Files:** create `src/types/abstract.jl`, `src/types/phase_data.jl`,
`src/types/frequency_data.jl`, `src/types/stability_result.jl`.
Source content lifted verbatim from
`lib/SigmaTauBase/src/SigmaTauBase.jl` (lines 5–67 of that file). Drop
the `module SigmaTauBase` wrapper, `using DocStringExtensions`, and
`export …` line — these move to the parent.

- [ ] **Step B2.1: Write `src/types/abstract.jl`.**

```julia
"""
    AbstractTimingData

Supertype of timing data records (`PhaseData`, `FrequencyData`).
"""
abstract type AbstractTimingData end
```

- [ ] **Step B2.2: Write `src/types/phase_data.jl`.**

```julia
"""
    PhaseData{T<:AbstractFloat}

Phase residuals `x(t)` sampled at uniform interval `tau0`.

$(TYPEDFIELDS)
"""
struct PhaseData{T<:AbstractFloat} <: AbstractTimingData
    "Phase samples in seconds."
    x::Vector{T}
    "Base sample interval τ₀ in seconds."
    tau0::Float64
end
```

- [ ] **Step B2.3: Write `src/types/frequency_data.jl`.**

```julia
"""
    FrequencyData{T<:AbstractFloat}

Fractional-frequency samples `y(t)` at uniform interval `tau0`.

$(TYPEDFIELDS)
"""
struct FrequencyData{T<:AbstractFloat} <: AbstractTimingData
    "Fractional-frequency samples (dimensionless)."
    y::Vector{T}
    "Base sample interval τ₀ in seconds."
    tau0::Float64
end
```

- [ ] **Step B2.4: Write `src/types/stability_result.jl`.**

```julia
"""
    StabilityResult

Unified return type for every stability calculation.

The `noise_type`, `ci_lower`, `ci_upper`, and `edf` vectors are empty when
the calculation was invoked with `calc_ci=false`.

$(TYPEDFIELDS)
"""
struct StabilityResult
    "Which deviation produced this result (e.g. `:adev`, `:mdev`)."
    deviation_type::Symbol
    "Analysis intervals τ in seconds."
    tau::Vector{Float64}
    "Stability deviation σ_y(τ) per interval."
    dev::Vector{Float64}
    "Noise-type symbol identified at each τ (empty unless `calc_ci=true`)."
    noise_type::Vector{Symbol}
    "Lower CI bound (empty unless `calc_ci=true`)."
    ci_lower::Vector{Float64}
    "Upper CI bound (empty unless `calc_ci=true`)."
    ci_upper::Vector{Float64}
    "Equivalent degrees of freedom (empty unless `calc_ci=true`)."
    edf::Vector{Float64}
end
```

### Task B3 — Rewrite `src/SigmaTau.jl`

**Files:** overwrite `src/SigmaTau.jl` (currently 8 lines of
`@reexport using SigmaTau{Base,Stability,Ensemble}`).

The new file owns:
- `using DocStringExtensions` (was per-subpackage)
- `using Statistics`, `using Distributions`, `using AbstractFFTs` —
  brought into the `Stab` submodule, NOT the parent
- `using LinearAlgebra`, `using StaticArrays` — into `Est`
- `module Stab` with the include order from
  `lib/SigmaTauStability/src/SigmaTauStability.jl` lines 16–28
  (utils.jl is included AFTER the cores in the existing module — keep
  that exact ordering, do not "tidy" it)
- `module Est` with the include order from
  `lib/SigmaTauEnsemble/src/SigmaTauEnsemble.jl` lines 7–8
- `using .Stab; using .Est` to flatten exports
- Re-export every symbol in the union of the two children's export lists

- [ ] **Step B3.1: Write the new `src/SigmaTau.jl`.**

```julia
module SigmaTau

using DocStringExtensions

# ── Shared types ────────────────────────────────────────────────────────
include("types/abstract.jl")
include("types/phase_data.jl")
include("types/frequency_data.jl")
include("types/stability_result.jl")

export AbstractTimingData, PhaseData, FrequencyData, StabilityResult

# ── Stab: clock-stability analysis ──────────────────────────────────────
module Stab
    using ..SigmaTau: AbstractTimingData, PhaseData, FrequencyData,
                      StabilityResult
    using Statistics
    using Distributions
    using DocStringExtensions

    """
    Package-wide default confidence factor used by every public deviation API
    (`adev`, `mdev`, `hdev`, `tdev`, `mhdev`, `htdev`, `totdev`, `mtotdev`,
    `htotdev`, `mhtotdev`) when `confidence` is not supplied.

    Set to 0.683 (1-sigma) — the time-and-frequency stability convention used
    by Stable32, AllanLab, allantools' published error bars, and the
    Greenhall–Riley uncertainty papers. Override per call by passing
    `confidence=0.95` (or any other level) explicitly.
    """
    const DEFAULT_CONFIDENCE = 0.683

    include("stab/core/allan.jl")
    include("stab/core/hadamard.jl")
    include("stab/core/total.jl")

    include("stab/noise/lag1.jl")
    include("stab/noise/synth.jl")
    include("stab/stats/edf.jl")

    include("stab/utils.jl")

    include("stab/api/allan.jl")
    include("stab/api/hadamard.jl")
    include("stab/api/total.jl")

    export _adev_core, _mdev_core, _tdev_core
    export _hdev_core, _mhdev_core
    export _totdev_core, _mtotdev_core, _htotdev_core, _mhtotdev_core

    export identify_noise, calculate_edf, confidence_intervals, bias_correction
    export DEFAULT_CONFIDENCE

    export adev, mdev, tdev
    export hdev, mhdev, htdev
    export ldev   # deprecated alias for htdev — remove in a future release
    export totdev, mtotdev, htotdev, mhtotdev
end

# ── Est: clock estimation ───────────────────────────────────────────────
module Est
    using ..SigmaTau: AbstractTimingData, PhaseData, FrequencyData,
                      StabilityResult
    using LinearAlgebra
    using StaticArrays
    using DocStringExtensions

    include("est/models/clocks.jl")
    include("est/estimators/filters.jl")

    export AbstractClockModel, TwoStateClock, ThreeStateClock, RelativisticClock
    export nstates, state_transition, process_noise, measurement_matrix, measurement_noise
    export AbstractEstimator, StandardKalmanFilter, UDFactorizedFilter, KuramotoOscillator
    export predict!, update!
    export PIDController, step!, steer_to_correction
end

# ── Flatten submodule exports onto the umbrella ─────────────────────────
using .Stab
using .Est

# Make submodules themselves accessible without qualification
export Stab, Est

# Plot recipes for `StabilityResult` live in the `SigmaTauRecipesBaseExt`
# package extension and load automatically when `RecipesBase` (or `Plots`) is.

end # module SigmaTau
```

> **Notes on this rewrite:**
> - `Stab` re-uses `_freq_to_phase` (defined in `stab/utils.jl`), and
>   `utils.jl` is included AFTER cores and noise/edf — same as the
>   current `SigmaTauStability.jl` ordering. Do not reorder.
> - Re-using `noise/synth.jl` requires `AbstractFFTs` — but the original
>   `Stability` Project.toml had `AbstractFFTs` in `[deps]` only, with
>   the actual FFT backend (`FFTW`) supplied by the test target. We
>   keep that contract: `synth.jl`'s `using AbstractFFTs` is its first
>   non-comment line; the parent `Stab` submodule does NOT itself
>   `using AbstractFFTs`.
> - `using .Stab; using .Est` brings every `export`ed name into
>   `SigmaTau`'s namespace; the explicit `export Stab, Est` then exposes
>   the submodules themselves so `SigmaTau.Stab.adev(...)` works.

### Task B4 — Update import statements in copied source files

Most copied files do NOT need edits — `_adev_core`, `_freq_to_phase`,
clock-model methods, etc. were unqualified inside their subpackage
modules and remain unqualified inside their new `Stab` / `Est`
submodules. The grep done in Phase A confirmed the only intra-package
type references are `PhaseData` / `FrequencyData` / `StabilityResult`
in the API wrappers — and those are now imported at the top of each
submodule via `using ..SigmaTau: …`, so they resolve unqualified.

Files that DO need edits:

- [ ] **Step B4.1: `src/stab/noise/synth.jl`** — keep `using AbstractFFTs`
  on line 12 unchanged. (Verifies copy.)

- [ ] **Step B4.2: `src/stab/noise/lag1.jl`** — keep `using Statistics`
  on line 2 unchanged. (Verifies copy.)

- [ ] **Step B4.3: `src/stab/stats/edf.jl`** — keep
  `using Statistics; using Distributions` on lines 2–3 unchanged.

  > These are top-level inside their old subpackage. Inside the new
  > `Stab` submodule they're harmless duplicates of the submodule-level
  > `using Statistics` / `using Distributions`. Leave them — removing
  > them is a churn risk and Julia tolerates redundant `using`.

- [ ] **Step B4.4: `ext/SigmaTauRecipesBaseExt.jl`** — change line 3:

  ```julia
  using SigmaTauBase: StabilityResult
  ```

  to:

  ```julia
  using SigmaTau: StabilityResult
  ```

  Nothing else in that file changes.

- [ ] **Step B4.5: Verify no stragglers.**

```bash
grep -rn "SigmaTauBase\|SigmaTauStability\|SigmaTauEnsemble" src/ ext/
```

Expected: zero matches. If any remain, edit per case.

### Task B5 — Rewrite root `Project.toml`

**Files:** overwrite `Project.toml`.

Merging policy:
- Keep the existing umbrella UUID `1434c4cc-a461-4fa3-9ce6-cb17a59e9a11`
  (it has not been registered yet, so we don't need a fresh UUID; the
  user's procedure said "generate a fresh UUIDv4" but reusing the
  existing one preserves the dev-path environments in his side projects).
  → flag this as a per-user choice in Phase A approval.
- Drop all three `SigmaTau{Base,Stability,Ensemble}` lines from `[deps]`
  and the matching `[sources.X]` blocks.
- Drop `[workspace]`.
- Union the three child `[deps]` sections into the single `[deps]`.
- Promote each subpackage's `[compat]` entry into the merged `[compat]`.
  Bump `julia` to `"1.11"` (was `"1.9"` at the umbrella; CLAUDE.md says
  the workspace targets 1.11).
- Move `Plots` from `[deps]` to `[weakdeps]`. *Reason:* the existing
  `SigmaTauStability/Project.toml` carries `Plots` as a hard `[deps]`
  entry, but plot recipes live in the package extension, so a hard
  Plots dep at the umbrella level would defeat the extension. Plots
  was never used inside `SigmaTauStability` source, only in the
  extension. → keep `RecipesBase` as the weakdep that triggers the
  extension; no `Plots` weakdep needed.
- `[extras]`/`[targets]`: union of the two test targets — `FFTW`,
  `Random`, `Statistics`, `Test`.

- [ ] **Step B5.1: Write the new `Project.toml`.**

```toml
name = "SigmaTau"
uuid = "1434c4cc-a461-4fa3-9ce6-cb17a59e9a11"
authors = ["Ian Lapinski <ianlapinski01@gmail.com>"]
version = "0.1.0"

[deps]
AbstractFFTs = "621f4979-c628-5d54-868e-fcf4e3e8185c"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
DocStringExtensions = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[weakdeps]
RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[extensions]
SigmaTauRecipesBaseExt = "RecipesBase"

[compat]
AbstractFFTs = "1"
Distributions = "0.25.125"
DocStringExtensions = "0.9"
RecipesBase = "1"
StaticArrays = "1"
julia = "1.11"

[extras]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["FFTW", "Random", "Test"]
```

> Notes:
> - `Statistics` is a stdlib (in 1.11) — no `[compat]` line needed for
>   it. `Test`, `Random`, `LinearAlgebra` likewise.
> - `Statistics` was in `[targets].test` for `SigmaTauEnsemble` because
>   that subpkg didn't list it as a hard dep; now that the merged
>   package depends on it (via `Stab`), it doesn't need to be in extras.

- [ ] **Step B5.2: Delete the obsolete `Manifest.toml` at repo root.**

```bash
rm -f Manifest.toml
```

Manifest.toml is gitignored (per `.gitignore` line `Manifest.toml`) and
will be regenerated on `Pkg.instantiate` / `Pkg.test`.

### Task B6 — Reorganize tests

- [ ] **Step B6.1: Update `test/stab/runtests.jl` imports.**

In the copied file, change the header:

```julia
using Test
using FFTW                                  # AbstractFFTs backend for noise synth
using SigmaTauBase
using SigmaTauStability
using SigmaTauStability: NEFF_RELIABLE, _gen_powerlaw_phase
```

to:

```julia
using Test
using FFTW                                  # AbstractFFTs backend for noise synth
using SigmaTau
using SigmaTau.Stab: NEFF_RELIABLE, _gen_powerlaw_phase
```

(`using SigmaTau` brings in `PhaseData`, `FrequencyData`, every `_*_core`,
`adev`, etc. via the umbrella's `using .Stab`.)

The `include("legacy_kernels.jl")` line stays — relative path is
unchanged because both files moved together into `test/stab/`.

- [ ] **Step B6.2: Update `test/est/runtests.jl` imports + legacy path.**

Change:

```julia
using Test
using LinearAlgebra
using Statistics
using StaticArrays
using SigmaTauEnsemble

const LEGACY_DIR     = joinpath(@__DIR__, "..", "..", "..", "legacy", "julia", "src")
```

to:

```julia
using Test
using LinearAlgebra
using Statistics
using StaticArrays
using SigmaTau

const LEGACY_DIR     = joinpath(@__DIR__, "..", "..", "legacy", "julia", "src")
```

(One `..` fewer because `test/est/` is two ancestors deep from repo root,
not three like `lib/SigmaTauEnsemble/test/`.)

If the test body references unqualified Est names like `predict!`,
`StandardKalmanFilter`, etc., they continue to resolve via
`using SigmaTau`. If any test references non-exported Est internals,
add a `using SigmaTau.Est: name1, name2` line — verify by grep:

```bash
grep -nE "SigmaTauEnsemble:" test/est/runtests.jl
```

If matches exist, rewrite each to `SigmaTau.Est:`.

- [ ] **Step B6.3: Author `test/types/runtests.jl` from scratch.**

```julia
using Test
using SigmaTau

@testset "Shared types" begin
    @testset "PhaseData" begin
        p = PhaseData([1.0, 2.0, 3.0], 1.0)
        @test p.x == [1.0, 2.0, 3.0]
        @test p.tau0 == 1.0
        @test p isa AbstractTimingData
    end

    @testset "FrequencyData" begin
        f = FrequencyData([0.1, 0.2], 0.5)
        @test f.y == [0.1, 0.2]
        @test f.tau0 == 0.5
        @test f isa AbstractTimingData
    end

    @testset "StabilityResult fields" begin
        r = StabilityResult(:adev, [1.0], [0.5], Symbol[], Float64[], Float64[], Float64[])
        @test r.deviation_type === :adev
        @test r.tau == [1.0]
        @test r.dev == [0.5]
        @test isempty(r.noise_type)
        @test isempty(r.edf)
    end
end
```

- [ ] **Step B6.4: Author the top-level `test/runtests.jl`.**

```julia
using Test

@testset "SigmaTau" begin
    @testset "types"  begin include("types/runtests.jl") end
    @testset "stab"   begin include("stab/runtests.jl")  end
    @testset "est"    begin include("est/runtests.jl")   end
end
```

### Task B7 — Update docs

- [ ] **Step B7.1: Rewrite `docs/Project.toml`.**

Old:

```toml
[deps]
Documenter = "..."
DocumenterCitations = "..."
SigmaTau = "1434c4cc-..."
SigmaTauBase = "..."
SigmaTauEnsemble = "..."
SigmaTauStability = "..."

[sources.SigmaTau]         path = ".."
[sources.SigmaTauBase]     path = "../lib/SigmaTauBase"
[sources.SigmaTauEnsemble] path = "../lib/SigmaTauEnsemble"
[sources.SigmaTauStability] path = "../lib/SigmaTauStability"
```

New:

```toml
[deps]
Documenter = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
DocumenterCitations = "daee34ce-89f3-4625-b898-19384cb65244"
SigmaTau = "1434c4cc-a461-4fa3-9ce6-cb17a59e9a11"

[sources.SigmaTau]
path = ".."

[compat]
Documenter = "1"
DocumenterCitations = "1"
julia = "1.11"
```

- [ ] **Step B7.2: Rewrite `docs/make.jl`.**

Replace the three `using` lines and three `DocMeta.setdocmeta!` lines:

```julia
using Documenter
using DocumenterCitations
using SigmaTau, SigmaTauBase, SigmaTauStability, SigmaTauEnsemble
…
DocMeta.setdocmeta!(SigmaTauBase,      :DocTestSetup, :(using SigmaTau); recursive=true)
DocMeta.setdocmeta!(SigmaTauStability, :DocTestSetup, :(using SigmaTau); recursive=true)
DocMeta.setdocmeta!(SigmaTauEnsemble,  :DocTestSetup, :(using SigmaTau); recursive=true)
```

with:

```julia
using Documenter
using DocumenterCitations
using SigmaTau
…
DocMeta.setdocmeta!(SigmaTau, :DocTestSetup, :(using SigmaTau); recursive=true)
```

And update the `modules` keyword to `makedocs`:

```julia
modules = [SigmaTau, SigmaTau.Stab, SigmaTau.Est],
```

(Documenter walks submodules with `recursive=true`, but listing them
explicitly catches `@docs`-block resolution failures earlier.)

The `pages` tree is left structurally alone, but the API Reference
filenames are renamed — see B7.3.

- [ ] **Step B7.3: Rename and rewrite the three reference pages.**

```bash
git mv docs/src/reference/base.md      docs/src/reference/types.md
git mv docs/src/reference/stability.md docs/src/reference/stab.md
git mv docs/src/reference/ensemble.md  docs/src/reference/est.md
```

Then in `docs/make.jl`, update the `pages` block:

```julia
"API Reference"   => [
    "reference/types.md",
    "reference/stab.md",
    "reference/est.md",
],
```

Inside each renamed page, the H1 title and `@docs` symbol lists stay
unqualified — Documenter resolves unqualified names against any module
in `modules=…`. Edit only the H1:

- `docs/src/reference/types.md` H1: `# Shared types`
- `docs/src/reference/stab.md` H1: `# SigmaTau.Stab — Stability analysis`
- `docs/src/reference/est.md` H1: `# SigmaTau.Est — Estimation`

Symbol lists in `@docs` blocks remain unqualified (`adev`, `mdev`, …).
Documenter accepts unqualified names when the symbol resolves uniquely
across the listed modules — verify after the docs build.

If any symbol fails to resolve (Documenter will warn:
`"could not find docs for adev"`), qualify it explicitly:

```
@docs
SigmaTau.Stab.adev
```

But try unqualified first; it usually works.

- [ ] **Step B7.4: Append "Restructure (2026-05-09)" entry to `docs/_pass4_apply_log.md`.**

Append at end of file:

```markdown

---

## Restructure — 2026-05-09

- Workspace `lib/{SigmaTauBase, SigmaTauStability, SigmaTauEnsemble}`
  collapsed into single package `SigmaTau` with submodules
  `SigmaTau.Stab` (was `SigmaTauStability`) and `SigmaTau.Est` (was
  `SigmaTauEnsemble`). Shared types live at the top level.
- `docs/src/reference/{base,stability,ensemble}.md` renamed to
  `{types,stab,est}.md`; `@docs` symbol lists left unqualified.
- `docs/make.jl` reduced to a single `using SigmaTau` and a single
  `DocMeta.setdocmeta!` call with `recursive=true`.
- All Code-link YAML metadata in `legdocs/.../Code links/` updated to
  reflect new module names and source paths. Note bodies untouched.
- `lib/` renamed to `lib.bak/` and gitignored as a recovery point.
```

### Task B8 — Update vault Code-link YAML

47 files in `legdocs/obsidian/SigmaTauVault/Code links/{base,stability,ensemble}/`
have YAML headers like:

```yaml
---
function: "adev"
module: "SigmaTauStability"
file: "lib/SigmaTauStability/src/api/allan.jl"
exported: true
subpackage: "stability"
type: code-link
tags: [code-link]
---
```

Each header has three rewrites:
1. `module:` value
2. `file:` value (path now under `src/`)
3. `subpackage:` value

Also rename the three folders.

- [ ] **Step B8.1: Rename the three folders.**

```bash
cd legdocs/obsidian/SigmaTauVault
git mv "Code links/base"      "Code links/types"
git mv "Code links/stability" "Code links/stab"
git mv "Code links/ensemble"  "Code links/est"
cd -
```

- [ ] **Step B8.2: Mass-update YAML in the renamed folders.**

For files in `Code links/types/`:
- `module: "SigmaTauBase"` → `module: "SigmaTau"`
- `file: "lib/SigmaTauBase/src/SigmaTauBase.jl"` →
  one of: `"src/types/abstract.jl"`, `"src/types/phase_data.jl"`,
  `"src/types/frequency_data.jl"`, `"src/types/stability_result.jl"`
  (per the symbol — match by filename of the note,
  e.g. `PhaseData.md` → `phase_data.jl`).
- `subpackage: "base"` → `subpackage: "types"`

For files in `Code links/stab/`:
- `module: "SigmaTauStability"` → `module: "SigmaTau.Stab"`
- `file:` path: replace `"lib/SigmaTauStability/src/"` prefix with
  `"src/stab/"`.
- `subpackage: "stability"` → `subpackage: "stab"`

For files in `Code links/est/`:
- `module: "SigmaTauEnsemble"` → `module: "SigmaTau.Est"`
- `file:` path: replace `"lib/SigmaTauEnsemble/src/"` prefix with
  `"src/est/"`.
- `subpackage: "ensemble"` → `subpackage: "est"`

Also fix any prose mention of the old module name in **headers and
"Source" lines** (these were grepped at plan time):

- "exported from module `SigmaTauStability`" → "exported from module `SigmaTau.Stab`"
- "exported from module `SigmaTauEnsemble`" → "exported from module `SigmaTau.Est`"
- "exported from module `SigmaTauBase`, re-exported by the umbrella `SigmaTau`" →
  "exported from `SigmaTau`"

Use a per-file edit (Edit tool with `replace_all=true`) rather than
sed across all files at once — the prose changes are context-dependent
and easier to verify per-file.

Suggested execution order:
1. Run `find "legdocs/obsidian/SigmaTauVault/Code links/types"  -name '*.md' | xargs grep -l SigmaTauBase` (and likewise for stab/est) to enumerate.
2. For each file, apply the three YAML rewrites + any prose mentions in the "## Source" line near the top.
3. After all files: re-grep to confirm zero matches remain (Phase C step C4).

### Task B9 — Update README.md

- [ ] **Step B9.1: Replace the "umbrella package that re-exports three…"
  paragraph and the subpackage table.**

Old (top of README):

```markdown
# SigmaTau.jl

Frequency-stability analysis and clock-ensemble estimation in Julia.

`SigmaTau` is an umbrella package that re-exports three focused subpackages:

| Subpackage | Purpose |
|---|---|
| [`SigmaTauBase`](lib/SigmaTauBase) | Core types … |
| [`SigmaTauStability`](lib/SigmaTauStability) | Stability deviations … |
| [`SigmaTauEnsemble`](lib/SigmaTauEnsemble) | Clock state-space models … |
```

New:

```markdown
# SigmaTau.jl

Frequency-stability analysis and clock-ensemble estimation in Julia.

`SigmaTau` ships shared types at the top level (`PhaseData`,
`FrequencyData`, `StabilityResult`) and two thematic submodules:

| Submodule | Purpose |
|---|---|
| `SigmaTau.Stab` | Allan / Hadamard / Total deviations, lag-1 noise ID, EDF + χ² confidence intervals, bias correction. |
| `SigmaTau.Est`  | Clock state-space models (`TwoStateClock`, `ThreeStateClock`, `RelativisticClock`) and Kalman estimators. AD-friendly out-of-place math via `StaticArrays`. |

All public symbols are re-exported from the umbrella, so casual code
(`using SigmaTau; adev(data, [1, 2, 4])`) is unchanged. Power-user code
can import the submodules directly:

\`\`\`julia
using SigmaTau                # adev, mdev, predict!, …
SigmaTau.Stab.adev(data, m)   # qualified
using SigmaTau.Stab           # bring stability symbols in directly
\`\`\`
```

- [ ] **Step B9.2: Update Install / Quickstart sections.**

Drop the "structured as a Julia 1.11 workspace" sentence — there is no
workspace anymore. The `pkg> add https://…` line stays. The Quickstart
code block is unchanged (`using SigmaTau; PhaseData(x, 1.0); adev(…)`).

- [ ] **Step B9.3: Update the "Running tests" section.**

Replace:

```markdown
Each subpackage has its own `test/runtests.jl` and can be tested in
isolation, e.g. `julia --project=lib/SigmaTauStability -e 'using Pkg; Pkg.test()'`.
```

with:

```markdown
Tests are organized into three subdirectories under `test/`
(`types/`, `stab/`, `est/`); `julia --project=. -e 'using Pkg; Pkg.test()'`
runs all three.
```

### Task B10 — Update CHANGELOG.md

- [ ] **Step B10.1: Add a Changed (BREAKING) entry under `## [Unreleased]`.**

```markdown
### Changed (BREAKING)

- Restructured from a three-subpackage workspace
  (`SigmaTauBase` / `SigmaTauStability` / `SigmaTauEnsemble`) into a
  single registerable package with two submodules (`SigmaTau.Stab`
  and `SigmaTau.Est`). Shared types (`PhaseData`, `FrequencyData`,
  `StabilityResult`) now live at the top level. All previous exports
  continue to be re-exported from `SigmaTau`, so user code that imported
  via `using SigmaTau` keeps working unchanged. Code that explicitly
  imported `using SigmaTauBase`, `using SigmaTauStability`, or
  `using SigmaTauEnsemble` must switch to `using SigmaTau` (or
  `using SigmaTau.Stab` / `using SigmaTau.Est` for the submodules).
- The repository no longer contains a `lib/` workspace tree.
- `docs/src/reference/{base,stability,ensemble}.md` renamed to
  `{types,stab,est}.md`.
```

Per CLAUDE.md's TODO ↔ CHANGELOG workflow:
- Audit `TODO.md` for any item that names the workspace structure.
  None today (verified before plan-write); add a note in the commit
  body if no TODO removal is needed.

- [ ] **Step B10.2: Refresh `project_overview.md`.**

The current overview lists three subpackages; switch the structural
section to describe the single-package + two-submodule layout. Keep
the per-component status matrix; relabel rows that name old module
names (`SigmaTauBase` row, etc.). Defer any wholesale rewrite —
target a minimal diff that just updates the structural language.

### Task B11 — Rename `lib/` → `lib.bak/` and gitignore it

- [ ] **Step B11.1: Move with `git mv`.**

```bash
git mv lib lib.bak
```

(Using `git mv` so git tracks it as a rename rather than delete-add;
this preserves blame history if Phase C surfaces something we need to
diff back to.)

- [ ] **Step B11.2: Append to `.gitignore`.**

Append at the end of `.gitignore`:

```text

# lib.bak/ — pre-restructure recovery point; safe to delete after
# the new structure has been validated for at least one development cycle.
lib.bak/
```

- [ ] **Step B11.3: Stage the deletion of the formerly-tracked files.**

```bash
git rm -r --cached lib.bak
```

(Removes the renamed paths from the index now that they're gitignored;
working-tree files stay on disk.)

---

## Phase C — Verify

### Task C1 — Tests pass

- [ ] **Step C1.1: Instantiate the new project.**

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Expected: resolve + precompile success. If it fails because of a
removed `[sources]` entry, recheck Project.toml.

- [ ] **Step C1.2: Run all tests.**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: every testset passes — including the legacy_kernels parity
tests at `rtol=1e-12` and the Stable32 parity (`rtol≥1e-4`). If a
testset fails because of an unqualified-import collision (likely
candidate: `step!` clash if Stab also defines one), qualify the import
in `test/est/runtests.jl` to `using SigmaTau.Est: step!`.

### Task C2 — Docs build clean

- [ ] **Step C2.1: Build docs.**

```bash
cd docs && julia --project=. make.jl
```

Expected: ≤ 2 warnings (matching the prior clean build), zero errors.
If `@docs` blocks fail with "could not find docs for X", qualify that
symbol in the reference page:

```markdown
@docs
SigmaTau.Stab.adev
```

### Task C3 — REPL smoke test

- [ ] **Step C3.1: Run a one-liner.**

```bash
julia --project=. -e '
  using SigmaTau
  @assert adev isa Function
  @assert SigmaTau.Stab isa Module
  @assert SigmaTau.Est  isa Module
  p = PhaseData(collect(1.0:100.0).^2, 1.0)
  r = adev(p, [1, 2, 4]; calc_ci=false)
  @assert r.deviation_type === :adev
  @assert length(r.dev) == 3
  println("Smoke test passed.")
'
```

Expected: `Smoke test passed.`

### Task C4 — No subpkg references leak

- [ ] **Step C4.1: Grep src, tests, docs, ext.**

```bash
grep -rn "SigmaTauStability\|SigmaTauEnsemble\|SigmaTauBase" \
  src/ ext/ test/ docs/Project.toml docs/make.jl docs/src/
```

Expected: zero matches.

- [ ] **Step C4.2: Grep vault Code-links.**

```bash
grep -rn "SigmaTauStability\|SigmaTauEnsemble\|SigmaTauBase" \
  "legdocs/obsidian/SigmaTauVault/Code links/"
```

Expected: zero matches.

### Task C5 — Append "Restructure complete" log entry

- [ ] **Step C5.1: Append to `docs/_pass4_apply_log.md`.**

```markdown
### Restructure complete — 2026-05-09

- All tests pass (`Pkg.test()` green).
- Docs build clean (warnings: <COUNT>; errors: 0).
- REPL smoke test confirmed `using SigmaTau`, `SigmaTau.Stab.adev`,
  `SigmaTau.Est`, top-level `PhaseData` all resolve.
- Vault Code-link grep returns zero subpackage references.
- `lib.bak/` retained as recovery point; gitignored.
```

(Fill in `<COUNT>` from the docs build output.)

### Task C6 — Single restructure commit

- [ ] **Step C6.1: Stage and commit.**

```bash
git add -A   # explicitly review staged files first; should NOT
             # include lib.bak/ (gitignored) or Manifest.toml
git status   # confirm
git commit -m "refactor: collapse workspace into single package with Stab/Est submodules

The three lib/ subpackages (SigmaTauBase, SigmaTauStability,
SigmaTauEnsemble) merge into a single registerable SigmaTau package.
Shared types live at the top level; clock-stability and
clock-estimation code becomes SigmaTau.Stab and SigmaTau.Est. The old
re-export contract is preserved so existing user code compiles
unchanged.

Phase B/C of docs/_restructure_plan.md."
```

> Per CLAUDE.md: do NOT add Co-authored-by trailers, "Generated with"
> footers, or any AI attribution to the commit message.

---

## What this plan does NOT do

- Does not register the package — Ian triggers `@JuliaRegistrator`
  manually after his own review of the restructured tree.
- Does not touch `SigmaTauVault/Concepts/` or `SigmaTauVault/Sources/`
  (concepts and sources are independent of Julia package structure).
- Does not modify the prose body of Code-link notes (only YAML headers
  and the "## Source" pointer line).
- Does not rebuild `docs/build/` artifacts beyond what the verification
  step requires. Goal is restructural parity, not a docs refresh.
- Does not delete `lib.bak/` — Ian deletes it manually after
  validating the new structure for one development cycle.
- Does not edit relativistic concept notes or "Planned implementation"
  callouts deferred from prior sessions.
- Does not change rtol parity contracts, Q-matrix math, or any Kalman
  semantics — purely a structural refactor.

---

## Risks and call-outs (for Phase A approval)

1. **UUID reuse vs fresh UUID.** Procedure said "generate a fresh
   UUIDv4". Plan keeps existing umbrella UUID. **Risk:** if Ian has
   ever published a 0.0.x of `SigmaTau.jl` to the General registry
   that ships the workspace layout, registering the new layout under
   the same UUID would be a backward-incompatible silent change. The
   package has not been registered (verified by absence of registry
   compat manifests), so reuse is fine — but flag for explicit OK.

2. **Sub-module name shortness.** Plan uses `Stab` and `Est` per spec.
   Some users may find `SigmaTau.Stab.adev` opaque vs. e.g.
   `SigmaTau.Stability.adev`. Spec is explicit; staying with `Stab`
   / `Est`.

3. **Vault Code-link refactor scale.** 47 files × 3 YAML rewrites + a
   prose line each = ~200 small edits. The plan applies these per-file
   to keep the diff reviewable. If any file's "## Source" prose has
   drifted from the canonical pattern, the find/replace might miss it
   — fall back to per-file edits and the C4.2 grep is the safety net.

4. **`Plots` from Stability `[deps]`.** `lib/SigmaTauStability/Project.toml`
   listed `Plots` in `[deps]` historically (compat 1.41.6). The merged
   `Project.toml` drops it because:
   - No file under `lib/SigmaTauStability/src/` does `using Plots`.
   - The recipes extension only needs `RecipesBase`.
   - Hard-Plots-dep on the umbrella defeats the recipes extension.
   If Ian relied on Plots loading transitively from `using SigmaTau`
   anywhere, that breaks. **Recommendation:** drop. Documented in the
   CHANGELOG as a side-effect of the restructure.

5. **`step!` name collision.** `Stab` does NOT export `step!` today
   (verified — it's an `Est` symbol). But `Base` defines `step` (no
   bang) used elsewhere; no conflict expected. Flag for Phase C.

6. **legacy_kernels.jl `LegacyKernels` module.** Stays as-is; it's
   namespaced inside the test file, no rename needed.

---

## Self-review

- **Spec coverage check:** every numbered procedure step in the user's
  Phase A/B/C maps to a Task above. Phase B steps 1–10 → B1–B11
  (split B5 into Project rewrite + Manifest delete; split B11 into
  rename + .gitignore + cache-purge). Phase C steps 1–5 → C1–C5
  (added C4 grep as separate task, C6 explicit commit step).
- **Placeholder scan:** no TBDs, no "implement later", every code block
  contains the actual content the executor needs to paste/run.
- **Type/symbol consistency:** submodule names (`Stab`, `Est`),
  reference-page filenames (`types.md`, `stab.md`, `est.md`), vault
  folder names (`types/`, `stab/`, `est/`), and `subpackage:` YAML
  values match across every task. Internal symbol names (`adev`,
  `_adev_core`, `_freq_to_phase`, `DEFAULT_CONFIDENCE`) preserved
  verbatim from the live source files.
- **Issues found and fixed inline:** dropped one redundant rim around
  `using Statistics` in `stab/stats/edf.jl` — kept the duplicate
  rather than adding a removal step (B4.3) because removing it risked
  semantic churn. The duplicate `using` is harmless.

---

## Execution Handoff

Plan complete. After Ian approves, two execution options:

1. **Subagent-Driven** — dispatch a fresh subagent per Phase B/C task,
   review between. Recommended only if conflicts/surprises in any task
   warrant fresh context.
2. **Inline Execution** — run all of Phase B in this session in
   order, run Phase C, report at the end. Faster — recommended given
   tasks are mostly file moves with no algorithmic risk.

**Default if not specified: Inline.**

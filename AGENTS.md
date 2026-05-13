# SigmaTau.jl — agent context

You are working on a Julia 1.11 package authored by Ian Lapinski.

## Authorship and attribution rules

- All changes are authored by Ian. Do not add yourself as a co-author or
  attribute work to "Codex", "ChatGPT", "OpenAI", or "AI" anywhere.
- Do not add "Co-authored-by: Codex" or similar trailers to commits.
- Do not add "Generated with Codex" footers, signatures, or comments to
  code, commit messages, PRs, changelogs, or docs.
- Commit messages, CHANGELOG entries, and code comments should be written
  in Ian's voice as if he wrote them directly. No first-person from you.
- Do not add comments like "// added by AI" or "# Codex: refactored this".

## Package architecture

Single Julia 1.11 package. Two submodules under `src/`, plus umbrella-level
IO that produces top-level types:

- `SigmaTau`      — umbrella module; re-exports shared types and owns
  IO (file readers, detrend, gap fill, result round-trip)
- `SigmaTau.Stab` — deviation cores + API + noise ID + EDF/CI
- `SigmaTau.Est`  — clock state-space models + Kalman filter + PID steering

Types (`PhaseData`, `FrequencyData`, `StabilityResult`) live in `src/types/`
and are re-exported from the umbrella `SigmaTau` module via `@reexport`.
Both submodules' exports are also flattened onto the umbrella, so callers
write `using SigmaTau` and get everything.

### File map (most-touched paths)

- `src/SigmaTau.jl` — umbrella module. Includes `types/`, `io/`, then
  `Stab` and `Est` submodules; `@reexport`s both. IO functions live at the
  umbrella level (not inside `Stab`) because they return top-level types.
- `src/io/` — umbrella-level IO. All files use `DelimitedFiles`, `FFTW`,
  `Statistics` imported once in `SigmaTau.jl`.
  - `read.jl`     — `read_phase`, `read_frequency` (with optional scaling,
    detrend, gap fill in one call)
  - `detrend.jl`  — `detrend(::PhaseData / ::FrequencyData)`, modes
    `:none | :mean | :endpoint | :linear`
  - `fillgaps.jl` — `fillgaps` (Howe & Schlossberger PTTI-2009 imputation)
  - `results.jl`  — `save_result`, `load_result` (self-describing
    tab-delimited round-trip, stdlib only)
- `src/stab/core/` — `_adev_core`, `_mdev_core`, etc., split by deviation family:
  - `allan.jl`    — `_adev_core`, `_mdev_core`, `_tdev_core`
  - `hadamard.jl` — `_hdev_core`, `_mhdev_core`
  - `total.jl`    — `_totdev_core`, `_mtotdev_core`, `_htotdev_core`, `_mhtotdev_core`
  - `mtie.jl`     — `_mtie_core`
  - `pdev.jl`     — `_pdev_core`
  Pure `Vector{Float64}` → array kernels.
- `src/stab/api/` — public API entry points, split by deviation family:
  - `allan.jl`, `hadamard.jl`, `total.jl`, `mtie.jl`, `pdev.jl`
  Each wraps `PhaseData`/`FrequencyData` → `StabilityResult`.
  `api/hadamard.jl` also exposes `htdev` and the deprecated alias `ldev`
  (slated for removal in a future release).
  New deviations need a `PhaseData` *and* `FrequencyData` method here.
- `src/stab/stats/edf.jl` — EDF/CI math (chi-squared, Greenhall–Riley fallbacks).
- `src/stab/noise/` — noise identification (`lag1.jl`, `synth.jl`).
- `src/stab/utils.jl` — shared helpers including `_freq_to_phase`.
- `DEFAULT_CONFIDENCE = 0.683` is defined at the top of the `Stab` submodule
  in `src/SigmaTau.jl` (not in `utils.jl`); it is the default `confidence`
  argument across every public deviation API.
- `src/est/models/clocks.jl` — `AbstractClockModel` and concrete
  `TwoStateClock`, `ThreeStateClock` (kwdef structs with `tau`, `q0`, `q1`,
  `q2`, optional `q3`); `RelativisticClock` is a stub that throws on
  method dispatch.
- `src/est/estimators/filters.jl` — `StandardKalmanFilter`
  (out-of-place, `StaticArrays`-based, AD-friendly), plus `UDFactorizedFilter`
  and `KuramotoOscillator` stubs, and `PIDController` + `step!` /
  `steer_to_correction` for closed-loop steering.
- `src/types/` — `abstract.jl`, `phase_data.jl`, `frequency_data.jl`,
  `stability_result.jl`.
- `ext/SigmaTauRecipesBaseExt.jl` — all plot recipes (loaded only when
  `RecipesBase` is available; declared in `[weakdeps]`).
- `reference/validation/` — parity fixtures. **Read-only.**
  - `stable32gen.DAT` — input data
  - `stable32out/stable32_data_full.csv` — Stable32 reference outputs (~5 sig figs, rtol ≥ 1e-4)
  - `allantools_out/allantools_data_full.csv` — allantools reference (full Float64, rtol ≈ 1e-11)
  - **Coverage gaps — no external reference exists:**
    - `mhtotdev` — not implemented in Stable32 or allantools. SigmaTau is the only library that computes it.
    - `htdev` (Hadamard time deviation) — not implemented in Stable32 or allantools. SigmaTau is the only library that computes it. The deprecated alias `ldev` resolves to the same function.
    - `mhdev` — defined in NIST SP1065 but not implemented in Stable32 or allantools. SigmaTau is (to our knowledge) the only library that actually computes it.
    Validate these three via `test/stab/legacy_kernels.jl` (MATLAB-era parity, rtol=1e-12) and internal consistency only.
- `test/runtests.jl` — root test entry point. Drives five sub-suites:
  - `test/types/runtests.jl`
  - `test/stab/runtests.jl` — includes `legacy_kernels.jl` (rtol=1e-12
    parity contract) and `allantools_cross_validation.jl` (allantools cross-checks)
  - `test/est/runtests.jl`
  - `test/io/runtests.jl` — `detrend.jl`, `fillgaps.jl`, `read.jl`
  - `test/umbrella_smoke.jl` — sanity check that re-exports resolve

### Agent-context pair

`CLAUDE.md` and `AGENTS.md` are both checked into the repo. When you change
one, mirror the change in the other. They diverge only on the agent-name in
attribution rules.

## Critical conventions — do not violate

- Core kernels (`_adev_core`, etc.) take `Vector{Float64}` and return raw
  arrays. Public API (`adev`, etc.) takes `PhaseData`/`FrequencyData` and
  returns `StabilityResult`. Never collapse these into one function.
- `StabilityResult` fields are non-parametric `Vector{Float64}`. Do not
  parameterize.
- `edf` is empty when `calc_ci=false`, populated when `calc_ci=true`.
  Preserve this contract.
- Kalman filter math is out-of-place via StaticArrays for AD-friendliness.
  In-place mutation is opt-in only.
- `legacy_compat=true` reproduces MATLAB-era `safe_sqrt` diagonal clamping
  bit-for-bit. Do not "fix" the clamp — it is a parity contract.
- Plot recipes live ONLY in `ext/SigmaTauRecipesBaseExt.jl`. Do not add
  plotting code to `src/SigmaTau.jl`.

## Verification standards

- Deviation kernels are cross-checked against Stable32 fixtures in
  `reference/validation/`. After any core kernel change, run those.
- Reference math: NIST SP1065 (Riley & Howe), Greenhall & Riley 2003,
  IEEE 1139-2022. Do not invent new χ² formulas or EDF expressions —
  cite the source.

## Development workflow — use Revise.jl

A persistent Julia REPL is available with Revise.jl loaded. Do not spawn
fresh `julia -e` invocations for verification when a REPL is available —
those pay the full JIT compilation cost (~30-60s) every time. Revise
hot-patches changes in ~100ms.

After editing a file:
1. Save the file. Revise picks up the change automatically.
2. Re-run the relevant function in the REPL to verify.
3. If you have sandboxed/cloud execution and no persistent REPL,
   batch your verification: make all related edits, then run a single
   test command at the end.

Revise CANNOT hot-patch the following — when you make these changes,
explicitly tell Ian to restart the Julia REPL:

- Adding/removing/reordering fields in a struct
- Changing a struct's type parameters
- Changes to any `Project.toml` or `Manifest.toml`
- New `@eval`'d definitions or some macro changes

## Testing

- Run all tests:              `julia --project=. -e 'using Pkg; Pkg.test()'`
- Inside the persistent REPL: `pkg> test` (faster, reuses session)
- `Random` is in `[extras]` — do not remove.

## Quality checks (run periodically, not every change)

- `using Aqua; Aqua.test_all(SigmaTau)` — method ambiguities, stale deps
- `using JET; report_package(SigmaTau)` — type instabilities across the full package

## TODO ↔ CHANGELOG workflow

Every shipped code change must, in the same commit:

1. Remove the matching item from `TODO.md`.
2. Add a Keep-a-Changelog entry under `## [Unreleased]` in
   `CHANGELOG.md` (terse, past tense, no marketing voice, no emoji).
3. Refresh `project_overview.md` if the change alters surface area
   (new exported function, new submodule role, new estimator).

If a change does not warrant a TODO/CHANGELOG entry (pure docs,
typo fixes), say so in the commit body so it's intentional.

## Specialised review delegation

- Kalman, Q-matrix, EDF, or filter-divergence design questions →
  delegate to the `kalman-filter-expert` subagent rather than
  reasoning end-to-end inline.
- After any `_*_core` deviation kernel edit, run the
  `stable32-parity-checker` subagent before declaring the change
  done — it filters Stable32 fixture results to just the rtol delta.

## Editing rules

- Never edit `Manifest.toml` files; they are gitignored locally.
- When adding a new deviation, follow the existing core/api split and
  add both a `PhaseData` and `FrequencyData` entry point (the latter
  delegates via `_freq_to_phase`).
- Match existing surrounding code style. 4-space indent, no trailing
  whitespace, docstrings on all exported functions.

## Style for prose deliverables

CHANGELOG entries follow Keep-a-Changelog. Use the same terse, factual
voice as existing entries. Past tense, no marketing language, no emoji.
Commit messages: imperative mood, ≤72 char subject, body explains *why*.

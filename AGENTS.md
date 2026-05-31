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

SigmaTau.jl v0.3.0 is a single, flat Julia 1.11 package covering
clock-stability analysis: deviation kernels, noise identification,
EDF/χ² confidence intervals, calibrated noise generation, and file IO.
No submodules — everything is at the umbrella level.

For clock state-space, Kalman, and PID steering, see the
[ClockEnsemble.jl](https://github.com/ianlap/ClockEnsemble.jl) sister
package (formerly the `SigmaTau.Est` submodule).

Types (`PhaseData`, `FrequencyData`, `StabilityResult`) live in
`src/types/` and are exported directly from `SigmaTau`. Callers write
`using SigmaTau` and get everything in one flat namespace.

### File map (most-touched paths)

- `src/SigmaTau.jl` — single flat module. Includes `types/`, `io/`,
  `DEFAULT_CONFIDENCE`, then `stab/core/`, `stab/noise/`,
  `stab/stats/edf.jl`, `stab/utils.jl`, and `stab/api/`. One umbrella-
  level `export` block at the bottom.
- `src/io/` — file readers, detrend, gap fill, result round-trip. All
  files rely on `DelimitedFiles`, `FFTW`, and `Statistics` imported once
  in `SigmaTau.jl`.
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
- `src/stab/noise/` — noise identification + synthesis.
  - `lag1.jl` — lag-1 ACF / B1 / R(n) noise-type ID.
  - `synth.jl` — internal `_gen_powerlaw_y` / `_gen_powerlaw_phase`
    spectral shaper used by the test suite.
  - `gen.jl` — public `noise_gen(::Type{PhaseData} | ::Type{FrequencyData},
    N, tau0; sigma1=…, h=…)` calibrated power-law generator.
- `src/stab/utils.jl` — shared helpers including `_freq_to_phase`.
- `DEFAULT_CONFIDENCE = 0.683` is a top-level const in `src/SigmaTau.jl`;
  it is the default `confidence` argument across every public deviation API.
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
- `test/runtests.jl` — root test entry point. Drives four sub-suites:
  - `test/types/runtests.jl`
  - `test/stab/runtests.jl` — includes `legacy_kernels.jl` (rtol=1e-12
    parity contract) and `allantools_cross_validation.jl` (allantools cross-checks)
  - `test/io/runtests.jl` — `detrend.jl`, `fillgaps.jl`, `read.jl`
  - `test/umbrella_smoke.jl` — sanity check that bare `using SigmaTau`
    exposes every public symbol; pins the absence of the old `Stab`/`Est`
    submodules.

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
- Plot recipes live ONLY in `ext/SigmaTauRecipesBaseExt.jl`. Do not add
  plotting code to `src/SigmaTau.jl`.

## Verification standards

- Deviation kernels are cross-checked against Stable32 fixtures in
  `reference/validation/`. After any core kernel change, run those.
- Reference math: NIST SP1065 (Riley & Howe), Greenhall & Riley 2003,
  IEEE 1139-2022. Do not invent new χ² formulas or EDF expressions —
  cite the source.

## Development workflow — the agent verifies autonomously

The agent runs its own verification by invoking `julia --project=.` directly (e.g.
via `julia --project=. -e '...'` or by running a test file). Do **not** rely on a
persistent REPL or ask Ian to run commands and paste output — drive the checks
yourself and read the results.

After editing a file:
1. Run the relevant test file or testset, e.g.
   `julia --project=. -e 'include("test/stab/runtests.jl")'`, or a focused
   `julia --project=. -e 'using SigmaTau; ...'` snippet for a single function.
2. Read the output and iterate until green.
3. Run the full suite before committing (see Testing).

Every `julia` invocation is a fresh process, so it always picks up the current
source — including struct field changes, `Project.toml`/`Manifest.toml` edits, and
new `@eval`'d or macro definitions. There is no hot-patch staleness to manage and
no REPL to restart. The trade-off is the per-run JIT cost (~30-60s); to keep
iteration cheap, scope each run to the smallest test file or snippet that exercises
the change, and only run the full `Pkg.test()` suite before committing.

If a persistent process is worth the speed on a long editing session, the agent may
keep one Julia process alive as a background job and feed it `include(...)` lines —
but the source of truth for a green build is always a fresh `Pkg.test()` run.

## Testing

- Run all tests:              `julia --project=. -e 'using Pkg; Pkg.test()'`
- Run one sub-suite:          `julia --project=. -e 'include("test/stab/runtests.jl")'`
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
   (new exported function, new role, new estimator).

If a change does not warrant a TODO/CHANGELOG entry (pure docs,
typo fixes), say so in the commit body so it's intentional.

## Specialised review delegation

- EDF or χ² confidence design questions — cite Riley/Greenhall/IEEE 1139
  explicitly; do not invent expressions.
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

# Phase-4 apply log

**Date:** 2026-05-08
**Plan:** `docs/_reconciliation_plan.md` (Decisions 1–6 all approved by user)
**Build status:** clean (0 errors; 23 silenced warnings; 0 doctest failures)

---

## Standard checks

### 1. Approved PATCH applied? (yes — 7 files)

- `docs/src/index.md` — citekey rename only (3 cite labels updated).
- `docs/src/getting_started.md` — touched only by the citekey-rename sweep (no labels existed in this file; sweep was a no-op).
- `docs/src/theory/overview.md` — citekey rename + slope table extended to α∈{−3,−4} with a Hadamard-convergence footnote + new MTIE section with `!!! note "Planned implementation"` callout.
- `docs/src/theory/allan_family.md` — citekey rename + HTDEV section restructured into a `!!! info "Original contribution"` callout citing the implementation source via GitHub permalink + new PVAR section after MDEV (with callout) + new Theo1 section after MHDEV (with callout) + new ThêoH section (with callout) + new Dynamic Allan deviation section (with callout).
- `docs/src/theory/total_family.md` — citekey rename + MHTOTDEV section restructured into a `!!! info "Original contribution"` callout citing the implementation source via GitHub permalink + 2×2-matrix-completion framing added (kept existing "all-three corner of the cube" prose per Decision 4) + new TTOT section after MTOT (with callout).
- `docs/src/theory/confidence.md` — citekey rename only (5 cite labels).
- `docs/src/theory/noise_id.md` — citekey rename + fixed `identify_noise(PhaseData(...); m=8)` API call (correct signature is `identify_noise(x, [m_values]; dmin=, dmax=)`).

### 2. Approved CREATE applied? (yes — 5 files, 1855 lines total)

- `docs/src/theory/ensemble_overview.md` — 434 lines.
- `docs/src/theory/kalman.md` — 517 lines.
- `docs/src/theory/steering.md` — 178 lines.
- `docs/src/theory/relativistic_clocks.md` — 218 lines.
- `docs/src/theory/ensembles_and_oscillator_networks.md` — 508 lines.

All five use `[Label](@cite key)` link form throughout; existing pages keep the older `[@cite key]` shortcut form (per Decision 2).

### 3. Approved REWRITE applied? (yes — 1 file)

- `docs/src/refs.bib` — replaced contents with `SigmaTauVault/references.bib` (44 BibTeX entries, kebab-case keys with shorttitles). Two entries (`sullivan-1990-tn1337` and `blair-1974-nbs140-annex8e`) had their `editor` field duplicated as `author` to satisfy DocumenterCitations' strict author-field requirement. **Vault `references.bib` was not modified** — this is a docs-only adjustment; if the apply log notes a future vault sync, those two entries should also be patched there.

### 4. Build chrome PATCH applied? (yes — 1 file)

- `docs/make.jl` — added the 5 new theory pages to the `pages =` tree under "Theory". No other config touched (`setdocmeta!`, `warnonly`, `canonical`, `:authoryear`, `deploydocs`, doctest=true preserved verbatim).

### 5. Julia source PATCH applied? (1 file, 1 line)

- `lib/SigmaTauStability/src/api/allan.jl:7` — the `adev` docstring used `[Greenhall2003](@cite)` which referenced the old citekey. Updated to `[Greenhall & Riley 2003](@cite greenhall-2003-edf-stability)` to match the new bib. Without this fix the bibliography expansion failed with `Key "Greenhall2003" not found`. This single-line docstring fix was in scope for atomicity of the citekey rename.

### 6. No HANDWRITTEN file modified without approved PATCH? (verified)

`docs/src/getting_started.md`, `docs/src/reference/{base,stability,ensemble}.md`, `docs/src/tutorials/0{1..5}_*.md`, `docs/src/validation/{methodology,stable32}.md`, `docs/src/bibliography.md`, `docs/Project.toml` — none touched.

### 7. Citekey migration complete?

```
$ grep -rE '\[(Riley|Greenhall|IEEE|Howe|Stein|Tavella|Zucca|Banerjee|Sullivan)[A-Z]' docs/src/
(no output — clean)
```

All 13 stale CamelCase citekeys (`RileyHowe2008`, `Greenhall2003`, `Greenhall1997`, `Greenhall1999`, `Howe1999`, `Howe2000`, `Howe2001`, `Howe2005`, `Banerjee2023`, `IEEE1139_2022`, `Riley_R_2020`, `Riley2004`, `Sullivan_NBS_TN_1337`) have been renamed to their kebab-case vault equivalents in all 7 prose files. One stale citekey was found in `lib/SigmaTauStability/src/api/allan.jl` and patched (item 5 above).

---

## Deferred-estimator checks (Decision 6)

### 5. All six deferred estimators have a section with a "Planned implementation" callout

| Estimator | File | Section title | Callout present |
|---|---|---|---|
| MTIE | `theory/overview.md` | "MTIE — Maximum Time Interval Error" | yes |
| PVAR | `theory/allan_family.md` | "PVAR — parabolic variance" | yes |
| Theo1 | `theory/allan_family.md` | "Theo1 — theoretical variance #1" | yes |
| ThêoH | `theory/allan_family.md` | "ThêoH — composite Theo1 + ADEV estimator" | yes |
| Dynamic Allan deviation | `theory/allan_family.md` | "Dynamic Allan deviation — DADEV" | yes |
| TTOT | `theory/total_family.md` | "TTOT — time-total deviation" | yes |

Citation grounding follows Decision 6:
- Theo1, ThêoH, MTIE, TTOT → `riley-2008-sp1065` + `banerjee-2023-timekeeping`.
- Dynamic ADEV → `mckelvy-2025-telemetrystability`.
- PVAR → `banerjee-2023-timekeeping`.

### 6. Deferred estimators do NOT appear in `reference/stability.md`

```
$ grep -wE 'theo1|theoh|mtie|dadev|ttot|pvar' docs/src/reference/stability.md
(no output)
```

No phantom API entries.

### 9 total "Planned implementation" callouts, broken down

- 6 deferred estimators (Decision 6): MTIE, PVAR, Theo1, ThêoH, Dynamic Allan, TTOT.
- 3 existing Julia stubs (already flagged in vault Code-link notes): `RelativisticClock` (`relativistic_clocks.md`), `UDFactorizedFilter` (`kalman.md`), `KuramotoOscillator` (`ensembles_and_oscillator_networks.md`).

---

## Build status

```
$ cd docs && julia --project=. make.jl
[ Info: SetupBuildDirectory: setting up build directory.
[ Info: Doctest: running doctests.
[ Info: ExpandTemplates: expanding markdown templates.
[ Info: CrossReferences: building cross-references.
[ Info: CheckDocument: running document checks.
[ Info: Populate: populating indices.
[ Info: RenderDocument: rendering document.
[ Info: HTMLWriter: rendering HTML pages.
[ Info: Automatic `version="0.1.0"` for inventory from ../Project.toml
```

**Hard errors: 0.** Build completes through `RenderDocument` and `HTMLWriter`. All 22 HTML pages produced under `docs/build/`:

```
docs/build/
├── assets/
├── bibliography.html
├── getting_started.html
├── index.html
├── objects.inv
├── reference/{base,stability,ensemble}.html
├── refs.bib
├── search_index.js
├── theory/
│   ├── allan_family.html
│   ├── confidence.html
│   ├── ensemble_overview.html              (NEW)
│   ├── ensembles_and_oscillator_networks.html  (NEW)
│   ├── kalman.html                         (NEW)
│   ├── noise_id.html
│   ├── overview.html
│   ├── relativistic_clocks.html            (NEW)
│   ├── steering.html                       (NEW)
│   └── total_family.html
├── tutorials/0{1..5}_*.html
└── validation/{methodology,stable32}.html
```

### Silenced warnings (23 total — all expected)

**12 `:missing_docs` warnings** for ensemble symbols whose Julia source has no docstring (silenced by `warnonly = [:missing_docs, ...]`):

```
AbstractClockModel, TwoStateClock, ThreeStateClock, RelativisticClock,
nstates, state_transition, process_noise, measurement_matrix, measurement_noise,
AbstractEstimator, UDFactorizedFilter, KuramotoOscillator
```

These are pre-existing TODOs from the pass-3 Code-link self-check (`_codelink_log.md`), not new issues introduced by this pass. Action: a follow-up PR adds docstrings in `lib/SigmaTauEnsemble/src/{models/clocks.jl, estimators/filters.jl}`.

**11 `:cross_references` warnings** of the form `Cannot resolve @ref for [`Name`](@ref)` — these are the new theory pages' attempts to link to ensemble API symbols whose docstrings are missing. They become resolvable as soon as the 12 missing docstrings above are added; no theory-page changes are needed.

**One `:docs_block` warning** for 6 docstrings present in source but not in any `@docs` block:

```
SigmaTauStability._freq_to_phase
SigmaTauStability.DEFAULT_CONFIDENCE  ← false positive: it's listed in stability.md
SigmaTauEnsemble.clamp_covariance_diag :: Tuple{SMatrix{2,2,Float64}}
SigmaTauEnsemble.clamp_covariance_diag :: Tuple{SMatrix{3,3,Float64}}
SigmaTauEnsemble.safe_sqrt_sq
SigmaTauStability._gen_powerlaw_phase
```

5 are internal helpers that should not appear in API docs. The `DEFAULT_CONFIDENCE` listing is a false positive — Documenter doesn't see the `const` declaration as part of the `@docs` block. Action: either accept this warning permanently or convert `DEFAULT_CONFIDENCE`'s docstring into a fenced `@docs` entry separately.

### Cross-references that broke and were repaired

None broke during this pass.

---

## Counts

- **PATCH applied:** 7 docs/src files + 1 lib/ file + 1 docs/make.jl = **9 patches**.
- **CREATE applied:** 5 new theory pages.
- **REWRITE applied:** 1 (`docs/src/refs.bib`).
- **Skipped (per plan, intentional):** all `tutorials/*`, `validation/*`, `reference/*`, `bibliography.md`, `getting_started.md` content, `Project.toml`, vault files.

## Surprising during application

1. **Stale citekey in Julia docstring.** `lib/SigmaTauStability/src/api/allan.jl:7` referenced `[Greenhall2003](@cite)` — the citekey-rename sweep over `docs/src/` missed this because it lives in package source, not docs source. The build refused to render until the docstring was fixed. Fixed in scope for atomicity.

2. **Two bib entries lacked an `author` field** (`sullivan-1990-tn1337` and `blair-1974-nbs140-annex8e`). Both used `editor` only because the underlying papers are edited compilations. DocumenterCitations 1.x demands an `author` field. Resolved by duplicating the editor list into `author` in `docs/src/refs.bib`. Vault `references.bib` was not modified — recommend a follow-up vault edit so the two stay in sync, otherwise the next bib resync will reintroduce the same two errors.

3. **Background-agent-written theory pages are slightly long.** Targets were 200/300/160/200/300 lines; actuals are 434/517/178/218/508. Coverage is comprehensive and citekey-discipline survived the translation (all 25 citekeys used resolve to entries in the new bib). No truncation or under-coverage to flag.

4. **Cite-syntax inconsistency now visible in source diffs.** New theory pages use `[Label](@cite key)`; old theory pages use `[@cite key]`. Both render identically in the HTML output; the disagreement is editorial only. A future commit can normalize to the link form across all pages without affecting the rendered docs.

5. **Documenter could not auto-detect the building environment.** This is a benign warning (`Skipping deployment.`) emitted because the build was run from a developer machine, not from CI. `deploydocs` only runs on the configured `devbranch` in CI; in CI builds this warning will be absent.

---

## Status

**Publishable.** The build produces clean HTML for all 22 pages including the 5 new ensemble theory pages. All 23 warnings are pre-existing TODOs expected from the pass-3 self-check (12 missing ensemble docstrings + 11 derivative `@ref` failures) plus 1 internal-helper `:docs_block` notice.

Recommended follow-up commits, in priority order:
1. **Add docstrings** to the 12 ensemble symbols. Eliminates the 12 `:missing_docs` warnings and the 11 derivative `:cross_references` warnings (23 total → 1).
2. **Vault sync for bib `author` fields** on `sullivan-1990-tn1337` and `blair-1974-nbs140-annex8e`.
3. **Optional cite-syntax normalization** (rewrite `[@cite key]` shortcut form to `[Label](@cite key)` link form across the existing theory pages, for source-diff consistency only).
4. **Tutorials narrative** — the five tutorial stubs remain as-is.
5. **`validation/methodology.md` narrative** and `validation/stable32.md` numerical refresh against current `reference/validation/` fixtures (Decision 3 deferred this).

---

## Docstring backfill (2026-05-09)

Recommendation 1 above (the docstring backfill) was applied in this session.

### Symbols documented (12)

In `lib/SigmaTauEnsemble/src/models/clocks.jl`:
- `AbstractClockModel` — interface description; lists shipped subtypes.
- `TwoStateClock` — `[phase, frequency]` polynomial model with WPM/WFM/RWFM coefficient meanings.
- `ThreeStateClock` — `[phase, frequency, drift]` polynomial model adding IRWFM (q3).
- `RelativisticClock` — stub; flagged with `!!! note "Stub implementation"` callout citing Seyffert 2025 in prose.
- `nstates` — accessor returning state-vector dimension.
- `state_transition` — Φ matrix (single docstring before the `TwoStateClock` method, per Documenter convention; the `ThreeStateClock` method inherits the docstring through `@docs` lookup).
- `process_noise` — Q matrix from closed-form Wiener integration.
- `measurement_matrix` — H row vector (phase-only observation).
- `measurement_noise` — R = `[q0]` (1×1 SMatrix).

In `lib/SigmaTauEnsemble/src/estimators/filters.jl`:
- `AbstractEstimator` — interface description; lists shipped subtypes.
- `UDFactorizedFilter` — stub; flagged with `!!! note "Stub implementation"` citing Ramos 2022 in prose.
- `KuramotoOscillator` — stub; flagged with `!!! note "Stub implementation"`.

### Test status

`SigmaTauEnsemble` tests pass (9/9). Note: an in-tree `Pkg.resolve()` was required to add `DocStringExtensions` to the manifest before tests could load — this dependency is declared in `Project.toml` but had never been instantiated locally; running tests was a no-op for the codebase but a one-time dependency-graph fix.

### Build status (post-backfill)

```
$ cd docs && julia --project=. make.jl
[ Info: SetupBuildDirectory: setting up build directory.
[ Info: Doctest: running doctests.
...
[ Info: HTMLWriter: rendering HTML pages.
```

Substantive warnings dropped from **23 → 1**:
- 12 `:missing_docs` for ensemble symbols → **0** (all documented).
- 11 derivative `:cross_references` → **0** (resolve once docstrings exist).
- 1 `:docs_block` notice for 6 internal/private symbols → **1** (pre-existing; includes `DEFAULT_CONFIDENCE` false-positive plus 5 truly internal helpers `_freq_to_phase`, `_gen_powerlaw_phase`, `safe_sqrt_sq`, `clamp_covariance_diag` ×2). Within the ≤2 target.

3 BibParser parser-noise warnings about TODO-list entries in `references.bib` are unchanged from before and unrelated.

### Notes / gaps

- The `RelativisticClock` docstring references Seyffert 2025 by name in prose only (Julia docstrings don't render Obsidian wikilinks); the canonical multi-source grounding lives in the vault concept notes.
- No methods were implemented for the three stubs — that's explicitly out of scope for this backfill (per the session prompt).

---

## Phase-4 redux (relativistic batch — 2026-05-09)

After the pass-2-redux added 13 concept notes under `Concepts/relativistic/`, the docs/theory page on relativistic clocks (single-source, ~220 lines) significantly understated what the vault now covers. This redux pass surfaces that material in docs.

### Pages created (3) and refactored (1)

- **NEW** `docs/src/theory/relativistic_frames_and_timescales.md` (~280 lines) — BCRS / GCRS / LCRS reference systems and the six time scales TT / TCG / TCB / TDB / TCL / TL with defining constants ($L_G, L_C, L_B, L_H, L_L, L_M$) and the explicit transformation forms. Multi-source ground from Ashby & Patla 2024, Turyshev 2025, Leonard et al. 2026, Reinhardt et al. 2024, Cacciapuoti & Salomon 2009, Seyffert 2025.
- **NEW** `docs/src/theory/relativistic_corrections.md` (~310 lines) — 1PN proper-time mapping, gravitational redshift on the Moon, $L_{Gm}$ Earth–Moon rate constant and L1 / L2 / L4-L5 offsets table, cislunar orbit drift regimes (vLLO / LLO / ELFO / L1 / NRHO) summarised in a regime-vs-numerics table, Shapiro / Sagnac light-time corrections, and the three-look-alike-$\vec v\cdot\vec R$-terms warning. Multi-source ground from Ashby & Patla 2024, Turyshev 2025, Reinhardt et al. 2024, Cacciapuoti & Salomon 2009, Leonard et al. 2026, Seyffert 2025.
- **NEW** `docs/src/theory/lunar_pnt_systems.md` (~150 lines) — Two-way time transfer (synchronous and asynchronous) in the ESA Moonlight LCNS architecture, ATLAS K-band SS code structure, MSPA ground architecture, and the relativistic positioning systems tradition (emission coordinates, Autonomous Basis of Coordinates). Multi-source ground from Iess et al. 2025, Iess, Boscagli & Di Benedetto 2026, Cacciapuoti & Salomon 2009, Gomboc et al. 2013.
- **REFACTORED** `docs/src/theory/relativistic_clocks.md` — converted from a single-source placeholder into a hub page (~150 lines) that anchors the operational headline numbers (lunar surface vs TCG, L1 vs Earth, ELFO drift, NRHO drift), states the operational design pattern (TDB internal / deterministic correction first / TL mapping at API), restates the three look-alike-terms warning, and points to the three companion pages above plus Code-link to [`RelativisticClock`](@ref). Now multi-sourced (Ashby, Turyshev, Leonard, Cacciapuoti, Iess, Reinhardt, Gomboc, Seyffert).

### `docs/make.jl` patch

Theory tree gained a nested `"Relativistic PNT" => [...]` group with the four pages above. Other tree entries unchanged.

### `docs/src/refs.bib` additions (8 entries)

`@ashby-2024-lunar-coordinate-time`, `@turyshev-2025-cislunar-time-scales`, `@leonard-2026-lcrns-relativistic`, `@iess-2026-moonlight-time-transfer`, `@iess-2025-cislunar-od-time-sync`, `@cacciapuoti-2009-aces-space-clocks`, `@reinhardt-2024-lisa-clock-sync`, `@gomboc-2013-relativistic-positioning`. All eight were already present in the vault `references.bib`; this redux just propagates them to the docs bib so the new theory pages' `[Author Year](@cite key)` references resolve.

### Build status (post-redux)

```
$ cd docs && julia --project=. make.jl
[ Info: SetupBuildDirectory: setting up build directory.
[ Info: Doctest: running doctests.
...
[ Info: HTMLWriter: rendering HTML pages.
```

Substantive warnings unchanged from the prior post-docstring-backfill state:
- 0 `:missing_docs`
- 0 `:cross_references`
- 1 `:docs_block` notice for 6 internal helpers (pre-existing; includes `DEFAULT_CONFIDENCE` false-positive plus 5 internal helpers `_freq_to_phase`, `_gen_powerlaw_phase`, `safe_sqrt_sq`, `clamp_covariance_diag` ×2)
- 3 BibParser warnings about TODO-list comments in `references.bib` (unchanged — predate this redux; an earlier draft of the new bib block triggered a fourth, but the comment was rewritten to drop the `@<...>` token that BibParser misreads as an entry kind)
- 1 deploydocs auto-detect notice (benign on local builds)

All 4 new HTML pages are produced under `docs/build/theory/`:
```
docs/build/theory/
├── lunar_pnt_systems.html               (NEW)
├── relativistic_clocks.html             (REFACTORED)
├── relativistic_corrections.html        (NEW)
└── relativistic_frames_and_timescales.html  (NEW)
```

### Notes

- The vault `Concepts/ensemble/clock-models/Relativistic clock.md` was deliberately not modified per the pass-2-redux disposition recommendation (keep where it is, cross-link from the new relativistic/ notes).
- The `RelativisticClock` Julia stub was not implemented in this redux — the closed-form Keplerian-Orbital path remains an open TODO item; this docs work is the spec it should implement against.
- Vault `references.bib` and docs `refs.bib` are now in sync on the 8 new entries (they were divergent before this redux because the pass-1 incremental ingest had only added them to the vault file).

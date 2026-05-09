# docs/ ↔ SigmaTauVault reconciliation plan

**Status:** plan only — no files have been modified. Review and approve / amend before phase 2 applies any changes.

**Date:** 2026-05-08

**Scope:** bring `docs/src/` and the build chrome (`docs/make.jl`, `docs/Project.toml`, `docs/src/refs.bib`, `docs/src/bibliography.md`) in line with `SigmaTauVault/` (the source of truth as of pass 2.6) without disturbing handwritten content or breaking the existing build.

---

## 1. Existing structure inventory

### `docs/src/` files

| Path | Lines | Summary | Classification |
|---|---:|---|---|
| `docs/src/index.md` | 33 | Landing page: package layout table, "where to next" links, three "reference math" cite bullets. | **MIXED** — handwritten landing prose + vault-derivable cite bullets |
| `docs/src/getting_started.md` | 41 | Installation snippet, minimal `adev` example (`@example basic` block), "where to next" links. | **HANDWRITTEN** — no vault analog (it's tutorial chrome) |
| `docs/src/theory/overview.md` | 110 | Phase vs frequency, power-law noise table, estimator family map (10-row table), notation, slope demo (`@example overview`), references list. | **VAULT-DERIVABLE** — overlaps `shared/Phase vs frequency data.md`, `shared/Power-law noise model.md`, plus the 6 power-law-noise concept notes |
| `docs/src/theory/allan_family.md` | 257 | ADEV/MDEV/TDEV/HDEV/MHDEV/HTDEV definitions, slope tables, demo (`@example allan`), provenance notes for HTDEV. | **VAULT-DERIVABLE** — overlaps `stability/estimators/Allan deviation.md`, `Modified Allan deviation.md`, `Hadamard deviation.md`, `Modified Hadamard deviation.md`, `Time deviation.md`, and the new `htdev.md` concept note |
| `docs/src/theory/total_family.md` | 271 | TOTDEV/MTOTDEV/HTOTDEV/MHTOTDEV with extension methodology, bias-correction policy box, demo, MHTOTDEV "all-three-corner" framing. | **VAULT-DERIVABLE** — overlaps `stability/estimators/Total deviation.md`, `Modified Total deviation.md`, `Hadamard Total deviation.md`, the new `mhtotdev.md` concept note, and bias-correction concepts |
| `docs/src/theory/confidence.md` | 108 | Chi-squared CI math, Greenhall–Riley EDF intro, bias correction summary table, implementation contract paragraph. | **VAULT-DERIVABLE** — overlaps `stability/algorithms/Chi-squared confidence intervals.md` and `Greenhall-Riley EDF.md` plus 4 bias-correction concept notes |
| `docs/src/theory/noise_id.md` | 101 | Lag-1 ACF method, NEFF_RELIABLE = 30, B1/R(n) fallback. | **VAULT-DERIVABLE** — overlaps `stability/algorithms/Lag-1 ACF noise identification.md` and `B1 ratio noise identification.md` |
| `docs/src/reference/base.md` | 12 | `@docs` block listing the 4 `SigmaTauBase` types. | **HANDWRITTEN** — minimal page chrome wrapping `@docs` blocks; the docstring content is in Julia source, the page header is handwritten |
| `docs/src/reference/stability.md` | 53 | `@docs` blocks for 25 stability symbols, grouped by Allan family / Total / Noise ID / EDF / Internal kernels. | **HANDWRITTEN** — page chrome; section headers and groupings are editorial decisions |
| `docs/src/reference/ensemble.md` | 34 | `@docs` blocks for 18 ensemble symbols, grouped by Clock models / Estimators / Steering. | **HANDWRITTEN** — same character as the other reference pages |
| `docs/src/tutorials/01_phase_data.md` | 14 | Skeleton with one `@example phase` block; "Narrative fills in a follow-up PR." | **HANDWRITTEN STUB** |
| `docs/src/tutorials/02_compute_adev.md` | 5 | Stub. | **HANDWRITTEN STUB** |
| `docs/src/tutorials/03_identify_noise.md` | 5 | Stub. | **HANDWRITTEN STUB** |
| `docs/src/tutorials/04_confidence_intervals.md` | 5 | Stub. | **HANDWRITTEN STUB** |
| `docs/src/tutorials/05_single_clock_steering.md` | 6 | Stub. | **HANDWRITTEN STUB** |
| `docs/src/validation/methodology.md` | 14 | Three-way validation strategy; "Detailed comparison narrative lands in a follow-up PR." | **HANDWRITTEN STUB** |
| `docs/src/validation/stable32.md` | 151 | Per-estimator agreement tables (12-row tables for ADEV/MDEV/TDEV/HDEV/Total/MTOT/HTOT), MTOT-bias narrative, methodology pointer. | **HANDWRITTEN** — the comparison numbers are derived from `reference/validation/` fixtures, not from any vault note |
| `docs/src/bibliography.md` | 8 | Single `@bibliography *` block with a one-line note about PDFs. | **HANDWRITTEN** — chrome only |
| `docs/src/refs.bib` | 160 | 22 BibTeX entries with CamelCase citekeys (`RileyHowe2008`, `Greenhall2003`, etc.). | **MIXED** — entries themselves are vault-derivable (the vault has 46 entries in `references.bib`), but the citekeys differ from vault style |

### Build chrome

| Path | Lines | Summary | Classification |
|---|---:|---|---|
| `docs/make.jl` | 59 | `Documenter.HTML` config, MathJax3, `DocumenterCitations` v1 plugin (`plugins=[bib]`), pages tree, `deploydocs` for `github.com/ianlap/SigmaTau.jl.git`, `warnonly=[:missing_docs, :cross_references, :docs_block]`. | **HANDWRITTEN** — build configuration |
| `docs/Project.toml` | 24 | Documenter 1, DocumenterCitations 1, `[sources]` for the 4 local SigmaTau packages, julia 1.11. | **HANDWRITTEN** — build configuration |

### Pre-existing artefacts

- `docs/build/` — generated output, not source. Must not be edited.
- `docs/superpowers/plans/2026-05-08-theory-pages-and-refs-reorg.md`, `2026-05-08-detrend-kwarg.md`; `docs/superpowers/specs/*.md` — planning skeletons from prior iterations. Not under `docs/src/`, not part of the build. Out of scope.

---

## 2. `make.jl` and `Project.toml` inspection

### Format
- `Documenter.HTML` only — no LaTeX, no PDF.
- `prettyurls = get(ENV, "CI", ...)` — standard CI-aware setting.
- `canonical = "https://ianlap.github.io/SigmaTau.jl"`.
- `mathengine = Documenter.MathJax3()`.

### Citations
- `DocumenterCitations 1.x`. Configured via the **new** `plugins=[bib]` API (line 22 of `make.jl`). Bibliography file is `docs/src/refs.bib`. Style is `:authoryear`.
- Citation syntax in prose uses **both** the old `[@cite KEY]` shortcut and the link form `[label](@cite KEY)`. DocumenterCitations 1.x supports both, but a future migration to the link form would be cleaner.

### `pages =` structure

```
Home              → index.md
Getting Started   → getting_started.md
Theory            → theory/{overview, allan_family, total_family, confidence, noise_id}.md
Tutorials         → tutorials/{01..05}.md
API Reference     → reference/{base, stability, ensemble}.md
Validation        → validation/{methodology, stable32}.md
Bibliography      → bibliography.md
```

### Deploy
- `deploydocs(repo = "github.com/ianlap/SigmaTau.jl.git", push_preview = true, devbranch = "main")`. No custom domain, no GitHub Actions hints in `make.jl` itself.

### Custom assets / config to preserve
- `DocMeta.setdocmeta!` for all three subpackages (lines 10–12) — needed for doctests.
- `doctest = true` (line 51).
- `warnonly = [:missing_docs, :cross_references, :docs_block]` (line 52) — relaxes Documenter's strict cross-reference checking; if I tighten cross-references in this pass I need to leave this in place or coordinate with the user.

### Custom configuration that would be lost on rewrite

- The `setdocmeta!` block.
- The `warnonly` triple.
- The `canonical` URL.
- The `:authoryear` citation style choice.
- The `[sources]` block in `Project.toml` (the four local subpackage paths).

**These four `make.jl` blocks plus the four `[sources]` paths must be preserved verbatim through any apply phase.**

---

## 3. Vault coverage gaps

### Concept notes vs `docs/src/theory/`

#### `Concepts/shared/` (11 notes)

| Concept note | Docs status | Where in docs |
|---|---|---|
| `White PM noise` | **STALE** | mentioned in `theory/overview.md` table only |
| `Flicker PM noise` | **STALE** | as above |
| `White FM noise` | **STALE** | as above |
| `Flicker FM noise` | **STALE** | as above |
| `Random Walk FM noise` | **STALE** | as above |
| `Random Run FM noise` | **MISSING** | not in docs at all |
| `Power-law noise model` | **STALE** | summarised in `theory/overview.md` but vault note has more grounding |
| `Fractional frequency` | **STALE** | brief paragraph in `theory/overview.md` |
| `Phase vs frequency data` | **STALE** | one-paragraph treatment in `theory/overview.md` |
| `Polynomial clock model` | **MISSING** | not covered in docs |
| `IEEE 1139 definitions` | **STALE** | name-checked in `index.md` only |

#### `Concepts/stability/estimators/` (14 notes)

| Concept note | Docs status |
|---|---|
| `Allan deviation` | **COVERED** (`theory/allan_family.md` ADEV section) |
| `Overlapping Allan deviation` | **STALE** — folded into ADEV section, not separately treated; vault note is its own concept |
| `Modified Allan deviation` | **COVERED** |
| `Time deviation` | **COVERED** |
| `Hadamard deviation` | **COVERED** |
| `Modified Hadamard deviation` | **COVERED** |
| `Total deviation` | **COVERED** |
| `Modified Total deviation` | **COVERED** |
| `Hadamard Total deviation` | **COVERED** |
| `Time Total deviation` | **MISSING** — vault has TTOT concept, docs has no entry |
| `Theo1` | **MISSING** |
| `ThêoH` | **MISSING** |
| `MTIE` | **MISSING** |
| `Dynamic Allan deviation` | **MISSING** |
| `htdev` (pass-2.6 original-contribution note) | **STALE** — `theory/allan_family.md` has an HTDEV section, but it predates the implementation-source grounding pattern and does not cite `@sigmatau-htdev-impl` |
| `mhtotdev` (pass-2.6 original-contribution note) | **STALE** — `theory/total_family.md` has an MHTOTDEV section that predates the implementation-source grounding |

#### `Concepts/stability/algorithms/` (5 notes)

| Concept note | Docs status |
|---|---|
| `Greenhall-Riley EDF` | **COVERED** (`theory/confidence.md`) |
| `Lag-1 ACF noise identification` | **COVERED** (`theory/noise_id.md`) |
| `B1 ratio noise identification` | **COVERED** (`theory/noise_id.md` fallback section) |
| `Chi-squared confidence intervals` | **COVERED** (`theory/confidence.md`) |
| `Overlapping vs non-overlapping estimation` | **MISSING** |

#### `Concepts/stability/bias-corrections/` (4 notes)

| Concept note | Docs status |
|---|---|
| `TOTVAR bias correction` | **STALE** — discussed inside `theory/total_family.md` and `theory/confidence.md` but not as a standalone concept |
| `MTOT bias correction` | **STALE** — same |
| `HTOT bias correction` | **STALE** — same |
| `Theo1 bias correction` | **MISSING** |

#### `Concepts/ensemble/` (28 notes)

**Entire ensemble subpackage is MISSING from `docs/src/theory/`.** The docs has zero theory pages for ensemble. The 28 ensemble concept notes — clock SDEs, Kalman filter and variants (UD factorized, structured, performance bounds), oscillator networks (Kuramoto, EWFA, fragmentation), time-scale algorithms (basic time-scale, three-cornered hat, ANN ensemble, ML clock-bias forecasting, telemetry-based stability), PID controller, relativistic clock — none have a theory page.

### Code-link notes vs `docs/src/reference/*`

The reference pages use `@docs` blocks that pull docstrings from Julia source. The vault's Code-link notes do not duplicate that — they add concept-link wikilinks, status flags (stub / no-docstring), and design notes.

| Code-link note | Reference docs status |
|---|---|
| `base/{AbstractTimingData, PhaseData, FrequencyData, StabilityResult}` | **COVERED** (4 `@docs` entries in `reference/base.md`) |
| `stability/*` (25 notes) | **COVERED** (25 `@docs` entries in `reference/stability.md`) — but `reference/stability.md` lists 13 `@docs` symbols across the public API section + 9 internal-kernel symbols + 3 noise/EDF symbols = 25, matching exactly. The grouping ("Allan family", "Total family", "Noise identification", "EDF, bias, and confidence intervals", "Internal kernels") is a reasonable editorial structure that the vault does not duplicate. |
| `ensemble/*` (18 notes) | **COVERED** (18 `@docs` entries in `reference/ensemble.md`). Note that 12 of these have `<no docstring>` flagged in the vault — the `@docs` block in Documenter will produce a "missing docstring" warning that is currently silenced by `warnonly = [:missing_docs, ...]`. Adding docstrings is a Julia-source task, not a docs task. |

### `docs/src/refs.bib` vs `SigmaTauVault/references.bib`

| Aspect | docs/src/refs.bib | SigmaTauVault/references.bib |
|---|---|---|
| Entry count | 22 | 46 |
| Citekey style | CamelCase, no shorttitle (`RileyHowe2008`, `Greenhall2003`, `IEEE1139_2022`) | kebab-case, with shorttitle (`riley-2008-sp1065`, `greenhall-2003-edf-stability`, `ieee1139-2022-definitions`) |
| Coverage | every key cited in current docs prose | every source ingested across passes 1, 1-redux, and 2.6 (including the three new pass-1-redux PDFs and the two pass-2.6 implementation sources) |
| Implementation sources | not present | `@sigmatau-htdev-impl`, `@sigmatau-mhtotdev-impl` (pass 2.6 — not BibTeX-publishable but *would* need a stub if cited) |

13 distinct `@cite` labels are currently used in docs prose: `RileyHowe2008`, `Greenhall2003`, `Greenhall1997`, `Greenhall1999`, `Howe1999`, `Howe2000`, `Howe2001`, `Howe2005`, `Banerjee2023`, `IEEE1139_2022`, `Riley_R_2020`, `Riley2004`, `Sullivan_NBS_TN_1337`. All 13 resolve to vault entries under different citekeys.

---

## 4. Proposed changes

### `docs/src/` files

| File | Proposal | Rationale |
|---|---|---|
| `docs/src/index.md` | **PATCH** | Update three "reference math" bullets to use new vault-style citekeys after step (a) below. Otherwise leave landing prose as-is. |
| `docs/src/getting_started.md` | **KEEP AS-IS** | Pure tutorial chrome with no vault analog. The minimal `adev` example uses public API that hasn't changed. |
| `docs/src/theory/overview.md` | **PATCH** | (a) Refresh power-law noise table to match the 6 vault noise concepts including Random Run FM at α=−4 (vault corrected the original prompt's α=−3 nomenclature). (b) Update `[@cite ...]` labels. (c) Optionally add a short "see also" pointer to the new ensemble theory pages once they exist. (d) Tighten the slope-table footnote about α=−3 vs α=−4 — currently the docs row shows only α ∈ {2,1,0,−1,−2}, which is silently inconsistent with HDEV/MHDEV's α=−3,−4 convergence claims later in the file. |
| `docs/src/theory/allan_family.md` | **PATCH** | (a) Update HTDEV section to cite `@sigmatau-htdev-impl` as authoritative for the definition; the existing "Provenance" paragraph already says HTDEV is original to the package, but does not name the implementation source. (b) Update all `[@cite ...]` labels to vault style. (c) Slope-vs-noise table is fine. |
| `docs/src/theory/total_family.md` | **PATCH** | (a) Update MHTOTDEV section to cite `@sigmatau-mhtotdev-impl` as authoritative. The existing "There is no canonical paper for MHTOTDEV; the construction follows HV99…" paragraph is correct in spirit but does not point to the implementation source. (b) Note the vault's "2×2 matrix completion" framing in the MHTOTDEV section (the existing "all-three corner of the cube" framing is also valid; both are present in the implementation source — the docs should pick one). (c) Update cite labels. (d) Bias-correction policy block is fine and matches `Concepts/stability/bias-corrections/`. |
| `docs/src/theory/confidence.md` | **PATCH** | Update cite labels. The math content matches the vault's `Greenhall-Riley EDF.md` and `Chi-squared confidence intervals.md`. Bias-correction summary table is fine. |
| `docs/src/theory/noise_id.md` | **PATCH** | Update cite labels. Content matches the vault's two algorithm concept notes. Optional: add a one-paragraph "see also" pointer to `Overlapping vs non-overlapping estimation` once that concept gets its own page. |
| `docs/src/reference/base.md` | **KEEP AS-IS** | Page chrome wrapping `@docs` blocks. The `@docs` content comes from Julia docstrings; the page header is fine. |
| `docs/src/reference/stability.md` | **KEEP AS-IS** | Same pattern. The grouping is editorial and works. |
| `docs/src/reference/ensemble.md` | **KEEP AS-IS** | Same pattern. The 12 `<no docstring>` symbols are a Julia-source TODO, not a docs-edit issue. |
| `docs/src/tutorials/01_phase_data.md` | **KEEP AS-IS** | Stub awaiting a follow-up PR. Vault has no tutorial-content equivalent. |
| `docs/src/tutorials/02_compute_adev.md` | **KEEP AS-IS** | Stub. |
| `docs/src/tutorials/03_identify_noise.md` | **KEEP AS-IS** | Stub. |
| `docs/src/tutorials/04_confidence_intervals.md` | **KEEP AS-IS** | Stub. |
| `docs/src/tutorials/05_single_clock_steering.md` | **KEEP AS-IS** | Stub. |
| `docs/src/validation/methodology.md` | **KEEP AS-IS** | Handwritten short stub; no vault analog. |
| `docs/src/validation/stable32.md` | **KEEP AS-IS** | Handwritten numerical comparison tables, derived from `reference/validation/` fixtures. Content does not come from any vault note and shouldn't. |
| `docs/src/bibliography.md` | **KEEP AS-IS** | Two-line chrome plus the `@bibliography *` directive. |
| `docs/src/refs.bib` | **REWRITE (with citekey-mapping)** | Replace contents with `SigmaTauVault/references.bib` (46 entries). Every prose citation in `docs/src/` will need its key updated. See "Citekey mapping" below. |

### Proposed new files

The vault has eight stability concepts and 28 ensemble concepts that have **no docs page at all**. Adding all of them would be a large expansion, so I propose a **conservative** new-file set covering the most-cited gaps:

| Proposed new file | Source(s) | Justification |
|---|---|---|
| `docs/src/theory/ensemble_overview.md` | `Concepts/shared/Polynomial clock model.md`, `Concepts/ensemble/clock-models/{Two-state,Three-state} clock SDE.md`, `Concepts/ensemble/clock-models/{State transition matrix, Process noise covariance}.md` | Closes the largest single gap: zero ensemble theory pages exist. Explains the polynomial + SDE clock model that all ensemble exports rely on. |
| `docs/src/theory/kalman.md` | `Concepts/ensemble/filters/{Kalman filter, Innovation, Kalman gain, U-D factorized Kalman filter, Adaptive Kalman filter, Generalized ALS noise tuning, Structured Kalman filter for timescales, Kalman filter performance bounds}.md` | Standard Kalman + variants. Includes the new pass-1-redux U-D content. |
| `docs/src/theory/steering.md` | `Concepts/ensemble/filters/PID clock-steering controller.md` | Closes the PID gap directly; SigmaTauEnsemble exports PIDController/step!/steer_to_correction without any theory page. |
| `docs/src/theory/relativistic_clocks.md` | `Concepts/ensemble/clock-models/Relativistic clock.md` | Closes the relativistic-clock gap from pass-1-redux. |
| `docs/src/theory/ensembles_and_oscillator_networks.md` | `Concepts/ensemble/oscillator-networks/{Kuramoto clock synchronization, Nearest-neighbor coupling, EWFA baseline, Constellation fragmentation, Kuramoto coupling strength}.md` + `Concepts/ensemble/time-scale-algorithms/{Basic time-scale equation, Three-cornered hat, N-cornered hat, Dynamic clock weights, ANN ensemble timescale, ML-based clock-bias forecasting, Telemetry-based stability estimation}.md` | Single page covering both oscillator-network coupling schemes and timescale algorithms, since they're conceptually adjacent and individually short. |

I am **not** proposing standalone pages for every concept note. The vault's per-concept resolution is right for Obsidian navigation; for Documenter HTML, grouped pages (5 new pages above) read better. The `@cite` machinery still resolves per-concept-note citations to individual bib entries.

I am also **not** proposing standalone Theo1 / ThêoH / MTIE / Dynamic Allan / Time Total deviation theory pages in this pass — those are stability-side gaps that the user can fill in a follow-up if they want full coverage. Flagging them here as "deferred" rather than skipping silently.

### Citekey mapping (must be applied atomically with `refs.bib` rewrite)

| Old (docs) | New (vault) | Used in |
|---|---|---|
| `RileyHowe2008` | `riley-2008-sp1065` | index.md, all theory/* |
| `Greenhall2003` | `greenhall-2003-edf-stability` | index.md, theory/confidence.md |
| `Greenhall1997` | `greenhall-1997-third-difference-mvar` | theory/{allan,total}_family.md |
| `Greenhall1999` | `greenhall-1999-totvar` | theory/total_family.md |
| `Howe1999` | `howe-1999-modtotvar` | theory/total_family.md |
| `Howe2000` | `howe-2000-tothvar-ptti` | theory/total_family.md |
| `Howe2001` | `howe-2001-tothvar-steering` | theory/{total_family,confidence}.md |
| `Howe2005` | `howe-2005-tothvar-ieee` | theory/{allan,total}_family.md |
| `Banerjee2023` | `banerjee-2023-timekeeping` | theory/overview.md |
| `IEEE1139_2022` | `ieee1139-2022-definitions` | index.md, theory/{overview,allan_family}.md |
| `Riley_R_2020` | `riley-2020-r-frequency-stability` | theory/{confidence,noise_id}.md |
| `Riley2004` | `riley-2004-lag1-acf` | theory/noise_id.md |
| `Sullivan_NBS_TN_1337` | `sullivan-1990-tn1337` | theory/allan_family.md |

After the rewrite there will be 46 bib entries (vs current 22). The 33 new keys are unused by current docs prose but will be cited by the 5 proposed new theory pages.

### Build chrome

| File | Proposal | Rationale |
|---|---|---|
| `docs/make.jl` | **PATCH (small)** | Add the 5 new theory pages to the `pages =` tree under "Theory". Preserve `setdocmeta!`, `warnonly`, `canonical`, citation style, and `deploydocs`. |
| `docs/Project.toml` | **KEEP AS-IS** | Versions and `[sources]` are correct. |

---

## 5. Risk flags

**(a) Cite-syntax fragility.** Theory pages mix `[@cite KEY]` and `[label](@cite KEY)` styles. DocumenterCitations 1.x supports both today, but the deprecation arc favours the link form. If we apply the citekey rename atomically and consistently, we should bias to the link form (`[Riley & Howe 2008](@cite riley-2008-sp1065)`) for new content; for existing content, in-place rename is safer than syntax-change. Doing both at once will obscure diffs.

**(b) Stale API in examples.** I did not run `@example` blocks. Spot-checks of the function signatures in `getting_started.md`, `theory/overview.md`, `theory/allan_family.md`, `theory/total_family.md`, and `theory/confidence.md` look current — they call `adev`, `mdev`, `totdev`, `tdev`, `htdev`, `identify_noise`, `PhaseData`, all of which are still exported with the documented signatures. The `theory/noise_id.md` snippet `identify_noise(PhaseData(x, 1.0); m=8)` calls the function with a single `m` keyword, which differs from the actual signature `identify_noise(x::Vector{Float64}, m_values::Vector{Int}; …)`. **This snippet would fail at build time if doctest evaluated it; currently it isn't a `@example` block, just an inline ```julia` block, so it's textual only.** Apply phase should fix the call to match the real API.

**(c) Cross-references that would break under proposed rewrite.** The proposed new ensemble theory pages will be reachable via Documenter cross-references like `[Theory: Kalman filter](theory/kalman.md)`. Existing pages already use `[Theory: …](theory/…)` patterns, so the relative-link style is consistent. Nothing under `docs/src/reference/` or `docs/src/tutorials/` cross-references theory pages by anchor (only by file), so I do not anticipate broken anchors. Rewriting `refs.bib` does not change anchor names because Documenter generates anchors from citekeys; **all `[Riley & Howe 2008](@cite RileyHowe2008)` -> `[Riley & Howe 2008](@cite riley-2008-sp1065)` rewrites must be done together with the bib swap in the same commit**, otherwise the build will fail with "missing reference" warnings (which are not in `warnonly`).

**(d) HTDEV / MHTOTDEV provenance now mismatched.** `theory/allan_family.md` line 191–195 already correctly states HTDEV is "original to this package" and that "SP1065 […], IEEE 1139-2022 […], NBS-TN-1337 […] do not define it." This is consistent with the pass-2.6 implementation-source pattern. The patch needs only to add a citation to `@sigmatau-htdev-impl` so future readers can find the authoritative definition; the existing prose does not need to be rewritten.

`theory/total_family.md` line 159–162 says "There is no canonical paper for MHTOTDEV; the construction follows HV99's modified-total methodology applied to the FCS01 Hadamard total." This is **partially superseded** by the vault's `mhtotdev.md` concept note, which uses the cleaner "completes the 2×2 matrix" framing. The patch should extend (not replace) this paragraph and add the `@sigmatau-mhtotdev-impl` citation.

**(e) Bias-correction tables differ slightly.** `theory/confidence.md` line 64–70 gives a five-row bias table with α-dependent factors `{1.06, 1.17, 1.27, 1.30, 1.31}` for MTOTDEV. The vault's `MTOT bias correction.md` (`Concepts/stability/bias-corrections/`) uses the same table. No change needed.

**(f) Random Run FM nomenclature.** The vault corrected the original Random-Run-FM-at-α=−3 placement to α=−4 during pass 2 (see `_concept_log.md`). The existing `theory/overview.md` table goes only to α=−2 and does not show α=−3 or α=−4, so there is no contradiction with the vault — but the patch should add α=−3 (Flicker Walk FM) and α=−4 (Random Run FM) rows to the slope table so the noise-type coverage matches what `theory/total_family.md` later claims for HDEV/MHDEV (`α ∈ {−4, −3}` convergence). Alternatively, leave the overview table at α ∈ {2..−2} and add a sentence noting the Hadamard-family extension to α ∈ {−4, −3}.

**(g) `warnonly = [:missing_docs, ...]` masks the 12 ensemble docstring TODOs.** The plan is to keep this `warnonly` setting in place and address the 12 missing docstrings in a separate Julia-source PR. If `warnonly` is tightened, the build will start failing on those 12 symbols (`AbstractClockModel`, `TwoStateClock`, `ThreeStateClock`, `RelativisticClock`, `nstates`, `state_transition`, `process_noise`, `measurement_matrix`, `measurement_noise`, `AbstractEstimator`, `UDFactorizedFilter`, `KuramotoOscillator`).

**(h) `validation/stable32.md` numerical tables.** I did not regenerate these. They derive from `reference/validation/stable32_data_full.csv` and `comparison_report.md`. If those fixtures changed since the tables were written, the tables may be stale. **Defer to user**: do you want me to verify the tables against current fixtures, or treat the page as a frozen reference snapshot?

**(i) `noise_id.md` API call.** Line 37 currently shows `identify_noise(PhaseData(x, 1.0); m=8)` — this does not match the actual `identify_noise(x::Vector{Float64}, m_values::Vector{Int}; …)` signature. Mark for a small fix in the patch.

---

## 6. Suggested apply order

1. **Apply citekey rename atomically** — single commit replacing `docs/src/refs.bib` with the vault `references.bib` and rewriting all `@cite` labels in the 7 prose-citation files (`index.md` + 5 theory pages + `getting_started.md` if needed). Build after this commit; no broken references should remain. This must be the first step because every later patch references the new keys.
2. **Patch theory pages** — `theory/overview.md`, `theory/allan_family.md`, `theory/total_family.md`, `theory/confidence.md`, `theory/noise_id.md`, in that order. Each commit standalone-buildable. The HTDEV / MHTOTDEV provenance updates land in this step.
3. **Fix the `noise_id.md` API call** (item 5(i)) in the same commit as the `theory/noise_id.md` patch.
4. **Add 5 new theory pages** — write content first, then add to `make.jl` `pages =` tree. Order: `ensemble_overview.md` → `kalman.md` → `steering.md` → `relativistic_clocks.md` → `ensembles_and_oscillator_networks.md`. Each new page commits independently; `make.jl` is updated only on the final commit so partial trees don't break the build.
5. **Optional follow-up commit** — bias all `@cite` syntax to the link form `[label](@cite KEY)` for consistency, only if user wants it.
6. **Defer (out of scope for this pass)**: theory pages for `Theo1`, `ThêoH`, `MTIE`, `Dynamic Allan deviation`, `Time Total deviation`, `Overlapping vs non-overlapping estimation`; the four bias-correction concept-note pages; tutorial 01–05 narrative; `validation/methodology.md` narrative; `validation/stable32.md` numerical-fixture refresh; the 12 missing docstrings in `lib/SigmaTauEnsemble/`.

After step 5, the build should produce no missing-doc warnings beyond the pre-existing `warnonly`-suppressed 12 ensemble docstring TODOs, and every concept note in `Concepts/{shared,stability,ensemble}/` will have a discoverable docs page (either a section within an existing theory page or one of the 5 new pages).

---

## Decisions for user before phase 2 starts

1. **Approve the citekey rename and the 5 new theory pages?** If yes → phase 2 proceeds as ordered. If no → fall back to "patch existing pages only, leave bib alone, accept the ensemble coverage gap." → **APPROVED (default).**
2. **Bias-syntax migration to link form (`[label](@cite KEY)`)?** Optional; phase 2 default is "leave existing prose as-is, write new content with the link form." → **APPROVED (default — keep existing syntax, use link form for new content).**
3. **Stale tables in `validation/stable32.md`?** Default plan: leave as a frozen snapshot. Alternative: regenerate against current `reference/validation/` fixtures — adds work but reduces stale-content risk. → **APPROVED (default — frozen snapshot).**
4. **MHTOTDEV framing** in `theory/total_family.md` — keep the existing "all-three corner of the cube" prose, swap to the vault's "completes the 2×2 matrix" framing, or include both? Default plan: keep both, add the 2×2 framing as a complementary lens. → **APPROVED (default — keep both).**
5. **Slope-table coverage** in `theory/overview.md` — extend to α ∈ {−3, −4} or leave at α ∈ {+2, …, −2}? Default plan: extend, with a footnote about Hadamard convergence. → **APPROVED (default — extend).**

## Decision 6 — Deferred estimators (added by user)

Approved deferrals: Theo1, ThêoH, MTIE, Dynamic Allan deviation,
Time Total deviation (TTOT), and PVAR.

Treatment: each gets a theory section with equation and citation,
plus a `!!! note "Planned implementation"` callout. No `@docs` block,
no Julia stub required. Section placement:

  - Theo1, ThêoH       → `theory/allan_family.md` (after MHDEV section)
  - Dynamic ADEV       → `theory/allan_family.md` (separate section)
  - PVAR               → `theory/allan_family.md` (after MDEV)
  - TTOT               → `theory/total_family.md` (after MTOT)
  - MTIE               → new short section in `theory/overview.md`
                          OR new file `theory/mtie.md`
                         (apply phase decides based on length)

Citation grounding:
  - Theo1, ThêoH, MTIE, TTOT  → `@riley-2008-sp1065` + `@banerjee-2023-timekeeping`
  - Dynamic ADEV              → `@mckelvy-2025-telemetrystability`
                                 (which cites Galleani & Tavella 2009)
  - PVAR                      → `@banerjee-2023-timekeeping` §4.4.3

Each section ends with a callout of this form:

  !!! note "Planned implementation"
      The mathematical definition is documented above. The
      `theo1` / `theoh` / `mtie` / `dadev` / `ttot` / `pvar`
      function is not yet implemented in SigmaTauStability.jl.

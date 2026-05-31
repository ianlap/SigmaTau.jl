# SigmaTau.jl — Roadmap

Working list of outstanding engineering work. Items move from this file
to `CHANGELOG.md` as soon as they land — every shipped change should
remove the matching entry here and add one under `## [Unreleased]` in
the changelog in the same commit.

> **Audit date**: 2026-05-21 (post-v0.3.0 cut).

---

## 🟡 Correctness / completeness

- [ ] **MHTOTDEV bias / EDF Monte Carlo (paper-grade).** Synthesize
  known-α noise via `_gen_powerlaw_y`; compute MHTOTVAR and MHVAR on the
  same realizations; the ratio `B(α) = E[MHTOTVAR]/E[MHVAR]` measures the
  bias (replacing the current "unbiased by policy" stance), and
  `2·E[V̂]²/Var[V̂]` measures the EDF. Re-fit `_coeff_mhtot` and wire a
  measured `bias_correction(:mhtot, …)`. Only the single canonical
  Greenhall form exists now (the `:linear` recipe was removed), so no
  per-recipe split is needed. Reproducible harness in `tools/`, heavy
  sweep on the workstation. (Wave 4 Phase B.)
- [ ] **TOTDEV m=512 Stable32 quirk follow-up.** The Stable32
  cross-validation testset skips the m=512 row because Stable32 reports
  a value ~1.5 % larger than the raw SP1065 result (allantools agrees
  with our `:howe`). Stable32 identifies that row as FLFM (α=−1) and
  appears to apply an α-aware correction. Either confirm via
  `legdocs/vendor/` how Stable32 derives the reported σ at FLFM-tagged
  rows and apply the matching correction in the test, or document the
  divergence as an irreducible Stable32-vs-SP1065 policy gap.
- [ ] **MTOTDEV EDF coefficients for α=−1 and α=−2 verification.**
  Current `_coeff_mtot` values for those two α are single-point fits
  against Stable32's `s32_5_12_26` fixture (with `c` assumed from
  SP1065). To pin both `b` and `c` independently, capture Stable32
  EDF dumps at one additional AF per α — AF=200 or AF=400 for α=−1,
  and any record dominated by RWFM at a second AF for α=−2 (the
  current fixture only has α=−2 at AF=1000). Replace the single-point
  fits with two-point fits in `src/stab/stats/edf.jl::_coeff_mtot`.
- [ ] **Modified-total kernel parity rtol under multi-thread runs.**
  Per the threading note in the v0.2.0 changelog, the inner `@threads`
  reduction reorders summation, so kernel-vs-legacy parity tests at
  `rtol = 1e-12`/`1e-11` may need to drop to `~1e-9` for the
  modified-total kernels on `--threads auto` CI runners. Verify on a
  multi-thread CI run and loosen only the testsets that actually drift.

---

## 🟢 Polish

- [ ] **More `examples/`** — Literate pipeline currently ships
  `00_julia_for_metrologists`, `01_phase_data`, `02_compute_adev`,
  `06_three_cornered_hat`. Next batch candidates:
  - A noise-ID walkthrough using `noise_gen` to synthesize a known
    composite α-mixture, then `identify_noise` to recover the
    dominant power law at each τ.
  - A spectral / Welch-PSD example exercising the `Sy` / `Sx` / `L`
    estimators (synthesize a known `h_α` mixture, recover the PSD slope).
- [ ] **Compat upper bounds** in the root `Project.toml`. The merged
  manifest already pins `Distributions = "0.25.125"` and lists `compat`
  for AbstractFFTs, DocStringExtensions, RecipesBase, StaticArrays,
  julia. Tighten upper bounds for `Distributions` and `StaticArrays`
  once the dep matrix has been exercised on the General registry.
- [ ] **Remove `ldev` alias** — now marked `@deprecate`; delete after
  v0.3.0 is tagged.

---

## 🟢 Housekeeping

- [ ] **Agent-context briefs** (`CLAUDE.md`, `AGENTS.md`) are tracked
  in the repo. When the package matures and conventions stabilise,
  decide whether to keep them as-is or to land a sanitised version
  (without the authorship rules) under `docs/contributing/`.

---

## Docs follow-ups

- Fill tutorial narrative — start with `01_phase_data.md` and
  `02_compute_adev.md`.
- Tighten `warnonly = []` in `docs/make.jl` once all public API has
  docstrings.
- Refine `docs/src/refs.bib` with DOIs and page numbers from the PDFs
  in `legdocs/papers/`.
- Convert remaining kernel docstrings (hdev, mhdev, totdev, etc.) to
  use `$(SIGNATURES)` + `jldoctest` blocks following the adev/mdev
  pattern.
- **Three-cornered-hat theory page** beyond the tutorial. Material
  exists at `legdocs/papers/deviations/three_cornered_hat_*` (web
  articles + Riley papers).
- **Preprocessing / uncertainty theory page.** Gaps + outliers material
  lives at `legdocs/papers/preprocessing/`.

---

## ✅ Recently shipped

See [CHANGELOG.md](CHANGELOG.md) for the annotated `## [0.3.0]` and
`## [0.2.0]` blocks.

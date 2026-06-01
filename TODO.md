# SigmaTau.jl — Roadmap

Working list of outstanding engineering work. Items move from this file
to `CHANGELOG.md` as soon as they land — every shipped change should
remove the matching entry here and add one under `## [Unreleased]` in
the changelog in the same commit.

> **Audit date**: 2026-05-21 (post-v0.3.0 cut).

---

## 🟡 Correctness / completeness

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

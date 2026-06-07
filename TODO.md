# SigmaTau.jl — Roadmap

Working list of outstanding engineering work. Items move from this file
to `CHANGELOG.md` as soon as they land — every shipped change should
remove the matching entry here and add one under `## [Unreleased]` in
the changelog in the same commit.

> **Audit date**: 2026-06-07 (post-v0.3.0 cut).

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

# SigmaTau.jl — Roadmap

Working list of outstanding engineering work, sorted by priority. Items move
from this file to `CHANGELOG.md` as soon as they land — every shipped change
should remove the matching entry here and add one under `## [Unreleased]` in
the changelog in the same commit.

> **Audit date**: 2026-05-07 (after legacy-parity + EDF-fallback batch).
> Last verified end-to-end: `132/132` Stability + `15/15` Ensemble tests
> passing locally; `using SigmaTau` precompiles in ~6 s.

---

## 🔴 Critical

_None._ The repository builds, all subpackages test green, and the legacy
numerical reference is locked in.

---

## 🟡 High (correctness / completeness)

- [ ] **FFT-based FLPM/FLFM noise synthesis** — extend the multi-noise
  validation testset (`MTOTDEV across noise regimes`) to cover the two
  power-law types that need a `1/f` filter. Reference generator lives in
  `legacy/julia/src/noise_gen.jl`.
- [ ] **PID steering port** — legacy `kalman_filter` includes a PID
  controller around the predict/update loop; the new `predict!`/`update!`
  is steering-free. Needed before any clock-steering examples can land.

---

## 🟡 Medium (new estimators / models)

- [ ] **`RelativisticClock`** — empty struct in
  [`clocks.jl`](lib/SigmaTauEnsemble/src/models/clocks.jl). Lunar-PNT
  future work; design the relativistic correction terms and `state_transition` /
  `process_noise` overrides.
- [ ] **`UDFactorizedFilter`** — empty struct in
  [`filters.jl`](lib/SigmaTauEnsemble/src/estimators/filters.jl). U–D
  Bierman/Thornton factorisation for low-observability lunar-distance
  measurements (numerical stability when `S` is near-singular).
- [ ] **`KuramotoOscillator`** — empty struct in same file. Nearest-neighbor
  phase coupling estimator targeted at SWaP-constrained pLEO ensembles.

---

## 🟢 Low (polish)

- [ ] **`Documenter.jl` site** under `docs/` — auto-generate API reference
  from docstrings, port the equation pages from the legacy
  `legacy/docs/equations/` tree, publish to GitHub Pages.
- [ ] **More `examples/`** — currently only `quickstart.jl`. Add:
  - Kalman-only ensemble walkthrough (no stability-deviation overlap)
  - `FrequencyData` ↔ `PhaseData` round-trip demo
  - Multi-clock ensemble scenario (once a multi-clock model lands)
- [ ] **Compat bounds** in subpackage `Project.toml`s — currently only
  `Distributions = "0.25.125"` is pinned in `SigmaTauStability`. Add upper
  bounds for `StaticArrays`, `Reexport`, `RecipesBase` to support
  General-registry registration cleanly.
- [ ] **`@reexport` of `tdev`, `FrequencyData`** — verified working but
  not separately smoke-tested at the umbrella level. Add a short
  `examples/`-driven smoke test or a top-level `using SigmaTau` test.

---

## ✅ Recently completed

See [CHANGELOG.md](CHANGELOG.md) for the full list. Headlines from this
session:

- Stable32 cross-validation (85 row checks against the
  `reference/validation/stable32gen.DAT` fixture; tight rtol on
  ADEV/MDEV/HDEV/MHDEV/TDEV, documented looser rtol on TOTDEV/HTOTDEV/MTOTDEV
  due to bias-correction policy and boundary-reflection differences with
  Stable32)
- Strict legacy-kernel parity (52 assertions, rtol=1e-12)
- TOTDEV/HTOTDEV α=2,1 EDF fallback
- MTOTDEV multi-noise validation (WPM/WHFM/RWFM)
- `tdev` API wrapper, `FrequencyData` dispatches, `StabilityResult.edf`
- `RecipesBase` package extension (replaces the `PlotRecipes.jl` stub)
- `NEFF_RELIABLE` 50 → 30 + boundary test
- GitHub Actions CI matrix
- README, CHANGELOG, project_overview, this TODO

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

_None._

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
- [ ] **More `examples/`** — currently `quickstart.jl` and
  `clock_steering.jl`. Add:
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

## 🟢 Housekeeping

- [ ] **Agent-context briefs** (`CLAUDE.md`, `AGENTS.md`) are gitignored
  while in active use. When the package matures and conventions stabilise,
  decide whether to delete them or to land a sanitised version (without the
  authorship rules) under `docs/contributing/`.

---

## ✅ Recently completed

See [CHANGELOG.md](CHANGELOG.md) for the full list. Headlines from this
session:

- PID steering controller (`PIDController`, `step!`, `steer_to_correction`)
  ported into `SigmaTauEnsemble`; `predict!` gained an optional `steering`
  kwarg that adds a correction vector to the propagated state. Verified
  end-to-end on a drifting-clock fixture (`examples/clock_steering.jl`,
  ~600× phase-residual reduction).
- Power-law phase-noise synthesizer (`_gen_powerlaw_phase`, AbstractFFTs-based);
  MTOTDEV multi-noise validation extended from 3 to all 5 SP1065 noise types
  (WPM, FLPM, WHFM, FLFM, RWFM); ADEV/MDEV/HDEV/MHDEV kernel-parity
  cross-checked across the same alpha range
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

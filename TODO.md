# SigmaTau.jl — Roadmap

Working list of outstanding engineering work, sorted by priority. Items move
from this file to `CHANGELOG.md` as soon as they land — every shipped change
should remove the matching entry here and add one under `## [Unreleased]` in
the changelog in the same commit.

> **Audit date**: 2026-05-09 (post single-package restructure + theory-pages
> + relativistic-PNT cluster + total-family threading + EDF stride fix).
> Top-level test layout: `test/runtests.jl` aggregates `test/types/`,
> `test/stab/` and `test/est/`. `using SigmaTau` precompiles in ~6 s.

---

## 🔴 Critical

_None._ The repository builds, every testset passes, and the legacy
numerical reference is locked in.

---

## 🟡 High (correctness / completeness)

- [ ] **MHTOTDEV bias / EDF Monte Carlo per detrend recipe.** Synthesize
  known-noise via `_gen_powerlaw_phase` for each α; compute MHTOT and
  MHDEV; the ratio yields the bias factor `B(α)`. Re-fit `_coeff_mhtot`
  empirically per detrend recipe (`:greenhall`, `:linear`). Track
  per-recipe — EDF is recipe-specific. Spec:
  `docs/superpowers/specs/2026-05-07-detrend-kwarg-design.md`
  → "Out-of-scope / future work".
- [ ] **TOTDEV `:howe` allantools cross-val tightening.** Update
  `test/stab/allantools_cross_validation.jl` so the `"Total"` branch
  uses `_totdev_core(...; detrend=:howe)` at `rtol=1e-7` (allantools'
  raw `totdev` already matches our `:howe` to ~7 sig figs on the shared
  fixture; no need for the `m=512` skip there since allantools doesn't
  apply Stable32's alpha-aware correction).
- [ ] **TOTDEV m=512 Stable32 quirk follow-up.** The Stable32
  cross-validation testset skips the m=512 row because Stable32 reports
  a value ~1.5% larger than the raw SP1065 result (allantools agrees
  with our `:howe`). Stable32 identifies that row as FLFM (alpha=-1)
  and appears to apply an alpha-aware correction. Either confirm via
  `legdocs/vendor/` how Stable32 derives the reported sigma at
  FLFM-tagged rows and apply the matching correction in the test, or
  document the divergence as an irreducible Stable32-vs-SP1065 policy
  gap.
- [ ] **EDF stride factor `S` for the totdev/htotdev WPM/FLPM fallback
  path** in `src/stab/stats/edf.jl`. The four overlapped variants
  (`:adev`, `:mdev`, `:hdev`, `:mhdev`) were corrected from `S=1` to
  `S=m` in the recent EDF fix. The Greenhall–Riley fallback that
  totdev/htotdev hit at α∈{1,2} still passes `S = 1`. Apply the same
  overlapped convention to the fallback or document why it stays
  non-overlapped.
- [ ] **HTOTDEV EDF off-by-one investigation** (`R-MED-6` from the
  theory-pages spec). Reconcile the EDF count returned for HTOTDEV
  against the FCS01 / GR03 expectation; fix or document the offset.
- [ ] **LDEV / HTDEV CI scaling formal verification** (`R-MED-5`).
  Confirm the `τ / √(10/3)` scaling propagates correctly through the
  χ² interval mapping in `confidence_intervals` for the `:htdev`
  branch.
- [ ] **Modified-total kernel parity rtol** under multi-thread runs.
  Per the threading note in `CHANGELOG.md`, the inner `@threads`
  reduction reorders summation, so kernel-vs-legacy parity tests at
  `rtol = 1e-12`/`1e-11` may need to drop to `~1e-9` for the
  modified-total kernels on `--threads auto` CI runners. Verify on a
  multi-thread CI run and loosen only the testsets that actually drift.
- [ ] **`prop!` — measurement-less covariance propagation.** Add a
  public `prop!(est::AbstractEstimator, model, dt; steering=nothing)`
  that advances `est.x ← Φ x` and `est.P ← Φ P Φ' + Q` unconditionally,
  i.e. **without** the `est.k > 0` gate that `predict!` uses today
  (`src/est/estimators/filters.jl:202`). Use case: forward holdover-
  budget prediction with no measurements (see
  `examples/05_holdover_comparison.jl`, which currently drives Φ/Q
  manually because `predict!` no-ops on a fresh `est`). Should land
  alongside a TwoStateClock + ThreeStateClock testset that locks in
  `prop!`-vs-analytical Q-integration parity at `rtol=1e-14`. After
  this lands, refactor `examples/05_holdover_comparison.jl` to
  call `prop!` instead of the inline Φ·P·Φ' + Q loop.

---

## 🟡 Medium (new estimators / models)

- [ ] **`RelativisticClock`** — empty struct in
  [`src/est/models/clocks.jl`](src/est/models/clocks.jl). Lunar-PNT
  future work; design the relativistic correction terms and
  `state_transition` / `process_noise` overrides. Docstring
  (`!!! note "Stub implementation"`) is in place per the recent
  Est-docstring batch.
- [ ] **`UDFactorizedFilter`** — empty struct in
  [`src/est/estimators/filters.jl`](src/est/estimators/filters.jl).
  U–D Bierman/Thornton factorisation for low-observability
  lunar-distance measurements (numerical stability when `S` is
  near-singular). Docstring stub callout in place.
- [ ] **`KuramotoOscillator`** — empty struct in the same file.
  Nearest-neighbor phase coupling estimator targeted at SWaP-constrained
  pLEO ensembles. Docstring stub callout in place.

---

## 🟢 Low (polish)

- [ ] **More `examples/`** — Literate pipeline now ships
  `01_phase_data`, `02_compute_adev`, `03_kalman_single_clock`,
  `04_kalman_pid_steering`, `05_holdover_comparison`. Next batch
  candidates:
  - Multi-clock ensemble scenario (once a multi-clock model lands).
  - `RelativisticClock` walk-through (depends on the stub being
    fleshed out per the Medium-priority entry).
  - Three-cornered-hat noise separation on three independent
    `PhaseData` records.
- [ ] **Compat upper bounds** in the root `Project.toml`. The merged
  manifest already pins `Distributions = "0.25.125"` and lists
  `compat` for AbstractFFTs, DocStringExtensions, RecipesBase,
  Reexport, StaticArrays, julia. Tighten upper bounds for
  `Distributions` and `StaticArrays` once the dep matrix has been
  exercised on the General registry.
- [ ] **Top-level `using SigmaTau` smoke test** — verify the umbrella
  re-exports `tdev`, `htdev`, the `FrequencyData` dispatches, and the
  `Stab` / `Est` submodule names. Currently exercised only indirectly
  via the per-submodule testsets.
- [ ] **Remove `ldev` deprecated alias.** `htdev` is now the canonical
  Hadamard time-deviation export; `ldev` is kept as a forwarding alias
  "for one release." Schedule deletion post-0.2.

---

## 🟢 Housekeeping

- [ ] **Agent-context briefs** (`CLAUDE.md`, `AGENTS.md`) are gitignored
  while in active use. When the package matures and conventions
  stabilise, decide whether to delete them or to land a sanitised
  version (without the authorship rules) under `docs/contributing/`.
- [ ] **Refresh `project_overview.md`** for the single-package layout.
  Every `lib/SigmaTauBase/...`, `lib/SigmaTauStability/...`,
  `lib/SigmaTauEnsemble/...` link in the per-component status tables
  and the §6 file-inventory tree is stale; rewrite to point into
  `src/{types,stab,est}/` and `test/{types,stab,est}/`. The top-of-file
  "Restructure note" promises this has happened — close the gap.

---

## Docs follow-ups (post-v1)

- Fill tutorial narrative — start with `01_phase_data.md` and
  `02_compute_adev.md`.
- Tighten `warnonly = []` in `docs/make.jl` once all public API has
  docstrings (substantive warnings dropped from 23 → 1 in the recent
  Est-docstring batch).
- Add `tutorials/06_masterclock.md` once C6 lands from the parallel
  implementation track.
- Add `tutorials/07_ensemble.md` once D1 lands.
- Refine `docs/src/refs.bib` with DOIs and page numbers from the PDFs
  in `legdocs/papers/`.
- Convert remaining kernel docstrings (hdev, mhdev, totdev, etc.) to
  use `$(SIGNATURES)` + `jldoctest` blocks following the adev/mdev
  pattern.
- **Three-cornered-hat theory page.** Material exists at
  `legdocs/papers/deviations/three_cornered_hat_*` (web articles +
  Riley papers); deferred from the May-08 theory pages design as a
  separate future page.
- **Preprocessing / uncertainty theory page.** Gaps + outliers
  material lives at `legdocs/papers/preprocessing/`; same deferral.
- Update the validation page `docs/src/validation/methodology.md` to
  cross-link the perf benches under `benchmarks/bench/` once
  long-record runtime numbers are stable.

---

## ✅ Recently completed

See [CHANGELOG.md](CHANGELOG.md) for the full annotated list. Headlines
from the active `## [Unreleased]` block, grouped thematically:

**Restructure**

- Three-subpackage workspace collapsed into a single registerable
  package; `Stab` and `Est` submodules at the top level. `lib/`
  workspace tree gone; `Project.toml` is a single-package manifest;
  package UUID regenerated.
- `docs/src/reference/{base,stability,ensemble}.md` renamed to
  `{types,stab,est}.md`; `docs/make.jl` collapsed to a single
  `using SigmaTau` and one recursive `DocMeta.setdocmeta!` call.

**Theory documentation**

- Relativistic-PNT cluster expanded from a Seyffert-only placeholder
  into four pages: `relativistic_clocks` (hub),
  `relativistic_frames_and_timescales` (BCRS / GCRS / LCRS,
  TT/TCG/TCB/TDB/TCL/TL with conversions), `relativistic_corrections`
  (1PN proper-time, redshift, $L_{Gm}$, Earth–Moon Lagrange offsets,
  Shapiro / Sagnac), and `lunar_pnt_systems` (TWSTFT, RPS / emission
  coordinates / ABC). Eight new bib entries.
- Theory section filled out: `overview`, `allan_family`, `total_family`,
  `confidence`, and `noise_id` from SP1065 / IEEE 1139-2022 / GR03 /
  FCS01 / RG04 source material; `Sullivan_NBS_TN_1337` and
  `Riley_R_2020` bib entries added.
- Twelve previously-undocumented `SigmaTauEnsemble` exports gained
  docstrings (`AbstractClockModel`, the three Clock structs, six
  contract methods, the three estimator stubs); substantive docs
  warnings 23 → 1.
- Theory page misattribution of HTDEV (formerly `ldev`) to
  IEEE 1139-2022 Annex C corrected — the construction is original.
- Holdover tutorial added; theory/validation and theory/publications
  pages added.

**Reference-material consolidation**

- `legdocs/` Option-C reorganisation. `stable32docs/` retired:
  Stable32 product PDFs into `legdocs/vendor/`; preprocessing papers
  under `legdocs/papers/preprocessing/`; HTML web articles renamed
  alongside their topical PDFs; one byte-identical PDF dupe removed.

**API surface**

- `correct_bias::Bool=true` kwarg on `totdev`, `mtotdev`, `htotdev`,
  `mhtotdev`. `mhtotdev` documents the kwarg as a no-op (FCS01 /
  SP1065 publish no MHTOT bias model). `calc_ci=false` fast path now
  populates `noise_type` when `correct_bias=true` (since bias
  correction needs α).
- `identify_noise(x, m_values; …, detrend::Bool=true)` kwarg —
  `true` retains the per-m quadratic detrend matching allantools'
  `autocorr_noise_id`; `false` skips the polynomial fit so α matches
  Stable32 point-for-point.
- `detrend::Symbol` kwarg on `totdev` / `mtotdev` / `htotdev` /
  `mhtotdev` (each kernel exposes its canonical recipe plus
  `:linear` and `:legacy` aliases).
- `_totdev_howe` helper implementing SP1065 eqn 25 verbatim;
  `totdev` default detrend is now `:howe`.
- `mhtotdev` default detrend is now `:greenhall`.
- Public `ldev` renamed to `htdev`; legacy `ldev` kept as a
  deprecated alias for one release.

**Numerics / performance**

- Total-family kernels (`_totdev_legacy/_howe`, `_mtotdev_*`,
  `_htotdev_*`, `_mhtotdev_*`) now multithreaded — `_totdev_*`
  parallelises the outer m-loop directly; the six modified-total
  kernels parallelise the inner subsequence loop with chunk-private
  accumulators. Sliding-window inner reduction replaces the
  per-subsequence cumulative-sum buffer; reused per-chunk scratch
  buffer pool eliminates per-m, per-chunk allocation churn.
- Default confidence factor lowered from `0.95` to `0.683` (1-sigma)
  and exposed as `SigmaTauStability.DEFAULT_CONFIDENCE`.
- EDF stride factor `S` corrected from `1` to `m` for the four
  overlapped variants (`:adev`, `:mdev`, `:hdev`, `:mhdev`) in
  `src/stab/stats/edf.jl` (`:tdev`/`:htdev` inherit transitively).
- Lag-1 ACF `dmax` for the four Hadamard-family kernels (`hdev`,
  `mhdev`, `htotdev`, `mhtotdev`) bumped from `2` to `3` per RG04 §6,
  unblocking α ∈ {−3, −4} reporting.
- Pre-record linear detrend removed from the noise-ID `_preprocess`
  step (neither Stable32 nor allantools applies one).

**Validation infrastructure**

- `tools/regen_allantools_fixtures.py` now writes `%.17e`
  (round-trip-exact Float64). The
  `reference/validation/allantools_out/allantools_data_full.csv`
  fixture has been regenerated at full machine precision.
- `test/stab/allantools_cross_validation.jl` ADEV/MDEV/HDEV/TDEV
  cross-check tightened from `rtol=1e-4` to `rtol=1e-11` against
  allantools.

**CI**

- `lib/SigmaTauBase` test-matrix entry dropped (it's a types-only
  package).
- Umbrella `using SigmaTau` job uses `Pkg.instantiate()` instead of
  `Pkg.resolve()` so registries register on a fresh runner.
- `test/est/runtests.jl` guards the `legacy/julia/src` includes
  behind `isfile`; CI checkouts skip the legacy-parity testsets and
  log an `@info` instead of failing at load time.
- CI restored on the single-package layout; switched to
  `julia-actions/docdeploy`.

**Project workflow**

- TODO ↔ CHANGELOG rule formalised: every shipped change moves the
  TODO entry to `## [Unreleased]` in the same commit.

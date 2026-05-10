# SigmaTau.jl — Roadmap

Working list of outstanding engineering work, sorted by priority. Items move
from this file to `CHANGELOG.md` as soon as they land — every shipped change
should remove the matching entry here and add one under `## [Unreleased]` in
the changelog in the same commit.

> **Audit date**: 2026-05-09 (post single-package restructure + theory-pages
> + relativistic-PNT cluster + total-family threading + EDF stride fix;
> Julia-ecosystem research pass added flicker-Markov, MTIE, IMM,
> calibration, and JSMD-lunar-stack entries).
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
- [ ] **TOTDEV m=512 Stable32 quirk follow-up.** The Stable32
  cross-validation testset skips the m=512 row because Stable32 reports
  a value ~1.5% larger than the raw SP1065 result (allantools agrees
  with our `:howe`). Stable32 identifies that row as FLFM (alpha=-1)
  and appears to apply an alpha-aware correction. Either confirm via
  `legdocs/vendor/` how Stable32 derives the reported sigma at
  FLFM-tagged rows and apply the matching correction in the test, or
  document the divergence as an irreducible Stable32-vs-SP1065 policy
  gap.
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
---

## 🟡 Medium (new metrics, estimators, models)

- [ ] **PDEV EDF / χ² confidence model.** Vernotte 2015 / 2020 derive
  the parabolic-variance EDF in closed form for the five canonical
  power-law noises; port the table into `_coeff_pdev` and wire CI
  into the `pdev` API (currently returns empty `noise_type` / CI /
  EDF vectors).
- [ ] **Flicker-noise Markov approximation in the KF.** `TwoStateClock` /
  `ThreeStateClock` cover only integer-α processes (WPM, WFM, RWFM);
  flicker (FPM α=1, especially FFM α=−1) is the dominant regime for
  Cs / Rb / H-maser clocks over the τ range we care about and is
  currently unmodelled. Approximate 1/f as a sum of N first-order
  Gauss–Markov processes with log-spaced poles `βᵢ` over the band of
  interest, append N states to the clock vector, build Φ/Q
  block-diagonally — stays AD-clean and StaticArrays-friendly with N
  as a type parameter. New `FlickerClock{N}` (or a
  `flicker::FlickerMarkov{N}` block on `ThreeStateClock`) in
  [`src/est/models/clocks.jl`](src/est/models/clocks.jl). References:
  Davis–Greenhall–Stacey 2005 "A Kalman filter clock algorithm for use
  in the presence of flicker frequency modulation noise";
  Galleani–Tavella 2010 "Time and the Kalman filter"; Zucca–Tavella
  2005. Validation: synthesise FFM via `_gen_powerlaw_phase`, fit the
  augmented filter, check residual whitening and recovered σ_y(τ)
  flicker floor.
- [ ] **`fit_clock_params(::PhaseData)` calibration helper.** Inverse
  problem: given a measured phase record, recover diffusion
  coefficients (`q0..q3` for `ThreeStateClock`, plus `(βᵢ, σᵢ)` for
  the flicker block when present). Backend:
  `EnsembleKalmanProcesses.jl` (Apache-2.0) — EKI/UKI is fast on this
  low-dimensional, derivative-free inverse fit and composes cleanly
  with whatever clock model is in play. Add as a weakdep + extension
  to keep `[deps]` slim.
- [ ] **IMM (Interacting Multiple Models) estimator** for clock fault
  detection — switching dynamics between healthy / glitch / aging
  modes. Reference implementation:
  `LowLevelParticleFilters.InteractingMultipleModels`. Pairs naturally
  with the existing noise-ID anomaly detector and a future
  "anomaly demo" example.
- [ ] **`RelativisticClock`** — empty struct in
  [`src/est/models/clocks.jl`](src/est/models/clocks.jl). Lunar-PNT
  future work; design the relativistic correction terms and
  `state_transition` / `process_noise` overrides. Docstring
  (`!!! note "Stub implementation"`) is in place per the recent
  Est-docstring batch. When this lands, adopt the
  JuliaSpaceMissionDesign stack as deps: `Tempo.jl`
  (TAI/UTC/TT/TDB/TCG/TCB epoch transformations, allocation-free,
  leap seconds), `Ephemerides.jl` (pure-Julia DE440/PA440 SPK/PCK
  reader, ForwardDiff-clean, faster than SPICE.jl), and
  `FrameTransformations.jl` (lunar frames via PA440 + user-defined
  frames). All MIT, all pure Julia, all AD-clean.
- [ ] **`UDFactorizedFilter`** — empty struct in
  [`src/est/estimators/filters.jl`](src/est/estimators/filters.jl).
  U–D Bierman/Thornton factorisation for low-observability
  lunar-distance measurements (numerical stability when `S` is
  near-singular). Reference: `LowLevelParticleFilters.SqKalmanFilter`
  (square-root form, same numerical-stability goal — read its source
  before implementing). Docstring stub callout in place.
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
  - **GP-based holdover prediction** via `TemporalGPs.jl` with a
    sum-of-Matérn-1/2 kernel — mathematically equivalent to the
    flicker Markov approximation (Matérn-1/2 ≡ OU); sits next to the
    existing TDEV / HTDEV / KF curves in
    `examples/05_holdover_comparison.jl` as a fourth reference line.
  - **UDE drift learning** via `DiffEqFlux.jl` — show how a
    `TwoStateClock`'s `predict!` composes with a NeuralODE that
    learns the unmodelled drift residual. Don't take SciML as a dep;
    keep it confined to the example.
- [ ] **`taus` enum API** (`AllTaus`, `Octave`, `HalfOctave`,
  `QuarterOctave`, `Decade`, `HalfDecade`) à la `AllanDeviations.jl`
  — cleaner than passing arrays. Backwards-compatible if added as an
  alternative to the current array form.
- [ ] **Compat upper bounds** in the root `Project.toml`. The merged
  manifest already pins `Distributions = "0.25.125"` and lists
  `compat` for AbstractFFTs, DocStringExtensions, RecipesBase,
  Reexport, StaticArrays, julia. Tighten upper bounds for
  `Distributions` and `StaticArrays` once the dep matrix has been
  exercised on the General registry.
- [ ] **Remove `ldev` deprecated alias.** `htdev` is now the canonical
  Hadamard time-deviation export; `ldev` is kept as a forwarding alias
  "for one release." Schedule deletion post-0.2.

---

## 🟢 Housekeeping

- [ ] **Agent-context briefs** (`CLAUDE.md`, `AGENTS.md`) are gitignored
  while in active use. When the package matures and conventions
  stabilise, decide whether to delete them or to land a sanitised
  version (without the authorship rules) under `docs/contributing/`.

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

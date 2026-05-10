# Changelog

All notable changes to **SigmaTau.jl** are tracked here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- **IO expansion: file readers, detrend, Howe gap imputation.** New
  top-level `src/io/` directory consolidating all data-side IO.
  `read_phase(path; …)` and `read_frequency(path; …)` parse 2-column
  (or N-column with `time_col=0`) text files via stdlib
  `DelimitedFiles.readdlm` and return `PhaseData` / `FrequencyData` with
  optional preprocessing (`scaling`, `detrend`, `fillgaps`).
  `detrend(::PhaseData; method=:linear)` and the matching
  `FrequencyData` method support `:linear`, `:endpoint`, `:mean`, and
  `:none` modes; multiple dispatch on the timing-data types keeps the
  bare name `detrend` collision-free with other packages' vector-only
  definitions. `fillgaps(::PhaseData)` / `fillgaps(::FrequencyData)`
  port the Howe & Schlossberger reflect-and-FFT-filter algorithm
  (PTTI 2009) so imputed samples preserve local noise character and
  AVAR/MDEV shape. `FFTW.jl` promoted from test-extra to core dep (also
  needed by upcoming spectral plots). 34 new tests under `test/io/`.
- `save_result(path, r)` and `load_result(path)` — round-trip a
  `StabilityResult` to/from a tab-separated text file using stdlib I/O only
  (no new dependencies). Handles both `calc_ci=true` and `calc_ci=false`; CI
  columns are written as `NaN` when absent and reconstructed as empty vectors
  on load. Exported at the top-level umbrella alongside `read_phase` /
  `read_frequency` (relocated from `SigmaTau.Stab`); the old call path
  still works through the umbrella re-export.
- `examples/00_julia_for_metrologists.jl` — Literate tutorial targeted at
  Stable32 users migrating to Julia. Covers `juliaup` installation, REPL /
  script / Pluto.jl usage modes, `.DAT` file loading idiom, first `adev` call,
  `StabilityResult` field guide, overlaying multiple deviations on one plot,
  and `save_result`. Wired into the docs Tutorials nav as tutorial 00.

- **B1(N, μ) and R(n)(af, b) closed-form regression coverage** in
  `test/stab/runtests.jl`. Cross-checks the constants in
  `_b1_theory` and `_rn_theory` against the canonical allantools
  reference (`ci.py`, Wallin 2018 citing Howe 2000) at machine
  precision for every (N, μ) and (af, b) consumed by the B1/R(n)
  noise-ID fallback. Audit confirmed the existing implementation
  is bit-identical to allantools — these tests now lock the
  contract so silent drift in either constant fails CI.
- **Pedagogical Kalman-filter reference** added to `docs/src/refs.bib`
  and cited in `docs/src/theory/kalman.md` References block:
  Chaudhari 2022, *Kalman Filter and its Variants* (UVA Link Lab,
  Learning in Robotics course notes, Chapter 3) — useful as a
  from-first-principles refresher of the KF/EKF/UKF/PF recursions
  this module implements.

### Fixed

- **Noise-ID was silently returning `:unknown` for any phase record
  with std below ~√eps ≈ 1.5e-8.** `_lag1_acf` (`src/stab/noise/lag1.jl`)
  guarded against 0/0 with `ssx < eps(Float64) * length(x)`, but that
  threshold mixes dimensions: for length-N data with std σ,
  `ssx ≈ σ²·N`, and the rule reduces to `σ² < eps`, i.e. tripping on
  any data with std below ~1.5e-8 regardless of N. Phase records in
  seconds (typically 1e-9..1e-12) hit it; they round-tripped through
  `_lag1_acf` → NaN → `:unknown`, and the rendered docs example
  (`tutorials/02_compute_adev.md`) showed 8 of 10 m values as
  `:unknown`. Replaced with a scale-invariant relative guard
  (`ssx ≤ eps · ‖x‖²`) that only catches genuinely-constant input.
  Cross-validated against `allantools.autocorr_noise_id 2024.06` on
  `reference/validation/stable32gen.DAT`: 8/9 m values match exactly,
  with one borderline disagreement at m=4 (documented). Two new
  regression testsets in `test/stab/runtests.jl` lock no-`:unknown`
  on the compute_adev fixture, scale invariance under ×1e6
  rescaling, and the allantools cross-check (8 m values).

### Changed

- `ldev` forwarding alias replaced with `Base.@deprecate ldev htdev`; callers
  now receive a deprecation warning and should migrate to `htdev`.
- Stub types `RelativisticClock`, `UDFactorizedFilter`, and `KuramotoOscillator`
  now throw `ArgumentError` with a descriptive message on any method call,
  replacing the previous silent `MethodError`.
- `docs/src/reference/stab.md` kernel section renamed to "Advanced / research
  kernels" with a usage note; `_mtie_core` and `_pdev_core` added to the
  `@docs` block (were exported but not documented).

### Added

- **Three-cornered-hat tutorial** at
  `examples/06_three_cornered_hat.jl`. Synthesises three independent
  free-running clocks via `_gen_powerlaw_phase`, builds the three
  pairwise difference records, runs `adev` on each, and solves the
  classical TCH linear system to recover each clock's individual
  σ_y(τ). End-to-end runs cleanly; recovered σ tracks ground truth
  to within ~5 % at small τ, illustrates "TCH break-points"
  (negative-variance recoveries clamped to zero) at long τ.
  Includes prose callout for the two real-world failure modes —
  correlated noises and one clock dominating. Wired into the docs
  Tutorials nav as `tutorials/06_three_cornered_hat.md`.
  `examples/Project.toml` gains `FFTW` for the noise-synthesis FFT
  backend.

- **KF 1σ holdover band via `prop!`** in
  `examples/05_holdover_comparison.jl`. New §6.5 captures the
  filter's converged `P_mature`, re-seeds a side-channel estimator
  per horizon, calls `prop!(side, clock, h·τ₀)` and reads
  `√P[1,1]` as the theoretical 1σ phase-error envelope. Plots
  alongside the existing TDEV / HTDEV / KF-RMS curves on the same
  log-log axis; the prose comparison is updated since `prop!` is
  no longer a "future feature".

- **`prop!(est, model, dt; steering=nothing)`** — unconditional
  covariance propagation alongside the existing `predict!` /
  `update!` loop. Advances `est.x ← Φ(dt) x` and
  `est.P ← Φ(dt) P Φ(dt)' + Q(dt)` regardless of `est.k`, never
  bumping the step counter. Powers the shaded ±1σ holdover-band
  pattern (`examples/05_holdover_comparison.jl` follow-up) and
  side-channel "what-if" projections without disturbing live filter
  sequencing.
- **`state_transition(model, dt)`** / **`process_noise(model, dt)`**
  overloads in `src/est/models/clocks.jl`. Existing single-arg
  methods (`state_transition(model)`, `process_noise(model)`)
  delegate to the new dt-aware forms with `model.tau`, so legacy
  callers see no behaviour change. Tested at `rtol=1e-14`
  Q-integration parity against hand-derived closed-form expressions,
  Φ group property + Q additivity composition (two prop!s ≡ one
  prop! at any (dt₁, dt₂)), and full prop!-vs-predict! parity once
  past the `est.k > 0` gate.
- **Top-level umbrella smoke test** at `test/umbrella_smoke.jl`.
  Verifies a bare `using SigmaTau` exposes every public symbol
  (Stab + Est + types), confirms `SigmaTau.Stab` / `SigmaTau.Est`
  are reachable as modules, exercises the `FrequencyData` dispatch
  on every deviation including the new MTIE / PDEV, and locks the
  `ldev` ≡ `htdev` deprecated-alias contract. Closes the
  "exercised only indirectly" gap from the previous TODO.
- **MTIE** (Maximum Time Interval Error) per ITU-T G.810. Public
  `mtie(::PhaseData, m_values)` + `mtie(::FrequencyData, …)` plus
  `_mtie_core` kernel in `src/stab/core/mtie.jl`. Sliding-window
  peak-to-peak phase excursion; reported as a deterministic envelope
  (no published EDF / χ² model, so CI fields are empty). Closes the
  only stability metric SigmaTau lacked vs `AllanDeviations.jl`.
- **PDEV** (parabolic deviation) per Vernotte–Lenczner–Bourgeois–
  Rubiola, IEEE T-UFFC 63(4) 2016 and Vernotte 2020. Public
  `pdev(::PhaseData, m_values)` + `pdev(::FrequencyData, …)` plus
  `_pdev_core` kernel in `src/stab/core/pdev.jl`. Built from a
  least-squares parabolic fit; recommended for ω-averaged frequency
  uncertainty. PDEV(τ₀) ≡ overlapping ADEV(τ₀) by construction at
  `m=1`. CI fields empty pending future EDF port (TODO).
- New testsets cross-checking MTIE on a hand fixture, monotonic ramp,
  constant phase, and naive-double-loop parity at `rtol=1e-15`; PDEV
  against the allantools reference formula at `rtol=1e-12`, plus the
  m=1 ≡ ADEV identity, linear-trend annihilation, and constant-phase
  invariants.

### Changed

- **HTDEV CI scaling formally verified** (closes R-MED-5). New
  assertions in `test/stab/runtests.jl` lock the identity that
  `htdev = mhdev × τ/√(10/3)` propagates correctly through the χ²
  CI mapping: dev, ci_lower, and ci_upper all scale by the same
  multiplicative factor, so the relative-CI ratios
  `ci_lower / dev` and `ci_upper / dev` are identical to MHDEV's at
  every τ. Same pattern as the existing TDEV / MDEV verification.
- **`predict!` now actually uses its `dt` argument.** Previously the
  signature accepted `dt` and silently ignored it in favour of
  `model.tau` (latent issue R-MED-8). With the dt-aware
  `state_transition(model, dt)` / `process_noise(model, dt)`
  overloads in place, `predict!` now re-derives Φ and Q for the
  caller-supplied `dt`. Backwards-compatible at every existing call
  site — they all pass `dt = model.tau` — but `predict!(est, model, h)`
  with `h ≠ model.tau` is now usable as a finer-than-tau or
  coarser-than-tau propagator instead of a no-op. Signature also
  loosened from `dt::Float64` to `dt::Real` for symmetry with `prop!`.
- **MTIE kernel reimplemented as a monotonic-deque sliding window**
  in `src/stab/core/mtie.jl`. Each m runs in O(N) total work via two
  parallel pre-allocated index deques (one for the running window
  max, one for the min). Bench at `N = 50 000` with 20 log-spaced m
  values (max m = N/2) shows a 138× speedup over the previous
  O(N·m) kernel; output is bit-identical (the existing
  naive-double-loop parity test at `rtol = 1e-15` continues to
  pass). No external dependency added.
- **TOTDEV `:howe` allantools cross-validation tightened** from
  `rtol = 0.15` (legacy-kernel boundary-policy floor) to `rtol = 1e-7`
  (`test/stab/allantools_cross_validation.jl`). Switched the `"Total"`
  branch from `LK.totdev_var` (legacy detrend) to
  `_totdev_core(...; detrend=:howe)` — apples-to-apples against
  allantools' SP1065-verbatim raw totdev. All 13 Total rows pass
  including m=512 (allantools doesn't apply Stable32's α-aware
  correction, so no skip needed there).
- **EDF stride factor `S` corrected from `1` to `m`** in the
  WPM/FLPM fallback path for `:totdev` and `:htotdev` in
  `src/stab/stats/edf.jl`. Matches the overlapped convention applied
  to the four overlapped variants (ADEV/MDEV/HDEV/MHDEV) — both
  TOTDEV and HTOTDEV operate on a stride-1 phase record (Howe's
  reflected-boundary extension preserves overlap), so the overlapped
  EDF is the consistent choice. Existing CI bounds shift slightly at
  α ∈ {1, 2} on long records.
- Repository housekeeping. `lib.bak/` (pre-restructure recovery
  snapshot) and `rough_changelog/` (six superseded
  implementation_plan/walkthrough drafts from May 7) deleted from the
  working tree; both were already gitignored. The three internal
  `docs/_pass4_apply_log.md`, `docs/_reconciliation_plan.md`, and
  `docs/_restructure_plan.md` planning logs moved out of Documenter's
  source root into `docs/superpowers/plans/`. The orphan
  `docs/MatsakisPIDControllersPTTI20-0022 (1).pdf` moved into
  `legdocs/papers/deviations/2020_matsakis_pid_controllers_ptti.pdf`
  (gitignored). Top-level `validation/` renamed to `benchmarks/` via
  `git mv` to disambiguate from `reference/validation/` (numerical
  fixtures) and `docs/src/validation/` (doc pages); every internal
  path reference inside the bench scripts and result reports updated
  in lockstep. `scratch.jl` and `scratch.py` added to `.gitignore`
  under a new "Scratchpad files" section. `.gitignore` updated to
  drop the `lib.bak/` and `rough_changelog/` stanzas and to rename
  `validation/plots/` → `benchmarks/plots/`.

### Changed (BREAKING)

- Restructured from a three-subpackage workspace
  (`SigmaTauBase` / `SigmaTauStability` / `SigmaTauEnsemble`) into a
  single registerable package with two submodules (`SigmaTau.Stab`
  and `SigmaTau.Est`). Shared types (`PhaseData`, `FrequencyData`,
  `StabilityResult`) now live at the top level. All previous exports
  continue to be re-exported from `SigmaTau`, so user code that imported
  via `using SigmaTau` keeps working unchanged. Code that explicitly
  imported `using SigmaTauBase`, `using SigmaTauStability`, or
  `using SigmaTauEnsemble` must switch to `using SigmaTau` (or
  `using SigmaTau.Stab` / `using SigmaTau.Est` for the submodules).
- Repository no longer contains the `lib/` workspace tree. `Project.toml`
  is a single-package manifest with no `[workspace]` or `[sources]`
  blocks. Package UUID regenerated; `julia` compat bumped to `1.11`.
- `Plots` removed from the merged `[deps]` (it was only declared on the
  legacy `SigmaTauStability/Project.toml`, never imported in source);
  plot recipes still load via the existing `RecipesBase` weakdep
  extension.
- `docs/src/reference/{base,stability,ensemble}.md` renamed to
  `{types,stab,est}.md`; `docs/make.jl` collapsed to a single
  `using SigmaTau` and one recursive `DocMeta.setdocmeta!` call.

### Added

- Theory: relativistic-PNT cluster expanded from a Seyffert-only
  placeholder into four pages — `relativistic_clocks` (hub), new
  `relativistic_frames_and_timescales` (BCRS / GCRS / LCRS, TT / TCG /
  TCB / TDB / TCL / TL with defining constants and conversions), new
  `relativistic_corrections` (1PN proper-time mapping, gravitational
  redshift on the Moon, $L_{Gm}$ Earth–Moon rate constant and
  Earth–Moon Lagrange-point offsets, cislunar drift regimes for
  vLLO / LLO / ELFO / L1 / NRHO, Shapiro and Sagnac), and new
  `lunar_pnt_systems` (TWSTFT synchronous / asynchronous in the ESA
  Moonlight architecture, relativistic positioning systems with
  emission coordinates / ABC). Multi-source grounding from Ashby &
  Patla 2024, Turyshev 2025, Leonard et al. 2026, Iess et al. 2025
  and 2026, Cacciapuoti & Salomon 2009, Reinhardt et al. 2024, and
  Gomboc et al. 2013, in addition to the existing Seyffert 2025.
  Eight bib entries added to `docs/src/refs.bib` for the new sources;
  pages nested under a "Relativistic PNT" sub-tree in `docs/make.jl`.
- Docstrings for twelve previously-undocumented `SigmaTauEnsemble`
  exports (`AbstractClockModel`, `TwoStateClock`, `ThreeStateClock`,
  `RelativisticClock`, `nstates`, `state_transition`, `process_noise`,
  `measurement_matrix`, `measurement_noise`, `AbstractEstimator`,
  `UDFactorizedFilter`, `KuramotoOscillator`); the three stub types
  carry explicit `!!! note "Stub implementation"` callouts. Substantive
  docs warnings fall from 23 to 1.
- Theory section: filled out `overview`, `allan_family`, `total_family`,
  `confidence`, and `noise_id` from SP1065 / IEEE 1139-2022 / GR03 /
  FCS01 / RG04 source material. Hybrid narrative/reference structure
  with small `@example` blocks for slope-vs-noise demonstrations and a
  bias-policy admonition on the total-family page that cross-links the
  validation page. Bib entries `Sullivan_NBS_TN_1337` (NBS-TN-1337,
  Sullivan/Allan/Howe/Walls 1990) and `Riley_R_2020` (Riley, *Frequency
  Stability Analysis Using R*, Hamilton Tech, Rev. C 2020) added to
  `docs/src/refs.bib`.
- Reference-material consolidation: `stable32docs/` merged into
  `legdocs/`. Stable32 product PDFs moved to `legdocs/vendor/`;
  preprocessing papers (gaps, outliers) under
  `legdocs/papers/preprocessing/`; HTML web articles renamed and placed
  alongside their topical PDFs; one byte-identical PDF dupe removed.
  Both folders remain gitignored; `.gitignore` updated to drop the
  now-defunct `stable32docs/` line.
- `correct_bias::Bool=true` keyword on `totdev`, `mtotdev`, `htotdev`,
  and `mhtotdev`. Default `true` preserves prior behavior (apply the
  SP1065 / FCS 2001 noise-type-dependent bias factor `B(α)` to the raw
  kernel output). Pass `correct_bias=false` to return the raw kernel
  value, which matches Stable32 and allantools (neither tool applies
  the bias correction by default). On `mhtotdev` the kwarg is a
  documented no-op — FCS 2001 and SP1065 publish no `B(α)` for MHTOT
  and the estimator is treated as unbiased regardless. Side benefit:
  the `calc_ci=false` fast path now also populates `noise_type` when
  `correct_bias=true` (previously left empty), since bias correction
  needs α; the field stays empty when both flags are false.
- `identify_noise(x, m_values; …, detrend::Bool=true)` keyword. The new
  default (`true`) keeps the per-m quadratic detrend that allantools
  applies in `autocorr_noise_id`; passing `detrend=false` skips the
  polynomial fit so α matches a Stable32-generated reference fixture
  point-for-point. The 5σ outlier filter on the full record is
  unchanged and unconditional.

### Changed

- Public function `ldev` renamed to `htdev` (Hadamard time deviation).
  The estimator and its formula are unchanged: `htdev` wraps `mhdev`
  and applies the `τ / √(10/3)` scaling, exactly as the previous `ldev`
  did. `StabilityResult.deviation_type` is now `:htdev`. The legacy
  `ldev` is kept as a deprecated alias that forwards to `htdev` for one
  release; downstream code should migrate. `htdev` is the documented
  canonical name throughout the API reference and theory pages.
- Total-family kernels (`_totdev_legacy`, `_totdev_howe`,
  `_mtotdev_greenhall`, `_mtotdev_linear`, `_htotdev_greenhall`,
  `_htotdev_linear`, `_mhtotdev_greenhall`, `_mhtotdev_linear`) now
  run multithreaded. The two `_totdev_*` kernels parallelize the
  outer m-loop directly (no inner reduction to thread). The six
  modified-total kernels parallelize the **inner subsequence loop**:
  the m-loop runs sequentially, and within each m the per-window
  `for n in 1:nsubs` reduction is partitioned into
  `min(nthreads, nsubs)` chunks; each chunk runs on a thread with
  a private `local_sum` accumulator. Chunks reduce at the end via
  `outer_sum = sum(chunk_sums)`.

  Per-m threading was tried on the modified kernels and lost.
  m-work is severely non-uniform (heavy m's dwarf light ones), so
  per-m `@threads :dynamic` left cores idle once the light m's
  finished while one core ground through the tail. Stacking inner
  threading on top (nested `@threads`) made it worse — task
  scheduling overhead with no spare cores to redistribute to.
  Subsequence-list parallelism is uniform (every chunk does the
  same FFT-size work for a given m), so it gets near-linear scaling
  per m without nesting.

  **Sliding-window inner reduction.** The cumulative-sum buffers (`cs`
  for `_mtotdev_*` / `_htotdev_*`, `S` for `_mhtotdev_*`) are gone.
  Their windowed differences `cs[j+km+1] − cs[j+(k−1)m+1]` are now
  maintained as running sums `a1`, `a2`, `a3` updated by
  `a += ext[j+km] − ext[j+(k−1)m]` per slide step. Saves a full
  pass over the data per subsequence (≈9m memory ops on the modified-
  total kernels) and one chunk-private allocation. Carry chain on
  the running sums prevents `@simd` vectorisation of the slide loop,
  but the saved memory traffic outweighs the lost SIMD on the heavy m's.

  **Buffer pool.** The remaining per-chunk scratch (`ext`, plus `d3_vec`
  on the `_mhtotdev_*` kernels) is allocated once at the top of the
  kernel call, sized for the largest `m` in `m_values`, and reused
  across every m. Cheap m's only touch the leading `3·seg_len` /
  `L3` slice. Eliminates the per-m, per-chunk allocation churn
  (was ~150 MB on a 2.6M-sample file, ~600 MB on a 3M-sample file)
  and the GC pauses it triggered. Net effect on the benchmarked
  `mtotdev(N=4·10⁵)` workload: ~415 s without pool → ~290–340 s
  with pool on 8 threads.

  **FP determinism caveat:** the inner reduction reorders summation, so
  results may drift by ~few ULPs across runs with different thread
  counts. Existing parity tests at `rtol = 1e-12`/`1e-11` may need to
  loosen to `~1e-9` for the modified-total kernels.

  To see the speedup, start Julia with `--threads auto` (or `-t N`);
  with `-t 1` the inner `@threads` is a no-op and behavior matches
  the sequential version exactly.
- Default confidence factor for every public deviation API (`adev`,
  `mdev`, `hdev`, `tdev`, `mhdev`, `htdev`, `totdev`, `mtotdev`,
  `htotdev`, `mhtotdev`) lowered from `0.95` to `0.683` (1-sigma).
  Now exposed as the exported constant
  `SigmaTauStability.DEFAULT_CONFIDENCE`. The new value matches the
  Stable32, AllanLab, allantools, and Greenhall–Riley convention; the
  prior `0.95` made cross-tool CI overlays disagree by a factor of
  ~1.96 even when the underlying EDFs matched. Pass
  `confidence=0.95` explicitly if you want the prior default.

### Removed

- Pre-record linear detrend in the noise-ID `_preprocess` step. Neither
  Stable32 nor allantools' `autocorr_noise_id` apply a full-record
  polynomial fit before the per-m loop, so the prior SigmaTau-only
  step shifted α away from both external references on records with
  any polynomial structure. The 5σ outlier filter that lived in the
  same helper is unchanged. The unused `_detrend_linear` helper was
  deleted in the same commit.

### Fixed

- Theory page: removed the misattribution of the Hadamard time
  deviation (HTDEV, the function previously named `ldev`) to
  IEEE 1139-2022 Annex C. IEEE 1139-2022 does not define HTDEV; the
  construction is original to this package. The Annex-C citation has
  been replaced with a provenance note in `theory/allan_family.md`.
- Lag-1 ACF `dmax` for the four Hadamard-family kernels (`hdev`,
  `mhdev`, `htotdev`, `mhtotdev`) bumped from `2` to `3`. Riley &
  Greenhall, *Power Law Noise Identification Using the Lag 1
  Autocorrelation* (EFTF 2004, paper #125), §6: "The dmax parameter
  should be set to 2 or 3 for an Allan or Hadamard (2 or 3-sample)
  variance analysis respectively." Hadamard-family estimators exist
  to resolve α ∈ {−3, −4} (FWFM / RRFM); capping the noise-ID
  iteration at `dmax = 2` prevents the algorithm from ever reporting
  those types and silently floors the identified α at `−2` (RWFM)
  even when the data has lower power-law content. `:htdev` inherits
  the fix transitively via `:mhdev`. `:adev`, `:mdev`, `:totdev`,
  `:mtotdev` correctly stay at `dmax = 2`.
- EDF stride factor `S` for the four overlapped variants
  (`:adev`, `:mdev`, `:hdev`, `:mhdev`) in
  `lib/SigmaTauStability/src/stats/edf.jl` corrected from `1`
  (non-overlapped convention) to `m` (overlapped convention). The
  Greenhall–Riley `M = 1 + ⌊S·(N − L) / m⌋` interpretation
  documented in AllanLab's `calculate_edf.m` has `S = 1` for
  non-overlapped and `S = m` for overlapped; we were computing
  EDFs under the non-overlapped convention and applying them to
  the overlapping deviation, producing artificially small EDFs
  (e.g. 127 instead of 8064 at `m = 64`, `N = 8192`) and the
  correspondingly wider χ² confidence bounds. `:tdev` and `:htdev`
  inherit the fix transitively via the `mdev` / `mhdev` factor
  scaling. The `:totdev` / `:htotdev` WPM/FLPM fallback paths in
  the same file still pass `S = 1`; left for a follow-up since
  they exercise a different code path.
- CI: `lib/SigmaTauBase` dropped from the test matrix. It is a
  types-only package with no `test/runtests.jl`; the previous
  workflow's `Pkg.develop(path="lib/SigmaTauBase")` line was
  papering over this by failing earlier with "same name or UUID
  as the active project". `lib/SigmaTauBase` is exercised
  indirectly by Stability/Ensemble (both `[deps]` it via
  `[sources.SigmaTauBase] path = "../SigmaTauBase"`) and by the
  umbrella `using SigmaTau` job, so the matrix entry is redundant.
  The redundant `Pkg.develop` line is also gone.
- CI: umbrella `using SigmaTau` job now uses `Pkg.instantiate()`
  instead of `Pkg.resolve()`. With no checked-in `Manifest.toml`,
  `resolve` in workspace mode failed to register the General
  registry on a fresh CI runner ("expected package Reexport
  [189a3867] to be registered"). `instantiate` triggers registry
  init.
- CI: `lib/SigmaTauEnsemble` test runner now guards the
  `legacy/julia/src` includes behind an `isfile` check. The legacy
  reference codebase is gitignored under `legacy/`, so CI checkouts
  don't have it; the legacy-parity testsets (`Phi/Q matrix parity`,
  `StandardKalmanFilter Parity (legacy_compat)`, `AD-clean`) are
  conditionally skipped when the directory is absent and a
  `@info "legacy/julia/src not present, skipping legacy-KF parity testsets"`
  message replaces the previous load-time `SystemError`. The
  PID / steering / TwoStateClock testsets still run unconditionally.
- `lib/SigmaTauStability/test/runtests.jl` "ADEV/HDEV across all 5
  power-law noise types" testset rtol relaxed from `1e-12` to
  `1e-11` for `_mdev_core` / `_mhdev_core` parity vs the legacy
  reference. On macOS x86_64 the two implementations agree
  bit-exactly (Δ_ol = 0 ULP for nearly all rows), but Linux x86_64
  LLVM picks a different SIMD reduction order and the values drift
  by ~10,000 ULPs (~4e-12 absolute) on this synthesised input —
  irreducible cross-platform codegen variance, not a math bug.
  3-way verification on the Stable32 fixture (ours / legacy /
  allantools) confirms the kernels are correct to ≤ 8.5e-14
  worst-case (full-precision allantools reference); 1e-11 is
  comfortable headroom.

### Changed

- `tools/regen_allantools_fixtures.py` now writes the CSV at `%.17e`
  (round-trip-exact Float64) instead of `%.6e` (~7 sig figs). The
  fixture at `reference/validation/allantools_out/allantools_data_full.csv`
  has been regenerated with full machine precision.
- `lib/SigmaTauStability/test/allantools_cross_validation.jl`
  ADEV/MDEV/HDEV/TDEV cross-check tightened from `rtol=1e-4` to
  `rtol=1e-11` against allantools, now that the fixture preserves
  full Float64. TOTDEV/HTOTDEV/MTOTDEV stay at their original
  boundary-policy floors (0.15 / 0.10 / 0.05). The B1 testset is now
  a machine-precision parity contract against an independent
  external implementation, not a smoke test.

### Added

- `detrend::Symbol` kwarg on the four total-family kernels and API
  wrappers (`totdev`, `mtotdev`, `htotdev`, `mhtotdev`). Each kernel
  exposes its canonical recipe and a `:linear` alternative that swaps
  the kernel's natural detrending for a per-window full LS slope+intercept
  while preserving the kernel's natural extension shape. `:legacy`
  aliases the pre-1.0 default of each kernel for backward compat.
- `_totdev_howe` helper implementing NIST SP1065 eqn 25 verbatim
  (Greenhall–Howe–Percival 1998 eq 3): no detrend, mean-flip endpoint
  reflection, sum over centers `n=2..N-1`. Cross-checked against
  Stable32's TOTDEV at `rtol=1e-4` (12 of 13 m values) and against
  allantools' `totdev` to ~7 significant figures on the same fixture.
  The single (m=512) Stable32 row that misses the rtol=1e-4 floor is
  identified by Stable32 as FLFM and carries an alpha-aware correction
  that diverges from the raw SP1065 value; documented in the test
  comment and `TODO.md`.
- allantools (Anders Wallin's Python library) wired in as a second
  external numerical reference alongside Stable32. New
  `tools/regen_allantools_fixtures.py` regenerates a row-aligned CSV
  at `reference/validation/allantools_out/allantools_data_full.csv`
  (one Sigma column per Stable32 (Type, AF) pair); new
  `lib/SigmaTauStability/test/allantools_cross_validation.jl` testset
  asserts the SigmaTau raw kernels against allantools output at
  `rtol=1e-4` for ADEV/MDEV/HDEV/TDEV and `rtol=0.05–0.10` for the
  total-family kernels (boundary-extension policy differs). The
  testset short-circuits silently when the fixture is absent so CI
  and pre-regen runs stay green. The Python script flushes per row
  and prints `[N/total] Type m=… → σ=…` progress lines because the
  total-family kernels are O(N²) on the 8192-sample fixture.

### Changed

- **Breaking:** `totdev` default detrend recipe is now `:howe` (SP1065
  eqn 25: no detrend, mean-flip endpoint reflection). Previous behavior
  (global LS detrend on top of the same reflection) is available via
  `totdev(...; detrend=:legacy)` or the `:linear` alias. Output values
  change for all τ; the new default matches allantools' raw `totdev` to
  ~7 significant figures and Stable32's `Total` column at `rtol=1e-4`
  (m=512 outlier excepted, see Added section). Bias correction
  `bias_correction(:totvar, ...)` is now correctly calibrated against
  the new `:howe` default — the previous `:legacy` path silently
  removed drift before the SP1065 bias factor was applied.
- **Breaking:** `mhtotdev` default detrend recipe is now `:greenhall`
  (per-window half-mean slope removal), aligning with the convention of
  MTOT and HTOT in the Hadamard-modified family. Previous behavior
  (per-window full LS) is available via `mhtotdev(...; detrend=:linear)`
  or `:legacy`. MHTOTDEV is novel to SigmaTau; no external numerical
  reference exists, so the recipe choice is a methodology decision
  rather than a parity contract.
- MTOTDEV and HTOTDEV default outputs are unchanged. Both keep
  `:greenhall` (per-window half-mean) as the default; the new `:linear`
  recipe is opt-in.
- MHTOTDEV bias-correction policy made explicit: `bias_correction`
  now short-circuits to `B = 1` for `var_type = :mhtot`, with the
  rationale promoted into the function docstring (FCS 2001 and
  NIST SP1065 publish no bias model for MHTOTDEV; Stable32 and
  AllanLab also treat it as unbiased). `mhtotdev` docstring
  rewritten to state this contract instead of the inline
  comment it carried before.
- `calculate_edf` docstring expanded with the WPM/FLPM (α ∈ {1, 2})
  fallback rationale for TOTDEV and HTOTDEV — the substitute
  ADEV/HDEV-style Greenhall–Riley formula is a documented policy
  choice, not a derivation, since SP1065 Table 9 and FCS 2001
  cover α ∈ {0, -1, -2} only.
- `Random` is no longer a hard dep of `SigmaTauStability`. The
  `_gen_powerlaw_phase` helper now uses the global RNG (callers seed via
  `Random.seed!`) so the convention "Random in [extras] for both
  subpackages" is preserved. `Random` stays in `[extras]` and the test
  target.
- Test runner under `lib/SigmaTauEnsemble/test/runtests.jl` now wraps the
  legacy `clock_model.jl` / `filter.jl` includes in a `LegacyKF` module
  and selectively imports the symbols it needs. This stops the legacy
  `step!`/`PIDController` from shadowing the new ones at test scope.
- `.gitignore` now excludes the `CLAUDE.md` and `AGENTS.md` agent-context
  briefs. They are local-only working documents and should not appear on
  the public repo.
- **Project workflow**: every shipped change must remove its TODO entry and
  add a CHANGELOG line under `## [Unreleased]` in the same commit. TODO.md
  rewritten to drop the now-stale "Author README", "Maintain CHANGELOG"
  bullets and re-prioritise remaining work into Critical / High / Medium /
  Low / Documentation buckets.

### Added

- PID steering controller in `SigmaTauEnsemble`. New `PIDController`
  struct (g_p / g_i / g_d gains + integral state) plus `step!(pid, x)` and
  `steer_to_correction(steer, ns, dt)`. `predict!(est, model, dt;
  steering=…)` now accepts an optional steering vector that's added to the
  propagated state mean — matches the legacy `filter_step!` semantics where
  a PID's last steer is folded in after Φ propagation.
- `examples/clock_steering.jl` — PID + Kalman walkthrough on a drifting
  clock; demonstrates ~600× residual-phase reduction vs the unsteered case.
- Power-law phase-noise synthesizer in `SigmaTauStability` at
  `src/noise/synth.jl`. `_gen_powerlaw_phase(α, N; tau0, rng)` generates an
  N-sample phase residual whose fractional-frequency PSD ∝ `f^α` via
  `f^(α/2)` shaping of white Gaussian noise (DC zeroed, integrated to phase).
  Non-exported helper. Requires an `AbstractFFTs` backend (e.g. `using
  FFTW`) loaded by the caller; `AbstractFFTs` is now a hard dep, `FFTW`
  is in test extras.
- MTOTDEV multi-noise validation extended from 3 → 5 power-law types
  (added FLPM and FLFM via the new synthesizer). New
  "ADEV/HDEV across all 5 power-law noise types" testset extends
  `_adev_core` / `_mdev_core` / `_hdev_core` / `_mhdev_core` kernel parity
  to the same alpha range. 120 new assertions; full Stability test count
  is now `339/339` passing.
- Stable32 cross-validation testset against the
  `reference/validation/stable32gen.DAT` fixture (8192 phase samples) and
  Stable32's published outputs (`stable32_data_full.csv`). Verifies the new
  kernels match Stable32's reported sigmas at:
  - `rtol=1e-4` for OADEV / Modified Allan / Overlapping Hadamard / Time
  - `rtol=0.05` for raw MTOTDEV (Stable32 reports unbiased; our API applies
    SP1065 B≈1.27 — kernel match is ~3% per the legacy comparison report)
  - `rtol=0.10` for raw HTOTDEV (~0.5% bias offset + boundary effects)
  - `rtol=0.15` for TOTDEV (close at short τ; documented boundary-reflection
    discrepancy at long τ)

  86 new assertions (85 row-driven + 1 sanity-count).
- Documenter.jl docs subproject under `docs/`. Ships skeleton v1 with
  theory / tutorials / reference / validation page tree. Reference
  pages auto-generated via `@docs` blocks. Stable32 comparison report
  migrated from `reference/validation/stable32out/` into
  `docs/src/validation/stable32.md`.
- DocStringExtensions integration across all three subpackages;
  struct docstrings on `PhaseData`, `FrequencyData`, `StabilityResult`
  with `$(TYPEDFIELDS)`.
- DocumenterCitations bibliography (`docs/src/refs.bib`) covering 20
  references; `@cite` markers wired into `adev` docstring and the
  homepage as a smoke test.
- `jldoctest` examples on `adev` and `mdev` (build-time asserted).
- `.github/workflows/Documentation.yml` deploying to GitHub Pages on
  push to `main` and on tags.
- `legdocs/` directory (gitignored) for source material lifted from
  the cross-language predecessor.
- MathJax3 math engine (replaces KaTeX default for full LaTeX support).
- Strict numerical parity testset for the eight stability kernels
  (`_adev_core`, `_mdev_core`, `_hdev_core`, `_mhdev_core`, `_totdev_core`,
  `_mtotdev_core`, `_htotdev_core`, `_mhtotdev_core`) against the legacy
  SigmaTau Julia reference. The legacy kernels are inlined verbatim under
  `lib/SigmaTauStability/test/legacy_kernels.jl` so the tests run on CI
  without depending on the gitignored `legacy/` tree. 52 new assertions
  at `rtol=1e-12`.
- MTOTDEV multi-noise validation: kernel + end-to-end pipeline checks
  against WPM, WHFM, and RWFM synthetic fixtures.
- TOTDEV/HTOTDEV WPM/FLPM EDF fallback: when `_coeff_totvar` /
  `_coeff_htot` return `(NaN, NaN)` for `α∈{1,2}` (since SP1065 Table 9 /
  FCS 2001 only cover `α∈{0,-1,-2}`), `calculate_edf` now falls back to
  the ADEV-style (`d=2`) or HDEV-style (`d=3`) Greenhall–Riley formula,
  yielding finite EDFs for every noise type instead of NaN.
- GitHub Actions CI: matrix on Julia 1.11 × {ubuntu, macOS} × all three
  subpackages, plus an umbrella `using SigmaTau` smoke job.
- `examples/quickstart.jl`: end-to-end walkthrough exercising
  `adev`/`mdev`/`tdev`/`hdev`/`totdev`, `FrequencyData` input, and a
  `ThreeStateClock` Kalman filter. Verified to run.

- `tdev(::PhaseData, m_values; …)` API wrapper. Wraps `mdev` and scales by
  `τ/√3`; CI bounds inherit MDEV's χ²/Gaussian limits scaled by the same
  factor. Now exported from `SigmaTauStability`.
- `FrequencyData` entry points for every stability deviation (`adev`, `mdev`,
  `tdev`, `hdev`, `mhdev`, `ldev`, `totdev`, `mtotdev`, `htotdev`,
  `mhtotdev`). Each converts via `_freq_to_phase` (`cumsum(y) · τ₀`) and
  delegates to the existing `PhaseData` method. Lives in
  `lib/SigmaTauStability/src/utils.jl`.
- `edf::Vector{Float64}` field on `StabilityResult`. Populated when
  `calc_ci=true`; empty when `calc_ci=false`. Lets users compute custom
  confidence intervals without re-running noise identification.
- `SigmaTauRecipesBaseExt` package extension at
  `ext/SigmaTauRecipesBaseExt.jl`. Provides a `RecipesBase.@recipe` for
  `StabilityResult` that draws a log-log τ–σ plot with optional
  `ci_lower`/`ci_upper` error bars; activates automatically when
  `RecipesBase` (or `Plots`) is loaded alongside `SigmaTau`.
- Boundary test for `identify_noise` at `N_eff ∈ {29, 31}` to lock in the
  `NEFF_RELIABLE` threshold.

### Changed

- `NEFF_RELIABLE` lowered from 50 → 30 in
  `lib/SigmaTauStability/src/noise/lag1.jl` per the legacy GEMINI.md §2
  mandate. `identify_noise` now switches from the B1-ratio fallback to the
  lag-1 ACF path at `N_eff = 30` instead of `N_eff = 50`.
- Plot recipes are no longer compiled into the umbrella module. The previous
  `src/PlotRecipes.jl` stub has been deleted; recipe code lives entirely in
  the `RecipesBase` package extension. The umbrella package no longer pulls
  any plotting dependency unless one is explicitly loaded.
- `lib/SigmaTau{Base,Stability,Ensemble}/Project.toml` declare a new
  `DocStringExtensions` dependency.
- `adev` and `mdev` docstrings rewritten with `$(SIGNATURES)`,
  `jldoctest` examples, and `[Greenhall2003](@cite)`.

### Removed

- `reference/validation/stable32out/comparison_report.md` and
  `comprehensive_comparison.md` (content migrated into
  `docs/src/validation/stable32.md`).

### Fixed / Verified

- Kalman filter parity now confirmed end-to-end: `Pkg.test()` against
  `SigmaTauEnsemble` reports **15/15 pass**, including the 4 `legacy_compat`
  parity assertions (phase, frequency, drift, P-history). The stale
  `test_output.log` from the prior session predated the
  `clamp_covariance_diag` wiring and is no longer representative.
- `Pkg.test()` against `SigmaTauStability` reports **46/46 pass**, covering
  the new `tdev` wrapper, `FrequencyData ↔ PhaseData` equivalence, the
  `edf` field, and the `NEFF_RELIABLE = 30` boundary at `N_eff ∈ {29, 31}`.
- `Random` was missing from both subpackages' `[extras]`; added so
  `Pkg.test()` resolves the standard-library dependency. Without this,
  `Random.seed!` calls inside test files raised `Package Random not found`
  before any assertion could run.

## [0.1.0] — 2026-05-07

Initial commit. Modularised legacy SigmaTau into three subpackages
(`SigmaTauBase`, `SigmaTauStability`, `SigmaTauEnsemble`) plus an umbrella
`SigmaTau` package. Workspace + path-source wiring per Julia 1.11 monorepo
semantics. MIT licensed.

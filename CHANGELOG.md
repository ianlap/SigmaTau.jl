# Changelog

All notable changes to **SigmaTau.jl** are tracked here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.2.0] — 2026-05-18

### Added

- **`noise_gen` — calibrated power-law clock-noise generator.** Public
  wrapper around the `_gen_powerlaw_y` spectral shaper that lets callers
  specify the noise mixture as a `Dict{Int, Float64}` keyed by SP1065
  α exponent (∈ {-2, -1, 0, 1, 2}). Two input modes, choose one:
  `sigma1[α] = σ_y(τ=τ₀)` (clock-spec-sheet style) or `h[α] = h_α`
  (PSD coefficient in `S_y(f) = h_α f^α`). Composite mixtures are drawn
  independently per α and summed. For WPM (α=2) and FPM (α=1) the
  `σ ↔ h` relation uses the Nyquist convention `f_h = 1/(2τ₀)`.
  Per-realization rescaling pegs each component's empirical
  `σ_y(τ=τ₀)` to the requested value exactly; the natural power-law
  slope carries to larger τ. Dispatches on the first type argument to
  return either `PhaseData` or `FrequencyData`
  (`noise_gen(PhaseData, …)` / `noise_gen(FrequencyData, …)`). Exported
  at the umbrella level. Internal `_gen_powerlaw_phase` is now a thin
  wrapper around the new `_gen_powerlaw_y` primitive — no change in
  RNG consumption order, so legacy tests that pin synthesized phase
  to a seed continue to match bit-for-bit.
- **`ClockEnsemble` model for joint Kalman time-scale estimation.** Stein
  2003 §V / Galleani–Tavella stacked-state formulation: the joint state
  concatenates per-clock states, Φ and Q are block-diagonal, and H selects
  the `N − 1` phase differences against a reference clock. The ensemble
  is itself an `AbstractClockModel`, so the existing `predict!` /
  `update!` / `prop!` loop on `StandardKalmanFilter` consumes it
  unmodified. `Base.:+` overloaded on `AbstractClockModel` for ergonomic
  construction (`ensemble = clockA + clockB + clockC`); homogeneous-only
  (mixing `TwoStateClock` and `ThreeStateClock` throws `ArgumentError`).
  Auto-derives Stein §VI–VII inverse-noise weights from each member's
  diffusion coefficients (`a_i ∝ 1/q1_i`, `b_i ∝ 1/q2_i`,
  `c_i ∝ 1/q3_i`); explicit `weights=` override supported. References:
  Stein 2003 *"Time Scales Demystified"* IFCS; Galleani–Tavella 2010
  *"Time and the Kalman filter"* IEEE CSM.
- **`prop!(est, model, dt; steering=nothing)`** — unconditional
  covariance propagation alongside `predict!` / `update!`. Advances
  `est.x ← Φ(dt) x` and `est.P ← Φ(dt) P Φ(dt)' + Q(dt)` regardless of
  `est.k`, never bumping the step counter. Powers shaded ±1σ
  holdover-band projections without disturbing live filter sequencing.
- **`state_transition(model, dt)`** / **`process_noise(model, dt)`**
  overloads on `AbstractClockModel`. The single-arg methods defer to the
  new dt-aware forms with `model.tau`, so legacy callers see no
  behaviour change. `predict!` and `prop!` now re-derive Φ and Q for the
  caller-supplied `dt`, covering arbitrary horizons rather than only the
  discretisation step.
- **IO at the umbrella level: file readers, detrend, gap fill, result
  round-trip.** New `src/io/` directory consolidating all data-side IO.
  `read_phase(path; …)` and `read_frequency(path; …)` parse 2-column (or
  N-column with `time_col=0`) text files via stdlib
  `DelimitedFiles.readdlm` and return `PhaseData` / `FrequencyData` with
  optional preprocessing (`scaling`, `detrend`, `fillgaps`).
  `detrend(::PhaseData/::FrequencyData; method=:linear)` covers
  `:linear`, `:endpoint`, `:mean`, and `:none`; multiple dispatch keeps
  the bare name `detrend` collision-free. `fillgaps(::PhaseData)` /
  `fillgaps(::FrequencyData)` port the Howe & Schlossberger
  reflect-and-FFT-filter imputation (PTTI 2009). `save_result` /
  `load_result` round-trip a `StabilityResult` to/from a tab-separated
  text file using stdlib I/O only — handles both `calc_ci=true` and
  `calc_ci=false`; CI columns are written as `NaN` when absent and
  reconstructed as empty vectors on load. `FFTW.jl` promoted from
  test-extra to core dep.
- **Three new stability metrics.**
  - **`tdev`** wraps `mdev` and scales by `τ/√3`; CI bounds inherit
    MDEV's χ²/Gaussian limits scaled by the same factor.
  - **`ttotdev`** wraps `mtotdev` with the same `τ/√3` rescaling for
    both `PhaseData` and `FrequencyData` inputs.
  - **`htdev`** (Hadamard time deviation) wraps `mhdev` with the
    `τ/√(10/3)` factor; the construction is original to this package
    and is not defined in Stable32 / allantools.
- **MTIE** (Maximum Time Interval Error) per ITU-T G.810. Public
  `mtie(::PhaseData, m_values)` + `mtie(::FrequencyData, …)` plus the
  `_mtie_core` kernel in `src/stab/core/mtie.jl`. Sliding-window
  peak-to-peak phase excursion implemented as an O(N) monotonic-deque
  scan (138× speedup over the naive double-loop at N=50000); output is
  bit-identical to the naive form. No CI fields — no published EDF /
  χ² model for the peak statistic.
- **PDEV** (parabolic deviation) per Vernotte–Lenczner–Bourgeois–Rubiola
  IEEE T-UFFC 63(4) 2016 and Vernotte 2020. Public
  `pdev(::PhaseData, m_values)` + `pdev(::FrequencyData, …)` plus
  `_pdev_core` kernel in `src/stab/core/pdev.jl`. Least-squares parabolic
  fit; PDEV(τ₀) ≡ overlapping ADEV(τ₀) at m=1 by construction. CI fields
  empty pending future EDF port.
- **Zero-arg convenience methods on every deviation API.** Calls like
  `adev(pd)`, `mdev(fd)`, `hdev(pd)`, …, `mtie(pd)`, `pdev(pd)` now
  resolve without a `m_values` argument. The default is an octave-
  spaced grid `1, 2, 4, …, 2^k` capped at the kernel's algorithmic
  m-max — derived from each `_*_core`'s `L`/`Ne` guard:

  | Deviation                                          | m_max          |
  |----------------------------------------------------|----------------|
  | `adev`, `totdev`, `pdev`                           | `(N − 1) ÷ 2`  |
  | `mdev`, `tdev`, `mtotdev`, `ttotdev`               | `N ÷ 3`        |
  | `hdev`, `htotdev`                                  | `(N − 1) ÷ 3`  |
  | `mhdev`, `htdev`, `mhtotdev`                       | `N ÷ 4`        |
  | `mtie`                                             | `N − 1`        |

  (HTOTDEV's general branch operates on `y = diff(x)` of length
  `N−1`, so its constraint matches HDEV's even though MTOTDEV — which
  runs on phase directly — uses `N ÷ 3`.)

  Shared via `SigmaTau.Stab._default_m_values(N, kernel::Symbol)`.
  All kwargs (`calc_ci`, `confidence`, `detrend`, `correct_bias`)
  pass through unchanged. The richer `taus` enum API
  (`AllTaus`/`Octave`/`Decade`/…) listed under
  [TODO.md](TODO.md) remains a future expansion; this is the
  one-call quality-of-life knob.
- `SigmaTau.Stab._phase_to_freq(::PhaseData) → FrequencyData` —
  internal companion to `_freq_to_phase`. Implements the canonical
  first difference `y[k] = (x[k+1] − x[k]) / τ₀` (length `N → N−1`).
  Test coverage verifies the round-trip identity (
  `_phase_to_freq ∘ _freq_to_phase` returns `y[2:end]`;
  `_freq_to_phase ∘ _phase_to_freq` returns `x[2:end] .- x[1]`) and
  that ADEV/MDEV/HDEV agree on the original phase record vs its
  differenced frequency form at matched `m_values`.
- **`FrequencyData` entry points for every deviation.** `adev`, `mdev`,
  `tdev`, `hdev`, `mhdev`, `htdev`, `totdev`, `mtotdev`, `ttotdev`,
  `htotdev`, `mhtotdev`, `mtie`, `pdev` all accept `FrequencyData` and
  convert via `_freq_to_phase` (`cumsum(y)·τ₀`) before delegating to the
  existing `PhaseData` method. Companion `_phase_to_freq(::PhaseData) →
  FrequencyData` implements the canonical first difference `y[k] =
  (x[k+1] − x[k]) / τ₀`.
- **PID steering controller.** `PIDController` struct (g_p / g_i / g_d
  gains + integral state) plus `step!(pid, x)` and
  `steer_to_correction(steer, ns, dt)`. `predict!(est, model, dt;
  steering=…)` and `prop!(…; steering=…)` accept an optional steering
  vector folded into the propagated state mean.
- **Power-law phase-noise synthesizer** at `src/stab/noise/synth.jl`.
  `_gen_powerlaw_phase(α, N; tau0, rng)` generates an N-sample phase
  residual whose fractional-frequency PSD ∝ `f^α` via `f^(α/2)` shaping
  of white Gaussian noise (DC zeroed, integrated to phase). Non-exported;
  requires an `AbstractFFTs` backend (e.g. `FFTW`).
- **Total-family kernels are multithreaded.** `_totdev_legacy` /
  `_totdev_howe` parallelise the outer m-loop. The six modified-total
  kernels (`_mtotdev_*`, `_htotdev_*`, `_mhtotdev_*`) parallelise the
  inner subsequence loop with chunk-private accumulators and a
  sliding-window inner reduction (`a += ext[j+km] − ext[j+(k−1)m]`)
  replacing the per-subsequence cumulative-sum buffer. A reused per-chunk
  scratch buffer pool eliminates the per-m, per-chunk allocation churn
  (was ~150 MB on a 2.6M-sample file). Net effect on `mtotdev(N=4·10⁵)`:
  ~415 s → ~290–340 s on 8 threads. FP determinism caveat: inner
  reduction reorders summation, so results may drift by a few ULPs across
  runs with different thread counts.
- **`StabilityResult.edf` field.** `Vector{Float64}`, populated when
  `calc_ci=true`, empty when `calc_ci=false`. Lets users compute custom
  confidence intervals without re-running noise identification.
- **`correct_bias::Bool=true` kwarg** on `totdev`, `mtotdev`, `htotdev`,
  `mhtotdev`. Default `true` preserves prior behavior; `false` returns
  the raw kernel value matching Stable32 and allantools. On `mhtotdev`
  the kwarg is a documented no-op — FCS 2001 / SP1065 publish no `B(α)`
  for MHTOT. The `calc_ci=false` fast path now also populates
  `noise_type` when `correct_bias=true`.
- **`detrend::Symbol` kwarg** on `totdev`, `mtotdev`, `htotdev`,
  `mhtotdev`. Each kernel exposes its canonical recipe plus a `:linear`
  alternative (per-window full LS slope+intercept) and a `:legacy` alias
  for the pre-1.0 default.
- **`_totdev_howe` helper** implementing NIST SP1065 eqn 25 verbatim
  (Greenhall–Howe–Percival 1998 eq 3): no detrend, mean-flip endpoint
  reflection, sum over centers `n=2..N-1`.
- **`identify_noise(x, m_values; …, detrend::Bool=true)` kwarg.** New
  default (`true`) keeps the per-m quadratic detrend allantools applies
  in `autocorr_noise_id`; `false` skips the polynomial fit so α matches
  a Stable32 reference fixture point-for-point. The 5σ outlier filter
  remains unconditional.
- **`SigmaTauRecipesBaseExt` package extension** at
  `ext/SigmaTauRecipesBaseExt.jl`. Provides a `RecipesBase.@recipe` for
  `StabilityResult` that draws a log-log τ–σ plot with optional
  `ci_lower`/`ci_upper` error bars; activates automatically when
  `RecipesBase` (or `Plots`) is loaded alongside `SigmaTau`.
- **Documenter.jl docs subproject under `docs/`.** Theory / Tutorials /
  Reference / Validation page tree. Reference pages auto-generated via
  `@docs` blocks; DocStringExtensions integration with `$(TYPEDFIELDS)`
  on `PhaseData` / `FrequencyData` / `StabilityResult`;
  DocumenterCitations bibliography (`docs/src/refs.bib`); MathJax3 math
  engine; `jldoctest` examples on `adev` and `mdev` build-time asserted.
  Theory section fills out `overview`, `allan_family`, `total_family`,
  `confidence`, `noise_id`, `ensemble_overview`, `kalman`, `steering`,
  and `clock_ensembles` from SP1065 / IEEE 1139-2022 / GR03 / FCS01 /
  RG04 / Stein 2003 / Tryon–Jones 1983 / Breakiron 2001 source material.
  PGFPlotsX backend for LaTeX-quality vector PDFs in `@example` blocks.
- **GitHub Actions CI.** Julia 1.11 × {ubuntu, macOS} matrix on the
  single package, plus `julia-actions/docdeploy` to GitHub Pages on
  push to `main` and on tags.
- **Eight Literate-driven tutorials** under `examples/`:
  `00_julia_for_metrologists`, `01_phase_data`, `02_compute_adev`,
  `03_kalman_single_clock`, `04_kalman_pid_steering`,
  `05_holdover_comparison`, `06_three_cornered_hat`,
  `07_clock_ensemble`. Auto-rendered into `docs/src/tutorials/` by
  Literate.jl at docs-build time.
- **Validation fixtures.** Stable32 cross-check against
  `reference/validation/stable32gen.DAT` (8192 phase samples) at
  `rtol=1e-4` for OADEV / Modified Allan / Overlapping Hadamard / Time
  and looser rtol on the total family (boundary-policy differs from
  Stable32). allantools cross-check at full Float64 precision
  (`rtol=1e-11`) for ADEV / MDEV / HDEV / TDEV via
  `tools/regen_allantools_fixtures.py`. Composite-noise fixture at
  `reference/validation/s32_5_12_26/` (N=32768 phase / frequency record
  from Stable32 NoiseGen with FPM→WFM→FFM→RWFM crossovers at τ ≈ 10,
  100, 1000): σ matches Stable32 to 5 sig figs across all 8 deviations
  (max |Δ| 4.3e-5).
- **Scaling-fit benchmark trio** under `benchmarks/bench/`:
  `scaling.jl` (SigmaTau, full public-API defaults; per-octave m-grid via
  `_default_m_values`), `scaling_allantools.py` (allantools, bare-kernel
  calls; capped N for the O(N²) `mtotdev`/`htotdev`), and
  `render_scaling.py` (cross-library comparison renderer). Each writes a
  machine-readable JSON (`scaling_sigmatau.json`,
  `scaling_allantools.json`). `predict(:kernel, N)` reads cached fits to
  estimate wall-time for arbitrary N. Power-law fit in log-log space
  (`T ≈ a · N^b`), with `R²` reported per kernel. Sub-call:
  `BenchmarkTools.@belapsed` on the Julia side, `timeit.Timer.repeat`
  (best of N) on Python.
- **`BenchmarkTools.@btime` quick-look scripts** at
  `benchmarks/bench/btime_sigmatau.jl` and `btime_allantools.py` —
  minimal-overhead per-kernel timings on synthetic randn input for
  `adev`/`mdev`/`hdev` (and their allantools equivalents
  `oadev`/`mdev`/`ohdev`). Synthetic input lets the scripts run anywhere
  without the `reference/clock_data/` fixtures the wall-clock
  `bench_sigmatau.jl` / `bench_allantools.py` scripts need.
  `BenchmarkTools` added to `benchmarks/bench/Project.toml` deps.
- **Legacy-kernel parity testset** (`test/stab/legacy_kernels.jl`)
  inlining the pre-restructure SigmaTau Julia reference kernels for the
  eight stability cores at `rtol=1e-12`. The legacy `legacy/julia/src/`
  tree is gitignored; CI checkouts skip the corresponding Kalman parity
  testsets behind an `isfile` guard with an `@info` message.
- **B1(N, μ) and R(n)(af, b) closed-form regression coverage** in
  `test/stab/runtests.jl`, locking the constants in `_b1_theory` and
  `_rn_theory` against the canonical allantools reference (`ci.py`,
  Wallin 2018 citing Howe 2000) at machine precision.
- **Top-level umbrella smoke test** at `test/umbrella_smoke.jl`.
  Verifies a bare `using SigmaTau` exposes every public symbol (Stab +
  Est + types), confirms `SigmaTau.Stab` / `SigmaTau.Est` are reachable
  as modules, exercises `FrequencyData` dispatch on every deviation, and
  locks the `ldev` ≡ `htdev` deprecated-alias contract.

### Changed

- **Restructured from a three-subpackage workspace
  (`SigmaTauBase` / `SigmaTauStability` / `SigmaTauEnsemble`) into a
  single registerable package with two submodules (`SigmaTau.Stab` and
  `SigmaTau.Est`).** Shared types (`PhaseData`, `FrequencyData`,
  `StabilityResult`) live at the top level. All previous exports
  continue to be re-exported from `SigmaTau`, so user code that imported
  via `using SigmaTau` keeps working unchanged. Code that explicitly
  imported `using SigmaTauBase`, `using SigmaTauStability`, or
  `using SigmaTauEnsemble` must switch to `using SigmaTau` (or
  `using SigmaTau.Stab` / `using SigmaTau.Est`). The `lib/` workspace
  tree is gone; `Project.toml` is a single-package manifest with no
  `[workspace]` / `[sources]` blocks; package UUID regenerated; `julia`
  compat bumped to `1.11`.
- **Public function `ldev` renamed to `htdev`.** The estimator and its
  formula are unchanged: `htdev` wraps `mhdev` and applies the
  `τ / √(10/3)` scaling. `StabilityResult.deviation_type` is now
  `:htdev`. The legacy `ldev` is now marked `Base.@deprecate ldev htdev`;
  callers receive a deprecation warning and should migrate. Delete after
  v0.2.0 is tagged (tracked in TODO.md).
- **Default confidence factor lowered from `0.95` to `0.683` (1-sigma)**
  on every public deviation API (`adev`, `mdev`, `hdev`, `tdev`, `mhdev`,
  `htdev`, `totdev`, `mtotdev`, `htotdev`, `mhtotdev`). Now exposed as
  the constant `SigmaTau.Stab.DEFAULT_CONFIDENCE`. The new value matches
  the Stable32, AllanLab, allantools, and Greenhall–Riley convention;
  the prior `0.95` made cross-tool CI overlays disagree by ~1.96× even
  when EDFs matched. Pass `confidence=0.95` explicitly to recover the
  prior default.
- **`totdev` default detrend recipe is now `:howe`** (SP1065 eqn 25: no
  detrend, mean-flip endpoint reflection). Previous behaviour (global LS
  detrend on top of the same reflection) is available via
  `totdev(...; detrend=:legacy)` or the `:linear` alias. Output values
  change for all τ; the new default matches allantools' raw `totdev` to
  ~7 significant figures and Stable32's `Total` column at `rtol=1e-4`.
  Bias correction `bias_correction(:totvar, ...)` is now correctly
  calibrated against the new `:howe` default.
- **`mhtotdev` default detrend recipe is now `:greenhall`** (per-window
  half-mean slope removal), aligning with MTOT and HTOT in the
  Hadamard-modified family. Previous behaviour (per-window full LS) is
  available via `mhtotdev(...; detrend=:linear)` or `:legacy`. MHTOTDEV
  is novel to SigmaTau; no external numerical reference exists, so the
  recipe choice is a methodology decision rather than a parity contract.
- **MHTOTDEV bias-correction policy made explicit.** `bias_correction`
  short-circuits to `B = 1` for `var_type = :mhtot`, with the rationale
  promoted into the function docstring (FCS 2001 and NIST SP1065 publish
  no bias model for MHTOTDEV; Stable32 and AllanLab also treat it as
  unbiased).
- **MTIE kernel reimplemented as a monotonic-deque sliding window** in
  `src/stab/core/mtie.jl`. Each m runs in O(N) total work via two
  parallel pre-allocated index deques (one for the running window max,
  one for the min). 138× speedup over the previous O(N·m) kernel at
  `N = 50 000` with 20 log-spaced m values; output bit-identical at
  `rtol = 1e-15`.
- **`predict!` now actually uses its `dt` argument.** Previously the
  signature accepted `dt` and silently ignored it in favour of
  `model.tau`. With the dt-aware overloads in place, `predict!` re-derives
  Φ and Q for the caller-supplied `dt`. Backwards-compatible at every
  existing call site (they all pass `dt = model.tau`); signature loosened
  from `dt::Float64` to `dt::Real`.
- **HTDEV CI scaling formally verified** (closes R-MED-5). Assertions in
  `test/stab/runtests.jl` lock the identity that
  `htdev = mhdev × τ/√(10/3)` propagates correctly through the χ² CI
  mapping: dev, ci_lower, and ci_upper all scale by the same factor, so
  relative-CI ratios `ci_lower / dev` and `ci_upper / dev` are identical
  to MHDEV's at every τ.
- **TOTDEV `:howe` allantools cross-validation tightened** from
  `rtol = 0.15` (legacy-kernel boundary-policy floor) to `rtol = 1e-7`.
  All 13 Total rows pass including m=512.
- **`tools/regen_allantools_fixtures.py` writes `%.17e`** (round-trip-
  exact Float64) instead of `%.6e`; ADEV/MDEV/HDEV/TDEV cross-checks
  tightened from `rtol=1e-4` to `rtol=1e-11`. TOTDEV/HTOTDEV/MTOTDEV stay
  at their original boundary-policy floors.
- **`docs/src/reference/{base,stability,ensemble}.md` renamed to
  `{types,stab,est}.md`**; `docs/make.jl` collapsed to a single
  `using SigmaTau` and one recursive `DocMeta.setdocmeta!` call.
- **Performance pass on `identify_noise` and `calculate_edf` — full-API
  call ~3.5× faster on fast kernels at N = 65 536.** Cumulative wall-time
  on `adev(pd; calc_ci=true)`: 1.945 ms → 0.55 ms; same magnitude across
  `mdev`/`tdev`/`hdev`/`mhdev`/`htdev`/`totdev`. Stable32 parity
  bit-identical against `reference/validation/s32_5_12_26/` (max σ
  rel-diff unchanged at ≤ 4.3 × 10⁻⁵; all 914 tests pass). Breakdown of
  contributing changes below.
- `identify_noise` (`src/stab/noise/lag1.jl`) refactored to cut
  allocation and per-m work. `_preprocess`'s 5σ outlier filter is now
  two-pass count-then-fill instead of allocating
  `z = abs.((x .- μ) ./ σ)` plus a boolean mask. `_detrend_quadratic`'s
  3×3 LS solve uses `StaticArrays` (`SMatrix{3,3,Float64}` /
  `SVector{3,Float64}`) instead of building heap matrices. `_lag1_acf`
  is three single-pass loops instead of allocating a centred copy and a
  broadcasted lag-product. `_simple_mdev` replaces the
  `cumsum(pushfirst!(copy(x), 0.0))` chain with one length-(N+1)
  prefix-sum buffer, mirroring the `_mdev_core` fix. Net: noise-ID time
  on a typical fast-kernel call drops ~2.1×.
- `calculate_edf` (`src/stab/stats/edf.jl`) hot path inlined — ~5×
  faster. `_compute_sw`/`_compute_sx`/`_compute_sz` are now `@inline`'d
  and accept `Float64` (was `Real`, dispatch-unstable). The inner
  closure `sx(u) = _compute_sx(u, F, alpha)` inside `_compute_sz` was
  preventing Julia from inlining the stencil; replaced with direct
  calls. Integer powers expanded to literal multiplications
  (`ta^5` → `ta*ta*…*ta`) and the reciprocal `1/F` is computed once
  outside the j-loop as `inv_F`.
- `_totdev_howe` (`src/stab/core/total.jl`) streams the SP1065 eqn 2
  mean-flip endpoint reflection on the fly rather than materialising a
  `(3N − 4)`-Float64 buffer. Inner n-loop split into three ranges —
  short scalar tails for `n ≤ m` and `n ≥ N − m + 1` (which need
  reflection), fully SIMD-vectorised central region for the bulk.
  Kernel-only allocation: 0.577 MiB → 0.005 MiB (constant in N; ~24·N
  bytes saved per call, ~72 MiB at N = 3 × 10⁶). `:legacy`/`:linear`
  code path unchanged.
- `_mdev_core` and `_mhdev_core` (`src/stab/core/allan.jl`,
  `src/stab/core/hadamard.jl`) replace
  `cumsum(pushfirst!(copy(x), 0.0))` (≈ 3 length-N allocations) with a
  single length-(N+1) buffer filled by an explicit scalar prefix-sum
  loop. `@simd` dropped since the carry chain blocks vectorisation
  either way. Cumulative bytes per call drops ~75 % on the synth bench
  (819 200 → 200 632 at N = 25 000).
- `_mtotdev_greenhall` half-mean slope (`src/stab/core/total.jl`) now
  amortised across overlapping windows: each chunk seeds its first
  window's `s1`, `s2` half-sums (cost O(m)) and updates them in O(1)
  per subsequent window (`s1 += x[n+hi] − x[n]`,
  `s2 += x[n+seg_len] − x[n+hi]`). Algorithmically correct but the
  wall-time savings are in run-to-run noise — the original `@simd`
  half-mean sums were already memory-bandwidth-bound on the same
  cached x[] reads the ext-build pass needs. Kept for code-clarity
  reasons; the irreducible O(N²) cost remains the per-window detrended
  extension build + slide-reduction.
- **Repository housekeeping.** `lib.bak/` (pre-restructure recovery
  snapshot) and `rough_changelog/` (six superseded implementation_plan
  drafts) deleted from the working tree. Three internal planning logs
  moved out of Documenter's source root into `docs/superpowers/plans/`.
  Top-level `validation/` renamed to `benchmarks/` via `git mv` to
  disambiguate from `reference/validation/` (numerical fixtures) and
  `docs/src/validation/` (doc pages). `scratch.jl` / `scratch.py` added
  to `.gitignore` under a new "Scratchpad files" section.
  `CLAUDE.md` and `AGENTS.md` agent-context briefs are now intentionally
  tracked in the repo.

### Fixed

- **`_lag1_acf` no-id guard was scale-dependent.** The threshold
  `ssx < eps(Float64) * length(x)` reduced to `σ² < eps`, tripping on
  any data with std below ~1.5e-8 regardless of N — phase records in
  seconds (typically 1e-9..1e-12) hit it and round-tripped through
  `_lag1_acf` → NaN → `:unknown`. Replaced with a scale-invariant
  relative guard (`ssx ≤ eps · ‖x‖²`). Cross-validated against
  `allantools.autocorr_noise_id 2024.06` — 8/9 m values match exactly on
  `reference/validation/stable32gen.DAT`.
- **B1 noise-ID fallback (`_noise_id_b1rn`) was using the wrong `N`.**
  The `N` argument to `_b1_theory` was `length(x_dec) − 2` (post-
  decimation count), which collapses the dynamic range of the theoretical
  B1 values at small `N_eff` and misclassifies long-τ points as steeper
  noise than they are. Per Howe 2005 *"Enhancements to GPS Operation
  Using a Total Hadamard Deviation"*, `N` is the input record's
  frequency-sample count (`length(x) − 1`). Independently, the variance
  of m-averaged frequencies was computed with the biased `N`-divisor
  estimator; Barnes / Howe B1 theory assumes the unbiased `N−1`-divisor
  (Bessel) variance. Both fixes align SigmaTau's B1 noise-ID to Stable32
  empirically.
- **Bias correction for `totdev` / `mtotdev` / `htotdev` was wrong in
  both convention and direction.** `bias_correction` returned a
  variance-scale ratio `B` but the three API call sites divided the
  deviation by full `B` (over-correcting by another factor of `√B`);
  separately, the `:htot` table inverted the FCS 2001 / Howe & Tasset
  2005 convention by storing `1/(1+a)` instead of `1+a`. After the fix
  the SP1065 unbias is applied as `σ ← σ/√B`, the `:htot` table matches
  Howe 2005 Table I exactly, and α coverage extended to {−4, …, 0}. The
  `bias_correction` docstring was rewritten to make the variance-scale
  convention explicit.
- **`_coeff_htot` EDF table was wrong by ~20 % at multiple α.** Values
  swapped to the canonical FCS 2001 / Howe & Tasset 2005 Table I
  `(b₀, b₁)` pairs. After the fix, HTOTDEV's EDF formula
  `(T/τ) / (b₀ + b₁·τ/T)` reproduces Stable32's reported EDF to better
  than 0.01 % for every AF in the paper's stated validity range
  `τ ≥ 16τ₀`. Closes the long-standing "HTOTDEV EDF off-by-one suspected"
  item.
- **`_coeff_mtot` EDF coefficients refit against Stable32 output.** The
  SP1065 / Stable32 manual values produce EDFs ~5–20 % below what
  Stable32 actually computes; the new values match Stable32 to better
  than 0.1 %. α=0 is well-determined from two AFs at `(1.330, 1.890)`;
  α=−1 and α=−2 are single-point fits with `c` assumed from SP1065,
  documented inline as interim (tracked in `TODO.md`).
- **EDF stride factor `S` corrected from `1` to `m`** for the four
  overlapped variants (`:adev`, `:mdev`, `:hdev`, `:mhdev`) in
  `src/stab/stats/edf.jl`. AllanLab's `calculate_edf.m` documents the
  Greenhall–Riley `M = 1 + ⌊S·(N − L) / m⌋` interpretation with `S = 1`
  for non-overlapped and `S = m` for overlapped. The non-overlapped
  convention applied to overlapping deviations was producing artificially
  small EDFs (e.g. 127 instead of 8064 at `m = 64`, `N = 8192`) and
  correspondingly wider χ² CIs. `:tdev` / `:htdev` inherit transitively.
- **EDF stride factor `S` corrected from `1` to `m`** in the WPM/FLPM
  fallback path for `:totdev` and `:htotdev`. Both operate on a stride-1
  phase record (Howe's reflected-boundary extension preserves overlap),
  so the overlapped EDF is the consistent choice.
- **Lag-1 ACF `dmax` for the four Hadamard-family kernels bumped from
  `2` to `3`.** Riley & Greenhall, *Power Law Noise Identification Using
  the Lag 1 Autocorrelation* (EFTF 2004), §6: "The dmax parameter should
  be set to 2 or 3 for an Allan or Hadamard (2 or 3-sample) variance
  analysis respectively." `:htdev` inherits the fix transitively via
  `:mhdev`.
- **HTOTDEV default m-grid capped at `(N − 1) ÷ 3`** instead of
  `N ÷ 3`. The HTOTDEV general branch operates on `y = diff(x)` of length
  `N − 1`, so its constraint matches HDEV's even though MTOTDEV — which
  runs on phase directly — uses `N ÷ 3`.
- **Theory page misattribution of HTDEV (formerly `ldev`) to IEEE
  1139-2022 Annex C removed.** IEEE 1139-2022 does not define HTDEV; the
  construction is original to this package. The Annex-C citation has
  been replaced with a provenance note in `theory/allan_family.md`.

### Removed

- **Stub types `RelativisticClock`, `UDFactorizedFilter`, and
  `KuramotoOscillator`** along with the throwing methods that made them
  inert. Dropped from `src/est/models/clocks.jl`,
  `src/est/estimators/filters.jl`, the umbrella export list,
  `test/est/runtests.jl`, `test/umbrella_smoke.jl`, and
  `docs/src/reference/est.md`. The `AbstractClockModel` and
  `AbstractEstimator` docstrings dropped the stub subtype callouts.
- **Speculative theory pages.** `theory/relativistic_clocks.md`,
  `theory/relativistic_corrections.md`,
  `theory/relativistic_frames_and_timescales.md`,
  `theory/lunar_pnt_systems.md`,
  `theory/ensembles_and_oscillator_networks.md`, and
  `theory/publications.md` deleted. The Stein 2003 time-scale equation
  and three-cornered-hat material that backs `ClockEnsemble` and
  `tutorials/06_three_cornered_hat.md` is consolidated into a new
  shorter `theory/clock_ensembles.md` page. The `theory/kalman.md` page
  was trimmed to the standard predict / update recursion plus the
  innovation / Kalman-gain / `prop!` discussion, dropping the U-D,
  generalised-ALS, adaptive, structured-KF, and Wu LTI-performance-bound
  sections (none of which back shipped code). The
  `theory/ensemble_overview.md` page dropped the
  frequency-jump-detection and clock-error-jumps sections (likewise
  unshipped).
- **Pre-record linear detrend from the noise-ID `_preprocess` step.**
  Neither Stable32 nor allantools' `autocorr_noise_id` apply a
  full-record polynomial fit before the per-m loop, so the prior
  SigmaTau-only step shifted α away from both external references on
  records with any polynomial structure. The 5σ outlier filter in the
  same helper is unchanged. The unused `_detrend_linear` helper was
  deleted in the same commit.
- **30 orphan bibliography entries** dropped from `docs/src/refs.bib`
  after the theory-page trim — none remained cited from any rendered
  page or source file.

### Deprecated

- `ldev` (forwarding alias for `htdev`). Marked with `Base.@deprecate`;
  delete after v0.2.0 is tagged.

## [0.1.0] — 2026-05-07

Initial commit. Modularised legacy SigmaTau into three subpackages
(`SigmaTauBase`, `SigmaTauStability`, `SigmaTauEnsemble`) plus an umbrella
`SigmaTau` package. Workspace + path-source wiring per Julia 1.11 monorepo
semantics. MIT licensed.

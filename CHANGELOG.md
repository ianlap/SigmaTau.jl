# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.6.0] - 2026-11-08

### Fixed

- EDF for `adev`, `mdev`, `tdev`, `hdev`, `mhdev` (and every caller of the
  Greenhall EDF) now follows the full Greenhall & Riley 2003 algorithm:
  exact basic sum only while J ≤ 100, Tables 1–3 beyond that, the
  `m′ = J_max/r` substitution, F = ∞ once m(d+1) > 100, and the white-PM
  closed form. Summing the whole series moved the EDF up to 3.5 % (CI bounds
  up to 0.25 %) away from Stable32 and allantools wherever J > 100.
- `htotdev` at `m = 1` (plain overlapping HDEV) no longer applies the
  TotHvar bias factor or the FCS 2001 table EDF; σ was 0.25 % high there on
  white FM.
- `FrequencyData` input is integrated to `N + 1` phase points
  (`x = [0; cumsum(y)]·τ₀`) as in SP1065, Stable32 and allantools. Every
  estimator on `FrequencyData` used to drop the first frequency sample;
  so did `sigmatau --freq`.
- `mtotdev` / `ttotdev` at odd `m ≥ 3` remove the subsegment trend with
  Stable32's half averages (first and last ⌊3m/2⌋ points); values move by
  0.2–0.4 % at m = 3 and now match allantools.

### Added

- Reference cross-check suite (`test/validation/`): every estimator on
  power-law records with α = +2 … −4, plus EDF, χ² intervals, noise ID,
  Welch spectra and dynamic maps, against allantools, scipy, Stable32's
  PHASE.DAT runs and the NBS 1000-point published values. A `Cross-check`
  workflow regenerates the allantools tables with a live install
  (`tools/validation/gen_allantools.py`); `tools/validation/report.jl`
  writes the agreement and discrepancy tables (`tools/validation/REPORT.md`).
- `identify_noise(...; b1_length=:averaged)`: the B1 fallback (fewer than
  30 averaged points) evaluates B1 theory at the averaged-point count
  instead of the full record length. Random-walk FM at long τ is found far
  more often (Monte Carlo, N = 1024, m = 40–128: 67 % vs 0 %); white FM
  somewhat less often (63 % vs 82 %). The default `:record` keeps
  Stable32's convention.

### Deprecated

Renames for 0.6. The old names keep working, with a deprecation warning,
until 0.7.

- `psd_y`, `psd_x`, `ssb_phase_noise` replace `Sy`, `Sx`, `L`, and their
  `taper=` keyword replaces `window=` (which now only means the dynamic
  deviations' window). `SpectralResult.window` is now `.taper`; `r.window`
  still reads it, with a warning.
- `ncornered_hat` replaces `nch`.
- `save(path, x)` writes records, results and suites, and `load(path)` reads
  any of them back (records through `read_phase` / `read_frequency`).
  `save_result`, `save_suite`, `load_result` and `load_suite` are deprecated.
- `stability(data; m=...)` replaces `taus=`, the same name the estimators use.
- `identify_noise(...; remove_drift=false)` replaces `detrend=false`, matching
  `mhtotdev`'s keyword and freeing the name of the `detrend` function.

### Changed

- The spectra's default `nperseg` (`psd_y`, `psd_x`, `ssb_phase_noise`) is now the largest power of two up to
  `N ÷ 4` (at least 16) instead of `min(N, 256)`, matching the explorer. Long
  records get finer low-frequency resolution by default; pass `nperseg=256`
  for the old spectra.
- `dadev` / `dhdev` `window` is now optional: the largest power of two up to
  `N ÷ 8` (at least 16), with `step = window ÷ 2` as before.
- `stability` takes `correct_bias` and `remove_drift` (for `mhtotdev`) as its only extra keywords; any other
  keyword is now a `MethodError` instead of being passed to some estimators
  and silently dropped for the rest.
- `stability` passes `confidence` through to `mtie` and `tierms` like the
  other estimators (both still report no intervals), and asking it for
  `:dadev` / `:dhdev` now says to call those functions directly instead of
  failing inside the τ grid.
- `remove_outliers(::PhaseData)` now cleans in frequency space: it
  differences the phase, gaps and fills the outlying frequency samples, and
  integrates back from `x[1]`. A phase step is now removed instead of being
  carried through with only the sample it lands on refilled. Existing `NaN`
  gaps in the phase record are filled first. `FrequencyData` is unchanged.
- `read_phase` / `read_frequency` parse columns with a small typed reader
  instead of `DelimitedFiles.readdlm`, take the median time step with a local
  helper instead of `Statistics.median`, and no longer box locals in the
  comment-header scan. Results are unchanged; the readers are now type-stable,
  and everything except FFTW-backed paths (`noise_gen`, `fillgaps`) passes
  `juliac --trim=safe` verification.
- **Breaking:** `mhtotdev` values change for `m ≥ 2`. The subsegment phase
  length is now `4m` (one full MHVAR span) instead of `3m + 1`. Bias and EDF
  coefficients were re-measured by the full Monte Carlo sweep
  (`tools/mc_mhtotdev.jl`, seed 20260531, R = 3000, N ≤ 32768); provenance in
  `tools/artifacts/mhtotdev_mc_full.json`.
- Repository trimmed for the public release: the extended theory pages, FAQ,
  plotting cookbook and migration guides moved to the `dev` branch; the docs
  build no longer needs a TeX install.
- `noise_gen` (and the internal `_gen_powerlaw_y`) now synthesizes power-law
  noise with the Kasdin–Walter discrete fractional filter instead of shaping
  FFT bins by `f^(α/2)`. The old shaper left white-PM phase slightly blue
  (lag-1 autocorrelation −0.2 instead of 0), which the noise identifier read
  as α ≈ 3 at τ₀ and reported as `:unknown`. Realizations for a given seed
  change. `sigma1` components are still calibrated empirically to σ_y(τ₀);
  `h` components are now scaled on the filter's low-frequency PSD so the
  generated spectrum is `h_α f^α` (before, h was converted to σ_y(τ₀) with the
  continuous-time formulas, which do not hold for sampled discrete processes).
- MHTOTDEV bias and EDF coefficients re-measured with the Kasdin–Walter
  synthesis (full sweep at `e3ba9da`, seed 20260531, R = 3000, N ≤ 32768;
  `tools/artifacts/mhtotdev_mc_full.json`). White FM is unchanged; the
  largest shifts are the white-PM EDF (`b, c` 1.687, 4.317 → 1.634, 3.833)
  and bias slope, and the random-walk EDF intercept (2.449 → 2.415).
- `mdev`, `tdev`, `mhdev` and `htdev` remove the least-squares phase line
  before forming their prefix sums (an exact invariance of the modified
  estimators). On a record with an uncorrected 1e-9 frequency offset the
  kernels previously lost about seven significant digits to cancellation;
  they now agree with BigFloat arithmetic to 1e-12.
- Thêo1's ThêoBR bias correction (`theo1`, `theoh`) computes Thêo1 at 48
  points of the eq. 33 ladder and interpolates the rest in log–log, within
  0.5 % of the full average (0.25 % on the deviation); `theo1` at N = 65 536
  drops from minutes to about 1 s. Short ladders are still averaged in full.
- The total family, Thêo1 and ThêoBR run their averaging factors on all
  Julia threads (`julia -t auto`); `mtotdev` / `htotdev` / `mhtotdev` at
  N = 65 536 take about 1 s on four threads instead of 3 s. The compiled
  `sigmatau` core uses every CPU core.
- `plot(::DynamicStabilityResult)` draws the SP1065 waterfall (one σ(τ)
  curve per window, stacked along t) and no longer offers a heatmap.

### Added

- Python bindings, the `sigmatau` package in `python/`. Every call runs in an
  embedded Julia through juliacall, so results match the Julia API exactly.
  The names and keywords are the 0.6 ones (`psd_y`, `ncornered_hat`,
  `save`/`load`, `m=`), plus `deviation(name, ...)` and `estimators()`.
  Tested by a pytest suite in its own CI workflow.
- `sigmatau`, a compiled command line in `cli/` (package `SigmaTauCLI`, built
  with `juliac --trim=safe`): `analyze`, `check`, `spectrum`, `dynamic`,
  `completion` and `version`, CSV or versioned JSON output
  (`schema_version: 1`), documented exit codes and shell completion for bash,
  zsh, fish and PowerShell. It replaces the explorer's `sigmatau-cli`
  (`app/web/cli`), and the explorer now runs it with `--format json`. The `Web core` workflow now builds and tests it. See
  `docs/src/guide/cli.md`.
- `sigmatau --remove-drift none|linear|quadratic`: fits a line (frequency
  offset) or a parabola (offset and linear frequency drift) to the phase over
  the non-gap samples and removes it; the JSON reports `removed_offset` and
  `removed_drift_per_day`.
- Explorer: a noise-identification strip (one per clock) and power-law slope
  guides on the Stability plot, confidence bands by default with a plot toggle
  and a bias-correction note in the caption, drift-removal controls with the
  phase residual on the Data page, and a "Compare at τ" chart under the
  multi-clock Dynamic maps (each clock's deviation over time at one τ).
- `StabilityResult` carries the record's `tau0`, sample count `n` and the
  interval `confidence` (`NaN` without intervals), and result files store
  them. Files from earlier versions load with `NaN` / 0. The seven-argument
  constructor still works and leaves them unknown.
- `estimators()`: the registry of the 16 estimators (name, label, family,
  σ_y or σ_x axis, whether `correct_bias` applies, one-line description), for
  menus, command lines and docs.
- `deviation(name, data, m=Octave; ci, confidence, correct_bias)`: run an
  estimator chosen at run time; written as plain branches so it compiles under
  `juliac --trim`. `stability` now runs through it and accepts names as
  strings too.
- `tau_mode("half-octave")` and friends: parse a τ-grid name into a `TauMode`,
  ignoring case and separators.
- `summarize(data)`: samples, τ₀, span, `NaN` gaps, 5σ outlier indices and
  linear frequency drift per day, as the explorer's Data page shows them.
- `dadev(pd)` / `dhdev(pd)` work with no grid and no window: the octave grid
  and a default window (see below).
- Every estimator (and `dadev` / `dhdev`) accepts any integer vector of
  averaging factors: `adev(pd, 1:10)`, `adev(pd, 2 .^ (0:6))`, `Int32` or view
  vectors work like a `Vector{Int}`, as do `stability(pd; m=1:10)`.
- `PhaseData` / `FrequencyData` accept any real vector (views, ranges,
  integers), copied into a `Vector` of its float type.
- `length(::PhaseData)` / `length(::FrequencyData)`.
- `identify_noise(data::PhaseData, m)` and `identify_noise(data::FrequencyData, m)`,
  with `m` an integer vector or a `TauMode`.
- `bench/`: a reproducible SigmaTau-vs-allantools timing harness.
- `app/web/`: a web explorer (HTML, CSS, Plotly.js) over an analysis core
  compiled with `juliac --trim`, and `app/desktop/`, a Tauri window around it.

### Removed

- The GLMakie desktop explorer (`app/StabilityExplorer.jl`), its Bonito
  browser edition (`app/StabilityWeb.jl`) and their shared `app/ExplorerCore.jl`
  and environment, replaced by `app/web/`.
- The unused `DelimitedFiles` dependency (the readers have their own parser).

## [0.5.0] - 2026-06-09

Initial public release.

- Deviations: `adev`, `mdev`, `tdev`, `hdev`, `mhdev`, `htdev`, `totdev`,
  `mtotdev`, `ttotdev`, `htotdev`, `mhtotdev`, `theo1`, `theoh`, `mtie`,
  `tierms`, `pdev`; dynamic `dadev` / `dhdev`; `stability` for a suite in one
  call; `nch` N-cornered hat.
- Per-τ lag-1 ACF noise identification (B1 / R(n) fallback), Greenhall–Riley
  EDF, χ² confidence intervals at `DEFAULT_CONFIDENCE = 0.683`, total-family
  bias correction, and Monte-Carlo-calibrated MHTOTDEV bias and EDF.
- `StreamingStability`: O(1)-per-sample real-time `adev`/`mdev`/`hdev`/`mhdev`
  (Dobrogowski & Kasznia 2007), exact against the batch kernels.
- `noise_gen` calibrated power-law synthesis; `Sy`, `Sx`, `L` Welch spectra.
- IO: `read_phase` / `read_frequency` with header, delimiter and unit options;
  `detrend`, `fillgaps` (Howe–Schlossberger), `find_outliers` /
  `remove_outliers`; `save` / `load_result` / `load_suite` round-trip files.
- `PhaseData` / `FrequencyData` with provenance, compact REPL display,
  `TauMode` grids (`AllTaus`, `Octave`, `HalfOctave`, `QuarterOctave`,
  `Decade`, `HalfDecade`), Plots recipes and Tables.jl support as package
  extensions, PrecompileTools workload.
- Validation fixtures against Stable32 and allantools; `CITATION.cff`.

# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Deprecated

Renames for 0.6. The old names keep working, with a deprecation warning,
until 0.7.

- `psd_y`, `psd_x`, `ssb_phase_noise` replace `Sy`, `Sx`, `L`, and their
  `taper=` keyword replaces `window=` (which now only means the dynamic
  deviations' window). `SpectralResult.window` is now `.taper`; `r.window`
  still reads it.
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
- `plot(::DynamicStabilityResult)` draws the SP1065 waterfall (one σ(τ)
  curve per window, stacked along t) and no longer offers a heatmap.

### Added

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
  vectors work like a `Vector{Int}`. `stability(pd; taus=1:10)` already did.
- `PhaseData` / `FrequencyData` accept any real vector (views, ranges,
  integers), copied into a `Vector` of its float type.
- `length(::PhaseData)` / `length(::FrequencyData)`.
- `identify_noise(data::PhaseData, m)` and `identify_noise(data::FrequencyData, m)`,
  with `m` an integer vector or a `TauMode`.
- `bench/`: a reproducible SigmaTau-vs-allantools timing harness.
- `app/StabilityWeb.jl`: a browser edition of the explorer (Bonito + WGLMakie)
  with native HTML controls for screen readers and keyboard use; shares
  `app/ExplorerCore.jl` with the window.
- `app/StabilityExplorer.jl`: a GLMakie desktop application covering file
  import, synthesis, preprocessing, the deviation suite with χ² intervals,
  the per-τ table, spectra, dynamic deviation waterfalls, and figure, table
  and record export.

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

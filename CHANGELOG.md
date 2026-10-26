# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed

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
- `mdev`, `tdev`, `mhdev` and `htdev` remove the least-squares phase line
  before forming their prefix sums (an exact invariance of the modified
  estimators). On a record with an uncorrected 1e-9 frequency offset the
  kernels previously lost about seven significant digits to cancellation;
  they now agree with BigFloat arithmetic to 1e-12.
- `plot(::DynamicStabilityResult)` draws the SP1065 waterfall (one σ(τ)
  curve per window, stacked along t) and no longer offers a heatmap.

### Added

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

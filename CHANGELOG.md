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

### Added

- `bench/`: a reproducible SigmaTau-vs-allantools timing harness.
- `app/StabilityExplorer.jl`: a Pluto notebook GUI covering file import,
  preprocessing, the deviation suite, plotting and export.

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

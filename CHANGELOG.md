# Changelog

All notable public changes to **SigmaTau.jl** are tracked here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

No changes yet.

## [0.5.0] - 2026-06-09

Initial public release.

### Added

- Flat `SigmaTau` API with `PhaseData`, `FrequencyData`, `StabilityResult`,
  `StabilitySuite`, and `SpectralResult`.
- Stability estimators for Allan, Modified Allan, Hadamard, total-family,
  MTIE, and parabolic deviation workflows:
  `adev`, `mdev`, `tdev`, `hdev`, `mhdev`, `htdev`, `totdev`, `mtotdev`,
  `ttotdev`, `htotdev`, `mhtotdev`, `mtie`, and `pdev`.
- `stability` for computing an ordered suite of deviations in one call.
- Tau-grid selectors: `AllTaus`, `Octave`, `HalfOctave`, `QuarterOctave`,
  `Decade`, `HalfDecade`, and `tau_values`.
- Greenhall/Riley EDF and chi-squared confidence intervals where a published
  model exists, with `DEFAULT_CONFIDENCE = 0.683`.
- Lag-1 ACF / B1 / R(n) noise identification and calibrated power-law noise
  generation through `identify_noise` and `noise_gen`.
- Spectral estimators `Sy`, `Sx`, and `L` using a Welch PSD core.
- File IO helpers for reading phase/frequency records, detrending, gap filling,
  and round-tripping result/suite TSV files.
- Optional plotting recipes through `RecipesBase` and row-table support through
  `Tables.jl` package extensions.
- Validation fixtures and tests against Stable32 and allantools, plus
  legacy-kernel parity tests for estimators without an external implementation.
- Tracked MHTOTDEV Monte Carlo fit provenance at
  `tools/artifacts/mhtotdev_mc_full.json`.
- Migration guides for Stable32 and allantools users, a plotting cookbook, and
  an FAQ/troubleshooting page.
- Tutorials on characterizing a drifting clock with the Hadamard family
  (`hdev`/`mhdev`/`htdev`/`mhtotdev`) and on importing real-world data files.
- `CITATION.cff` and a README citation section.

### Changed

- Pre-public API cleanup settled on `ci` as the confidence-interval keyword and
  `htdev` as the canonical Hadamard time-deviation name.
- `PhaseData` and `FrequencyData` default `tau0` to `1.0`.
- Result types have compact REPL display methods instead of dumping full field
  arrays.
- Public docs, examples, validation notes, and GitHub templates were trimmed for
  the first public release.

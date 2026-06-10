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
- `mhtotdev` removes the record's least-squares frequency drift by default
  (`remove_drift=true`): the total kernel is not drift-immune (the boundary
  extension re-admits residual drift), per-window drift removal was measured
  to damage the statistic in both the phase and frequency domains, and the
  bias/EDF coefficients are calibrated with the global removal in the loop.
- `detrend` gains a `:quadratic` method (least-squares drift removal on
  phase data).
- `StabilityResult.neff`: the number of analysis windows the kernel averaged
  at each τ (the Stable32 "#" / allantools `ns` analog), derived from each
  kernel's loop bounds and populated even when `ci=false`; 0 where the kernel
  returns NaN. The result/suite TSV format gains a matching `neff` column
  (files written without it load with zeros), and the Tables.jl extension
  emits it as a column.
- Exported `ci_lower(r)` / `ci_upper(r)` accessors returning the confidence
  bounds as plain `Vector{Float64}`s (empty-in, empty-out).
- Migration guides for Stable32 and allantools users, a plotting cookbook, and
  an FAQ/troubleshooting page.
- Tutorials on characterizing a drifting clock with the Hadamard family
  (`hdev`/`mhdev`/`htdev`/`mhtotdev`) and on importing real-world data files.
- `CITATION.cff` and a README citation section.

### Changed

- Pre-public API cleanup settled on `ci` as the confidence-interval keyword and
  `htdev` as the canonical Hadamard time-deviation name.
- `StabilityResult` stores confidence bounds as a single `ci` field of
  `(lo, hi)` named tuples (`r.ci[i].lo` / `r.ci[i].hi`), replacing the
  separate `ci_lower` / `ci_upper` vector fields; the on-disk TSV format
  keeps flat numeric `ci_lower` / `ci_upper` columns, as does the Tables.jl
  extension.
- `PhaseData` and `FrequencyData` default `tau0` to `1.0`.
- Result types have compact REPL display methods instead of dumping full field
  arrays.
- Public docs, examples, validation notes, and GitHub templates were trimmed for
  the first public release.

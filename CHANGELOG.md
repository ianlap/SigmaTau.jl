# Changelog

All notable public changes to **SigmaTau.jl** are tracked here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

No changes yet.

## [0.5.0] - 2026-06-09

Initial public release.

### Added

- Thêo1 family per NIST SP1065 §5.2.15–5.2.16: `theo1` (eq. 30 with the
  Howe & Peppler 0.75 normalization, reported at the effective
  `τ = 0.75·m·τ0`; `correct_bias=true` default applies the ThêoBR
  automatic bias removal of eq. 33, `correct_bias=false` gives raw
  Thêo1) and `theoh` (the eq. 34 ADEV/ThêoBR hybrid with the 20 %-of-T
  crossover). EDF/CIs use the SP1065 Table 3 per-noise-type Thêo1
  formulas (ADEV EDF on the ThêoH ADEV segment). Even-m grid handling in
  `tau_values`/`_default_m_values` (`:theo1` grids round to even and
  dedupe; `:theoh` grids are τ-grid factors). Raw `theo1` is
  cross-validated against allantools 2024.06 to machine precision;
  allantools reports the same values at `τ = m·τ0` (resolved in favor of
  the SP1065/Stable32 effective-τ convention) and has no ThêoBR/ThêoH
  counterpart. Both integrate with `stability` (`devs=(:theo1, :theoh)`).
- `tierms` — RMS Time Interval Error per NIST SP1065 §5.2.18 eq. 37
  (`√(Σ(x[i+m]−x[i])²/(N−m))`), the estimator Stable32 and allantools
  (`tierms`) compute. A σ_x quantity reported at `τ = m·τ0` with the
  `N − m` span count in `neff`; like `mtie` it has no published EDF /
  χ² model, so `ci`/`edf`/`noise_type` stay empty. Cross-validated
  against allantools 2024.06 to ≤2e-15 relative on the Stable32
  fixture. Integrates with `stability` (`devs=(:tierms,)`).
- `StreamingStability` — real-time streaming accumulators for
  `adev`/`mdev`/`hdev`/`mhdev` per the Dobrogowski–Kasznia 2007 IEEE FCS
  scheme (running sums of squared second differences, eqs. 6–9; the
  inner-sum/overall-sum form for the modified family, eqs. 10–14;
  generalized to third differences for the Hadamard pair). `push!` feeds
  one phase sample in O(1) per averaging factor, `append!` feeds chunks,
  `snapshot` returns a `StabilityResult` (point estimates + `neff`; no
  CI — on-demand EDF/CI is future work), `nsamples` the stream length.
  Memory is a ring buffer of ≤ `4·maximum(m_values)+1` samples. At every
  sample count the streamed estimates match the batch kernels exactly
  (equivalence locked to rtol 1e-10 in test/stab/streaming.jl). Totals
  and Thêo cannot stream (whole-record extension/sampling).
- `nch` — N-cornered-hat noise separation: recovers each clock's
  individual deviation from the full set of pairwise comparisons
  (`σ̂²ᵢ = (Rᵢ − S/(N−1))/(N−2)`), with the classic Gray–Allan
  three-cornered hat as the `N = 3` case. Takes an upper-triangular
  matrix of pairwise `StabilityResult`s (or the three results directly
  in A−B, B−C, C−A order), keeps the input `deviation_type`, reports
  non-positive variance estimates as NaN, and propagates the
  elementwise-minimum `neff`; no CI is fabricated. Reproduces the
  manual solution of the three-cornered-hat tutorial
  (examples/06) to machine precision.
- `dadev` / `dhdev` — dynamic (time-resolved) Allan and Hadamard
  deviations (Galleani & Tavella 2009; McKelvy et al. 2025): a window of
  `window` phase samples slides across the record in `step`-sample
  increments (default `window ÷ 2`) and the overlapping ADEV / HDEV
  kernel is evaluated inside each window, giving a 2-D σ(t, τ) map that
  localizes non-stationary stability events a static deviation averages
  over. Returns the new `DynamicStabilityResult` (window-center times
  `t`, `tau`, a `t × τ` `dev` matrix with NaN at unsupported averaging
  factors, `window`, `tau0`); no CI machinery — the literature gives no
  EDF model for the time-resolved map. `TauMode` grids clamp to the
  window length, and a `window = N` map reproduces the static
  `adev`/`hdev` values exactly. The RecipesBase extension renders the
  map as a log10(σ) heatmap over (t, τ) with a log-scale τ axis.
- Flat `SigmaTau` API with `PhaseData`, `FrequencyData`, `StabilityResult`,
  `StabilitySuite`, and `SpectralResult`.
- Stability estimators for Allan, Modified Allan, Hadamard, total-family,
  Thêo, MTIE, TIE rms, and parabolic deviation workflows:
  `adev`, `mdev`, `tdev`, `hdev`, `mhdev`, `htdev`, `totdev`, `mtotdev`,
  `ttotdev`, `htotdev`, `mhtotdev`, `theo1`, `theoh`, `mtie`, `tierms`,
  and `pdev`.
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
- `PhaseData` / `FrequencyData` carry a `source` provenance field: `"user"`
  for directly constructed records (the default — existing constructors are
  unchanged), the originating file path for records from `read_phase` /
  `read_frequency`. `detrend`, `fillgaps`, and `remove_outliers` propagate
  it, and the compact data-record displays show it when it is not `"user"`.
- `save(path, data)` writes `PhaseData` / `FrequencyData` as a plain
  two-column text file (sample time `(i-1)·τ₀` in seconds, value) behind a
  `# SigmaTau <kind> data` / `# source:` / `# tau0:` comment header, e.g. to
  hand `noise_gen` output to software outside SigmaTau. `read_phase` /
  `read_frequency` skip leading `#` comment lines automatically and pick up
  `# tau0:` when the keyword is omitted, so saved files round-trip with no
  arguments. `save` is also the canonical generic for results and suites —
  `save(path, r)` / `save(path, suite)` delegate to `save_result` /
  `save_suite`, which remain available as aliases.
- `find_outliers` / `remove_outliers`: Stable32-style median-absolute-
  deviation outlier check (flag samples where
  `|yᵢ − median| > nsigma·MAD/0.6745`, default `nsigma=5`; phase records are
  tested on their first differences, flagging the sample each step lands
  on), with removal imputing the flagged samples through the Howe
  `fillgaps` machinery.
- PrecompileTools workload covering every public entry point: the first
  call in a fresh session runs in milliseconds instead of paying full JIT
  compilation.
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

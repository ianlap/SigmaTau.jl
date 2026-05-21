# SigmaTau.jl — Roadmap

Working list of outstanding engineering work. Items move from this file
to `CHANGELOG.md` as soon as they land — every shipped change should
remove the matching entry here and add one under `## [Unreleased]` in
the changelog in the same commit.

> **Audit date**: 2026-05-21 (post-v0.3.0 cut).

---

## 🟡 Correctness / completeness

- [ ] **PDEV EDF / χ² confidence model.** Vernotte 2015 / 2020 derive
  the parabolic-variance EDF in closed form for the five canonical
  power-law noises; port the table into `_coeff_pdev` and wire CI into
  the `pdev` API (currently returns empty `noise_type` / CI / EDF
  vectors).
- [ ] **MHTOTDEV bias / EDF Monte Carlo per detrend recipe.** Synthesize
  known-noise via `_gen_powerlaw_phase` for each α; compute MHTOT and
  MHDEV; the ratio yields the bias factor `B(α)`. Re-fit `_coeff_mhtot`
  empirically per detrend recipe (`:greenhall`, `:linear`). Track
  per-recipe — EDF is recipe-specific.
- [ ] **TOTDEV m=512 Stable32 quirk follow-up.** The Stable32
  cross-validation testset skips the m=512 row because Stable32 reports
  a value ~1.5 % larger than the raw SP1065 result (allantools agrees
  with our `:howe`). Stable32 identifies that row as FLFM (α=−1) and
  appears to apply an α-aware correction. Either confirm via
  `legdocs/vendor/` how Stable32 derives the reported σ at FLFM-tagged
  rows and apply the matching correction in the test, or document the
  divergence as an irreducible Stable32-vs-SP1065 policy gap.
- [ ] **MTOTDEV EDF coefficients for α=−1 and α=−2 verification.**
  Current `_coeff_mtot` values for those two α are single-point fits
  against Stable32's `s32_5_12_26` fixture (with `c` assumed from
  SP1065). To pin both `b` and `c` independently, capture Stable32
  EDF dumps at one additional AF per α — AF=200 or AF=400 for α=−1,
  and any record dominated by RWFM at a second AF for α=−2 (the
  current fixture only has α=−2 at AF=1000). Replace the single-point
  fits with two-point fits in `src/stab/stats/edf.jl::_coeff_mtot`.
- [ ] **Modified-total kernel parity rtol under multi-thread runs.**
  Per the threading note in the v0.2.0 changelog, the inner `@threads`
  reduction reorders summation, so kernel-vs-legacy parity tests at
  `rtol = 1e-12`/`1e-11` may need to drop to `~1e-9` for the
  modified-total kernels on `--threads auto` CI runners. Verify on a
  multi-thread CI run and loosen only the testsets that actually drift.

---

## 🟡 Post-v0.3.0 metrics (spectral)

- [ ] **`Sy(::PhaseData | ::FrequencyData)` — fractional-frequency
  PSD.** One-sided power spectral density of `y(t)` per IEEE 1139-2022
  §3.4, units `1/Hz`. Welch's method on the fractional-frequency
  sequence (after `_phase_to_freq` for `PhaseData` input) with a
  Hann/Hamming window, configurable segment length and overlap, and
  one-sided normalization so that `∫₀^{f_h} S_y(f) df` equals the
  variance of `y` over the analysis bandwidth. Returns a new
  `SpectralResult` (frequency vector + PSD vector + window/segment
  metadata) that mirrors `StabilityResult`'s self-describing layout.
  Validation: closes the loop with `noise_gen` — synthesize WPM /
  WFM / FFM / RWFM at known `h_α`, recover the slope and offset to
  within Welch's expected variance bound.
- [ ] **`Sx(::PhaseData | ::FrequencyData)` — phase PSD.** Same
  estimator pipeline applied to phase residuals (or `_freq_to_phase`
  of a frequency input). Units `s²/Hz`. Theoretical relation
  `S_x(f) = S_y(f) / (2πf)²` (IEEE 1139-2022 §3.3) → testset asserts
  the two estimators agree under that conversion to within
  ensemble-Welch tolerance.
- [ ] **`L(::PhaseData; f_carrier)` — single-sideband phase noise
  ℒ(f).** IEEE 1139-2022 §3.5: `ℒ(f) = (1/2) S_φ(f)` in dBc/Hz,
  with `S_φ(f) = (2π f_carrier)² S_x(f)` converting our
  seconds-valued phase record `x(t)` into radians at the user-
  supplied carrier frequency. Required kwarg `f_carrier` (Hz) — no
  meaningful default. Returns the same `SpectralResult` type with
  `units=:dBc_per_Hz`. Useful for direct comparison against
  spec-sheet plots from oscillator vendors. Cross-check via
  `noise_gen`: WPM at known σ_y(τ₀=1) yields a flat `ℒ(f)` whose
  level matches the analytic formula in SP1065 Table 3.

Common machinery for all three:

  - New file `src/stab/spectral.jl` (or `src/stab/spectral/welch.jl`
    + a thin api wrapper) with a Welch `_welch_psd(y, fs; nperseg,
    noverlap, window)` core that uses `FFTW` (already a dep).
  - New `SpectralResult` type alongside `StabilityResult` in
    `src/types/spectral_result.jl`, exported from the umbrella.
  - Plot recipe in `ext/SigmaTauRecipesBaseExt.jl` for log-log PSD
    overlays with α-slope reference lines.
  - Tests under `test/stab/spectral.jl`: round-trip with
    `noise_gen`, `S_x ↔ S_y` consistency, ℒ(f) carrier-scaling
    identity, segment/overlap edge cases.

---

## 🟢 Polish

- [ ] **More `examples/`** — Literate pipeline currently ships
  `00_julia_for_metrologists`, `01_phase_data`, `02_compute_adev`,
  `06_three_cornered_hat`. Next batch candidates:
  - A noise-ID walkthrough using `noise_gen` to synthesize a known
    composite α-mixture, then `identify_noise` to recover the
    dominant power law at each τ.
  - A spectral / Welch-PSD example once the `Sy` / `Sx` / `L`
    metrics ship (see TODO above).
- [ ] **`taus` enum API** (`AllTaus`, `Octave`, `HalfOctave`,
  `QuarterOctave`, `Decade`, `HalfDecade`) à la `AllanDeviations.jl` —
  cleaner than passing arrays. Backwards-compatible if added as an
  alternative to the current array form.
- [ ] **Compat upper bounds** in the root `Project.toml`. The merged
  manifest already pins `Distributions = "0.25.125"` and lists `compat`
  for AbstractFFTs, DocStringExtensions, RecipesBase, StaticArrays,
  julia. Tighten upper bounds for `Distributions` and `StaticArrays`
  once the dep matrix has been exercised on the General registry.
- [ ] **Remove `ldev` alias** — now marked `@deprecate`; delete after
  v0.3.0 is tagged.

---

## 🟢 Housekeeping

- [ ] **Agent-context briefs** (`CLAUDE.md`, `AGENTS.md`) are tracked
  in the repo. When the package matures and conventions stabilise,
  decide whether to keep them as-is or to land a sanitised version
  (without the authorship rules) under `docs/contributing/`.

---

## Docs follow-ups

- Fill tutorial narrative — start with `01_phase_data.md` and
  `02_compute_adev.md`.
- Tighten `warnonly = []` in `docs/make.jl` once all public API has
  docstrings.
- Refine `docs/src/refs.bib` with DOIs and page numbers from the PDFs
  in `legdocs/papers/`.
- Convert remaining kernel docstrings (hdev, mhdev, totdev, etc.) to
  use `$(SIGNATURES)` + `jldoctest` blocks following the adev/mdev
  pattern.
- **Three-cornered-hat theory page** beyond the tutorial. Material
  exists at `legdocs/papers/deviations/three_cornered_hat_*` (web
  articles + Riley papers).
- **Preprocessing / uncertainty theory page.** Gaps + outliers material
  lives at `legdocs/papers/preprocessing/`.
- Update the validation page `docs/src/validation/methodology.md` to
  cross-link the perf benches under `benchmarks/bench/` once long-record
  runtime numbers are stable.

---

## ✅ Recently shipped

See [CHANGELOG.md](CHANGELOG.md) for the annotated `## [0.3.0]` and
`## [0.2.0]` blocks.

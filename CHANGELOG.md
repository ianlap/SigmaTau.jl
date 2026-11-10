# Changelog

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

Clock models, the Kalman filter and timescales.

### Added

- `FrequencyData(p)` and `PhaseData(y)` convert a record to the other kind
  (first difference, or running sum from 0), keeping τ₀, `source` and
  `start`; `PhaseData` refuses a frequency record with gaps. In Python,
  `PhaseData.to_frequency()` and `FrequencyData.to_phase()`.
- `noise_id_nnls` and `powerlaw_fit` take the averaging factors as any
  integer vector, a range such as `1:64` included.
- Python: `telecom_mask`, `mask_check` (a `MaskCheck`; Julia's `pass` is
  `passed`), `TELECOM_MASKS` and `summarize`.
- `integrated_jitter`: integrated phase jitter (rms radians and seconds) and
  residual FM (rms hertz) over an offset band, from a `ssb_phase_noise`,
  `psd_x` or `psd_y` result or from spectrum points such as an analyzer's
  ℒ(f) (NIST SP1065 §5.2.19; IEEE Std 1139-1999 eq. A.16). Also in
  Python.
- `adev` and `hdev` take `overlapping=false` for the non-overlapped
  estimators (NIST SP1065 §5.2.2 and §5.2.8; Stable32's Normal Allan and
  Hadamard, allantools' `adev` and `hdev`), with the Greenhall–Riley EDF
  for stride τ. Also in Python.
- Explorer: Espresso dark theme; the header button follows the system
  setting or picks light or dark, and plots, the heatmap and exports follow
  the theme.
- Explorer import dialog: a typed scale factor (a number or an expression such
  as `1/(2*pi*518.295836590863e12)`) for readings that are not in the offered
  units, kept when the dialog is reopened; GNSS clock tables, whitespace
  `.csv` files, one-column files and comment-line names load as downloaded.
- `read_data` / `--input-format` read BIPM `UTC - UTC(k)` files (`:bipm_utc`)
  and RINEX 3 clock files (`:rinex_clk`, the IGS and CODE clock products: one
  satellite or station clock per call, picked with `column` / `--column`; old
  files with `D` exponents too). The explorer opens both as downloaded, unpacks
  `.gz` files, and lists a RINEX file's clocks to pick from, each on its own
  sample interval.
- Docs: "Coming from Stable32" and "Coming from allantools" guides, mapping
  their statistics, options and defaults onto SigmaTau.
- ITU-T wander masks for MTIE and TDEV: `TELECOM_MASKS` (G.811 PRC,
  G.813 SEC and G.8262 EEC options 1 and 2, and the G.823 network limits at
  PRC, SSU, SEC and PDH interfaces), `telecom_mask` for the limit at any τ
  and `mask_check` for pass or fail with the margin. `sigmatau analyze
  --mask NAME|all` adds them to the JSON, and the explorer's Stability page
  draws the chosen mask over MTIE and TDEV with pass or fail.
- Prebuilt downloads on each release: the `sigmatau` command line as a
  `.tar.gz` for Linux x86_64 and macOS (Apple silicon and Intel), next to
  the explorer's `.AppImage`, `.deb` and `.dmg`. The desktop release
  workflow now also runs on TagBot's `v…` tags and adds the files to that
  release; a manual run lists its artifacts in the run summary.
- `PhaseData` / `FrequencyData` carry the UTC time of their first sample
  (`start`, read with `start_time`), set by the readers from ISO 8601 and
  MJD timestamps, kept by `detrend`, `fillgaps`, `remove_outliers` and
  `sigmatau --range`, and written by `save` as a `# start:` line. The CLI's
  record JSON and the explorer's record check show it.
- `correct_bias` also takes the convention of a reference, on the total
  estimators, `theo1`, `deviation` and `stability`: `:published` keeps to
  the published bias tables (no factor on PM for TOTDEV and HTOTDEV, the
  τ ≥ 16τ₀ factor of MTOTDEV and HTOTDEV at every τ), `:stable32` does what
  Stable32 does for each estimator (published factors for TOTDEV and
  HTOTDEV, raw MTOTDEV, TTOTDEV and Thêo1), and `false` remains the raw
  statistic that allantools reports. `bias_correction` takes
  `published=true` for the same factors; the Python wrappers take
  `"published"` and `"stable32"`. MHTOTDEV has no published factor.
- Python: `noise_id_nnls`, `StreamingStability`, `start_time`, and a `start`
  field on `PhaseData` / `FrequencyData`.
- `noise` keyword on every estimator that reports intervals (`adev` …
  `theoh`), on `deviation` and on `stability`: `:auto` (identified at each
  τ, the default), one noise type for every τ, a vector with one per
  averaging factor, or a `PowerLawFit`, whose strongest term at each τ is
  taken. The stated type sets the EDF and, in the total family, the bias
  factor. `sigmatau analyze` and `report` take `--noise auto|fit|TYPE` and
  the JSON gains `analysis.noise`; the explorer has a Noise type menu under
  Options; the Python functions take the same `noise` argument.
- Cross-checked against AllanDeviations.jl 0.3.0 in the Cross-check CI job
  (every shared estimator within 1.5e-9 relative), and timed beside
  allantools in the benchmark.
- Explorer: end-to-end browser smoke test, run in the Web core CI job; no
  Plotly error when leaving the Data page.

- Interactive noise fit: `noise_fit` starts a power-law model on an ADEV,
  MDEV, HDEV or MHDEV result (or a total form), `add_component!` fits one
  noise type's level on a chosen τ range with the slope fixed and the
  confidence intervals as weights, `undo!` takes the last one back, and
  `noise_residual` shows what is left. Every change refits all the levels
  against each other, so the model does not depend on the order the
  components were added in. The model feeds `powerlaw_sigma`,
  `powerlaw_psd` and `clock_model`.
- Single-clock Kalman filter: `kalman_filter` runs a `ClockModel` over a
  phase record (gaps predicted across, Joseph-form update, least-squares
  start), `kalman_smooth` is the Rauch–Tung–Striebel smoother,
  `holdover_error` gives the 1σ time error after free-running for τ, and
  `innovation_check` tests the normalized innovations for unit variance and
  whiteness. `tune_q` retunes the white and random-walk FM levels for the
  smallest phase prediction error at chosen horizons (each level is resolved
  only to the grid step, a factor of 10^(2·decades/(points−1)), 3.2 by
  default), `prediction_error`
  reports that error, and `simulate_clock` draws a record and its true
  states from a model.
- Measured holdover test: `holdover(data; train, tau, drift, anchor, limit)`
  learns the frequency (and drift) from a training stretch, free-runs from
  the last reading and compares with the record from every start (Delporte et
  al., PTTI 2000). It gives the RMS, a quantile, the worst case and the mean of
  the time error at each holdover time, and the time at which the quantile
  reaches a limit.
- Predicted holdover: `holdover(model; train, tau, tau0, …)` predicts the
  error the measured test finds from a `PowerLawFit`, a `ClockModel` or a
  `Dict` of h_α, with `budget`, the variance each noise type adds
  (`HoldoverPrediction`; Python too). Exact for the sampled spectrum,
  flicker noises included; matches the closed forms of Delporte et al.
  (PTTI 2000) within 0.4 %.
- `sigmatau holdover FILE` runs the holdover test from the command line:
  `--train`, `--tau`, `--fit-drift`, `--no-anchor`, `--limit`,
  `--confidence`; CSV with one row per holdover time, or JSON, where `--fit`
  adds the error predicted from the power-law fit with its budget by noise
  type.
- Clock ensembles: `ClockComparisons` holds several clocks read against one
  reference (`cc[i]` one clock, `cc[i, j]` a difference), and
  `simulate_ensemble` draws one from a `ClockModel` per clock, with the true
  states.
- Timescales: `timescale(cc; algorithm=:ewfa | :at1)` forms a scale from
  the clocks' predicted phases, with equal weights (EWFA) or AT1's learned
  inverse-variance weights; missing readings are skipped and `timescale_phase`
  gives the scale against the reference. AT1 follows Weiss & Weissert 1991:
  weights from the prediction-error variance with the bias term `K` (eqs
  3 to 7), each clock's frequency constant from its own τ_min (eq. 10), the
  variance memory `N_τ` of 20 days and the starting variance τ²σ_y²(τ). The
  exact bias correction (`bias=:exact`), a weight cap (`wmax`) and an outlier
  test (`outlier`) are options.
- AT2: `timescale(cc; algorithm=:at2)` is AT1 with the paper's Kalman
  frequency update (eqs 13 to 18), also in the CLI (`--algorithm at2`) and the
  explorer. `steps=:best` or `:real` (`--steps best|real` in the CLI) adds the
  paper's frequency-step detection (19 to 22): a clock whose frequency jumps is
  left out for its τ_min and the scale rerun from the step, with `best`
  replacing the scale from the step on and `real` keeping the values already
  issued, as the paper does. The threshold of eq 20 is rewritten as a variance
  in consistent units (the printed form mixes time and frequency), a
  documented departure, and it defaults to 5 standard deviations rather than
  the paper's 4, which fired falsely about once per 3×10⁴ clock-samples on
  noisy clocks (`step_threshold=4` restores it). Off by default.
- Greenhall Kalman timescale: `timescale(cc; algorithm=:greenhall, models, pivot)`
  runs one Kalman filter over every clock with the covariance reduction
  after each step (the reduced Kalman timescale at `R = 0`, IEMx above), so
  the scale is quiet at short τ as well as long; without white PM its
  weights are the optimal 1/q1 ones. The scale is read at clock `pivot`.
  `algorithm=:jones_tryon` is the same filter without the reduction, the
  natural Jones–Tryon scale (JT83): a clock's reading less its estimated
  phase. In
  simulation `implicit_mean(r, truth)` gives the scale's own phase, the mean
  of the clocks' true phases, free of the pivot clock's reading noise.
- Stein's Kalman timescale: `timescale(cc; algorithm=:stein, models, pivot, wmax)`
  differences every clock against the pivot, runs a Kalman filter on each
  difference and gets the pivot's own shocks from the constraints that the
  weighted shocks of all clocks sum to zero, with separate weights for the
  phase, frequency and aging noise, each inversely proportional to the
  matching noise variance (Stein 2003; the KAS-2 algorithm, with `wmax` as its
  weight cap). On white FM clocks it matches `:greenhall` to a few percent.
  Also in the command line (`--algorithm stein`) and Python.
- `clock_models(cc)` fits a `ClockModel` per clock from the comparisons alone
  (Groslambert variances, then a power-law fit with white PM, flicker PM,
  white FM and random-walk FM, `alphas = [2, 1, 0, -2]`), so `:greenhall` and
  `:jones_tryon` need no hand-built `models`. `clock_model` folds flicker PM
  into the reading noise `R` as the PM the filter does not follow (Van
  Dierendonck, McGraw & Brown 1984); left as white FM, flicker PM inflated
  `q1`. `powerlaw_fit` accepts a `StabilityResult`.
- EWFA follows section 3.1 of Kettering's 2023 thesis: equal weights, and
  the scale's frequency is the mean of the clocks' frequency steps, so
  without gaps it is the plain mean of the phases. The `tau_y` option (an
  exponential average of each clock's frequency, which the thesis does not
  have) and `--tau-y` are gone; a clock without a reading still keeps its
  prediction rather than stepping the scale.
- Explorer: the Ensemble page offers the Greenhall and Jones–Tryon Kalman
  timescales, with a pivot clock; Jones–Tryon shows no weights (#166).
- Explorer: AT2's frequency-step detection (best or real) on the Ensemble page (#174).
- Explorer: the Ensemble page offers Stein's Kalman timescale (KAS-2), with the
  pivot clock and the largest weight (#200).
- Explorer: a × on each clock in the Clocks list removes it; the Ensemble page
  forms the scale from the clocks that are left (#182); the Ensemble page
  shows the core's error instead of waiting forever (#185).
- Command line: `sigmatau analyze --component NOISE:A:B` builds the `--fit`
  model by hand, and the JSON `powerlaw` gains `mode`, `components`,
  `residual` and `clock_model`. `sigmatau timescale FILE --value-cols 2,3,4`
  forms a timescale (`--algorithm`, `--wmax`; `--synth N`
  simulates four clocks).
- Explorer: **Edit** in the noise model panel builds the model by hand, with
  the MHDEV it leaves drawn dashed, and is saved with the session. The new
  **Ensemble** page forms a timescale from the clocks of one loaded file, or
  simulated ones, and plots each clock's weight and the ADEV of the scale
  and of each clock against it.
- Python: `ClockModel`, `clock_model`, `noise_fit`, `kalman_filter`,
  `kalman_smooth`, `holdover_error`, `innovation_check`, `tune_q`,
  `simulate_clock`, `ClockComparisons`, `simulate_ensemble`, `timescale`,
  `timescale_phase`, `clock_models`, `implicit_mean`, `holdover` and
  `HoldoverResult` wrap the Julia calls of the same names.
- Tutorial "From one clock to an ensemble": noise fit, Kalman filter,
  holdover and an AT1 timescale on simulated clocks.
- `groslambert(...; interval=:klts)`: KLTS confidence intervals (Lantz et al.
  2019) for three clocks, valid for a clock quieter than the others; the χ²
  interval is no longer given below 1 EDF.
- Flicker FM and flicker PM states in the Kalman clock model:
  `clock_model(fit; flicker=true)` and `clock_models(cc; flicker=true)` give
  a `ClockModel` with `flicker_fm`, `fm_band`, `flicker_pm`, `pm_band`, and
  `kalman_filter`, `kalman_smooth`, `holdover_error`, `prediction_error`,
  `tune_q`, `simulate_clock`, `simulate_ensemble` and the `:greenhall` and
  `:jones_tryon` scales run on it. Flicker FM is the Markov frequency states
  of Davis, Greenhall & Stacey (Metrologia 42, 2005); flicker PM states are
  the package's own, held as measurement noise with memory. States sit only
  across the τ where the flicker holds a quarter of the model's Allan
  variance. Off by default; `:stein` refuses such models.

### Changed

- `hdev`, `mhdev`, `htdev` and `htotdev` identify the noise type on the
  record less its least-squares frequency drift, as `mhtotdev` already did,
  and take `remove_drift` (default `true`); `deviation`, `stability` and
  the Python wrappers pass it on. A drift used to turn every long-τ label
  into random-walk FM, which on white FM left `htotdev` 13 to 15 % high and
  the intervals of all four too wide. The values of `hdev`, `mhdev` and
  `htdev` do not change; intervals and the `htotdev` bias factor follow the
  label. Without drift, flicker FM near the longest τ is read as white FM
  more often. `remove_drift=false` identifies the noise on the record as
  given, which is Stable32's convention.
- `confidence_intervals(devs, edfs, confidence)` no longer takes the noise
  types and record length, which only the removed Gaussian stand-in used
  (see Fixed). The five-argument form still works in Julia, with a
  deprecation warning; the Python wrapper takes the three arguments.
- The readers stream files line by line and keep only the columns they
  return: a 2-million-line two-column file reads in a third of the time and
  a quarter of the memory, a named-column table in a quarter of the time
  and a twentieth of the memory.
- `fillgaps=true` places samples on the grid of the regular time steps
  instead of the smallest step, so jittered timestamps no longer scatter
  the record.
- `identify_noise(...; dmax=3)`, and with it the Hadamard family (`hdev`,
  `mhdev`, `htdev`, `htotdev`, `mhtotdev`), now labels flicker-walk FM
  (`:FWFM`, α = −3) and random-run FM (`:RRFM`, α = −4) instead of
  reporting both as `:RWFM`. HDEV, MHDEV and HTOTDEV take their EDF, and
  HTOTDEV its bias factor, at the identified α: read as random-walk FM,
  random-run FM had an HDEV EDF 26 % too high and an HTOTDEV 6 % too low.
  MHTOTDEV keeps random-walk FM's model for both labels. The explorer, the
  CLI report and the Python docs name the two labels. At long τ some
  random-walk FM records are now labelled flicker-walk FM (6 % with 64
  averaged points left, 11 % with 32).
- `identify_noise` and every estimator's noise type: with fewer than 30
  averaged points the B1 ratio is compared with B1(N, μ) at the number of
  averaged points (`b1_length=:averaged`, SP1065 §5.7), so random-walk FM
  is reported at long τ; below 6 averaged points the label is the one
  measured at 6. The labels, and with them the bias factor of the total
  estimators (TOTDEV among them), change at the longest τ, where this
  departs from Stable32; `b1_length=:record` keeps the old labels.
- MHTOTDEV bias and EDF models replaced by exact ones. The bias factor is
  the exact expectation of the statistic, a constant per noise type (0.900,
  0.822, 0.792, 0.783, 0.787, 0.800, 0.820 for white PM to random-run FM);
  the τ/T slope is gone. The EDF model is `b·T/τ − c + d·τ/T + e·(τ/T)³`,
  fitted to the exact EDF through τ = T/4, where the old two-term model
  gave about half the true value. Flicker-walk and random-run FM have their
  own rows. Corrected deviations move by at most 0.7 % (2.5 % on random-run
  FM); intervals in the last octave narrow. `tools/mhtot_models.jl`
  replaces `tools/mc_mhtotdev.jl`.
- TOTDEV and HTOTDEV on white and flicker PM now divide out their bias and
  take the EDF of the total statistic itself, both from the exact moments of
  each estimator (`tools/total_pm_models.jl`, docs "Total estimators on PM
  noise"). HTOTDEV read 25 to 38 % high in variance on white PM and 4 to
  28 % on flicker PM; TOTDEV reads high in proportion to τ/T (33 % on white
  PM at τ = T/4). Stable32 applies no factor on PM, so the defaults differ
  from it on rows labelled PM; `correct_bias=false` matches it.
  **Breaking:** `bias_correction` now takes `(noises, var_type, m_values, N)`
  in place of `(noises, var_type, taus, T)`, in Julia and Python.
- MTOTDEV, TTOTDEV, HTOTDEV and MHTOTDEV divide by the exact bias factor at
  each averaging factor from 2 to 15, where the published constants (which
  hold from τ = 16τ₀) were up to 6 % off in variance; the flicker rows are
  the mean of two flicker models. MHTOTDEV's white PM EDF below m = 16,
  which the model overstated by up to 54 %, takes the exact ratio. HTOTDEV
  on FM noise at τ = 2τ₀ to 15τ₀ now differs from Stable32, which keeps the
  published factor, by up to 3 % in σ.
- `sigmatau analyze` and `report` take `--bias exact|published|stable32|none`,
  the library's `correct_bias` conventions; `--no-bias-correction` remains
  as `--bias none`. **Breaking:** the JSON field `analysis.bias_correction`
  is now that word instead of `true`/`false`. The explorer's bias checkbox
  is a menu of the same four; saved sessions still open, as Exact (was on)
  or None (was off).
- The test suite runs the Aqua.jl package checks.
- Noise identification: from m = 16 a white or flicker PM label stands only
  if R(n) = MVAR/AVAR agrees (nearer its flicker-PM value than the 1/2 of
  white FM); otherwise the noise is reported as white FM. With few averaged
  points the lag-1 and B1 tests mistook FM for PM, which gave ADEV and HDEV
  intervals several times too narrow. Labels on the Stable32 cross-check
  record now agree at every τ.
- Default calls are faster: the drift fit inside noise identification is a
  closed-form fit that agrees with the general fit to rounding, no longer
  goes through filtered generator sums, and is shared between m = 1 and the
  R(n) check. The fast estimators' default call is about 15× faster than on
  0.6.0 (about 0.4 ms on 25 000 points).
- `theo1(...; correct_bias=:stable32)` now reproduces Stable32's Thêo1 run
  with automatic noise identification, which is ThêoBR over the ladder
  SP1065 prints (to N/6), on every `⌈N/3000⌉`-th phase for records above
  Stable32's 3000-point limit; checked on all 181 rows of its result file
  for PHASE.DAT and on its two bias-removed rows for the 8192-point record.
  It returned raw Thêo1, which Stable32 gives only with Alpha set to 0.
  `correct_bias=false` still gives raw Thêo1.
- `analyze --spectrum` on a record with gaps runs the estimators and reports
  `"spectrum": {"error": …}` instead of stopping; `sigmatau spectrum` still
  stops. The explorer analyses such a record and its Spectrum page says why
  it has no spectrum.

### Removed

- **Breaking:** the names deprecated in 0.6: `Sy`, `Sx`, `L`, `nch`,
  `save_result`, `save_suite`, `load_result`, `load_suite`, the five-argument
  `confidence_intervals`, `SpectralResult.window`, `stability(; taus=)` and
  `identify_noise(; detrend=)`; their replacements are listed under 0.6.0,
  Deprecated.

### Fixed

- Explorer: a core that fails to start shows its error message instead of
  the last line of Julia's output (#179).
- `mhtotdev` at `m = 1` uses the exact bias factor 124/225 and MHDEV's EDF.
  The fitted bias and EDF models do not hold at τ₀, and the 68.3 % interval
  never covered the true value there.
- `totdev` at `m = 1` takes ADEV's EDF for every noise type: at τ₀ TOTVAR is
  AVAR term for term, and SP1065 Table 7 was applied to the FM noises.
- `mtotdev` and `ttotdev` take their EDF from Howe and Peppler (2001),
  Table 2, from τ = 16τ₀, in place of SP1065 Table 8, which runs from 16 %
  high to 41 % low against simulation; Table 2 stays within 8 % for
  T/τ ≥ 16.
- ThêoBR (`theo1` by default, and the Thêo segment of `theoh`) averages its
  AVAR/Thêo1 ratio over `i = 0…⌊N/30 − 3⌋`, as Howe (2006), McGee and Howe
  (2007) and Taylor and Howe (2010) define it. SP1065 eq. 33 prints
  `⌊N/6 − 3⌋`, which runs AVAR out to half the record: the ratio scattered
  by 25 to 30 % and biased the corrected variance by −9 to +11 %, against
  8 to 11 % and −5 to +3 % now. Output moves by a few percent. Below
  N = 90, where the papers' formula has no rung, the first rung is used
  alone.
- Flicker FM keeps its label at medium and long τ. On phase decimated by m
  the lag-1 estimate of flicker FM settles at α = −1.45 and of random-walk
  FM at −2.4, so rounding to the nearest integer sent pure flicker FM to
  random-walk FM in a large share of records: it kept its own label in
  only 48 to 72 % for m = 8 to 64 (N = 4096). For m ≥ 2 the limits between
  white FM, flicker FM and random-walk FM now sit halfway between those
  readings, at −0.72 and −1.92; flicker FM then keeps its label in 87 to
  100 % of those records and random-walk FM in 98 % or more. Stable32
  rounds, so its labels differ there.
- Flicker PM is identified as flicker PM at long τ. The lag-1 test runs on
  phase decimated by m, which whitens flicker PM, so from m ≈ 16 a pure
  flicker-PM record was labelled white PM (in 64 % of records at m = 16 and
  93 % at m = 64, N = 4096) and `adev` reported about six times the true
  EDF there. A white-PM result at m ≥ 4 is now confirmed with the R(n)
  ratio of SP1065 §5.5.2. Labels of the other four noise types, and of
  every shipped Stable32 and allantools fixture, are unchanged.
- `timescale` no longer turns the whole scale into `NaN` when a clock has no
  reading at the first sample: a clock joins the scale, at zero weight, at
  its first reading.
- `timescale_phase` reads the scale through any clock with a reading at
  each sample, instead of returning `NaN` wherever the first clock has a
  gap.
- The default bias correction of `mtotdev` and `ttotdev` went the wrong way.
  Mod Totvar reads low against MVAR, by 0.5 % for white PM up to 32 % for
  random-walk FM (Howe and Peppler 2001, Table 2; Howe and Vernotte 1999),
  and the correction divided it by the factors of SP1065 Table 11, 1.06 to
  1.31, which lowered σ by a further 3 to 13 %. It now divides by the Howe
  and Peppler ratios (`bias_correction(…, :mtot, …)` returns 0.995, 0.851,
  0.771, 0.717, 0.679), so the default tracks MDEV, and σ is 3 % (white PM)
  to 39 % (random-walk FM) higher than before. At `m = 1`, where the
  estimator is exactly MDEV/√2, the correction is that factor.
  `correct_bias=false` is unchanged and still matches Stable32 and
  allantools.
- Confidence intervals where the EDF is below 1 are χ² intervals like every
  other, and so the widest on the plot. They used to be replaced by a
  Gaussian interval on the whole record length, which gave the least
  certain point the tightest bounds: ±2.4 % at the last τ of `mhtotdev` on a
  1024-sample record, where its neighbour had −29 % / +290 %. The same
  applied to `theo1` and `theoh` at τ beyond about 0.6 of the record. A
  point whose EDF model gives no positive value now has `NaN` bounds.
- An averaging factor below 1 is an `ArgumentError` in every estimator. A
  negative one used to return numbers read from outside the record, and 0
  gave `Inf`, `NaN` or a `DivideError`.
- `totdev`, `mtotdev`, `ttotdev`, `htotdev`, `mhtotdev`, `pdev`, `theo1`,
  `theoh`, `dadev` and `dhdev` on a record with `NaN` gaps throw an
  `ArgumentError` that points to `fillgaps`, instead of returning `NaN` at
  every τ with noise labels attached.
- `find_outliers` flags nothing when the MAD is zero, as on a coarsely
  quantized record where more than half the steps are equal. It used to
  flag every sample off the median, and `remove_outliers`, `summarize` and
  `sigmatau check` followed.
- `deviation(name, data)` on a `FrequencyData` record uses the same τ grid
  as `adev(data)` and `stability(data)`; it stopped one step short.
- The readers take a one-column file as its values (a `# data interval
  1.0 s` note gives τ₀), never use the value column as its own time axis,
  and split a blank-separated `.csv` on blanks.
- Timestamps rounded to a coarse resolution (MJD to six decimals, 0.0864 s)
  give the true τ₀ and grid: the step is the mean of the regular steps,
  rounded to the shortest value every timestamp fits, where the median gave
  1.0368 s on 1 s data and joined the gaps.
- A table column logged more sparsely than the table (5 min values in a
  30 s file) keeps its own step instead of failing to fill; a lone column
  beside the time column is picked as the clock.
- `fillgaps` anchors each fill to the filtered samples beside its gap
  instead of refiltering the whole record per gap: a month of 1 s data with
  13 500 gaps fills in 5 s instead of over 7 minutes.
- The Linux AppImage of the explorer could not run its analysis core: the
  AppImage launcher's `LD_LIBRARY_PATH` made the core load the AppImage's
  copy of libjulia, which then looked for `libgcc_s.so.1` beside it. The app
  now starts the core without that variable.
- A reading missing from a timestamped file is a `NaN` sample in its slot,
  which the gap-tolerant estimators skip. Without `fillgaps` the samples
  around it used to be joined end to end, shortening the record and moving
  every τ after the gap. `sigmatau --freq` asks for `--fillgaps` on a
  frequency record with missing samples instead of integrating `NaN` into
  the rest of the phase.
- `read_data` reads a time column of MJD values in days whatever its name
  (MJD 30000 to 90000, steps under 0.05 d, more than four decimals), a
  `ratio` column as frequency, and unnamed columns of median size under
  1e-11 as frequency; ROCIT files now read as downloaded.
- `push!` on a `StreamingStability` throws on a sample that is not finite
  instead of turning every later estimate into `NaN`.
- `ncornered_hat` no longer gives a χ² interval below 1 EDF, as
  `groslambert` now does.
- `totdev` with an explicit averaging factor above N − 1, where the
  reflected extension ends, read outside the record and returned an
  arbitrary value with an interval; it now gives a `NaN` row with `neff` 0,
  as the other estimators do past their limits.
- Readers: exact timestamps keep their own sample interval on short records
  (86400 s no longer reads as 86000 s below 55 rows); a unit in a table
  column's name sets phase or frequency (`ppb`, `ppm`, `ns/s` are frequency;
  `ns`, `ms` phase) and a time unit in the time column's name is applied
  (`time (ms)`); `read_data` keeps `delim` for plain columns, and a
  comma-separated file without column names reads whatever its extension.
- Explorer: page buttons act once per click again. After a page was
  revisited, "Pin run" pinned twice and the noise model's Edit could not be
  opened.

## [0.6.0] - 2026-11-08

### Fixed

- `fillgaps` fills a gap at the start or end of a record, joined to the
  data on its one side, instead of throwing `BoundsError`; a layout it
  cannot fill now throws its "gap topology" error instead of looping
  forever.
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

- `groslambert`: the Groslambert covariance for three or more clocks, each
  clock's Allan deviation from the cross Allan variance of the two
  difference records that share it, with χ² intervals from the covariance
  model of Vernotte & Lantz (2019): each record's Allan-variance EDF for its
  noise type, propagated in closed form. `ncornered_hat` now fills `edf` and
  `ci` the same way when its inputs carry EDF.
- Live measurement: `LiveStability` keeps ADEV, MDEV, HDEV and MHDEV up to
  date as readings arrive (phase, fractional frequency or counter Hz), equal
  to the batch estimators at every reading. Sources: `InstrumentSource` for
  any line-based stream (an LXI SCPI socket, a serial port, a log file) and
  `SimulatedSource`. `sigmatau live SOURCE` streams the running σ(τ) as JSON
  lines, and the explorer's Live page plots it and hands the readings to the
  full analysis.
- Analysis range: `--range A:B` (seconds from the first sample) and
  `--tau-min` / `--tau-max` in the CLI. In the explorer, **Select range** on
  the Data page drags out the part of a record to analyse, the options take
  τ limits, and **Pin run** on the Stability page keeps a run's curves to
  compare with the next one.
- Power-law noise model: `powerlaw_fit` fits the h_α of
  S_y(f) = Σ h_α f^α to a record's σ(τ) (ADEV, MDEV, HDEV or MHDEV, with
  each component's exact response for a sampled record, from m = 1) or to a
  spectrum, including a datasheet's ℒ(f) points; `powerlaw_sigma` and
  `powerlaw_psd` convert a model to σ(τ), S_y(f), S_x(f) or ℒ(f).
  `sigmatau analyze --fit` adds the model to the JSON, and the explorer's
  **Fit** toggle draws it as dotted lines with a panel of the h_α and the
  model's ADEV at round τ.
- Two-channel estimators: `cross_psd_y` / `cross_psd_x` (Welch
  cross-spectrum with the floor `√(S_aa S_bb / nseg)`) and `cross_adev`
  (cross Allan variance with its 1σ uncertainty from the same covariance
  model and each channel's noise-type EDF),
  for a device measured against two references or one signal split into two
  channels. `sigmatau cross FILE --with FILE2` on the command line.
- Reports: `sigmatau report FILE` writes one self-contained HTML file with
  the record summary, the settings, σ(τ) plots, per-τ tables and the
  command that reproduces it; a browser prints it to PDF. The explorer's
  **File** menu (formerly Export) writes the same kind of report for every
  loaded clock, with the SHA-256 of each input file.
- Explorer sessions: **File → Save session** stores the clocks, their files,
  import settings and analysis settings in one JSON file; **Open session**
  (or dropping the file on the window) restores it and reruns.
- `noise_id_nnls(data[, m]; deviation=:mhdev, alphas, weights=:edf)`: a
  weighted nonnegative least-squares fit of the measured σ̂²(τ) profile to
  power-law components. One fit over the whole τ range gives the amplitude of
  every α in the basis and the dominant noise type at each τ; it sits beside
  the per-τ lag-1 ACF classifier `identify_noise`, which is unchanged. The
  NNLS problem is solved exactly by enumerating active sets. Default bases:
  α ≥ -2 for `:mdev`, down to α = -4 for `:mhdev`, `:hdev` and `:mhtotdev`,
  one PM column for the unmodified `:hdev`.
- `noise_gen` covers flicker-walk FM (α = −3) and random-run FM (α = −4), and
  takes noise-type names as keys (`Dict(:WHFM => 1e-12)`). Its spectrum,
  `h_α (sin(πfτ₀)/(πτ₀))^α`, is checked against theory for every α from +2
  to −4 in the test suite.
- `read_data` and `detect_format` read instrument files by layout: TimeLab
  `.tim` (τ₀ and scale factor from the header), Stable32 `.phd`/`.frd`
  (MJD time tags, `Type:` header), SCPI counter readings (Hz converted to
  fractional frequency against a nominal frequency, `9.91E+37` as a missing
  reading, time-interval readings as phase) and named-column logs (unit
  suffixes such as `(ns)` or `ppb`, ISO 8601 or MJD timestamps). The CLI
  detects them too (`--input-format`, `--column`, `--nominal`), the
  explorer's import dialog recognises them, and Python has the same names.
- Gap-tolerant `adev`, `mdev`, `tdev`, `hdev`, `mhdev`, `htdev`, `tierms`
  and `mtie`: `NaN` samples are skipped instead of turning σ into `NaN` at
  every τ. On phase a term counts when every sample it reads was measured;
  on frequency, when no missing sample lies inside its span. `neff` counts
  the terms used and the EDF and confidence interval follow from it. ADEV
  matches allantools' `gradev` value and term count to 2e-15. The CLI and
  the explorer pick this up without a flag.
- `drift_fit(data; model)`: Stable32's frequency drift fits (linear,
  quadratic, MIL-PRF-55310 log and diffusion) with standard errors, the
  drift rate at the end of the record per day and its uncertainty, and the
  fitted curve. Gaps are skipped. Coefficients and standard errors match
  scipy least squares (polynomials to 1e-13, log to 2e-8, diffusion to 1e-10).
  The CLI adds `sigmatau drift FILE --model MODEL`, and the explorer's Data
  plot gets a Frequency drift view with the fit and its residuals.
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
- Python bindings, the `sigmatau` package in `python/`. Every call runs in an
  embedded Julia through juliacall, so results match the Julia API exactly.
  The names and keywords are the 0.6 ones (`psd_y`, `ncornered_hat`,
  `save`/`load`, `m=`), plus `deviation(name, ...)` and `estimators()`.
  Tested by a pytest suite in its own CI workflow. A `python-v*` tag
  publishes it to PyPI (`pip install sigmatau`), pinned to the same
  SigmaTau.jl version from the General registry.
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
- **Breaking:** `mhtotdev` values change for `m ≥ 2`. Each subsegment is now
  `4m` phase samples (one full MHVAR span) instead of `3m + 1`, and its
  least-squares quadratic is removed before the reflection instead of a line.
  MHTOTDEV is now immune to linear frequency drift on its own and converges
  for power laws down to α = −4, like HDEV and MHDEV; `remove_drift` now
  only affects noise identification. It reads low by a nearly constant factor
  (bias B ≈ 0.90 for white PM, ≈ 0.79 for FM noise, corrected by default).
  Bias and EDF coefficients were re-measured by the full Monte Carlo sweep
  (`tools/mc_mhtotdev.jl` at `9d7a408`, seed 20260531, R = 3000, N ≤ 32768);
  provenance in `tools/artifacts/mhtotdev_mc_full.json`.
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

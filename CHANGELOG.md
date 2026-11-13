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
  whiteness. `tune_q` retunes the white and random-walk FM levels by
  maximizing the likelihood of the innovations over `q1` and `q2` (Tryon
  and Jones 1983, Eq. 5.1) by Nelder–Mead, with `burnin` and `maxiter`, and
  returns `(model, cost, converged, evaluations)`; `prediction_error`
  reports the phase prediction error at chosen horizons, and
  `simulate_clock` draws a record and its true states from a model.
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
- Explorer: a Holdover page shows the time error after a training stretch
  against holdover time (RMS, a chosen share of the trials, the worst trial
  and the noise model's prediction) and how long each clock stays within a
  limit.
- Holdover with a stretch that grows with the holdover time:
  `holdover(…; grow=k)`, `sigmatau holdover --grow K`, a "2τ" stretch on the
  explorer's Holdover page. At `grow = 2` (3 with drift) it equals √6·TDEV
  (√20·HTDEV) at τ₀ and runs a constant factor above on FM noise; the page
  and the CLI's JSON (`reference`) show that curve.
- Holdover with a Kalman filter as the predictor:
  `holdover(data, model::ClockModel; …)` carries the filter's estimate ahead
  from every sample past the training stretch; `HoldoverResult.claimed` is
  the filter's own standard deviation for that error, to set beside the
  measured RMS (Python too).
- `sigmatau holdover --kalman` and a "Predicted by" choice on the explorer's
  Holdover page: the holdover test with a Kalman filter as the predictor, on
  a clock model fitted to the first training stretch alone, with the error
  the filter itself claims (`claimed`) beside the one measured. The CSV
  gains a trailing `claimed` column, empty without `--kalman`, so its header
  changes for existing readers.
- `sigmatau holdover --kalman --format json` gives `fit_tau_max`, the
  longest averaging time of the fit behind the predictor's clock model, and
  the explorer dots the claimed curve past it, where the claim extrapolates
  the model.
- Clock steering: `steering_gains` sets the gains of a steering loop by pole
  placement, from a time constant or from poles, with or without an integral
  gain (Matsakis 2019, 2020); `steering_poles` gives the poles, time
  constants and stability of any gains; `steer(data, model, gains)` replays
  a phase record through the loop on one Kalman filter pass, and
  `steer(model, gains)` predicts the steady-state RMS of phase, frequency
  and steers and the Allan deviation of the steered clock. Also in Python.
- `lqg_gains`: steering gains that minimize a quadratic cost in phase,
  frequency and steer, the linear quadratic Gaussian design (Koppang &
  Leland 1999); it reproduces the three designs that paper prints. Also in
  Python.
- `sigmatau steer`: replays a record through a steering loop and writes the
  Allan deviation free-running, steered and as the clock model predicts.
  Gains from `--time-constant` (with `--integral`) or `--cost P:F:U`; the
  clock model is fitted to the first `--train` stretch alone.
- `steer(model, gains; readings=true)`: the predicted time error and Allan
  deviation of the readings of a steered clock, with the noise each reading
  carries, to set beside a replayed record. `sigmatau steer`'s `predicted`
  column now uses it; before, it was for the clock without that noise and
  sat below the steered column at short averaging times. Also in Python.
- Explorer: a Steering page. It replays a steering loop on the record, with
  gains from a time constant or a cost, and shows the Allan deviation
  free-running, steered and as the clock model predicts, or the steered
  phase with the frequency correction in force. Keys: 7 Steering, 8
  Ensemble, 9 Live.
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
  phase. In simulation `implicit_mean(r, truth)` gives the scale's own
  phase, the mean of the corrected clocks, free of the pivot clock's reading
  noise.
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
- `sigmatau timescale --flicker` fits the Greenhall and Jones–Tryon scales'
  clock models with flicker states and reports the number each clock carries
  (`"flicker": {"fm", "pm"}`); the explorer's Ensemble page has the same
  choice. Python's `ClockModel` mirrors `flicker_fm`, `fm_band`,
  `flicker_pm` and `pm_band`, and `clock_model` and `clock_models` take
  `flicker=True`.
- Python: the keywords and call forms of the Julia functions that the
  wrappers lacked. `dadev` and `dhdev` without `window`; `identify_noise`
  with `b1`, `b1_length` and a `TauMode` grid; `read_data` with `delim` and
  `keep_zeros`, and `detect_format` with `delim`; `powerlaw_fit` of a
  `StabilityResult`; `powerlaw_sigma` and `powerlaw_psd` of a `ClockModel`
  or a `NoiseFit`; `holdover` of a `NoiseFit`; `tune_q` with `burnin`;
  `timescale` with `step_threshold`, `lmax` and `step_inflate`;
  `clock_model` with `fm_band` and `pm_band`.
- `groslambert` takes phase records with missing samples (`NaN`) under
  `interval=:chi2`: a term that reads one is left out (SP1065 §10.4), `neff`
  counts the terms kept and the EDF follows them. `clock_models` therefore
  fits models from comparisons with missing readings, at the octaves of each
  clock's interval between readings, and `sigmatau timescale` and the
  explorer form the Kalman scales on files with gaps.
- Telecom masks of G.8272 (`:G8272_PRTC_A`, `:G8272_PRTC_B`), G.8272.1
  (`:G8272_1`), G.811.1 (`:G811_1`), G.812 (`:G812_I`, `:G812_II_III`) and
  G.8262.1 (`:G8262_1`) in `telecom_mask`, `mask_check`,
  `sigmatau analyze --mask`, Python and the explorer's mask menu, which is
  grouped into clocks and network limits. Each mask names the edition and
  the tables its limits were read from: `source` in the JSON of `--mask`,
  and on the explorer's card.
- `drift_fit` takes `method` for the linear model: `:frequency_fit` (the
  default, as before), `:phase_fit` (a least-squares parabola through the
  phase), `:bisection` (the three-point fit), `:second_difference` (the mean
  second difference of the phase) and `:four_point` (Greenhall 1997). These
  are the drift methods of NIST SP1065 §5.21 and §5.22 and of Stable32's
  Drift function, and the four-point estimator. `DriftFit` has a `method`
  field, and a standard error the method does not give is `NaN`.
  `detrend(data, fit)` removes the model of a fit. In Python,
  `drift_fit(..., method=)`, `detrend(data, fit)` and `DRIFT_METHODS`.
- `sigmatau drift --method`; its CSV row has a `method` column after
  `model`. Explorer: an Estimator list on the Frequency drift view.
- `sigmatau hat`: the stability of each of three or more clocks measured
  against one another, by the N-cornered hat on ADEV, MDEV, HDEV, MHDEV,
  PDEV, TDEV or HTDEV, or by the Groslambert covariance (`--covariance`),
  with a χ² or the KLTS interval (`--klts`). The clocks are columns of a
  file read against one reference, which counts as one more clock
  (`--no-reference` leaves it out), or three columns A−B, B−C and C−A
  (`--pairs`).
- `ncornered_hat` and `groslambert` take a `ClockComparisons`: clocks read
  against one reference, with the reference as one more clock unless
  `reference=false`.
- `sample_statistics`, `sample_histogram` and `acf`: the statistics of a
  record's samples (with the standard deviation of the fractional frequency
  that Stable32 shows), their histogram with what a normal curve over it
  needs (Scott's rule for the bins), and their autocorrelation at an
  averaging factor with the lag-1 estimate of α, as Stable32's Statistics,
  Histogram and Autocorrelation functions show them. In Python,
  `st.sample_statistics`, `st.sample_histogram` and `st.acf`.
- `sigmatau check` reports the statistics of the samples, and with
  `--histogram` (`--bins`) and `--acf` (`--m`, `--lags`) their histogram and
  autocorrelation, of phase or frequency (`--quantity`). Its CSV row has
  nine more columns.
- Explorer: Histogram and Autocorrelation views on the Data page, with the
  statistics beside the histogram.
- `bias_b1`, `bias_b2` and `bias_b3`: Barnes's bias functions for frequency
  averages taken with dead time between them (NBS TN 375; NIST TN 1318), for
  an exponent μ from −2 to 2 or a noise type, `bias_b3` for any number of
  samples. In Python, `st.bias_b1`, `st.bias_b2` and `st.bias_b3`.
- `adev` takes `dead_time` on a `FrequencyData` record whose samples average
  over `tau0` and start every `tau0 + dead_time`: the deviation and its
  interval are divided by √(B2·B3) for the stated `noise` (SP1065 §5.15).
  Also `dead_time=` on `st.adev` and `st.deviation` in Python, and
  `sigmatau analyze --dead-time`. The "Dead time" page of the documentation
  gives the functions, the limits of the correction and the agreement with
  the published tables.
- Explorer, Ensemble page: an N-cornered hat view beside the time scale. It
  separates the stability of each clock by the hat on ADEV, MDEV, HDEV,
  MHDEV, PDEV, TDEV or HTDEV, by the Groslambert covariance, or by the
  covariance with the KLTS interval, from clocks read against one reference
  or from three pairwise columns, and plots each clock with its interval
  over the pairs as measured.
- `spectrum_sigma`: the Allan, modified Allan, Hadamard or modified Hadamard
  deviation from a table of ℒ(f), S_φ(f), S_x(f) or S_y(f), integrated
  through the transfer function of the estimator over the band of the table.
  The result states the band and what the noise below and above it could
  add, which has no limit where the deviation depends on the bandwidth of
  the measurement (white and flicker PM in the Allan and Hadamard
  deviations). Also in Python.
- `sigmatau phasenoise TABLE --carrier HZ`: σ(τ) from a table of phase
  noise, as ADEV, MDEV, TDEV, HDEV or MHDEV, with `--density`, `--tau` and
  `--band`, in CSV or JSON.
- `quantity=:Sphi`, S_φ(f) in rad²/Hz, in `powerlaw_fit`, `powerlaw_psd` and
  `integrated_jitter`, and `"Sphi"` in Python.
- Python: `ClockComparisons.from_records`, one phase record per clock.
- Explorer: **Phase noise…** on the Stability page draws the σ(τ) that a
  table of ℒ(f) gives, beside the deviations of the record.
- `groslambert(...; weighted = true)` weighs each pair of the other clocks
  by the inverse of the product of its records' Allan variances (Python:
  `weighted=True`).
- `sigmatau timescale`: the JSON output of a Kalman scale has `models`, the
  model fitted to each clock, and a message on stderr names a clock whose
  fit has no white FM; the explorer's Ensemble page shows the fitted white
  FM.
- `powerlaw_sigma(h, τ)` takes one τ as a number and returns a number
  (Python: a float).

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
  random-walk FM records are now labeled flicker-walk FM (6 % with 64
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
  from it on rows labeled PM; `correct_bias=false` matches it.
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
- `:greenhall` with some noiseless-reading clocks and some with reading
  noise: the noiseless clocks are treated as one in the implicit-mean
  weights, so the reported weights follow clock stability (before, the clock
  with reading noise was pulled toward 1/k) and the scale changes slightly
  with them (up to about 15 % in ADEV at short τ in the mixed case).
  Frequencies are unchanged.
- `theo1` and `theoh` sum Thêo1 in cache-sized blocks across all threads; a
  default `theo1` call on 557 000 points takes 48 s on 4 threads. Values
  move by rounding only.
- `psd_x`, `cross_psd_x` and `ssb_phase_noise` estimate S_x as the Welch S_y
  over the response of the first difference, 4 sin²(πfτ₀)/τ₀². A Welch
  estimate on the phase itself read high at the lowest bins on FM noise
  (random-walk FM about 100 times at the second bin) and was swamped by a
  frequency offset. The value at f = 0 is now `NaN`, and `nperseg` counts
  frequency samples, at most one fewer than the phase samples. The Welch
  estimators remove each segment's taper-weighted mean instead of its plain
  mean, which read white PM `1 + nperseg/(12π²)` times high at the first
  bin; `psd_y` changes at its first two values only and equals scipy from
  the second bin.
- `fillgaps` follows Howe and Schlossberger (2022) and their script where it
  had departed (the reflected samples are inserted unfiltered and the
  filter, now their 16-sample form, only sets the line that joins them),
  sets the series' trend aside while filling, and takes the series it fills
  (phase or frequency, `form=:auto`) and the filter's averaging span from
  the record's own two-sample deviation. A filled record read 8 % low at
  short τ, a phase record with a frequency offset several times high, and a
  frequency record of white PM 55 times high at long τ after a 30 % fill;
  with one to six gaps ADEV is now within 5 % of the gap-free record's on
  every power-law noise. The docstring states what a fill does not recover:
  many gaps in a frequency record, a strong periodic term (remove it before
  filling), MDEV at the longest τ on white PM. `form=:phase` or `:frequency`
  forces the series; Python takes `form=`.
- Noise labels where an averaging factor leaves fewer than 30 averaged
  points. A result is PM only if it is PM at the longest τ with 30 points
  and R(n) at the τ itself agrees; otherwise it is the FM type of the B1
  ratio. `hdev`, `mhdev`, `htotdev` and `mhtotdev` with `remove_drift=true`
  no longer use the B1 ratio on the drift-removed record: an FM type there
  is the one at the longest τ with 30 points
  (`identify_noise(...; b1=false)`). Default `htotdev` on flicker FM and
  redder was 7–20 % low in variance at 16 points and fewer; `totdev`,
  `mtotdev` and `htotdev` on white PM were up to 33, 12 and 19 % high.
  Flicker PM with fewer than 6 averaged points is now reported as white FM.
- `--series N` keeps each stretch's lowest and highest sample and its first
  gap, instead of an evenly strided subset.
- `adev` and `hdev`: the overlapping kernels vectorize again on Julia 1.12
  and newer, about 8× faster than since `overlapping=false` was added;
  values move by rounding at most.
- Every estimator's `NaN`-gap check is one vectorized pass; a call on a
  25 000-point record saves about 9 µs.
- Documentation: every docstring checked against the code and its cited
  source and rewritten in one style; bibliography verified against publisher
  records, with DOIs, and 21 sources added.
- Documentation site, Tutorials 1 to 3 and the READMEs rewritten from the
  code and the cited sources, in six parts: tutorials, how-to guides,
  explorer, background, reference, bibliography. New: a manual of the
  explorer with screenshots of simulated records, a command-line reference
  whose samples the docs build produces by running the program, a guide to
  gaps, outliers and drift, and a page on related software. The cross-check
  group "NBS 1000-point" is now "1000-point test suite" (SP1065 §12.4), and
  the benchmark table of a measured record is removed with its two result
  files.
- Documentation: a background page, "The Hadamard time deviation", on where
  the 10/3 of `htdev` comes from, and √6·TDEV and √20·HTDEV as the errors of
  predicting the next τ-average with a line and with a parabola.
- Documentation: the clock model, Kalman filter, clock ensembles, holdover
  and steering pages and Tutorial 4 are rewritten, and the docstrings of the
  clock-model, Kalman, time-scale, comparison and holdover functions follow
  the layout of the rest, each citation with its equation or section. Every
  figure given in simulation was measured again. Corrected: what
  `:jones_tryon` loses against `:greenhall`, that it carries a drift state,
  what leaving flicker FM out of a clock model costs, and where a parabola
  predicts better than a line.
- Explorer: the help is the manual. **?** opens a menu with a search through
  the manual and the sections on the page in view, F1 opens the help of that
  page, and the help page shows the manual by topic, written from
  `docs/src/explorer/` by `docs/explorer_help.jl`; the documentation build
  checks that the two agree. Links to the documentation open in the browser,
  from the desktop app too. The top bar fits the window again with nine
  tabs.
- `bench/README.md` re-timed on a Fedora workstation, with a Thêo1 timing on
  the public allantools Cs 5071A record.
- The Stable32 reader (`.phd`, `.frd`) reads a value of zero as a gap, as
  Stable32 does (all zeros in frequency data; all but the first and last
  point in phase data), and reads Stable32's older header of bare `File:` …
  `Tau:` lines. `read_data(path; keep_zeros=true)`, CLI `--keep-zeros` and
  the explorer's *Keep zeros as readings* box keep them as readings for
  quantized records.
- `clock_model` keeps the flicker FM and flicker PM of a fit as Markov
  states by default (`flicker = true`). `clock_model(powerlaw_fit(p))` used
  to leave the fitted flicker FM out, and its filter then claimed 1.2 to 8
  times less error than it made on a clock with a flicker floor.
  `flicker = false` gives the former model, which the `:stein` time scale
  needs.
- `timescale` with `:stein` starts the difference of each clock from its
  filter's estimate over the clock's first tenth of the record, smoothed
  back to its first reading, and starts the scale at the weighted mean of
  its clocks. A clock that joins late no longer shifts the rate of the scale
  by many times what the noise of the clocks gives, the rate and drift of
  the scale no longer depend on the pivot at the start, and models with
  aging no longer run the scale away. Every `:stein` result changes.
- `timescale` with `:at1` or `:at2`: the starting prediction-error variances
  are an N-cornered hat over the clocks and the reference, the default
  `tau_min` is taken from the octaves up to a tenth of the record and beyond
  them only while the Allan deviation still falls, and the scale starts at
  the mean rate of its clocks under their starting weights. Scales of clocks
  that differ much in noise change from their first sample; `:ewfa` and the
  Kalman scales do not change.
- `clock_models` names the clock whose curve leaves the fit too few points,
  in place of the `powerlaw_fit` error.
- `sigmatau analyze --mask` works with CSV output: each check of MTIE or
  TDEV is one line on standard error, silenced by `--quiet`.
- `noise_gen` draws the components of a mixture in order of α, not in the
  order of the Dict, which differs between Julia versions: a seeded record
  now agrees across versions to rounding error, and a seeded mixture differs
  from earlier releases.
- `timescale`: `x`, `y` and `xhat` of a clock are `NaN` before its first
  reading in every scale, and a record whose first rows hold no reading
  forms every scale from the first row that has one.
- `sigmatau timescale`: `--fillgaps` no longer hands a clock readings it did
  not have; the filled samples are left out of the scale.
- Explorer: the import dialog proposes every column named for a reading
  (`maser1 phase (s)`, `UTC-UTC(PTB) (ns)`) as a clock, of the first one's
  kind and stated unit, so a file of several clocks loads in one step
  ("Load 6 clocks"); a file without column names still proposes one. A
  clock's automatic name gives way to the column of a newly opened file; a
  typed name stays.
- Explorer: refusals of the core that named command-line flags and library
  functions (`--algorithm greenhall`, `wmax`, `powerlaw_fit … alphas`, the
  steady state of a steer every n samples) are said in the page's words. The
  Ensemble ADEV plot has no legend.
- Python: an argument that counts (an averaging factor, a number of samples,
  a window, a segment length, a clock index, a column) must be an `int` or a
  NumPy integer; a float is a `TypeError`, where it was cut to its whole
  part.
- Explorer: on the Data page the plot keeps a usable height in a small
  window (the page scrolls); the hat plot's legend sits above the curves,
  not over them.
- `clock_models` weighs the pairs by their noise, so the fitted models
  change for four clocks or more; a quiet clock beside noisy ones is now
  resolved.
- `timescale` refuses `steps` with a Kalman scale, `models` with an
  averaging scale, and a `tau_min`, `memory` or `step_threshold` that is out
  of range.
- `clock_model(...; flicker = false)` raises an error for a fit whose only
  frequency noise is flicker FM, which would leave a model with a frozen
  frequency; fit with `alphas = [2, 0, -2]` or keep the flicker.
  `sigmatau analyze --fit` writes `clock_model: null` for such a fit and
  keeps the rest.
- `kalman_filter`, `prediction_error`, `tune_q`, `holdover` and `steer`
  refuse a `FrequencyData` record with a message to convert it with
  `PhaseData(y)`, where they gave a `MethodError` or asked for `tau0`. In
  Python they raise a `TypeError` naming `data.to_phase()`.
- `clock_model` of a `NoiseFit`, or of a `PowerLawFit` of a spectrum, places
  the flicker states within the τ the fit saw (its longest τ, or 1/(2f) at
  the spectrum's lowest frequency), not out to 10⁶ τ₀.
  `clock_model(fit; tau0)` refuses a τ₀ that differs from the fit's, where
  it used it silently.
- `stability(data; devs = "adev")` runs one estimator, where a string was
  split into letters. `deviation` takes the name as a string and `tau_mode`
  as a symbol. A single averaging factor (`adev(p, 4)`), `m_values` or an
  unknown keyword to `stability` gets a message that says what is meant.
- Where a deviation is `NaN` (a τ the record cannot reach), its EDF is `NaN`
  and its noise type empty, in the library, the CLI and Python; `pdev`, the
  total estimators and `theo1` reported an EDF from a formula outside its
  range there.
- Reader errors say what is wrong with the file: an empty file holds no
  samples, a gzip file must be unpacked first, and a file without a time
  column or `# tau0:` header has no sample interval; each names the reader
  that was called.
- CLI errors name command-line options and not library keywords or functions
  (`--tau0`, `--fillgaps`, `--column`); an output file that cannot be
  written is named with the reason; `timescale --synth` with one clock is a
  usage error (exit 1).
- `sigmatau dynamic --window` below 4 (5 with `--dev dhdev`) is a usage
  error (exit 1) that names the option; it was an analysis error (exit 3)
  saying the record was too short for adev or hdev.

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
- ThêoBR (`theo1` by default, and the Thêo segment of `theoh`) keeps the
  ladder of SP1065 eq. 33, IEEE Std 1139-2022 eq. A.37 and Stable32,
  `i = 0…⌊N/6 − 3⌋`. `correct_bias=:published` takes the ladder of Howe
  (2006), McGee and Howe (2007) and Taylor and Howe (2010), `⌊N/30 − 3⌋`,
  which stops at a tenth of the record: where PM gives way to FM beyond that
  point it applies the PM ratio to FM noise and reads 0.18 to 0.77 of the
  true variance at long τ, though on a record of one FM noise its factor
  scatters 2.5 times less.
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
  flicker-PM record was labeled white PM (in 64 % of records at m = 16 and
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
  1024-sample record, where its neighbor had −29 % / +290 %. The same
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
- `remove_outliers` replaces only the outliers: samples already missing stay
  missing, where they were filled whenever the record held an outlier and
  left when it held none. Fill them with `fillgaps` (`--fillgaps`) when that
  is wanted.
- A record without a time column can have its `NaN` readings filled
  (`fillgaps=true`, `--fillgaps`, the explorer's Fill gaps), so a one-column
  frequency file with missing readings can be analyzed; the explorer turns
  Fill gaps on for such a file.
- Explorer: hover labels are readable in both themes; the Data page's phase
  plot breaks at a gap instead of dropping to zero and draws every outlier;
  τ min then τ max can be typed in one go; the run note counts the analysis
  range; PDEV's white PM slope guide is τ^−3/2; PNG and SVG export holds
  every plot of the page.
- The explorer's local server answers 400 to a request body that is not a
  JSON object or whose `file` is malformed, where it used to drop the
  connection.
- `--series N` keeps the record's first and last sample, so the phase plot
  spans the whole record.
- `sigmatau holdover --kalman` and the explorer's Kalman predictor keep the
  flicker FM and flicker PM that the fit of the training stretch finds, as
  Markov states; without them the filter of a flicker FM clock claimed a
  tenth of its error or stopped on NaN. `holdover(data, model)` refuses a
  model with no noise.
- Explorer: *Fill gaps* no longer stays on from the previous file, and its
  box stays on the Data page when a frequency record with gaps fails to read
  without it.
- AT2 step detector (`timescale(...; algorithm=:at2, steps=…)`): a clock
  kept out of the scale after a frequency step is still estimated, as in
  Weiss & Weissert 1991 §4, so it comes back with its new frequency. Before,
  it came back with the old one, was taken out again, and a large step moved
  the scale at the return.
- `implicit_mean(r, truth)` returns Greenhall's implicit ensemble mean,
  `Σ aᵢ (xᵢ − x̂ᵢ)` (Greenhall 2006, Eqs. (3.8), (3.13)). It returned
  `Σ aᵢ xᵢ`, which is the same only for clocks of one noise shape: on clocks
  that differ in their random-walk FM it read up to 2.6 times noisier than
  the scale at long τ. `TimescaleResult` has a new field `xhat` (Julia and
  Python).
- The interval of bias-removed Thêo1 (`theo1`, and `theoh` beyond 20 % of
  the record) carries the scatter of the bias ratio on FM noise: its EDF was
  that of Thêo1 alone, and the 68.3 % interval held the true value in 15 to
  62 % of simulated FM records; it now holds it in 63 to 84 %.
  `calculate_edf` accepts `:theobr`. The EDF of Thêo1 is no longer taken
  below 1, which gives `theo1` an interval on random-walk FM beyond 63 % of
  the record. Deviations are unchanged.
- `timescale(...; algorithm=:greenhall)`: a clock without a reading has no
  weight in the implicit ensemble mean and runs on its prediction (Greenhall
  2006, §10). It kept a weight through a gap, and the scale took a phase
  step when the clock returned, or when a clock or the pivot joined late (in
  simulation, a median of 440 times the scale's rms step after a 2000-sample
  gap on Greenhall's three clocks; now 1.5). Records with every clock read
  are unchanged.
- `timescale`: in `:ewfa`, `:at1`, `:at2` and `:stein` a clock read again
  after a gap joins as a new clock does and no longer moves the scale. It
  went straight back into the average with the prediction error of the whole
  gap (in simulation, a median step of 60 to 580 times the scale's rms step
  after a 2000-sample gap on four white-FM clocks; now about 1). The AT2
  step detector's windows do not reach behind the sample a clock joined at.
  `:jones_tryon` still steps at a return, and `:stein` needs its pivot read
  throughout; both are documented.
- `timescale` with `:ewfa`, `:at1` or `:at2` diverged to NaN when the clocks
  were read less often than the grid. The interval of the equations is now
  the time since the last sample with a reading (Weiss & Weissert 1991,
  §2.8), so a sparse record gives the scale of its readings alone; records
  read at every sample are unchanged.
- In `timescale` with `:ewfa`, `:at1` or `:at2`, a clock that joins late
  starts from its own first readings (Weiss & Weissert 1991, §2.7). It took
  the largest starting variance of the other clocks and a starting frequency
  of zero, so a good clock that joined late had almost no weight for the
  rest of the record and one with a rate offset moved the scale's rate.
- A clock whose first reading comes after the start of the record is started
  against the scale, so the scale no longer depends on the reference the
  readings were taken against: `:ewfa`, `:at1` and `:at2` take its frequency
  against the scale over its own first stretch, and `:greenhall` and
  `:jones_tryon` start it from a clock already in the scale, without a step.
  A clock first read while no other clock is read no longer turns the
  `:ewfa`, `:at1` and `:at2` scales to NaN.
- `timescale` with `:greenhall` or `:jones_tryon` and models with a drift
  state no longer carries the drift of a parabola through the first 16
  readings as a constant drift of the scale. The drift starts at zero (Jones
  and Tryon 1983, Fig. 1), with a variance wide enough that the drifts of
  the clocks against each other are still learned.
- `read_phase` and `read_frequency` take the unit of a time column from its
  name in the line that `header` skips, as `read_data` does, so a column
  named MJD is read in days at any step; `sigmatau timescale --value-cols`
  reads a table with named columns without `--header`, and its read errors
  name the CLI flags.
- `sigmatau model` writes to `--out`, which it accepted and ignored, and
  refuses a FILE and the options of the record, which it never used; its
  errors are answered in JSON. `sigmatau help` sets the long command
  signatures apart from their descriptions, `help model` names the three
  deviations it writes, the help of `live` no longer speaks of a FILE, and
  `--nominal` and `--tau0` state their defaults for `live`, `model` and
  `--synth`. A test checks that the command-line reference names every
  command and option.
- Shell completion: every command completes its options and every option
  with fixed values completes them, in bash, zsh, fish and PowerShell. The
  scripts are generated from the option table and tested in each shell.
  Before, `cross`, `steer`, `timescale` and `model` completed no options,
  and no shell completed the values of `--input-format`, `--bias`,
  `--noise`, `--mask`, `--algorithm` or `--steps`.
- The compiled `sigmatau` requires Julia 1.13 to build. A build on Julia
  1.12 compiled without a warning, but the total deviations and Thêo1 failed
  in it with "failed (CompositeException)", because a program that 1.12
  trims cannot start a task. The build now stops on 1.12, and an error
  inside a threaded loop is reported by what the loop threw.
- Explorer: every page fits its window from 960 px (the ninth tab had
  squeezed the top bar, and with several clocks the Dynamic maps lay over
  each other); an error on a record with gaps names the
  *Fill gaps by reflection* box in place of a command-line flag, and the
  Ensemble page says which scales a file with missing samples can form;
  labels that said something untrue are corrected (quadratic drift removal
  and SP1065, the log drift model, the slow estimators, EWFA's largest
  weight, what a run updates); the interface spells as the manual does.
- Python: `clock_model` of a `PowerLawFit` uses Julia's default `tau_max`,
  the longest τ of a fit of a deviation. It had used 10⁶ τ₀, which with the
  flicker states on left the model with Markov states far past the record.
- Python: `delim` of `read_phase` and `read_frequency` reaches Julia; every
  call with it had raised `JuliaError`. `StreamingStability` raises
  `ValueError` on a sample that is not finite. The docstrings are corrected
  against the Julia ones, among them the sign of `timescale_phase`, what
  `remove_outliers` does and the formats of `read_data`, and `help()` of an
  estimator lists the keywords it shares with the others.
- `identify_noise`, `calculate_edf` and `bias_correction` throw an
  `ArgumentError` for an averaging factor below 1, as the estimators do;
  before, 0 raised `DivideError` or `InexactError` and a negative factor
  returned a label, an EDF or a bias factor.
- `sigmatau model` refuses `--format csv`: it wrote JSON all the same, and
  an error later on that line came as plain text. Its error object names
  `"source":"model"`, as its result does.
- The TDEV mask of `:G811` ends below 10 000 s, as G.811 clause 6.1 states;
  τ = 10 000 s had a limit of 30 ns.
- `steer` and `sigmatau steer`: the predicted Allan deviation is right for
  loops with a time constant of 10⁶ sample intervals or more; it was wrong,
  or 0.0, from there on.
- `sigmatau spectrum` and `analyze --spectrum` accept `--nperseg` at or
  above the length of the record (the segment is cut to the record) and work
  on a record of 16 samples.
- `sigmatau live` stops when the program reading its output closes the pipe,
  and no command prints the runtime's `broken pipe` line on stderr any more.
- `holdover --kalman` and `steer` say that the clock model needs a training
  stretch of 65 samples, and `steer` a record of 81, in place of
  `powerlaw_fit … pass fewer alphas`; the reference had said 64 and 16
  samples.
- `timescale` refuses `--wmax` outside 0 to 1 and `--steps` with any
  algorithm but `at2` as usage errors (exit 1); `sigmatau help nosuch` is a
  usage error; `--delim` help states its default.
- `holdover(data, model)` leaves out a trial whose starting reading is
  missing, as the least-squares test does.
- `powerlaw_sigma` and `holdover(model)` accept a model with a drift state
  and no random-run FM.
- `read_data(file; column = n)` reads column `n` of a plain column file; it
  was ignored.
- `confidence` outside 0 to 1 is an error.
- Clearer errors: `ttotdev` on a gapped record, `calculate_edf` and
  `bias_correction` with vectors of different lengths, `load` of a file
  `save` did not write, `simulate_clock(m, 0)`.
- `ttotdev` with `ci = false` returns the noise types its bias correction
  used, as `mtotdev` does.
- Explorer: on the Data page in a narrow window a long file name wraps
  inside its card, and the plot's title and hint fold before the range
  control goes to a second line; they had pushed it out of the card and
  squeezed the hint to 41 px. An error in the status bar ends when the
  inputs change or a field is edited; it stayed until the next run.
- `timescale(...; algorithm = :stein)`: the first row of `weights` is the
  scale's starting weights, and a model with a flicker level but no flicker
  states is no longer refused.
- An unknown noise type is reported as "unknown noise type :XX; use an
  integer α or one of WHPM, FLPM, WHFM, FLFM, RWFM, FWFM, RRFM" by
  `add_component!`, `clock_model`, `powerlaw_sigma`, `powerlaw_psd` and
  `holdover`; the message began "noise_gen:" whoever called, and listed the
  types in no order.
- Explorer: Holdover and Steering name a clock whose computation fails, with
  the reason, on a line above the plot beside the clocks that worked; they
  dropped it without a word.
- Explorer: removing a clock on the Synthesis page no longer raises a
  TypeError; the plot stays until the new preview arrives. The report's
  SHA-256 of a Live record is that of the file Save readings writes. The
  Ensemble sub line names the typed reference; clocks of two files with one
  name say which file the later one is from; the record check's start time
  is rounded to the microsecond. On the Data page a message over the plot
  covers the plot row only, in all four views, so the view buttons stay
  reachable.
- Python: `timescale(pivot=...)` and `simulate_ensemble(reference=...)`
  check the index against the clocks and name the value passed, and
  `reference=-1` no longer means the ideal reference; a NumPy integer is a
  seed; `ClockComparisons` makes Julia's checks when built; an unknown noise
  type is named in the error; `holdover(record, model)` refuses `grow`,
  `drift` and `anchor`.
- Explorer: on the Stability page a pinned curve of a single averaging time
  is drawn as a point, where it drew nothing, and the τ axis labels every
  other decade or fewer once it spans more than six, where the labels
  collided in a narrow window.
- `sigmatau --column NAME` on a file without column names (a plain column
  file, TimeLab, Stable32, counter log or BIPM file) is a usage error (exit
  1); it was ignored, and the default column was read.

### Security

- `serve.py` runs a request without a file only when it starts with
  `--synth N`, and its API refuses another host name on loopback and a POST
  that is not JSON.

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

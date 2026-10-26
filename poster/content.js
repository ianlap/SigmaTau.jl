// Poster copy for the SigmaTau.jl annual-review poster. Every number here is
// taken verbatim from the repository (README, docs/src/validation.md,
// docs/src/mhtotdev.md, bench/README.md, tools/artifacts/mhtotdev_mc_full.json,
// src/edf.jl). Edit this file, then run `node build_poster.js`.
module.exports = {
  title: 'Frequency stability with quantified uncertainty: per-τ noise identification, χ² intervals, and a calibrated modified-Hadamard total deviation',
  subtitle: 'SigmaTau.jl, an open and validated time-and-frequency stability library in Julia',
  authors: 'Ian Lapinski',
  affiliation: 'The University of Alabama   ·   github.com/ianlap/SigmaTau.jl   ·   MIT license',

  columns: [
    // ───────────────────────── column 1: the measurement problem and the method
    [
      { type: 'panel', heading: 'A clock record is finite, so every σ(τ) point is an estimate with an interval',
        body: [
          'Clock noise is power-law, S_{y}(f) ∝ f^{α} with α ∈ {+2, +1, 0, −1, −2} (white PM to random-walk FM). From N phase samples at spacing τ₀, the number of independent samples at averaging time τ = mτ₀ is not N/m: it depends on the noise type, on τ/T and on the estimator.',
          'A usable stability point therefore needs four quantities computed together: the deviation, the noise type at that τ, the equivalent degrees of freedom (EDF) and the χ² interval.',
        ] },
      { type: 'panel', heading: 'One call returns deviation, noise type, EDF and interval at every τ',
        bullets: [
          ['Estimators: ', 'Allan, modified Allan, time; Hadamard, modified Hadamard, Hadamard time; five totals; Thêo1/ThêoBR/ThêoH; MTIE, TIE rms, PDEV; time-resolved and real-time variants.'],
          ['Noise type per τ: ', 'lag-1 autocorrelation of the m-decimated, quadratically detrended phase differences, with the B1 / R(n) ratio as fallback.'],
          ['EDF: ', 'Greenhall–Riley 2003 for the difference families; the published total-variance tables (SP1065; Howe et al.) for TOTDEV, MTOTDEV and HTOTDEV, with documented fallbacks outside their ranges; Monte Carlo for MHTOTDEV (column 2).'],
          ['Interval: ', 'χ² bounds on σ² with ν = EDF, 68.3 % by default. Total-family bias corrections are applied by default and switchable.'],
          ['Definitions: ', 'NIST SP1065 and IEEE Std 1139-2022; every departure from a reference is a named, documented switch.'],
        ],
        figure: 'F1_suite', figH: 6.2,
        caption: 'Fig. 1. Deviation suite on the shipped 8192-point validation record (τ₀ = 1 s) with 68.3 % χ² intervals. MDEV rolls white PM off as τ^{−3/2} and separates it from flicker PM; the totals reach τ = T/4 with tighter long-τ intervals.' },
      { type: 'panel', heading: 'Real-time and multi-clock use',
        bullets: [
          ['Streaming: ', 'ADEV, MDEV, HDEV and MHDEV accumulate in O(1) per sample per τ (Dobrogowski–Kasznia 2007, extended to third differences) from a ring buffer of at most 4·m_{max} + 1 samples; the running deviation equals the batch kernel on the samples seen so far.'],
          ['Time-resolved and ensemble: ', 'dynamic σ(t, τ) maps for non-stationary records; N-cornered hat for N ≥ 3 clocks with independent noises.'],
        ] },
    ],
    // ───────────────────────── column 2: MHTOTDEV and its Monte Carlo calibration
    [
      { type: 'panel', heading: 'MHTOTDEV: a modified-Hadamard total deviation with no published coefficients',
        body: [
          'The Hadamard (third-difference) family rejects linear frequency drift; the modified form separates white from flicker PM; the total form extends each window by reflection to recover confidence at long τ. MHTOTDEV combines all three. Every other total estimator inherits a published bias and EDF table; this one had none, so both were measured.',
        ],
        figure: 'F6_mhtot_construction', figH: 4.1,
        caption: 'Fig. 2. Construction (Greenhall scheme, as for MTOTDEV and HTOTDEV): each of the N − 4m + 1 windows of length L = 4m has its half-mean slope removed and is extended by two time-reversed copies; the third difference at lag m is averaged over m positions, squared, normalised by 6m²τ² and averaged over all windows. The public call first removes the record\'s least-squares frequency drift.' },
      { type: 'panel', heading: 'Monte Carlo against synthesized noise measures both the bias and the EDF',
        bullets: [
          ['Design: ', 'α ∈ {+2, +1, 0, −1, −2}; N ∈ {1024, …, 32768}; octave m grid; R = 3000 realizations per cell; one independent random stream per (seed, α, N, r). The statistic calibrated is the shipped pipeline, including its default drift removal.'],
          ['Reference: ', 'MHDEV on the same realization (its EDF follows from the Greenhall–Riley spectral sum). Bias B = E[MHTOTVAR] / E[MHVAR]; EDF = 2 mean(V)^{2} / var(V), bootstrap standard error ≈ √(2/R) = 2.6 %.'],
          ['Fits: ', 'over τ/τ₀ ≥ 16, the validity window of the total-variance literature. EDF ν = b T/τ − c (Mod-Totvar; the TotHvar form is better by at most ΔR² = 8 × 10⁻⁴, at α = ±2) and bias B = b₀ + b₁ τ/T, both by weighted least squares.'],
        ],
        figures: ['F4_mhtot_bias', 'F5_mhtot_edf'], figH: 3.85,
        caption: 'Fig. 3. Shipped models. Left: ≈ unbiased for the PM noises, reads high for red FM (B ≈ 2.2 for random-walk FM at small τ/T). Right: EDF, R² ≥ 0.9982 for every α; below T/τ ≈ 5–6 the interval is indicative only, while σ itself is valid to the reach limit τ = T/4. Table: shipped coefficients (seed 20260531, R = 3000, N ≤ 32768).',
        table: {
          header: ['α', 'noise', 'bias b₀', 'bias b₁', 'EDF b', 'EDF c', 'EDF R²'],
          rows: [
            ['+2', 'white PM',       '1.0308', '−0.1197', '1.6869', '4.3166', '0.9982'],
            ['+1', 'flicker PM',     '0.9883', '−0.0686', '1.1671', '3.2107', '0.9988'],
            ['0',  'white FM',       '1.0545', '−0.1706', '1.0861', '3.3451', '0.9991'],
            ['−1', 'flicker FM',     '1.2947', '−0.7193', '1.0095', '3.3254', '0.9991'],
            ['−2', 'random-walk FM', '2.1686', '−5.5884', '0.7631', '2.4490', '0.9985'],
          ] } },
    ],
    // ───────────────────────── column 3: validation, performance, real-time use, next steps
    [
      { type: 'panel', heading: 'Agreement with Stable32 to its printed precision and with allantools to 10⁻¹³',
        body: [
          'On the shipped 8192-point record, ADEV, MDEV, TDEV and HDEV match Stable32 to ≤ 0.005 % (its five-significant-figure output) and allantools to ≤ 1 × 10⁻¹³; the raw total kernels match allantools to ≤ 4.4 × 10⁻¹⁵. Confidence bounds agree with Stable32 to ≤ 0.05 % at 11 of 12 τ. The one recurring residual is a convention: at τ = 256 s Stable32 assigns α = +2 where the detrended identifier assigns α = 0; one keyword reproduces Stable32 at every τ. Estimators neither reference implements (MHDEV, HTDEV, MHTOTDEV, ThêoH, nch, dynamic, streaming) are checked against inlined textbook kernels and algebraic identities.',
        ],
        figure: 'F2_validation', figH: 4.7,
        caption: 'Fig. 4. Relative difference to Stable32 at every τ for seven estimators. Filled points lie at or below the 5 × 10⁻⁵ floor set by Stable32\'s printed precision; open points are the τ = 256 s convention. MTOTDEV is compared raw, as Stable32 reports the uncorrected statistic.' },
      { type: 'panel', heading: 'Modified-total analysis of a day-long record takes seconds, not half an hour',
        body: [
          'Same samples, same τ grid, same statistic; allantools is single-threaded as it has no parallel path. On a measured 406 763-sample phase record MTOTDEV takes 0.469 s against 1 803.6 s. The gap is algorithmic as well as compiled-versus-interpreted: the reflected extension is streamed through a per-thread buffer and the window sum is split across threads.',
        ],
        figure: 'F3_speed', figH: 4.2,
        caption: 'Fig. 5. Kernel-only speed-up (no intervals, raw totals) on an N = 25 000 white-PM record (bench/, mean of 30 realizations, 12 cores; diamonds: independent 4-core re-run). Ratios are the portable result; absolute times depend on the machine.' },
      { type: 'panel', heading: 'Next steps',
        bullets: [
          ['Laboratory records: ', 'apply the drift-tolerant Hadamard family and the real-time accumulators to the Rb-87 clock comparisons.'],
          ['Calibration: ', 'extend the Monte Carlo method to the Thêo family and gapped records.'],
        ] },
    ],
  ],

  references: [
    'W. J. Riley, Handbook of Frequency Stability Analysis, NIST SP 1065 (2008).',
    'IEEE Std 1139-2022, Standard Definitions of Physical Quantities for Fundamental Frequency and Time Metrology.',
    'C. A. Greenhall and W. J. Riley, Uncertainty of stability variances based on finite differences, Proc. 35th PTTI (2003).',
    'D. A. Howe, R. L. Beard, C. A. Greenhall, F. Vernotte and W. J. Riley, A total estimator of the Hadamard function used for GPS operations, Proc. 32nd PTTI (2000).',
    'F. Vernotte and D. A. Howe, Generalization of the total variance approach to the different classes of structure functions, Proc. EFTF (2000).',
    'A. Dobrogowski and M. Kasznia, Real-time assessment of Allan deviation and time deviation, Proc. IEEE FCS / EFTF (2007).',
    'Reference implementations compared: Stable32 (W. J. Riley) and allantools (A. E. Wallin).',
  ],
  acknowledgment: 'Acknowledgments: [funding / program / advisor placeholder]. Validation fixtures, benchmark scripts and the Monte Carlo harness are in the repository.',
};

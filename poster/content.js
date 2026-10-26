// Poster copy for the SigmaTau.jl annual-review poster (expert audience:
// industry and government partners in time and frequency). Every number is
// taken verbatim from the repository (README.md, docs/src/validation.md,
// docs/src/mhtotdev.md, bench/README.md, tools/artifacts/mhtotdev_mc_full.json,
// src/edf.jl, src/streaming.jl). Edit this file, then run `node build_poster.js`.
// Markup: _{x} subscript, ^{x} superscript.
module.exports = {
  title: 'SigmaTau.jl: unified, validated and fast time-and-frequency stability analysis in Julia',
  subtitle: 'Every estimator family with per-τ noise identification, EDF and χ² intervals in one call · matches Stable32 and allantools · modified-total kernels 1,600 to 4,200× faster · a new Monte-Carlo-calibrated MHTOTDEV',
  authors: 'Ian Lapinski',
  affiliation: 'The University of Alabama   ·   github.com/ianlap/SigmaTau.jl   ·   MIT license   ·   Julia ≥ 1.11',

  columns: [
    // ───────────────────────── column 1: what the library offers
    [
      { type: 'panel', heading: 'One API, every estimator, a χ² interval with every point',
        bullets: [
          ['Allan: ', 'ADEV, MDEV, TDEV.  Hadamard: HDEV, MHDEV, HTDEV.  Total: TOTDEV, MTOTDEV, TTOTDEV, HTOTDEV, MHTOTDEV.  Thêo1 / ThêoBR / ThêoH.  MTIE, TIE rms, PDEV.  Time-resolved DADEV / DHDEV.  Real-time ADEV / MDEV / HDEV / MHDEV.  N-cornered hat.  S_{y}, S_{x}, ℒ(f), calibrated power-law synthesis.'],
          ['Per τ, in the same call: ', 'the deviation, the noise type (lag-1 autocorrelation with B1 / R(n) fallback), the Greenhall–Riley EDF, and the χ² bounds at 68.3 % (switchable). Total-family bias corrections are on by default and switchable.'],
          ['Input: ', 'phase or frequency, any τ grid; read, detrend, gap-fill and outlier tools; plain-text round-trip of results.'],
          ['Real-time and multi-clock: ', 'ADEV, MDEV, HDEV and MHDEV accumulate in O(1) per sample per τ from a ring buffer of at most 4·m_{max} + 1 samples and reproduce the batch kernels on the samples seen so far to ≤ 2 × 10⁻¹¹ relative (Dobrogowski–Kasznia 2007, extended to third differences); dynamic σ(t, τ) maps; N-cornered hat for N ≥ 3 clocks.'],
        ],
        code: ['p = read_phase("clock.DAT"; tau0 = 1.0)', 'r = adev(p)            # τ grid, σ_y, noise type, EDF, CI', 'suite = stability(p; devs = (:adev, :mdev, :hdev, :mhtotdev))'],
        table: { name: 'Feature comparison', colW: [3.6, 2.3, 2.3, 3.2], leftCols: [0],
          header: ['', 'Stable32', 'allantools', 'SigmaTau'],
          rows: [
            ['Platform', 'Windows GUI, closed', 'Python library', 'Julia library + desktop app'],
            ['χ² intervals with per-τ noise ID', 'yes', 'manual', 'every call'],
            ['Total-family bias correction', 'TOTDEV, HTOTDEV', 'none', 'all five, switchable'],
            ['MDEV-class Hadamard (MHDEV, HTDEV, MHTOTDEV)', 'no', 'no', 'yes'],
            ['Thêo1 / ThêoBR / ThêoH', 'yes', 'Thêo1 only', 'yes'],
            ['Dynamic (time-resolved) deviation', 'yes', 'no', 'yes'],
            ['Real-time accumulators', 'no', 'ADEV, HDEV, TDEV', 'ADEV, MDEV, HDEV, MHDEV'],
            ['N-cornered hat', 'three clocks', 'three clocks', 'N ≥ 3'],
            ['Multithreaded kernels', 'n/a', 'no', 'yes'],
          ] },
        tableAfter: true },
      { type: 'panel', heading: 'Deviation, noise type and interval at every τ',
        figure: 'F1_suite', figH: 6.0,
        caption: 'Fig. 1. One call on the shipped 8192-point validation record (τ₀ = 1 s): five estimators with 68.3 % χ² intervals from the noise type identified at each τ (top strip). The totals reach further (TOTDEV to τ = T/2, MHTOTDEV to T/4) with tighter long-τ bounds; MDEV separates white from flicker PM.' },
    ],
    // ───────────────────────── column 2: evidence, validation and performance
    [
      { type: 'panel', heading: 'Matches Stable32 to its printed precision and allantools to machine precision',
        body: [
          'Shipped fixture: an 8192-point phase record with the Stable32 and allantools outputs beside it. ADEV, MDEV, TDEV and HDEV agree with Stable32 to ≤ 0.005 % (its five significant figures) and with allantools to ≤ 1 × 10⁻¹³; raw total kernels agree with allantools to ≤ 4.4 × 10⁻¹⁵; confidence bounds agree with Stable32 to ≤ 0.05 % at 11 of 12 τ. Where a default differs from a reference (MTOTDEV bias policy, noise-ID detrending), the per-τ tables name the keyword that reproduces it.',
        ],
        figure: 'F2_validation', figH: 4.3,
        caption: 'Fig. 2. Relative difference to Stable32 at every τ for seven estimators. Filled points lie at or below the 5 × 10⁻⁵ floor set by Stable32\'s printed precision; open points are the τ = 256 s noise-ID convention (Stable32 α = +2, detrended identifier α = 0). MTOTDEV compared raw, as Stable32 reports the uncorrected statistic.' },
      { type: 'panel', heading: 'Modified totals on a 4 × 10⁵-sample record: 0.47 s instead of 30 min',
        body: [
          'Same samples, same τ grid, same statistic, kernel-only timing; allantools is single-threaded as it has no parallel path. On a measured 406 763-sample phase record MTOTDEV takes 0.469 s against 1 803.6 s and HTOTDEV 0.455 s against 1 928.0 s (3,842× and 4,233×). The gap is algorithmic as well as compiled-versus-interpreted: the reflected extension is streamed through a per-thread buffer and the window sum is split across threads.',
        ],
        figure: 'F3_speed', figH: 4.2,
        caption: 'Fig. 3. Speed-up per kernel on an N = 25 000 white-PM record (bench/, 12-core workstation; diamonds: independent 4-core re-run). Ratios are the portable result; absolute times depend on the machine.' },
      { type: 'panel', heading: 'Reproducible by construction',
        bullets: [
          ['Shipped with the code: ', 'the validation fixture with Stable32 and allantools outputs, the per-τ comparison tables, the benchmark scripts for both libraries, and the Monte Carlo harness with its seed and tracked fit artifact.'],
          ['Next: ', 'real-time monitoring of the laboratory Rb-87 clock comparisons; Monte Carlo calibration of the Thêo family and of gapped records.'],
        ] },
    ],
    // ───────────────────────── column 3: new science and reproducibility
    [
      { type: 'panel', heading: 'MHTOTDEV: a new total estimator, calibrated by Monte Carlo',
        body: [
          'Modified Hadamard total deviation combines drift rejection (third differences), white-from-flicker PM separation (modified form) and long-τ confidence (reflected extension). No published bias or EDF existed, so both were measured against synthesized noise: α ∈ {+2 … −2}, N ∈ {1024 … 32768}, R = 3000 realizations per cell, one seeded stream per (α, N, r), fits over τ/τ₀ ≥ 16. Bias B = E[MHTOTVAR]/E[MHVAR] against MHDEV on the same realization; EDF from the χ² variance relation, bootstrap standard error 2.6 %.',
          'Kernel: the Greenhall scheme of MTOTDEV and HTOTDEV on windows of length 4m (half-mean slope removal, 3× time-reverse extension, modified third difference), after global least-squares drift removal.',
        ],
      },
      { type: 'panel', heading: 'Unbiased for PM noise, reads high for red FM; EDF fits with R² ≥ 0.9982',
        figures: ['F4_mhtot_bias', 'F5_mhtot_edf'], figH: 4.0,
        caption: 'Fig. 4. Shipped models B = b₀ + b₁ τ/T and ν = b T/τ − c (seed 20260531). Below T/τ ≈ 5–6 the interval is indicative only; σ itself is valid to the reach limit τ = T/4.',
        table: { name: 'MHTOTDEV coefficients',
          header: ['α', 'noise', 'bias b₀', 'bias b₁', 'EDF b', 'EDF c', 'R²'],
          rows: [
            ['+2', 'white PM',       '1.0308', '−0.1197', '1.6869', '4.3166', '0.9982'],
            ['+1', 'flicker PM',     '0.9883', '−0.0686', '1.1671', '3.2107', '0.9988'],
            ['0',  'white FM',       '1.0545', '−0.1706', '1.0861', '3.3451', '0.9991'],
            ['−1', 'flicker FM',     '1.2947', '−0.7193', '1.0095', '3.3254', '0.9991'],
            ['−2', 'random-walk FM', '2.1686', '−5.5884', '0.7631', '2.4490', '0.9985'],
          ] } },
      { type: 'panel', heading: 'Time-resolved stability: σ(t, τ) maps for non-stationary records',
        figure: 'F7_dynamic', figH: 5.0,
        caption: 'Fig. 5. DADEV on a synthesized 8191-sample record (τ₀ = 1 s): white FM at σ_{y}(1 s) = 1 × 10⁻¹² for the first half, then 3 × 10⁻¹² plus flicker FM at 1 × 10⁻¹²; window 1024 samples, step 128. Each cell is an ordinary overlapping ADEV of its window, so the change is located to within one window length; a 3-D view is in the desktop explorer.' },
    ],
  ],

  references: [
    'W. J. Riley, Handbook of Frequency Stability Analysis, NIST SP 1065 (2008).',
    'IEEE Std 1139-2022, Standard Definitions of Physical Quantities for Fundamental Frequency and Time Metrology.',
    'C. A. Greenhall and W. J. Riley, Uncertainty of stability variances based on finite differences, Proc. 35th PTTI (2003).',
    'D. A. Howe, R. L. Beard, C. A. Greenhall, F. Vernotte and W. J. Riley, A total estimator of the Hadamard function used for GPS operations, Proc. 32nd PTTI (2000).',
    'F. Vernotte and D. A. Howe, Generalization of the total variance approach to the different classes of structure functions, Proc. 14th EFTF (2000).',
    'A. Dobrogowski and M. Kasznia, Real-time assessment of Allan deviation and time deviation, Proc. IEEE FCS / EFTF (2007).',
    'Reference implementations compared: Stable32 (W. J. Riley) and allantools (A. E. Wallin).',
  ],
  acknowledgment: 'Acknowledgments: [funding / program / advisor placeholder]. Documentation: ianlap.github.io/SigmaTau.jl',
};

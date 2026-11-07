// Poster copy for the SigmaTau.jl annual-review poster (expert audience). Statements only;
// every number is verbatim from the repository (README.md, docs/src/validation.md,
// docs/src/mhtotdev.md, bench/README.md, tools/artifacts/mhtotdev_mc_full.json, src/edf.jl,
// src/streaming.jl, poster/data/*.csv). Markup: _{x} subscript, ^{x} superscript.
// Edit, then run `node build_poster.js`.
module.exports = {
  "title": "SigmaTau.jl: time-and-frequency stability analysis in Julia",
  "subtitle": "Every estimator family, per-τ noise ID, EDF and χ² intervals in one call · matches Stable32 and allantools · modified-total kernels 1,600 to 4,200× faster · new Monte-Carlo-calibrated MHTOTDEV",
  "authors": "Ian Lapinski",
  "affiliation": "The University of Alabama   ·   github.com/ianlap/SigmaTau.jl   ·   MIT license   ·   Julia ≥ 1.11",
  "columns": [
    [
      {
        "type": "panel",
        "heading": "One API: every estimator, noise ID, EDF, χ² bounds",
        "bullets": [
          [
            "Estimators: ",
            "ADEV/MDEV/TDEV; HDEV/MHDEV/HTDEV; TOTDEV, MTOTDEV, TTOTDEV, HTOTDEV, MHTOTDEV; Thêo1/ThêoBR/ThêoH; MTIE, TIE rms, PDEV; DADEV/DHDEV; S_{y}, S_{x}, ℒ(f); calibrated power-law synthesis."
          ],
          [
            "Per τ: ",
            "σ, noise type (lag-1; B1/R(n) fallback), Greenhall–Riley EDF, 68.3 % χ² bounds (switchable). Total-family bias correction default on, switchable."
          ],
          [
            "Input: ",
            "phase or frequency, any τ grid; read/detrend/gap-fill/outlier tools; plain-text round-trip."
          ],
          [
            "Real-time, multi-clock: ",
            "ADEV/MDEV/HDEV/MHDEV in O(1) per sample per τ; ring buffer ≤ 4·m_{max} + 1 samples; ≤ 2 × 10⁻¹¹ relative to batch kernels on samples seen so far (Dobrogowski–Kasznia 2007, extended to third differences). N-cornered hat, N ≥ 3."
          ]
        ],
        "code": [
          "p = read_phase(\"clock.DAT\"; tau0 = 1.0)",
          "r = adev(p)            # τ grid, σ_y, noise type, EDF, CI",
          "suite = stability(p; devs = (:adev, :mdev, :hdev, :mhtotdev))"
        ],
        "table": {
          "name": "Feature comparison",
          "colW": [
            3.6,
            2.3,
            2.3,
            3.2
          ],
          "leftCols": [
            0
          ],
          "header": [
            "",
            "Stable32",
            "allantools",
            "SigmaTau"
          ],
          "rows": [
            [
              "Platform",
              "Windows GUI, closed",
              "Python library",
              "Julia library + desktop app"
            ],
            [
              "χ² intervals with per-τ noise ID",
              "yes",
              "manual",
              "every call"
            ],
            [
              "Total-family bias correction",
              "TOTDEV, HTOTDEV",
              "none",
              "all five, switchable"
            ],
            [
              "Modified Hadamard Class (MHDEV, HTDEV, MHTOTDEV)",
              "no",
              "no",
              "yes"
            ],
            [
              "Thêo1 / ThêoBR / ThêoH",
              "yes",
              "Thêo1 only",
              "yes"
            ],
            [
              "Dynamic (time-resolved) deviation",
              "yes",
              "no",
              "yes"
            ],
            [
              "Real-time accumulators",
              "no",
              "ADEV, HDEV, TDEV",
              "ADEV, MDEV, HDEV, MHDEV"
            ],
            [
              "N-cornered hat",
              "three clocks",
              "three clocks",
              "N ≥ 3"
            ],
            [
              "Multithreaded kernels",
              "n/a",
              "no",
              "yes"
            ]
          ]
        },
        "tableAfter": true
      },
      {
        "type": "panel",
        "heading": "Deviation, noise type and χ² interval at every τ",
        "figure": "F1_suite",
        "figH": 7.4,
        "caption": "Fig. 1. Shipped 8192-point record, τ₀ = 1 s, one call: five estimators, per-τ noise type (top strip), 68.3 % χ² intervals. TOTDEV reaches T/2, MHTOTDEV T/4, tighter long-τ bounds."
      }
    ],
    [
      {
        "type": "panel",
        "heading": "Matches Stable32 to printed precision, allantools to machine precision",
        "body": [
          "Shipped fixture: 8192-point phase record, Stable32 and allantools outputs. ADEV/MDEV/TDEV/HDEV: ≤ 0.005 % vs Stable32 (five significant figures); ≤ 1 × 10⁻¹³ vs allantools. Raw total kernels: ≤ 4.4 × 10⁻¹⁵ vs allantools. Confidence bounds: ≤ 0.05 % vs Stable32 at 11 of 12 τ."
        ],
        "figure": "F2_validation",
        "figH": 4.3,
        "caption": "Fig. 2. Seven estimators vs Stable32, every τ. Filled: ≤ 5 × 10⁻⁵ printed-precision floor. Open: τ = 256 s, Stable32 α = +2, detrended α = 0. MTOTDEV raw."
      },
      {
        "type": "panel",
        "heading": "4 × 10⁵-sample modified totals: 0.46 s vs 30 min",
        "body": [
          "Measured 406 763-sample phase record, kernel-only timing, same samples, τ grid and statistic; allantools single-threaded (no parallel path). MTOTDEV 0.461 s vs 1 810 s (3,926×); HTOTDEV 0.461 s vs 1 926 s (4,178×). Reflected extension streamed per thread; window sum split across threads."
        ],
        "figure": "F3_speed",
        "figH": 4.1,
        "caption": "Fig. 3. Speed-up per kernel, N = 25 000 white-PM record (bench/), 12-thread workstation; diamonds: 4-core container run. Ratios portable; absolute times machine-dependent."
      },
      {
        "type": "panel",
        "heading": "Desktop explorer: the whole workflow without code",
        "figure": "F8_explorer",
        "figH": 4.4,
        "caption": "Fig. 4. GLMakie desktop app: load or synthesize a record; preprocessing, estimators, τ grid, confidence level; σ_{y} and σ_{x} curves with χ² intervals, per-τ table, spectra, dynamic map; PNG/PDF/SVG and TSV export."
      }
    ],
    [
      {
        "type": "panel",
        "heading": "MHTOTDEV: new total estimator, calibrated by Monte Carlo",
        "body": [
          "Monte Carlo: α ∈ {+2 … −2}, N ∈ {1024 … 32768}, R = 3000 per cell, one seeded stream per (α, N, r), fit window τ/τ₀ ≥ 16. Bias B = E[MHTOTVAR]/E[MHVAR] vs MHDEV, same realization; EDF from χ² variance relation; bootstrap SE 2.6 %.",
          "Kernel: Greenhall MTOTDEV/HTOTDEV scheme on 4m windows (half-mean slope removal, 3× time-reverse extension, modified third difference) after global least-squares drift removal."
        ]
      },
      {
        "type": "panel",
        "heading": "PM unbiased, red FM reads high; EDF R² ≥ 0.9988",
        "figures": [
          "F4_mhtot_bias",
          "F5_mhtot_edf"
        ],
        "figH": 3.9,
        "caption": "Fig. 5. Models B = b₀ + b₁ τ/T, ν = b T/τ − c (seed 20260531). Interval indicative only below T/τ ≈ 5–6; σ valid to τ = T/4.",
        "table": {
          "name": "MHTOTDEV coefficients",
          "header": [
            "α",
            "noise",
            "bias b₀",
            "bias b₁",
            "EDF b",
            "EDF c",
            "R²"
          ],
          "rows": [
            [
              "+2",
              "white PM",
              "1.0249",
              "−0.0078",
              "1.6342",
              "3.8334",
              "0.9989"
            ],
            [
              "+1",
              "flicker PM",
              "0.9879",
              "−0.0630",
              "1.1651",
              "3.1934",
              "0.9988"
            ],
            [
              "0",
              "white FM",
              "1.0545",
              "−0.1706",
              "1.0861",
              "3.3451",
              "0.9991"
            ],
            [
              "−1",
              "flicker FM",
              "1.2946",
              "−0.7752",
              "1.0096",
              "3.3662",
              "0.9992"
            ],
            [
              "−2",
              "random-walk FM",
              "2.1678",
              "−5.6243",
              "0.7641",
              "2.4152",
              "0.9988"
            ]
          ]
        }
      },
      {
        "type": "panel",
        "heading": "DADEV σ(t, τ) map: change located to within one window",
        "figure": "F7_dynamic",
        "figH": 4.0,
        "caption": "Fig. 6. Synthetic 8191-sample record, τ₀ = 1 s: white FM σ_{y}(1 s) = 1 × 10⁻¹² first half, then 3 × 10⁻¹² plus flicker FM 1 × 10⁻¹²; window 1024, step 128."
      },
      {
        "type": "panel",
        "heading": "Every result reproducible from the repository",
        "bullets": [
          [
            "Shipped: ",
            "validation fixture with Stable32/allantools outputs; per-τ comparison tables; benchmark scripts for both libraries; Monte Carlo harness, seed, tracked fit artifact."
          ],
          [
            "Next: ",
            "real-time monitoring of laboratory Rb-87 clock comparisons; Monte Carlo calibration of the Thêo family and gapped records."
          ]
        ]
      }
    ]
  ],
  "references": [
    "W. J. Riley, Handbook of Frequency Stability Analysis, NIST SP 1065 (2008).",
    "IEEE Std 1139-2022, Standard Definitions of Physical Quantities for Fundamental Frequency and Time Metrology.",
    "C. A. Greenhall and W. J. Riley, Uncertainty of stability variances based on finite differences, Proc. 35th PTTI (2003).",
    "D. A. Howe, R. L. Beard, C. A. Greenhall, F. Vernotte and W. J. Riley, A total estimator of the Hadamard function used for GPS operations, Proc. 32nd PTTI (2000).",
    "F. Vernotte and D. A. Howe, Generalization of the total variance approach to the different classes of structure functions, Proc. 14th EFTF (2000).",
    "A. Dobrogowski and M. Kasznia, Real-time assessment of Allan deviation and time deviation, Proc. IEEE FCS / EFTF (2007).",
    "Reference implementations compared: Stable32 (W. J. Riley) and allantools (A. E. Wallin)."
  ],
  "acknowledgment": "Acknowledgments: [funding / program / advisor placeholder]. Documentation: ianlap.github.io/SigmaTau.jl"
};

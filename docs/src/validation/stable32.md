# Validation: Stable32

This page documents agreement between SigmaTau.jl outputs and reference
fixtures generated from `test/fixtures/validation/stable32gen.DAT` — an
8192-point phase record at τ₀ = 1 s. Two references are compared at every
τ:

- **Stable32** (`stable32out/stable32_data_full.csv`): σ, the identified
  noise exponent α, and the 68.3 % confidence bounds, all printed to
  5 significant figures. That rounding sets the agreement floor of
  roughly `5e-5` for any comparison against this file.
- **allantools** (`allantools_out/allantools_data_full.csv`): σ printed
  at full Float64 precision, so kernel-level agreement can be checked to
  machine precision.

All SigmaTau values below were computed from the shipped fixture with the
current code:

```julia
p = read_phase("test/fixtures/validation/stable32gen.DAT";
               time_col=0, value_col=1, tau0=1.0, header=10)
```

Default settings throughout; the total-family tables additionally show
`correct_bias=false` ("raw") columns.

## Summary

| Estimator | σ vs Stable32 | σ (raw) vs allantools | CI bounds vs Stable32 |
|---|---|---|---|
| OADEV, MDEV, TDEV, OHDEV | ≤ 0.005 % (fixture precision) | ≤ 1e-13 | ≤ 0.05 % at 11 of 12 τ¹ |
| TOTDEV | ≤ 0.005 % at all 13 τ | ≤ 2.8e-15 | ≤ 0.07 % at 10 of 13 τ¹ ² |
| MTOTDEV (raw) | ≤ 0.004 % | ≤ 1.1e-15 | within 3.2 % (≤ 0.02 % at τ ≤ 8 s)³ |
| HTOTDEV | ≤ 0.003 % at 11 of 12 τ¹ | ≤ 4.4e-15 | ≤ 0.1 % at 9 of 12 τ¹ ² |

¹ The single recurring residual is a noise-identification convention, not
a numerical disagreement — see below.
² At the two flicker-PM τ (32 s, 64 s) the bounds differ by up to 0.8 %
(TOTDEV) / 1.1 % (HTOTDEV): the published TOTVAR/HTOT EDF tables stop at
α = 0, and the two programs substitute different EDF models above it.
³ Stable32's MTOT EDF model and SP1065 Table 8 drift apart at long τ;
the bounds agree to ≤ 0.02 % at τ ≤ 8 s and to within 3.2 % at the
longest τ.

Three per-estimator facts behind those rows:

- **TOTDEV** — both programs apply the Howe/Walter `B = 1 − a·τ/T`
  correction for FLFM/RWFM. SigmaTau's default (`correct_bias=true`)
  matches Stable32's σ to ≤ 0.005 % at every τ including the FLFM row at
  τ = 512 s; `correct_bias=false` reproduces allantools' raw SP1065
  eq. 25 kernel to ≤ 2.8e-15.
- **HTOTDEV** — Stable32 *applies* the Howe 2005 bias correction by
  default (the factors are exactly `1/√(1+a)`), and so does SigmaTau:
  `correct_bias=true` (the default) matches Stable32 to ≤ 0.003 % at
  every matching-α τ. `correct_bias=false` reproduces allantools' raw
  kernel to ≤ 4.4e-15.
- **MTOTDEV** — here the two programs deliberately report different
  statistics. Stable32 reports the raw estimator; SigmaTau's raw kernel
  matches it to ≤ 0.004 % (and allantools to ≤ 1.1e-15). SigmaTau's
  default divides by `√B(α)` (SP1065 Table 11) to unbias toward the
  parent modified-Allan statistic, which places the default output
  2.9 % (WHPM) to 12.3 % (FLFM) below Stable32 on this fixture — up to
  ≈ 14.5 % at RWFM per the table. Pass `correct_bias=false` to reproduce
  Stable32 and allantools.

### The noise-identification residual

At τ = 256 s (every estimator) Stable32 assigns α = +2 where SigmaTau's
default identification assigns α = 0. This is a convention difference in
the lag-1 ACF identifier: SigmaTau quadratically detrends each decimated
subseries before classifying (`detrend=true` default); Stable32 does not.
`identify_noise(x, m_values; detrend=false)` reproduces Stable32's α at
all 12 τ on this fixture. Because the EDF — and, for the total family,
the bias factor — depend on α, the flipped row carries different CI
bounds (and, for HTOTDEV, a 0.25 % σ difference). The rows are marked ¹
in the tables below.

### Confidence intervals

Stable32 has computed its confidence intervals with the Greenhall–Riley
2003 EDF algorithm since v1.41 (stated in the paper itself, p. 268), and
SigmaTau implements the same algorithm: wherever the two programs agree
on α, the 68.3 % bounds agree to ≤ 0.05 % for the core estimators. The
remaining total-family differences at α > 0 (footnote ² above) come from
EDF substitutions outside the published coefficient tables, not from a
different interval method.

## Detailed comparison

The tables below compare Stable32, allantools, and SigmaTau at each
averaging interval τ. Stable32 columns carry 5 significant figures (the
fixture's precision); SigmaTau and allantools columns are truncated to 6.
CI columns are the 68.3 % `[lower, upper]` bounds — Stable32's
MinSigma/MaxSigma against SigmaTau's `ci_lower`/`ci_upper`. Rows marked ¹
are the τ = 256 s noise-identification flip.

### Overlapping Allan

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|
| 1 | 1.0097e+00 | 1.00975e+00 | 1.00975e+00 | [9.9898e-01, 1.0209e+00] | [9.9892e-01, 1.0209e+00] |
| 2 | 5.0444e-01 | 5.04443e-01 | 5.04443e-01 | [4.9907e-01, 5.1000e-01] | [4.9903e-01, 5.1003e-01] |
| 4 | 2.5021e-01 | 2.50214e-01 | 2.50214e-01 | [2.4755e-01, 2.5297e-01] | [2.4753e-01, 2.5299e-01] |
| 8 | 1.3063e-01 | 1.30632e-01 | 1.30632e-01 | [1.2924e-01, 1.3207e-01] | [1.2923e-01, 1.3208e-01] |
| 16 | 6.2750e-02 | 6.27503e-02 | 6.27503e-02 | [6.2081e-02, 6.3442e-02] | [6.2076e-02, 6.3447e-02] |
| 32 | 3.1804e-02 | 3.18042e-02 | 3.18042e-02 | [3.1141e-02, 3.2512e-02] | [3.1136e-02, 3.2517e-02] |
| 64 | 1.6085e-02 | 1.60852e-02 | 1.60852e-02 | [1.5668e-02, 1.6537e-02] | [1.5668e-02, 1.6537e-02] |
| 128 | 8.1063e-03 | 8.10630e-03 | 8.10630e-03 | [7.5737e-03, 8.7695e-03] | [7.5735e-03, 8.7700e-03] |
| 256¹ | 4.5004e-03 | 4.50037e-03 | 4.50037e-03 | [4.4513e-03, 4.5511e-03] | [4.0951e-03, 5.0557e-03] |
| 512 | 3.1642e-03 | 3.16419e-03 | 3.16419e-03 | [2.7368e-03, 3.8843e-03] | [2.7365e-03, 3.8851e-03] |
| 1024 | 3.1373e-03 | 3.13732e-03 | 3.13732e-03 | [2.6177e-03, 4.1767e-03] | [2.6175e-03, 4.1769e-03] |
| 2048 | 4.6131e-03 | 4.61311e-03 | 4.61311e-03 | [3.5914e-03, 7.7529e-03] | [3.5909e-03, 7.7548e-03] |

### Modified Allan

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|
| 1 | 1.0097e+00 | 1.00975e+00 | 1.00975e+00 | [9.9898e-01, 1.0209e+00] | [9.9892e-01, 1.0209e+00] |
| 2 | 3.5650e-01 | 3.56499e-01 | 3.56499e-01 | [3.5252e-01, 3.6062e-01] | [3.5249e-01, 3.6065e-01] |
| 4 | 1.2611e-01 | 1.26111e-01 | 1.26111e-01 | [1.2434e-01, 1.2796e-01] | [1.2433e-01, 1.2797e-01] |
| 8 | 4.6674e-02 | 4.66743e-02 | 4.66743e-02 | [4.5785e-02, 4.7618e-02] | [4.5779e-02, 4.7624e-02] |
| 16 | 1.5818e-02 | 1.58178e-02 | 1.58178e-02 | [1.5400e-02, 1.6272e-02] | [1.5397e-02, 1.6275e-02] |
| 32 | 6.0448e-03 | 6.04484e-03 | 6.04484e-03 | [5.7952e-03, 6.3299e-03] | [5.7936e-03, 6.3318e-03] |
| 64 | 2.5945e-03 | 2.59450e-03 | 2.59450e-03 | [2.4463e-03, 2.7736e-03] | [2.4453e-03, 2.7747e-03] |
| 128 | 1.7347e-03 | 1.73473e-03 | 1.73473e-03 | [1.5955e-03, 1.9182e-03] | [1.5954e-03, 1.9183e-03] |
| 256¹ | 1.8372e-03 | 1.83717e-03 | 1.83717e-03 | [1.6585e-03, 2.0896e-03] | [1.6358e-03, 2.1375e-03] |
| 512 | 2.2115e-03 | 2.21151e-03 | 2.21151e-03 | [1.8801e-03, 2.8157e-03] | [1.8800e-03, 2.8161e-03] |
| 1024 | 2.7671e-03 | 2.76714e-03 | 2.76714e-03 | [2.2099e-03, 4.2001e-03] | [2.2094e-03, 4.2012e-03] |
| 2048 | 5.2906e-03 | 5.29055e-03 | 5.29055e-03 | [3.8697e-03, 1.3794e-02] | [3.8687e-03, 1.3800e-02] |

### Time

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|
| 1 | 5.8298e-01 | 5.82978e-01 | 5.82978e-01 | [5.7676e-01, 5.8940e-01] | [5.7672e-01, 5.8944e-01] |
| 2 | 4.1165e-01 | 4.11650e-01 | 4.11650e-01 | [4.0705e-01, 4.1641e-01] | [4.0702e-01, 4.1644e-01] |
| 4 | 2.9124e-01 | 2.91241e-01 | 2.91241e-01 | [2.8715e-01, 2.9551e-01] | [2.8713e-01, 2.9554e-01] |
| 8 | 2.1558e-01 | 2.15579e-01 | 2.15579e-01 | [2.1147e-01, 2.1994e-01] | [2.1144e-01, 2.1997e-01] |
| 16 | 1.4612e-01 | 1.46119e-01 | 1.46119e-01 | [1.4226e-01, 1.5032e-01] | [1.4223e-01, 1.5034e-01] |
| 32 | 1.1168e-01 | 1.11680e-01 | 1.11680e-01 | [1.0707e-01, 1.1695e-01] | [1.0704e-01, 1.1698e-01] |
| 64 | 9.5868e-02 | 9.58679e-02 | 9.58679e-02 | [9.0391e-02, 1.0249e-01] | [9.0356e-02, 1.0253e-01] |
| 128 | 1.2820e-01 | 1.28198e-01 | 1.28198e-01 | [1.1791e-01, 1.4175e-01] | [1.1790e-01, 1.4176e-01] |
| 256¹ | 2.7154e-01 | 2.71537e-01 | 2.71537e-01 | [2.4513e-01, 3.0884e-01] | [2.4177e-01, 3.1592e-01] |
| 512 | 6.5373e-01 | 6.53730e-01 | 6.53730e-01 | [5.5577e-01, 8.3233e-01] | [5.5573e-01, 8.3245e-01] |
| 1024 | 1.6359e+00 | 1.63595e+00 | 1.63595e+00 | [1.3065e+00, 2.4831e+00] | [1.3062e+00, 2.4838e+00] |
| 2048 | 6.2556e+00 | 6.25562e+00 | 6.25562e+00 | [4.5756e+00, 1.6310e+01] | [4.5744e+00, 1.6317e+01] |

### Overlapping Hadamard

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|
| 1 | 1.0646e+00 | 1.06459e+00 | 1.06459e+00 | [1.0522e+00, 1.0774e+00] | [1.0522e+00, 1.0775e+00] |
| 2 | 5.3205e-01 | 5.32049e-01 | 5.32049e-01 | [5.2587e-01, 5.3845e-01] | [5.2584e-01, 5.3849e-01] |
| 4 | 2.6192e-01 | 2.61924e-01 | 2.61924e-01 | [2.5888e-01, 2.6507e-01] | [2.5886e-01, 2.6509e-01] |
| 8 | 1.3854e-01 | 1.38538e-01 | 1.38538e-01 | [1.3693e-01, 1.4021e-01] | [1.3692e-01, 1.4022e-01] |
| 16 | 6.6018e-02 | 6.60177e-02 | 6.60177e-02 | [6.5250e-02, 6.6813e-02] | [6.5245e-02, 6.6818e-02] |
| 32 | 3.3500e-02 | 3.35003e-02 | 3.35003e-02 | [3.2736e-02, 3.4321e-02] | [3.2740e-02, 3.4316e-02] |
| 64 | 1.6951e-02 | 1.69515e-02 | 1.69515e-02 | [1.6477e-02, 1.7469e-02] | [1.6477e-02, 1.7469e-02] |
| 128 | 8.4500e-03 | 8.44998e-03 | 8.44998e-03 | [7.8517e-03, 9.2096e-03] | [7.8514e-03, 9.2102e-03] |
| 256¹ | 4.5672e-03 | 4.56718e-03 | 4.56718e-03 | [4.5122e-03, 4.6242e-03] | [4.1228e-03, 5.1951e-03] |
| 512 | 2.9510e-03 | 2.95105e-03 | 2.95105e-03 | [2.5185e-03, 3.7257e-03] | [2.5182e-03, 3.7265e-03] |
| 1024 | 2.4608e-03 | 2.46084e-03 | 2.46084e-03 | [2.0108e-03, 3.4682e-03] | [2.0106e-03, 3.4687e-03] |
| 2048 | 2.8169e-03 | 2.81691e-03 | 2.81691e-03 | [2.1331e-03, 5.4864e-03] | [2.1326e-03, 5.4890e-03] |

### Total

`SigmaTau corrected` is the default (`correct_bias=true`); it differs
from the raw column only at the FLFM row (τ = 512 s), where the
Howe/Walter `B = 1 − a·τ/T` correction applies — and that is the row
where Stable32 matches the corrected value, confirming Stable32 applies
the same correction. The SigmaTau CI column belongs to the default
(corrected) output.

| τ (s) | Stable32 | allantools (raw) | SigmaTau raw | SigmaTau corrected | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|:---|
| 1 | 1.0097e+00 | 1.00975e+00 | 1.00975e+00 | 1.00975e+00 | [9.9884e-01, 1.0210e+00] | [9.9892e-01, 1.0209e+00] |
| 2 | 5.0443e-01 | 5.04428e-01 | 5.04428e-01 | 5.04428e-01 | [4.9898e-01, 5.1006e-01] | [4.9902e-01, 5.1002e-01] |
| 4 | 2.5026e-01 | 2.50255e-01 | 2.50255e-01 | 2.50255e-01 | [2.4755e-01, 2.5305e-01] | [2.4757e-01, 2.5303e-01] |
| 8 | 1.3060e-01 | 1.30605e-01 | 1.30605e-01 | 1.30605e-01 | [1.2919e-01, 1.3206e-01] | [1.2920e-01, 1.3205e-01] |
| 16 | 6.2782e-02 | 6.27816e-02 | 6.27816e-02 | 6.27816e-02 | [6.2103e-02, 6.3483e-02] | [6.2107e-02, 6.3478e-02] |
| 32 | 3.1842e-02 | 3.18422e-02 | 3.18422e-02 | 3.18422e-02 | [3.1345e-02, 3.2364e-02] | [3.1174e-02, 3.2556e-02] |
| 64 | 1.6090e-02 | 1.60901e-02 | 1.60901e-02 | 1.60901e-02 | [1.5791e-02, 1.6407e-02] | [1.5673e-02, 1.6542e-02] |
| 128 | 8.1074e-03 | 8.10739e-03 | 8.10739e-03 | 8.10739e-03 | [7.5805e-03, 8.7619e-03] | [7.5801e-03, 8.7624e-03] |
| 256¹ | 4.4756e-03 | 4.47557e-03 | 4.47557e-03 | 4.47557e-03 | [4.4264e-03, 4.5264e-03] | [4.0807e-03, 5.0127e-03] |
| 512 | 3.1028e-03 | 3.05583e-03 | 3.05583e-03 | 3.10283e-03 | [2.6974e-03, 3.7709e-03] | [2.6974e-03, 3.7709e-03] |
| 1024 | 2.8694e-03 | 2.86942e-03 | 2.86942e-03 | 2.86942e-03 | [2.4276e-03, 3.6939e-03] | [2.4273e-03, 3.6948e-03] |
| 2048 | 3.4123e-03 | 3.41229e-03 | 3.41229e-03 | 3.41229e-03 | [2.7445e-03, 5.0548e-03] | [2.7439e-03, 5.0558e-03] |
| 4096 | 4.0832e-03 | 4.08324e-03 | 4.08324e-03 | 4.08324e-03 | [3.1056e-03, 7.7435e-03] | [3.1048e-03, 7.7483e-03] |

### Modified Total

Stable32 and allantools both report the raw estimator here; the SigmaTau
raw column matches both (≤ 0.004 % / ≤ 1.1e-15). The corrected column is
SigmaTau's default — the raw value divided by `√B(α)` from SP1065
Table 11, which on this fixture is 2.9 % (WHPM), 7.5 % (FLPM), 11.3 %
(WHFM), or 12.3 % (FLFM) below Stable32. The SigmaTau CI column belongs
to the raw output for like-for-like comparison with Stable32's bounds.

| τ (s) | Stable32 | allantools (raw) | SigmaTau raw | SigmaTau corrected | Stable32 CI | SigmaTau CI (raw) |
|---:|:---|:---|:---|:---|:---|:---|
| 1 | 7.1400e-01 | 7.13999e-01 | 7.13999e-01 | 6.93497e-01 | [7.0639e-01, 7.2186e-01] | [7.0634e-01, 7.2191e-01] |
| 2 | 3.5625e-01 | 3.56253e-01 | 3.56253e-01 | 3.46024e-01 | [3.5227e-01, 3.6037e-01] | [3.5225e-01, 3.6040e-01] |
| 4 | 1.2570e-01 | 1.25704e-01 | 1.25704e-01 | 1.22095e-01 | [1.2394e-01, 1.2755e-01] | [1.2393e-01, 1.2756e-01] |
| 8 | 4.6303e-02 | 4.63031e-02 | 4.63031e-02 | 4.49735e-02 | [4.5421e-02, 4.7239e-02] | [4.5415e-02, 4.7245e-02] |
| 16 | 1.5819e-02 | 1.58191e-02 | 1.58191e-02 | 1.53648e-02 | [1.5495e-02, 1.6165e-02] | [1.5472e-02, 1.6191e-02] |
| 32 | 5.8856e-03 | 5.88559e-03 | 5.88559e-03 | 5.44123e-03 | [5.6979e-03, 6.0933e-03] | [5.6612e-03, 6.1389e-03] |
| 64 | 2.4483e-03 | 2.44833e-03 | 2.44833e-03 | 2.26348e-03 | [2.3398e-03, 2.5737e-03] | [2.3192e-03, 2.6017e-03] |
| 128 | 1.5120e-03 | 1.51199e-03 | 1.51199e-03 | 1.34167e-03 | [1.4073e-03, 1.6442e-03] | [1.3982e-03, 1.6589e-03] |
| 256¹ | 1.5468e-03 | 1.54684e-03 | 1.54684e-03 | 1.37260e-03 | [1.4288e-03, 1.7001e-03] | [1.3889e-03, 1.7747e-03] |
| 512 | 1.8467e-03 | 1.84668e-03 | 1.84668e-03 | 1.61965e-03 | [1.5780e-03, 2.3251e-03] | [1.5713e-03, 2.3468e-03] |
| 1024 | 2.2545e-03 | 2.25454e-03 | 2.25454e-03 | 2.00058e-03 | [1.8657e-03, 3.0666e-03] | [1.8461e-03, 3.1572e-03] |
| 2048 | 4.0510e-03 | 4.05095e-03 | 4.05095e-03 | 3.59464e-03 | [3.1147e-03, 7.2325e-03] | [3.0964e-03, 7.4599e-03] |

### Hadamard Total

Stable32 applies the Howe 2005 bias correction by default and so does
SigmaTau: the corrected (default) column matches Stable32 at every
matching-α τ; the raw column matches allantools. The two columns differ
by exactly `√(1+a)` per Howe 2005 — for example 0.9975 at WHFM and 0.9225
at FLFM (τ = 512 s, where Stable32 = 3.3578e-03 matches the corrected
3.35777e-03, not the raw 3.09753e-03). The SigmaTau CI column belongs to
the default (corrected) output.

| τ (s) | Stable32 | allantools (raw) | SigmaTau raw | SigmaTau corrected | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|:---|
| 1 | 1.0646e+00 | 1.06459e+00 | 1.06459e+00 | 1.06459e+00 | [1.0531e+00, 1.0765e+00] | [1.0522e+00, 1.0775e+00] |
| 2 | 5.9406e-01 | 5.94056e-01 | 5.94056e-01 | 5.94056e-01 | [5.8764e-01, 6.0069e-01] | [5.8712e-01, 6.0124e-01] |
| 4 | 3.0296e-01 | 3.02955e-01 | 3.02955e-01 | 3.02955e-01 | [2.9968e-01, 3.0634e-01] | [2.9942e-01, 3.0662e-01] |
| 8 | 1.5785e-01 | 1.57846e-01 | 1.57846e-01 | 1.57846e-01 | [1.5614e-01, 1.5961e-01] | [1.5600e-01, 1.5976e-01] |
| 16 | 7.7624e-02 | 7.76241e-02 | 7.76241e-02 | 7.76241e-02 | [7.6785e-02, 7.8492e-02] | [7.6716e-02, 7.8565e-02] |
| 32 | 3.9227e-02 | 3.92275e-02 | 3.92275e-02 | 3.92275e-02 | [3.8615e-02, 3.9870e-02] | [3.8338e-02, 4.0182e-02] |
| 64 | 1.9757e-02 | 1.97566e-02 | 1.97566e-02 | 1.97566e-02 | [1.9389e-02, 2.0146e-02] | [1.9204e-02, 2.0360e-02] |
| 128 | 9.9148e-03 | 9.89001e-03 | 9.89001e-03 | 9.91483e-03 | [9.3157e-03, 1.0647e-02] | [9.3118e-03, 1.0652e-02] |
| 256¹ | 5.1795e-03 | 5.17950e-03 | 5.17950e-03 | 5.19250e-03 | [5.1226e-03, 5.2383e-03] | [4.7578e-03, 5.7731e-03] |
| 512 | 3.3578e-03 | 3.09753e-03 | 3.09753e-03 | 3.35777e-03 | [2.9050e-03, 4.1199e-03] | [2.9047e-03, 4.1207e-03] |
| 1024 | 2.4120e-03 | 2.40600e-03 | 2.40600e-03 | 2.41204e-03 | [2.0370e-03, 3.1179e-03] | [2.0368e-03, 3.1187e-03] |
| 2048 | 2.9382e-03 | 2.93086e-03 | 2.93086e-03 | 2.93822e-03 | [2.3269e-03, 4.6000e-03] | [2.3265e-03, 4.6009e-03] |

## Methodology

See [Validation Methodology](methodology.md) for the shipped fixture strategy.

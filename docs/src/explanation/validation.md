# Validation

SigmaTau is checked against two independent references:

1. **Stable32** (W. Riley), the de facto reference application for
   time-and-frequency stability analysis. Its outputs are printed to five
   significant figures, which sets an agreement floor of about `5e-5`
   relative.
2. **allantools** (A. Wallin), an open-source Python library whose outputs
   are full Float64, so kernel-level agreement can be checked to machine
   precision.

Estimators that neither reference implements (`mhdev`, `htdev`, `mhtotdev`,
`theoh`, `ncornered_hat`, the dynamic and streaming variants) are tested against
inlined textbook kernels and algebraic identities in `test/stab/`.

Where the references disagree with each other, or with SigmaTau's default
output, the difference is traced to a named definitional choice (a bias
correction policy or a noise-identification convention), and the option that
reproduces each reference is stated. The comparisons below use the shipped
fixture `test/fixtures/validation/stable32gen.DAT`, an 8192-point phase
record at τ₀ = 1 s, and `tools/validation/tables.jl` rewrites the per-τ
tables from it:

```julia
p = read_phase("test/fixtures/validation/stable32gen.DAT";
               time_col=0, value_col=1, tau0=1.0, header=10)
```

The Stable32 reference is `stable32out/stable32_data_full.csv` (σ, the
identified noise exponent α, and the 68.3 % confidence bounds); the
allantools reference is `allantools_out/allantools_data_full.csv`.

This page walks through one fixture in detail. The
[Cross-check report](crosscheck.md) summarizes every asserted comparison
(against Stable32, allantools, the NBS 1000-point set and scipy),
regenerated from the validation test group by
`tools/validation/report.jl`.

Default settings throughout, with two exceptions: the total-family tables
additionally show `correct_bias=false` ("raw") columns, and the two Hadamard
tables use `remove_drift=false`, Stable32's convention for the noise type.

## Summary

| Estimator | σ vs Stable32 | σ (raw) vs allantools | CI bounds vs Stable32 |
|---|---|---|---|
| OADEV, MDEV, TDEV, OHDEV | ≤ 0.005 % (fixture precision) | ≤ 1e-13 | ≤ 0.06 % at 9 of 12 τ¹ |
| TOTDEV | ≤ 0.005 % on FM rows with the same α¹; 0.02 to 0.34 % lower on the six PM rows² | ≤ 2.8e-15 | ≤ 0.07 % on FM rows with the same α¹; within 1.1 % on the PM rows² |
| MTOTDEV (raw) | ≤ 0.004 % | ≤ 1.1e-15 | within 1.3 % up to τ = 512 s¹ (≤ 0.02 % at τ ≤ 8 s)³ |
| HTOTDEV | ≤ 0.003 % on FM rows with the same α¹; 10 to 15 % lower on the six PM rows² | ≤ 4.4e-15 | ≤ 0.1 % on FM rows with the same α¹; 10 to 15 % lower on the PM rows² |

¹ The recurring residuals are two noise-identification conventions, at
τ = 256 s and at τ ≥ 1024 s, not numerical disagreements; see below.
² On the rows labelled white PM (τ = 2 s to 16 s) and flicker PM (32 s,
64 s) both estimators read high and the published tables stop at α = 0.
SigmaTau divides out the bias and takes the EDF of the total statistic
itself ([Total estimators on PM noise](total_pm.md)); Stable32 reports the
raw value. `correct_bias=false` matches it to fixture precision.
³ From τ = 16 s Stable32's MTOT EDF follows no published table and at
most τ runs 1.2 to 1.5 times the Howe and Peppler (2001, Table 2) value
SigmaTau uses, so its intervals are usually the narrower ones; the bounds
agree to ≤ 0.02 % at τ ≤ 8 s and to within 1.3 % wherever the two programs
agree on the noise type.

Three per-estimator facts behind those rows:

- **TOTDEV**: both programs apply the Howe/Walter `B = 1 − a·τ/T`
  correction for FLFM/RWFM. SigmaTau's default (`correct_bias=true`)
  matches Stable32's σ to ≤ 0.005 % on the FM rows where the noise types
  agree, including the FLFM row at τ = 512 s. On the PM rows it also divides
  out `B = 1 + a·τ/T`, which Stable32 does not, and reads 0.02 % (τ = 2 s)
  to 0.34 % (64 s) lower. At τ ≥ 1024 s SigmaTau reads
  flicker FM and random-walk FM where Stable32 reads white FM, and the
  correction puts its σ 3.1 %, 10.9 % and 26.5 % higher (¹).
  `correct_bias=false` reproduces allantools' raw SP1065 eq. 25 kernel to
  ≤ 2.8e-15.
- **HTOTDEV**: Stable32 *applies* the Howe 2005 bias correction by
  default (the factors are exactly `1/√(1+a)`), and so does SigmaTau:
  `correct_bias=true` (the default) matches Stable32 to ≤ 0.003 % at
  every matching-α τ on FM. On the PM rows Stable32 reports the raw value,
  which is 12 to 18 % above its own Hadamard deviation; SigmaTau's default
  is 10 to 15 % lower and within 6 % of HDEV (0.03 % at τ = 2 s).
  `correct_bias=false` reproduces allantools' raw kernel to ≤ 4.4e-15.
- **MTOTDEV**: here the two programs deliberately report different
  statistics. Stable32 reports the raw estimator; SigmaTau's raw kernel
  matches it to ≤ 0.004 % (and allantools to ≤ 1.1e-15). The raw
  estimator reads low against MDEV, and SigmaTau's default divides by
  `√B(α)` (Howe and Peppler 2001, Table 2, `B < 1`) to unbias it, which
  places the default output 0.25 % (WHPM) to 18 % (FLFM) and 21 % (RWFM)
  above Stable32 on this fixture, and 41 % at τ = τ₀, where
  the raw estimator is exactly MDEV/√2. Pass `correct_bias=false` to
  reproduce Stable32 and allantools.

### The noise-identification residuals

Because the EDF (and, for the total family, the bias factor) depend on α, a
row where the two programs identify different noise types carries different
CI bounds, and for TOTDEV and HTOTDEV a different σ. Those rows are marked ¹
in the tables below. Two conventions produce them.

At τ = 256 s (every estimator) Stable32 assigns α = +2 where SigmaTau's
default identification assigns α = 0. This is a convention difference in
the lag-1 ACF identifier: SigmaTau quadratically detrends each decimated
subseries before classifying (`remove_drift=true` default); Stable32 does
not. For HTOTDEV the row carries a 0.25 % σ difference.

From τ = 512 s fewer than 30 averaged points remain and both programs turn
to the B1 ratio. They agree on flicker FM at 512 s. At τ ≥ 1024 s Stable32
reads white FM; SigmaTau compares the ratio with B1(N, μ) at the number of
averaged points, the N that B1 is defined for (SP1065 §5.7), and reads
flicker FM at 1024 s and random-walk FM beyond, which is what σ(τ) does on
this record: ADEV is flat from 512 s to 1024 s and rises by a factor 1.47 to
2048 s. The labels at 2048 s and 4096 s are the one measured with six
averaged points (m = 1365), the fewest SigmaTau runs the ratio on. With
eight points or fewer either reading is weak evidence. `identify_noise(x, m_values; remove_drift=false,
b1_length=:record)` reproduces Stable32's α at every τ but 256 s, where the
R(n) confirmation of white PM gives +1, and the test suite checks TOTDEV's σ
with Stable32's α at all 13 τ.

The Hadamard estimators differ in one more way. By default they identify the
noise on the record less its least-squares frequency drift
(`remove_drift=true`), which Stable32 does not do, and on this record that
leaves white FM at τ = 1024 s and 2048 s, as in Stable32, where the record
as given reads flicker FM and random-walk FM. HDEV itself is unchanged; its
default interval at 2048 s is [2.1326e-03, 5.4895e-03]. HTOTDEV follows the
bias factor of the label, so its default values there are Stable32's,
2.4120e-03 and 2.9382e-03. The two Hadamard tables below pass
`remove_drift=false`, which identifies the noise on the record as given,
like the other tables.

### Confidence intervals

Stable32 has computed its confidence intervals with the Greenhall–Riley
2003 EDF algorithm since v1.41 (stated in the paper itself, p. 268), and
SigmaTau implements the same algorithm: wherever the two programs agree
on α, the 68.3 % bounds agree to ≤ 0.06 % for the core estimators. The
total-family differences at α > 0 (footnote ² above) come from the bias
factor and EDF SigmaTau derives where the published tables stop, not from a
different interval method.

## Per-τ comparison

The tables below compare Stable32, allantools, and SigmaTau at each
averaging interval τ. Stable32 columns carry 5 significant figures (the
fixture's precision); SigmaTau and allantools columns are truncated to 6.
CI columns are the 68.3 % `[lower, upper]` bounds, Stable32's
MinSigma/MaxSigma against SigmaTau's `ci_lower`/`ci_upper`. Rows marked ¹
are those where the identified noise type differs from Stable32's.

### Overlapping Allan

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|
| 1 | 1.0097e+00 | 1.00975e+00 | 1.00975e+00 | [9.9898e-01, 1.0209e+00] | [9.9892e-01, 1.0209e+00] |
| 2 | 5.0444e-01 | 5.04443e-01 | 5.04443e-01 | [4.9907e-01, 5.1000e-01] | [4.9903e-01, 5.1003e-01] |
| 4 | 2.5021e-01 | 2.50214e-01 | 2.50214e-01 | [2.4755e-01, 2.5297e-01] | [2.4753e-01, 2.5299e-01] |
| 8 | 1.3063e-01 | 1.30632e-01 | 1.30632e-01 | [1.2924e-01, 1.3207e-01] | [1.2923e-01, 1.3208e-01] |
| 16 | 6.2750e-02 | 6.27503e-02 | 6.27503e-02 | [6.2081e-02, 6.3442e-02] | [6.2076e-02, 6.3447e-02] |
| 32 | 3.1804e-02 | 3.18042e-02 | 3.18042e-02 | [3.1141e-02, 3.2512e-02] | [3.1136e-02, 3.2517e-02] |
| 64 | 1.6085e-02 | 1.60852e-02 | 1.60852e-02 | [1.5668e-02, 1.6537e-02] | [1.5666e-02, 1.6540e-02] |
| 128 | 8.1063e-03 | 8.10630e-03 | 8.10630e-03 | [7.5737e-03, 8.7695e-03] | [7.5735e-03, 8.7699e-03] |
| 256¹ | 4.5004e-03 | 4.50037e-03 | 4.50037e-03 | [4.4513e-03, 4.5511e-03] | [4.0951e-03, 5.0557e-03] |
| 512 | 3.1642e-03 | 3.16419e-03 | 3.16419e-03 | [2.7368e-03, 3.8843e-03] | [2.7366e-03, 3.8850e-03] |
| 1024¹ | 3.1373e-03 | 3.13732e-03 | 3.13732e-03 | [2.6177e-03, 4.1767e-03] | [2.5690e-03, 4.3934e-03] |
| 2048¹ | 4.6131e-03 | 4.61311e-03 | 4.61311e-03 | [3.5914e-03, 7.7529e-03] | [3.4274e-03, 1.0302e-02] |

### Modified Allan

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|
| 1 | 1.0097e+00 | 1.00975e+00 | 1.00975e+00 | [9.9898e-01, 1.0209e+00] | [9.9892e-01, 1.0209e+00] |
| 2 | 3.5650e-01 | 3.56499e-01 | 3.56499e-01 | [3.5252e-01, 3.6062e-01] | [3.5249e-01, 3.6065e-01] |
| 4 | 1.2611e-01 | 1.26111e-01 | 1.26111e-01 | [1.2434e-01, 1.2796e-01] | [1.2433e-01, 1.2797e-01] |
| 8 | 4.6674e-02 | 4.66743e-02 | 4.66743e-02 | [4.5785e-02, 4.7618e-02] | [4.5779e-02, 4.7624e-02] |
| 16 | 1.5818e-02 | 1.58178e-02 | 1.58178e-02 | [1.5400e-02, 1.6272e-02] | [1.5397e-02, 1.6275e-02] |
| 32 | 6.0448e-03 | 6.04484e-03 | 6.04484e-03 | [5.7952e-03, 6.3299e-03] | [5.7936e-03, 6.3318e-03] |
| 64 | 2.5945e-03 | 2.59450e-03 | 2.59450e-03 | [2.4463e-03, 2.7736e-03] | [2.4453e-03, 2.7748e-03] |
| 128 | 1.7347e-03 | 1.73473e-03 | 1.73473e-03 | [1.5955e-03, 1.9182e-03] | [1.5954e-03, 1.9183e-03] |
| 256¹ | 1.8372e-03 | 1.83717e-03 | 1.83717e-03 | [1.6585e-03, 2.0896e-03] | [1.6358e-03, 2.1374e-03] |
| 512 | 2.2115e-03 | 2.21151e-03 | 2.21151e-03 | [1.8801e-03, 2.8157e-03] | [1.8799e-03, 2.8163e-03] |
| 1024¹ | 2.7671e-03 | 2.76714e-03 | 2.76714e-03 | [2.2099e-03, 4.2001e-03] | [2.2039e-03, 4.2393e-03] |
| 2048¹ | 5.2906e-03 | 5.29055e-03 | 5.29055e-03 | [3.8697e-03, 1.3794e-02] | [3.7901e-03, 1.9043e-02] |

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
| 256¹ | 2.7154e-01 | 2.71537e-01 | 2.71537e-01 | [2.4513e-01, 3.0884e-01] | [2.4178e-01, 3.1591e-01] |
| 512 | 6.5373e-01 | 6.53730e-01 | 6.53730e-01 | [5.5577e-01, 8.3233e-01] | [5.5572e-01, 8.3250e-01] |
| 1024¹ | 1.6359e+00 | 1.63595e+00 | 1.63595e+00 | [1.3065e+00, 2.4831e+00] | [1.3030e+00, 2.5063e+00] |
| 2048¹ | 6.2556e+00 | 6.25562e+00 | 6.25562e+00 | [4.5756e+00, 1.6310e+01] | [4.4814e+00, 2.2517e+01] |

### Overlapping Hadamard

SigmaTau columns use `remove_drift=false` (see the noise-identification
residual above).

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|
| 1 | 1.0646e+00 | 1.06459e+00 | 1.06459e+00 | [1.0522e+00, 1.0774e+00] | [1.0522e+00, 1.0775e+00] |
| 2 | 5.3205e-01 | 5.32049e-01 | 5.32049e-01 | [5.2587e-01, 5.3845e-01] | [5.2584e-01, 5.3849e-01] |
| 4 | 2.6192e-01 | 2.61924e-01 | 2.61924e-01 | [2.5888e-01, 2.6507e-01] | [2.5886e-01, 2.6509e-01] |
| 8 | 1.3854e-01 | 1.38538e-01 | 1.38538e-01 | [1.3693e-01, 1.4021e-01] | [1.3692e-01, 1.4022e-01] |
| 16 | 6.6018e-02 | 6.60177e-02 | 6.60177e-02 | [6.5250e-02, 6.6813e-02] | [6.5245e-02, 6.6818e-02] |
| 32 | 3.3500e-02 | 3.35003e-02 | 3.35003e-02 | [3.2736e-02, 3.4321e-02] | [3.2732e-02, 3.4326e-02] |
| 64 | 1.6951e-02 | 1.69515e-02 | 1.69515e-02 | [1.6477e-02, 1.7469e-02] | [1.6474e-02, 1.7473e-02] |
| 128 | 8.4500e-03 | 8.44998e-03 | 8.44998e-03 | [7.8517e-03, 9.2096e-03] | [7.8514e-03, 9.2102e-03] |
| 256¹ | 4.5672e-03 | 4.56718e-03 | 4.56718e-03 | [4.5122e-03, 4.6242e-03] | [4.1228e-03, 5.1951e-03] |
| 512 | 2.9510e-03 | 2.95105e-03 | 2.95105e-03 | [2.5185e-03, 3.7257e-03] | [2.5182e-03, 3.7265e-03] |
| 1024¹ | 2.4608e-03 | 2.46084e-03 | 2.46084e-03 | [2.0108e-03, 3.4682e-03] | [1.9715e-03, 3.6922e-03] |
| 2048¹ | 2.8169e-03 | 2.81691e-03 | 2.81691e-03 | [2.1331e-03, 5.4864e-03] | [2.0598e-03, 7.3497e-03] |

### Total

`SigmaTau corrected` is the default (`correct_bias=true`); it differs
from the raw column from τ = 2 s to 64 s, the rows labelled white or flicker
PM, where TOTVAR reads high by a factor that grows with τ/T, and from
τ = 512 s, where SigmaTau reads flicker FM or
random-walk FM and the Howe/Walter `B = 1 − a·τ/T` correction applies. At
512 s Stable32 reads flicker FM too and matches the corrected value,
confirming Stable32 applies the same correction; beyond it Stable32 reads
white FM and matches the raw one. The SigmaTau CI column belongs to the
default (corrected) output.

| τ (s) | Stable32 | allantools (raw) | SigmaTau raw | SigmaTau corrected | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|:---|
| 1 | 1.0097e+00 | 1.00975e+00 | 1.00975e+00 | 1.00975e+00 | [9.9884e-01, 1.0210e+00] | [9.9892e-01, 1.0209e+00] |
| 2 | 5.0443e-01 | 5.04428e-01 | 5.04428e-01 | 5.04346e-01 | [4.9898e-01, 5.1006e-01] | [4.9894e-01, 5.0993e-01] |
| 4 | 2.5026e-01 | 2.50255e-01 | 2.50255e-01 | 2.50174e-01 | [2.4755e-01, 2.5305e-01] | [2.4749e-01, 2.5295e-01] |
| 8 | 1.3060e-01 | 1.30605e-01 | 1.30605e-01 | 1.30520e-01 | [1.2919e-01, 1.3206e-01] | [1.2912e-01, 1.3197e-01] |
| 16 | 6.2782e-02 | 6.27816e-02 | 6.27816e-02 | 6.27000e-02 | [6.2103e-02, 6.3483e-02] | [6.2024e-02, 6.3399e-02] |
| 32 | 3.1842e-02 | 3.18422e-02 | 3.18422e-02 | 3.17908e-02 | [3.1345e-02, 3.2364e-02] | [3.1123e-02, 3.2503e-02] |
| 64 | 1.6090e-02 | 1.60901e-02 | 1.60901e-02 | 1.60354e-02 | [1.5791e-02, 1.6407e-02] | [1.5618e-02, 1.6488e-02] |
| 128 | 8.1074e-03 | 8.10739e-03 | 8.10739e-03 | 8.10739e-03 | [7.5805e-03, 8.7619e-03] | [7.5801e-03, 8.7624e-03] |
| 256¹ | 4.4756e-03 | 4.47557e-03 | 4.47557e-03 | 4.47557e-03 | [4.4264e-03, 4.5264e-03] | [4.0807e-03, 5.0127e-03] |
| 512 | 3.1028e-03 | 3.05583e-03 | 3.05583e-03 | 3.10283e-03 | [2.6974e-03, 3.7709e-03] | [2.6974e-03, 3.7709e-03] |
| 1024¹ | 2.8694e-03 | 2.86942e-03 | 2.86942e-03 | 2.95976e-03 | [2.4276e-03, 3.6939e-03] | [2.4567e-03, 3.9928e-03] |
| 2048¹ | 3.4123e-03 | 3.41229e-03 | 3.41229e-03 | 3.78565e-03 | [2.7445e-03, 5.0548e-03] | [2.9051e-03, 6.8227e-03] |
| 4096¹ | 4.0832e-03 | 4.08324e-03 | 4.08324e-03 | 5.16513e-03 | [3.1056e-03, 7.7435e-03] | [3.7325e-03, 1.5796e-02] |

### Modified Total

Stable32 and allantools both report the raw estimator here; the SigmaTau
raw column matches both (≤ 0.004 % / ≤ 1.1e-15). The corrected column is
SigmaTau's default, the raw value divided by `√B(α)` from Howe and Peppler
2001, Table 2, which on this fixture is 0.25 % (WHPM), 8.4 % (FLPM), 13.9 %
(WHFM), 18.1 % (FLFM) or 21.4 % (RWFM) above Stable32, and √2 above it at
τ = 1 s. It is
the column to compare with the Modified Allan table above: from 128 s to
512 s, where the raw value is 13 to 17 % below MDEV, the corrected one is
within 5 %. The SigmaTau CI column belongs to the raw output for
like-for-like comparison with Stable32's bounds.

| τ (s) | Stable32 | allantools (raw) | SigmaTau raw | SigmaTau corrected | Stable32 CI | SigmaTau CI (raw) |
|---:|:---|:---|:---|:---|:---|:---|
| 1 | 7.1400e-01 | 7.13999e-01 | 7.13999e-01 | 1.00975e+00 | [7.0639e-01, 7.2186e-01] | [7.0634e-01, 7.2191e-01] |
| 2 | 3.5625e-01 | 3.56253e-01 | 3.56253e-01 | 3.57147e-01 | [3.5227e-01, 3.6037e-01] | [3.5225e-01, 3.6040e-01] |
| 4 | 1.2570e-01 | 1.25704e-01 | 1.25704e-01 | 1.26020e-01 | [1.2394e-01, 1.2755e-01] | [1.2393e-01, 1.2756e-01] |
| 8 | 4.6303e-02 | 4.63031e-02 | 4.63031e-02 | 4.64193e-02 | [4.5421e-02, 4.7239e-02] | [4.5415e-02, 4.7245e-02] |
| 16 | 1.5819e-02 | 1.58191e-02 | 1.58191e-02 | 1.58588e-02 | [1.5495e-02, 1.6165e-02] | [1.5461e-02, 1.6203e-02] |
| 32 | 5.8856e-03 | 5.88559e-03 | 5.88559e-03 | 6.38007e-03 | [5.6979e-03, 6.0933e-03] | [5.6568e-03, 6.1446e-03] |
| 64 | 2.4483e-03 | 2.44833e-03 | 2.44833e-03 | 2.65402e-03 | [2.3398e-03, 2.5737e-03] | [2.3167e-03, 2.6053e-03] |
| 128 | 1.5120e-03 | 1.51199e-03 | 1.51199e-03 | 1.72196e-03 | [1.4073e-03, 1.6442e-03] | [1.3961e-03, 1.6625e-03] |
| 256¹ | 1.5468e-03 | 1.54684e-03 | 1.54684e-03 | 1.76164e-03 | [1.4288e-03, 1.7001e-03] | [1.3855e-03, 1.7819e-03] |
| 512 | 1.8467e-03 | 1.84668e-03 | 1.84668e-03 | 2.18088e-03 | [1.5780e-03, 2.3251e-03] | [1.5789e-03, 2.3221e-03] |
| 1024¹ | 2.2545e-03 | 2.25454e-03 | 2.25454e-03 | 2.66255e-03 | [1.8657e-03, 3.0666e-03] | [1.8173e-03, 3.3142e-03] |
| 2048¹ | 4.0510e-03 | 4.05095e-03 | 4.05095e-03 | 4.91611e-03 | [3.1147e-03, 7.2325e-03] | [2.9780e-03, 9.9754e-03] |

### Hadamard Total

Stable32 applies the Howe 2005 bias correction by default and so does
SigmaTau: the corrected (default) column matches Stable32 at every
matching-α τ on FM; the raw column matches allantools. From τ = 2 s to 64 s,
the PM rows, Stable32 reports the raw value and the corrected column is the
one to compare with the Overlapping Hadamard table above. On FM the two columns differ
by exactly `√(1+a)` per Howe 2005, for example 0.9975 at WHFM and 0.9225
at FLFM (τ = 512 s, where Stable32 = 3.3578e-03 matches the corrected
3.35777e-03, not the raw 3.09753e-03). SigmaTau columns use
`remove_drift=false`, and the CI column belongs to the bias-corrected output.

| τ (s) | Stable32 | allantools (raw) | SigmaTau raw | SigmaTau corrected | Stable32 CI | SigmaTau CI |
|---:|:---|:---|:---|:---|:---|:---|
| 1 | 1.0646e+00 | 1.06459e+00 | 1.06459e+00 | 1.06459e+00 | [1.0531e+00, 1.0765e+00] | [1.0522e+00, 1.0775e+00] |
| 2 | 5.9406e-01 | 5.94056e-01 | 5.94056e-01 | 5.32224e-01 | [5.8764e-01, 6.0069e-01] | [5.2707e-01, 5.3753e-01] |
| 4 | 3.0296e-01 | 3.02955e-01 | 3.02955e-01 | 2.64207e-01 | [2.9968e-01, 3.0634e-01] | [2.6178e-01, 2.6670e-01] |
| 8 | 1.5785e-01 | 1.57846e-01 | 1.57846e-01 | 1.35887e-01 | [1.5614e-01, 1.5961e-01] | [1.3467e-01, 1.3714e-01] |
| 16 | 7.7624e-02 | 7.76241e-02 | 7.76241e-02 | 6.64022e-02 | [7.6785e-02, 7.8492e-02] | [6.5814e-02, 6.7006e-02] |
| 32 | 3.9227e-02 | 3.92275e-02 | 3.92275e-02 | 3.53632e-02 | [3.8615e-02, 3.9870e-02] | [3.4756e-02, 3.6003e-02] |
| 64 | 1.9757e-02 | 1.97566e-02 | 1.97566e-02 | 1.76829e-02 | [1.9389e-02, 2.0146e-02] | [1.7312e-02, 1.8078e-02] |
| 128 | 9.9148e-03 | 9.89001e-03 | 9.89001e-03 | 9.91483e-03 | [9.3157e-03, 1.0647e-02] | [9.3118e-03, 1.0652e-02] |
| 256¹ | 5.1795e-03 | 5.17950e-03 | 5.17950e-03 | 5.19250e-03 | [5.1226e-03, 5.2383e-03] | [4.7578e-03, 5.7731e-03] |
| 512 | 3.3578e-03 | 3.09753e-03 | 3.09753e-03 | 3.35777e-03 | [2.9050e-03, 4.1199e-03] | [2.9047e-03, 4.1207e-03] |
| 1024¹ | 2.4120e-03 | 2.40600e-03 | 2.40600e-03 | 2.60814e-03 | [2.0370e-03, 3.1179e-03] | [2.1422e-03, 3.6207e-03] |
| 2048¹ | 2.9382e-03 | 2.93086e-03 | 2.93086e-03 | 3.33786e-03 | [2.3269e-03, 4.6000e-03] | [2.5337e-03, 6.3995e-03] |

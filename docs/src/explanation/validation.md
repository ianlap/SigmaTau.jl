# Validation

The results of this package rest on five kinds of evidence, described
below. The second half of the page follows one record through a comparison
with Stable32 and allantools at each averaging time. The
[cross-check report](crosscheck.md) lists every comparison the test suite
asserts, with its tolerance and the largest difference found.

## Evidence

| Evidence | What it establishes | Where |
|:---|:---|:---|
| Other programs | The deviations, EDF, χ² limits, lag-1 estimates, spectra and drift fits agree with allantools, AllanDeviations.jl, scipy and Stable32 | `test/validation/` |
| Published values | The deviations of the 1000-point test suite of [SP1065 §12.4](@cite riley-2008-sp1065) are reproduced | `test/validation/` |
| Definitions | Each estimator equals a direct evaluation of its defining sum | `test/stab/` |
| Exact moments | The bias factors and EDF models derived in this package equal the exact mean and variance of the estimator on Gaussian noise | `test/stab/`, `tools/` |
| Simulation | How often the noise identification is right, and how often an interval holds the true value | the Background pages |

### Other programs

[Related software](related.md) describes the programs. allantools and
AllanDeviations.jl compute in double precision, and the suite requires the
deviations to agree with theirs within ``10^{-9}`` (``10^{-8}`` for the
modified Allan and time deviations). Stable32 prints five significant
digits, and the suite requires agreement with its result files within
``10^{-4}`` for the deviations and ``10^{-3}`` for the confidence limits.
The comparisons run on seven power-law records of 4096 points, from white
PM to random-run FM, on the test record that Stable32 distributes and on a
flicker PM record from its noise generator.

Where two programs follow different conventions the suite compares like
with like. Against allantools, which returns the total deviations
uncorrected, it calls [`mtotdev`](@ref) with `correct_bias=false`. The
report states each convention. The differences that remain, such as the
EDF each program gives the total deviations at short averaging times, are
listed in the report and not asserted.

A continuous-integration job installs allantools 2024.6 and
AllanDeviations.jl 0.3.0, regenerates their tables and repeats the
comparison on every pull request, so that the stored tables cannot drift
from the programs.

### Published values

[SP1065 §12.4](@cite riley-2008-sp1065) defines a 1000-point frequency
record by a random number generator (its Eq. (73)) and lists its
deviations in Table 31. The suite generates the record and reproduces the
overlapping Allan, modified Allan, time, overlapping Hadamard, total and
Hadamard total deviations of the table within ``3 × 10^{-7}``, the
resolution of the seven digits printed.

### Definitions

`test/stab/reference_kernels.jl` holds an implementation of each variance
written directly from its definition, sharing no code with the package.
The optimized kernels agree with it within ``10^{-12}`` on a 4096-point
record of white PM and random-walk FM. For the estimators that no other
program computes ([`mhdev`](@ref), [`htdev`](@ref), [`mhtotdev`](@ref))
this is the primary test.

### Exact moments

A squared deviation is a quadratic form in the phase samples, so its mean
and its variance on Gaussian noise of a given spectrum follow exactly,
without simulation. [Total estimators on PM noise](total_pm.md) and
[MHTOTDEV bias and EDF](mhtotdev.md) derive bias factors and EDF models
this way, with tools that do not call the package. The tests work in the
other direction: they extract the quadratic form from each kernel and hold
the bias factors to its exact mean and the EDF models to its exact EDF.

### Simulation

None of the four says how a default call behaves on noise of unknown type.
The Background pages therefore give, in simulation, how often the
noise identification is right
([Noise identification](noise_identification.md)), how much bias a default
call leaves ([Bias correction](bias.md)) and how often the interval holds
the true value ([Confidence intervals](confidence.md)). Each statement
gives the record length and the number of records. These runs are not part
of the test suite.

## One record in detail

`test/fixtures/validation/stable32gen.DAT` is a phase record of 8192
points at 1 s, in Stable32's file format. White PM dominates up to about
250 s, flicker FM and random-walk FM beyond. The fixtures hold the tables
that Stable32 and allantools computed from it
(`stable32out/stable32_data_full.csv`,
`allantools_out/allantools_data_full.csv`).

```@example validation
using SigmaTau
file = joinpath(pkgdir(SigmaTau), "test", "fixtures", "validation",
                "stable32gen.DAT")
p = read_phase(file)
adev(p, [1, 16, 256, 2048])
```

### Noise type

The EDF, and for the total deviations the bias factor, depend on the noise
type, so the comparison starts there. Stable32 prints the same noise type
for every overlapping estimator at a given averaging time:

| ``τ`` (s) | 1 to 16 | 32, 64 | 128 | 256 | 512 | 1024 | 2048, 4096 |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Stable32 | white PM | flicker PM | white FM | white PM | flicker FM | white FM | white FM |
| This package | white PM | flicker PM | white FM | white FM | flicker FM | flicker FM | random-walk FM |
| Hadamard family, default | white PM | flicker PM | white FM | white FM | flicker FM | flicker FM | flicker FM |

The programs differ at three places.

- **256 s.** With 32 averaged points this is the longest averaging time at
  which both use the lag-1 autocorrelation. Stable32 reads white PM and
  this package white FM. [`identify_noise`](@ref) with
  `remove_drift=false`, which leaves the averaged series as it is, reads
  flicker PM. The Allan deviation falls as ``τ^{-0.85}`` from 128 s to
  256 s and as ``τ^{-0.5}`` from 256 s to 512 s.
- **1024 s and beyond.** With fewer than 30 averaged points both programs
  turn to the ``B_1`` ratio. This package compares it with ``B_1(N, μ)``
  at the number of averaged points and reads flicker FM at 1024 s and
  random-walk FM beyond; the Allan deviation is flat from 512 s to 1024 s
  and rises by a factor 1.47 to 2048 s. `identify_noise` with
  `remove_drift=false` and `b1_length=:record` returns white FM there, as
  Stable32 does. The labels at 2048 s and 4096 s come from six averaged
  points, the fewest the ratio is computed on
  ([Noise identification](noise_identification.md)).
- **The Hadamard family.** By default [`hdev`](@ref),
  [`htotdev`](@ref) and their modified forms identify the noise on the
  record less its frequency drift, and there without the ``B_1`` ratio.
  From 512 s on, where fewer than 30 averaged points remain, their label
  is the one found at 273 s, the longest averaging time that has 30:
  flicker FM, to the end of their tables at 2048 s. With
  `remove_drift=false` they use the labels of the second row. The
  Hadamard deviation is the same either way. Its default interval at
  2048 s is ``[2.08, 6.64] × 10^{-3}``, and Stable32's, on white FM,
  ``[2.13, 5.49] × 10^{-3}``.

### Deviations

Without bias correction the deviations agree with those of allantools
within ``10^{-13}``, for the seven estimators tabulated below.

With Stable32's conventions stated, which are its noise type at each
averaging time (`noise=`) and its bias factors (`correct_bias=:stable32`),
the deviations agree with Stable32 within ``5 × 10^{-5}``, one unit of the
last digit it prints, at every averaging time and for the time total
deviation as well.

A default call differs from Stable32 in these rows:

| Estimator | Rows | Default call against Stable32 | Cause |
|:---|:---|:---|:---|
| Total | 2 s to 64 s | 0.02 to 0.34 % lower | the bias factor on PM noise ([Total estimators on PM noise](total_pm.md)) |
| Total | 1024 s, 2048 s, 4096 s | 3.1 %, 10.9 %, 26.5 % higher | the noise type, and with it the factor ``1 - a\,τ/T`` |
| Modified total, time total | every row | 41 % higher at ``τ_0``, 0.1 to 0.25 % on white PM, 8.4 % on flicker PM, 13.9 %, 18.1 % and 21.4 % on white, flicker and random-walk FM | Stable32 prints these two uncorrected |
| Hadamard total | 2 s to 64 s | 10 to 15 % lower | the bias factor on PM noise |
| Hadamard total | 256 s | 0.25 % higher | the noise type |
| Hadamard total | 1024 s, 2048 s | 8.1 % higher | the noise type |

At ``τ_0`` the modified total variance is half the modified Allan variance
on any record, which is the 41 % ([Bias correction](bias.md)). With
`remove_drift=false`, as in the table further down, the Hadamard total
deviation takes random-walk FM at 2048 s and is 13.6 % above Stable32
there.

### Confidence limits

Both programs compute the EDF of the Allan and Hadamard families by the
algorithm of
[Greenhall and Riley (2003)](@cite greenhall-2003-edf-stability), which
Stable32 has used since version 1.41 according to that paper, and both
take χ² limits. With Stable32's noise types stated, the 68.3 % limits
differ by at most:

| Estimator | Largest difference |
|:---|:---|
| Overlapping Allan, modified Allan, time, overlapping Hadamard | 0.06 % |
| Total | 0.06 %, except 0.5 to 1.2 % at 32 s, 64 s and 256 s, which Stable32 labels flicker PM and, at 256 s, white PM |
| Hadamard total | 0.26 % |
| Modified total, time total | 0.01 % up to 8 s, 0.1 to 1.2 % from 16 s to 512 s, 5.0 % at 1024 s and 7.4 % at 2048 s |

From ``16τ_0`` the EDF that Stable32 gives the modified total deviation
follows no published table
([cross-check report](crosscheck.md)). In a default call the limits also
move wherever the noise type differs: at 256 s, where this package reads
white FM, its interval for the Allan deviation is
``[4.10, 5.06] × 10^{-3}`` against Stable32's ``[4.45, 4.55] × 10^{-3}``.

### Thêo1 with the bias removed

Stable32's ThêoH table for this record is in the fixtures. Up to 1024 s
it gives the Allan deviation, as [`theoh`](@ref) does. Beyond, it gives
the bias-removed Thêo1 at 1536 s and 3072 s. [`theo1`](@ref) with
`correct_bias = :stable32` agrees with it to the five digits printed, and
`theo1` by default is 1.8 % lower:

| ``τ`` (s) | Stable32 | `theo1`, `:stable32` | `theo1` | `theo1`, `:published` | `theo1`, uncorrected |
|---:|:---|:---|:---|:---|:---|
| 1536 | 4.1340e-03 | 4.1340e-03 | 4.0596e-03 | 1.6896e-03 | 2.8344e-03 |
| 3072 | 5.1100e-03 | 5.1100e-03 | 5.0181e-03 | 2.0885e-03 | 3.5035e-03 |

The 1.8 % come from the samples used. Both take the ratio of the Allan
variance to Thêo1 over averaging times up to half the record
([SP1065 Eq. (33)](@cite riley-2008-sp1065);
[IEEE Std 1139-2022, Eq. (A.37)](@cite ieee1139-2022-definitions)). The
default takes it from all 8192 samples, and Stable32 from at most 3000.
For this record `correct_bias = :stable32` takes every third sample and
adds the term of
[Taylor and Howe (2010), Eq. (4)](@cite taylor-2010-fast-theobr) for it,
0.0011 in the ratio, which they derived for samples averaged in blocks.
That Stable32 takes every third sample is inferred from this agreement:
averages of three give values 6.5 % higher. On the 1001-point record of
the [cross-check report](crosscheck.md), which needs no such step, the
agreement holds at each of 181 averaging factors, and the default is
within 0.1 % of it.

`correct_bias = :published` is lower by a factor 2.4. The Allan deviation
is ``3.14 × 10^{-3}`` at 1024 s and ``4.61 × 10^{-3}`` at 2048 s. The
default values continue it, and those of `:published` do not: its 68.3 %
interval at 1536 s, ``[1.35, 2.59] × 10^{-3}``, lies below both.

The cause is the range over which the ratio is averaged. `:published`
follows [Howe (2006), Eq. (6)](@cite howe-2006-theoh) and averages it up
to a tenth of the record, 819 s here, where PM noise dominates and the
ratio is 0.2 to 0.7. Its mean is 0.355. Up to half the record FM noise
dominates and the ratio reaches 3.4, and the mean over that range is 2.05.
Either mean is then applied at every averaging time. The correction
assumes one noise type over the whole range, and this record does not have
one ([Bias correction](bias.md) gives the size of the error in simulation,
for both ranges).

## Comparison at each averaging time

[`tools/validation/tables.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/tools/validation/tables.jl)
writes the tables below from the record and the two reference tables.
Stable32's columns carry the five digits it prints. The intervals are the
68.3 % limits. Rows marked ¹ are those where the noise type of this
package differs from Stable32's.

The tables show default calls, with two exceptions. The two Hadamard
tables use `remove_drift=false`, so that their noise types are those of
the other tables. The total tables add the uncorrected deviation, which is
what allantools returns.

### Overlapping Allan

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 interval | SigmaTau interval |
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

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 interval | SigmaTau interval |
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

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 interval | SigmaTau interval |
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

Computed with `remove_drift=false`.

| τ (s) | Stable32 | allantools | SigmaTau | Stable32 interval | SigmaTau interval |
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

The interval is that of the default call. At 512 s, where both programs
read flicker FM, Stable32's value is the corrected one.

| τ (s) | Stable32 | allantools | Uncorrected | Default | Stable32 interval | SigmaTau interval |
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

Stable32 prints this deviation uncorrected, so the interval shown is that of
the uncorrected call.

| τ (s) | Stable32 | allantools | Uncorrected | Default | Stable32 interval | SigmaTau interval |
|---:|:---|:---|:---|:---|:---|:---|
| 1 | 7.1400e-01 | 7.13999e-01 | 7.13999e-01 | 1.00975e+00 | [7.0639e-01, 7.2186e-01] | [7.0634e-01, 7.2191e-01] |
| 2 | 3.5625e-01 | 3.56253e-01 | 3.56253e-01 | 3.56628e-01 | [3.5227e-01, 3.6037e-01] | [3.5225e-01, 3.6040e-01] |
| 4 | 1.2570e-01 | 1.25704e-01 | 1.25704e-01 | 1.25956e-01 | [1.2394e-01, 1.2755e-01] | [1.2393e-01, 1.2756e-01] |
| 8 | 4.6303e-02 | 4.63031e-02 | 4.63031e-02 | 4.64076e-02 | [4.5421e-02, 4.7239e-02] | [4.5415e-02, 4.7245e-02] |
| 16 | 1.5819e-02 | 1.58191e-02 | 1.58191e-02 | 1.58588e-02 | [1.5495e-02, 1.6165e-02] | [1.5461e-02, 1.6203e-02] |
| 32 | 5.8856e-03 | 5.88559e-03 | 5.88559e-03 | 6.38007e-03 | [5.6979e-03, 6.0933e-03] | [5.6568e-03, 6.1446e-03] |
| 64 | 2.4483e-03 | 2.44833e-03 | 2.44833e-03 | 2.65402e-03 | [2.3398e-03, 2.5737e-03] | [2.3167e-03, 2.6053e-03] |
| 128 | 1.5120e-03 | 1.51199e-03 | 1.51199e-03 | 1.72196e-03 | [1.4073e-03, 1.6442e-03] | [1.3961e-03, 1.6625e-03] |
| 256¹ | 1.5468e-03 | 1.54684e-03 | 1.54684e-03 | 1.76164e-03 | [1.4288e-03, 1.7001e-03] | [1.3855e-03, 1.7819e-03] |
| 512 | 1.8467e-03 | 1.84668e-03 | 1.84668e-03 | 2.18088e-03 | [1.5780e-03, 2.3251e-03] | [1.5789e-03, 2.3221e-03] |
| 1024¹ | 2.2545e-03 | 2.25454e-03 | 2.25454e-03 | 2.66255e-03 | [1.8657e-03, 3.0666e-03] | [1.8173e-03, 3.3142e-03] |
| 2048¹ | 4.0510e-03 | 4.05095e-03 | 4.05095e-03 | 4.91611e-03 | [3.1147e-03, 7.2325e-03] | [2.9780e-03, 9.9754e-03] |

### Hadamard Total

The interval is that of the default call with `remove_drift=false`. On the
FM rows Stable32's value is the corrected one.

| τ (s) | Stable32 | allantools | Uncorrected | Default | Stable32 interval | SigmaTau interval |
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

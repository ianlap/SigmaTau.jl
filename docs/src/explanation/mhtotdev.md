# MHTOTDEV bias and EDF

MHTOTDEV (modified Hadamard total deviation) is original to SigmaTau, and no
published reference gives its bias factor or its equivalent degrees of
freedom (EDF). Every other total estimator takes a published table: TOTDEV
from SP1065 [riley-2008-sp1065](@cite), MTOTDEV and HTOTDEV from Howe and
Peppler [howe-2001-total-estimators](@cite). For MHTOTDEV both are computed
here, exactly, from the definition of the statistic. The computation is
[`tools/mhtot_models.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/tools/mhtot_models.jl);
it uses no random numbers and runs in about five minutes.

## The statistic

For each window of 4m phase samples, MHTOTVAR removes the window's
least-squares quadratic, extends it by mirror reflection at both ends to 12m
samples, and averages the squared m-averaged third difference over the
8m + 1 positions of the extension. The record's value is the mean over all
N − 4m + 1 windows.

The quadratic removal matches the third difference, which annihilates phase
polynomials up to degree 2. It makes the statistic immune to a linear
frequency drift on its own and lets it converge down to random-run FM
(α = −4), like HDEV and MHDEV. The reference it is compared with is MHVAR,
the square of [`mhdev`](@ref), at the same averaging factor.

## Method

Every step above is linear in the phase, so a window's statistic is a
quadratic form w′Qw, and the record's is x′Q̄x with Q̄ the mean of the
window forms. For Gaussian noise with phase covariance Σ,

```math
E[\mathrm{MHTOTVAR}] = \operatorname{tr}(\bar Q\,\Sigma), \qquad
\mathrm{edf} = \frac{2\,E^2}{\mathrm{Var}}
             = \frac{\operatorname{tr}(\bar Q\,\Sigma)^2}{\operatorname{tr}\big((\bar Q\,\Sigma)^2\big)} ,
```

the same definition of the EDF as Greenhall and Riley use for the unextended
estimators [greenhall-2003-edf-stability](@cite). Q annihilates quadratics,
so any valid covariance of the samples gives the same value for every noise
down to α = −4. The tool builds Σ = FF′ with F the response to white noise:
the identity for white PM, Kasdin's half-order filter for flicker PM, and
one cumulative sum for each step of 2 down in α.

The tool builds its matrices from the definition and never calls the package
kernel. The test suite closes the loop the other way: it extracts F′Q̄F from
the kernel itself and holds the shipped models to its exact moments.

## Bias

B = E[MHTOTVAR] / E[MHVAR] depends on the averaging factor m and the noise
type only, not on the record length or on τ/T.

| Noise | m = 2 | 3 | 4 | 8 | 16 | 64 | 256 | from m = 16 |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| white PM | 0.9000 | 0.8584 | 0.9013 | 0.9011 | 0.9007 | 0.9002 | 0.9000 | 0.9000 |
| flicker PM | 0.8475 | 0.8233 | 0.8340 | 0.8269 | 0.8239 | 0.8221 | 0.8218 | 0.8218 |
| white FM | 0.8281 | 0.8059 | 0.8064 | 0.7972 | 0.7938 | 0.7921 | 0.7917 | 0.7917 |
| flicker FM | 0.8052 | 0.7948 | 0.7927 | 0.7869 | 0.7846 | 0.7831 | 0.7828 | 0.7828 |
| random-walk FM | 0.8009 | 0.7964 | 0.7946 | 0.7903 | 0.7883 | 0.7870 | 0.7867 | 0.7867 |
| flicker-walk FM | 0.8025 | 0.8063 | 0.8057 | 0.8028 | 0.8011 | 0.7999 | 0.7996 | 0.7996 |
| random-run FM | 0.8145 | 0.8239 | 0.8243 | 0.8223 | 0.8209 | 0.8199 | 0.8196 | 0.8196 |

MHTOTDEV reads low on every noise: removing a quadratic from each window
absorbs part of that window's variance before the reflection.
[`bias_correction`](@ref) divides by the exact value at each m from 2 to 15
and by the large-m constant in the last column from m = 16, where it is
within 0.4 %. White PM dips at odd m, to 0.858 at m = 3. The flicker rows are
the mean of the two flicker models of
[Total estimators on PM noise](total_pm.md), which are within 1.3 % of it at
m = 2 and 0.5 % from m = 3. At m = 1 the detrended four-point window leaves
a single degree of freedom and B = 124/225 on every noise.

## EDF

The exact EDF, computed for records of 512, 1024 and 2048 points, is a
function of s = τ/T to two digits. The model is

```math
\mathrm{edf} = \frac{b}{s} - c + d\,s + e\,s^3, \qquad s = \tau/T,\quad T = (N-1)\,\tau_0 ,
```

fitted by least squares on the relative residual over 16 ≤ m ≤ N/4, the
three record lengths together:

| Noise | b | c | d | e | largest residual |
|:---|---:|---:|---:|---:|---:|
| white PM | 1.5949 | 6.535 | 20.18 | −52.1 | 4.5 % |
| flicker PM | 1.0656 | 4.223 | 10.19 | 7.8 | 2.6 % |
| white FM | 0.9420 | 3.506 | 6.50 | 23.3 | 2.0 % |
| flicker FM | 0.8938 | 3.124 | 3.94 | 31.2 | 2.5 % |
| random-walk FM | 0.8596 | 2.676 | 0.43 | 48.8 | 2.8 % |
| flicker-walk FM | 0.8129 | 2.115 | −4.27 | 83.8 | 2.9 % |
| random-run FM | 0.6508 | 1.991 | −2.63 | 78.3 | 1.0 % |

The first two terms are the form Howe and Peppler use for the modified total
variance [howe-2001-total-estimators](@cite). They are not enough here: the
exact EDF turns upward in the last octave, where the reflection of the few
remaining windows does the most work, and two terms that fit the range below
T/16 give about half the exact EDF at τ = T/4. The terms in s and s³ carry
that rise.

At the reach limit τ = T/4 a single window is left. MHVAR has one degree of
freedom there; MHTOTVAR has 4.0 on white PM, 2.7 on flicker PM, 2.2 on white
FM, 1.9 on flicker FM, 1.6 on random-walk FM, 1.4 on flicker-walk FM and 1.2
on random-run FM.

!!! note "Below m = 16"
    The model is fitted from m = 16 and is used down to m = 2. On the FM
    noises it is within 6 % of the exact EDF from m = 3 and reads 5 to 13 %
    low, the safe side, at m = 2. On flicker PM it reads 4 to 10 % high at
    m = 2 and 3, depending on the flicker model, and is within 4 % from
    m = 4. On white PM it would read high, by 5 % at m = 8, 15 %
    at m = 4, 35 % at m = 3 and 54 % at m = 2, so there the model is
    multiplied by the exact ratio at each m from 2 to 15 (0.65 at m = 2,
    0.75 at m = 3, 0.88 at m = 4, 0.97 at m = 8).

## Degrees of freedom gained over MHDEV

The exact EDF of MHTOTVAR over that of MHVAR at the same τ:

| Noise | m = 16 | T/16 | T/8 | T/4 |
|:---|---:|---:|---:|---:|
| white PM | 1.38 | 1.41 | 1.50 | 4.0 |
| flicker PM | 1.19 | 1.21 | 1.27 | 2.7 |
| white FM | 1.10 | 1.12 | 1.17 | 2.2 |
| flicker FM | 1.05 | 1.06 | 1.10 | 1.9 |
| random-walk FM | 1.01 | 1.02 | 1.06 | 1.6 |
| flicker-walk FM | 0.99 | 1.00 | 1.03 | 1.4 |
| random-run FM | 0.98 | 0.98 | 1.00 | 1.2 |

The extension buys the most on the phase noises and in the last octave. On
the red FM noises below T/8, MHTOTDEV's advantage over MHDEV is drift
immunity inside each window and reach, not precision.

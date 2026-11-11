# MHTOTDEV bias and EDF

The modified Hadamard total deviation, [`mhtotdev`](@ref), is the total
form of the modified Hadamard deviation. It is constructed in this
package. [Vernotte and Howe (2000)](@cite vernotte-2000-total-structure)
outline a total form of the same third-difference variance, with a linear
trend removed from each subsequence and a quadratic from the whole record,
and give no bias or EDF for it. No published table covers the estimator
defined here, so its bias factor and its EDF model are derived on this
page, exactly, from the definition.

The computation is
[`tools/mhtot_models.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/tools/mhtot_models.jl).
It uses no random numbers and runs in about five minutes.

## The statistic

For each window of ``4m`` phase samples the estimator

1. removes the least-squares quadratic of the window,
2. extends the window at both ends by even reflection, to ``12m`` samples,
3. averages the squared third difference of ``m``-sample phase means over
   the ``8m + 1`` positions of the extension.

The variance is the mean over the ``N - 4m + 1`` windows, divided by
``6τ^2``.

The third difference annihilates a phase polynomial of degree 2, and the
quadratic removed in step 1 is that polynomial. The estimate is therefore
insensitive to frequency offset and to linear frequency drift, and it
converges for every noise type down to random-run FM, like the Hadamard
and modified Hadamard deviations. Its parent, the variance it is compared
with, is the modified Hadamard variance, the square of [`mhdev`](@ref), at
the same averaging factor.

## Method

Each step is linear in the phase. The statistic of one window is a
quadratic form in its samples, and that of the record is
``x^\mathsf{T}\bar{Q}x``, with ``\bar{Q}`` the mean of the window forms.
For Gaussian noise whose phase has covariance ``Σ``,

```math
E = \operatorname{tr}(\bar{Q}Σ), \qquad
ν = \frac{2E^2}{\mathrm{Var}}
  = \frac{\operatorname{tr}(\bar{Q}Σ)^2}{\operatorname{tr}\big((\bar{Q}Σ)^2\big)} ,
```

with the EDF as defined by
[Greenhall and Riley (2003), Eq. (1)](@cite greenhall-2003-edf-stability).
Because ``\bar{Q}`` annihilates quadratics, any valid covariance of the
samples gives the same value, for every noise type down to random-run FM.
The tool takes ``Σ = F F^\mathsf{T}``, where ``F`` is the response to white
noise: the identity for white PM, the half-order filter of
[Kasdin and Walter (1992)](@cite kasdin-1992-discrete-simulation) for
flicker PM, and one cumulative sum for each step of 2 down in ``α``.

The tool builds its matrices from the definition and does not call the
package. The test suite works in the other direction: it extracts the form
from the package kernel and holds the models to its exact moments. Applied
to the modified Hadamard variance itself, the same computation reproduces
the EDF of the Greenhall and Riley algorithm within 0.3 % on the four
noise types with even ``α``, and within 0.7 % on the flicker types, where
the algorithm assumes the other flicker model.

## Bias

The bias factor ``B = E[\text{MHTOTVAR}]/E[\text{MHVAR}]`` depends on the
averaging factor and on the noise type. It does not depend on the record
length or on ``τ/T``.

| Noise | ``m = 2`` | 3 | 4 | 8 | 16 | 64 | 256 | used from ``m = 16`` |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| White PM | 0.9000 | 0.8584 | 0.9013 | 0.9011 | 0.9007 | 0.9002 | 0.9000 | 0.9000 |
| Flicker PM | 0.8475 | 0.8233 | 0.8340 | 0.8269 | 0.8239 | 0.8221 | 0.8218 | 0.8218 |
| White FM | 0.8281 | 0.8059 | 0.8064 | 0.7972 | 0.7938 | 0.7921 | 0.7917 | 0.7917 |
| Flicker FM | 0.8052 | 0.7948 | 0.7927 | 0.7869 | 0.7846 | 0.7831 | 0.7828 | 0.7828 |
| Random-walk FM | 0.8009 | 0.7964 | 0.7946 | 0.7903 | 0.7883 | 0.7870 | 0.7867 | 0.7867 |
| Flicker-walk FM | 0.8025 | 0.8063 | 0.8057 | 0.8028 | 0.8011 | 0.7999 | 0.7996 | 0.7996 |
| Random-run FM | 0.8145 | 0.8239 | 0.8243 | 0.8223 | 0.8209 | 0.8199 | 0.8196 | 0.8196 |

The estimator reads low on every noise type, by 10 to 22 % in variance:
the quadratic removed from each window takes part of the window's variance
with it. [`bias_correction`](@ref) returns the exact value at each ``m``
from 2 to 15 and the constant of the last column from ``m = 16``, where
the two differ by less than 0.3 %. White PM dips at odd ``m``, to 0.858
at ``m = 3``. The flicker rows are the mean of the two flicker models of
[Total estimators on PM noise](total_pm.md), which lie within 1.3 % of it
at ``m = 2``.

At ``m = 1`` the window has four samples, and the quadratic leaves one
degree of freedom. The variance is then 124/225 of the modified Hadamard
variance term for term, on any record.

## EDF

The exact EDF depends on the record length and the averaging factor
through ``s = τ/T`` alone, within 1 % for records of 512, 1024 and 2048
points. The model is

```math
ν = \frac{b}{s} - c + d\,s + e\,s^3, \qquad s = \frac{τ}{T}, \quad T = (N-1)τ_0 ,
```

fitted by least squares on the relative residual over ``16 ≤ m ≤ N/4``,
the three record lengths together:

| Noise | ``b`` | ``c`` | ``d`` | ``e`` | Largest residual |
|:---|---:|---:|---:|---:|---:|
| White PM | 1.5949 | 6.535 | 20.18 | −52.1 | 4.5 % |
| Flicker PM | 1.0656 | 4.223 | 10.19 | 7.8 | 2.6 % |
| White FM | 0.9420 | 3.506 | 6.50 | 23.3 | 2.0 % |
| Flicker FM | 0.8938 | 3.124 | 3.94 | 31.2 | 2.5 % |
| Random-walk FM | 0.8596 | 2.676 | 0.43 | 48.8 | 2.8 % |
| Flicker-walk FM | 0.8129 | 2.115 | −4.27 | 83.8 | 2.9 % |
| Random-run FM | 0.6508 | 1.991 | −2.63 | 78.3 | 1.0 % |

The first two terms are the form of the total variance EDF
([Greenhall et al. (1999), Eq. (26)](@cite greenhall-1999-totvar)). They
are not enough here. The exact EDF turns upward in the last octave, where
the reflection does the most for the few windows that remain, and two
terms fitted below ``T/16`` give a third to a half of the exact EDF at
``τ = T/4``. The terms in ``s`` and ``s^3`` carry the rise.

At ``τ = T/4``, the longest averaging time, one window is left. The
modified Hadamard variance has one degree of freedom there. The modified
Hadamard total variance has 4.0 on white PM, 2.7 on flicker PM, 2.2 on
white FM, 1.9 on flicker FM, 1.6 on random-walk FM, 1.4 on flicker-walk FM
and 1.2 on random-run FM.

!!! note "Below m = 16"
    The model is fitted from ``m = 16`` and used from ``m = 2``. On the FM
    noises it is within 6 % of the exact EDF from ``m = 3`` and 4 to 13 %
    low at ``m = 2``, which widens the interval. On flicker PM it is 2 to
    10 % high at ``m = 2`` and 3, depending on the flicker model, and
    within 4 % from ``m = 4``. On white PM the fitted form would be high by
    53 % at ``m = 2``, 34 % at ``m = 3``, 14 % at ``m = 4`` and 3 % at
    ``m = 8``, so it is multiplied by the exact ratio at each ``m`` from 2
    to 15, which brings it within 1.4 %. At ``m = 1`` the EDF is that of
    [`mhdev`](@ref).

## Degrees of freedom gained

The exact EDF of the modified Hadamard total variance over that of the
modified Hadamard variance at the same ``τ``, for a record of 1024 points:

| Noise | ``m = 16`` | ``T/16`` | ``T/8`` | ``T/4`` |
|:---|---:|---:|---:|---:|
| White PM | 1.38 | 1.41 | 1.50 | 4.0 |
| Flicker PM | 1.19 | 1.21 | 1.27 | 2.7 |
| White FM | 1.10 | 1.12 | 1.17 | 2.2 |
| Flicker FM | 1.05 | 1.06 | 1.10 | 1.9 |
| Random-walk FM | 1.01 | 1.02 | 1.06 | 1.6 |
| Flicker-walk FM | 0.99 | 1.00 | 1.03 | 1.4 |
| Random-run FM | 0.98 | 0.98 | 1.00 | 1.2 |

The extension gains the most on PM noise and in the last octave. On
random-walk FM and redder noise below ``T/8`` it gains nothing in
precision. What the estimator offers there is its insensitivity to drift
within each window.

## Accuracy in simulation

In simulation (1024 points, 2000 records per noise type) with the noise
type stated, the corrected variance averages within 2 % of the modified
Hadamard variance up to ``τ = T/8`` and within 4 % at ``T/4``, on all
seven noise types, and the 68.3 % interval holds the true value in 65 to
77 % of records ([Bias correction](bias.md),
[Confidence intervals](confidence.md)).

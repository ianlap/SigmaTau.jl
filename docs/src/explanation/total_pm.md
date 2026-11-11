# Total estimators on PM noise

The published bias factors and EDF of the total deviation cover white,
flicker and random-walk FM
([Greenhall et al. (1999), Table I](@cite greenhall-1999-totvar);
[SP1065 §5.11](@cite riley-2008-sp1065)), and those of the Hadamard total
deviation cover the five FM noise types from ``τ = 16τ_0``
([Howe et al. (2005), Table I](@cite howe-2005-tothvar)). On white and
flicker PM both estimators read high, and neither their own FM tables nor
the EDF of their parents describe them. This page derives the factors and
the EDF models that the package uses on PM noise, and the exact factors of
the modified and Hadamard totals below ``16τ_0``.

Every number in the derivation is computed without random numbers by
[`tools/total_pm_models.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/tools/total_pm_models.jl),
which runs in under a minute.

## Method

A total variance is a quadratic form ``x^\mathsf{T} Q x`` in the phase
record ``x``: the extension, the differences and the average are all
linear. For Gaussian noise whose phase has covariance
``Σ = F F^\mathsf{T}``, its mean and its EDF are

```math
E = \operatorname{tr}(G), \qquad
ν = \frac{2E^2}{\mathrm{Var}} = \frac{\operatorname{tr}(G)^2}{\operatorname{tr}(G^2)},
\qquad G = F^\mathsf{T} Q F ,
```

with the EDF as defined by
[Greenhall and Riley (2003), Eq. (1)](@cite greenhall-2003-edf-stability).
The bias factor ``B`` is ``E`` over the expected value of the parent, the
Allan or the Hadamard variance, on the same noise.

The tool builds ``Q`` from the definitions of the estimators and does not
call the package. The test suite works in the other direction: it extracts
``Q`` from the package kernels and holds the models to the same moments.

## Two models of flicker PM

White PM has ``Σ = I``. Flicker PM has no unique discrete form, and two are
used.

- The **sampled** process is continuous ``1/f`` phase noise averaged over
  each sampling interval. Its generalized autocovariance at lag ``k`` is
  ``-Δ^2[t^2 \ln t](k)``, where ``Δ^2`` is the second central difference.
  It is the model of the EDF algorithm of
  [Greenhall and Riley (2003)](@cite greenhall-2003-edf-stability), which
  the Allan and Hadamard families use.
- The **filtered** process is white noise through the half-order filter of
  [Kasdin and Walter (1992)](@cite kasdin-1992-discrete-simulation), the
  model of [`noise_gen`](@ref).

On the bias the two agree within 2 % from ``m = 16``, and the factors below
are fitted to their mean. On the EDF they do not agree: at the same
low-frequency level the filtered process has more power near the Nyquist
frequency, and 2 to 22 % more degrees of freedom. The EDF models follow the sampled process, like the
rest of the package. On records from `noise_gen` they are therefore up to
17 % low, and the interval is up to 9 % wider than it needs to be.

## Bias

With ``s = τ/T`` and ``T = (N-1)τ_0``:

| | White PM | Flicker PM |
|:---|:---|:---|
| Total | ``1 + \frac{8}{3}\,\frac{⌊m/2⌋}{N-2}`` | ``1 + (0.593 + 0.068 \ln m)\,s`` |
| Hadamard total, even ``m`` | ``1.3838 - 0.2759/m`` | ``\frac{5.148 + 1.440 \ln m}{4.75 + \ln m} - \frac{0.114}{m}`` |
| Hadamard total, odd ``m`` | ``1.3697 - 0.8198/m`` | ``\frac{5.148 + 1.440 \ln m}{4.75 + \ln m} - \frac{0.396}{m}`` |

**Total variance.** The bias grows in proportion to ``τ/T``, as it does
with the opposite sign on flicker and random-walk FM. At ``τ = T/4`` the
variance reads 33 % high on white PM and, for a record of 1024 points,
24 % high on flicker PM. The cause is the extension: the odd reflection
about each end puts the end sample, at twice its weight, into every term
that reaches past the record. The white PM form is exact. The flicker PM
form is within 1.0 % of both flicker models up to ``τ = T/2``.

**Hadamard total variance.** The bias depends on ``m`` and not on the
record length, as it does on FM noise, and it is large at every ``τ``. On
white PM the variance reads 25 % high at ``m = 2``, 37 % at ``m = 16`` and
38 % beyond. On flicker PM it reads 4 % or 11 % high at ``m = 2``,
depending on the flicker model, 21 % at ``m = 16`` and 28 % at
``m = 256``. The cause is
the drift removal: the slope taken from each subsequence of ``3m``
frequency samples is the difference of two half-averages, which three or
four phase samples determine, and on PM noise its error is as large as the
noise being measured. Odd and even ``m`` follow different curves because
the half-averages share the middle phase sample only when ``3m`` is even.

The white PM forms are within 0.1 % of the exact value at every ``m``.
The flicker PM form is within 1.0 % of both flicker models from
``m = 16``. Below that the models part, and at ``m = 2`` the form is 3.3 %
above one and 3.7 % below the other.

## Averaging factors below 16

The published factors of the modified total and Hadamard total variances
hold from ``τ = 16τ_0``
([Howe and Peppler (2001), Table 2](@cite howe-2001-total-estimators);
[Howe et al. (2005), Table I](@cite howe-2005-tothvar)). Below it the
expectation depends on ``m``, and [`bias_correction`](@ref) uses the exact
ratio at each ``m`` from 2 to 15, computed as above.

The Hadamard total variance is the modified total computation applied to
the frequency samples in place of the phase samples. One table therefore
serves both estimators: its rows are white PM to random-walk FM for the
modified total and white FM to random-run FM for the Hadamard total. The
ratio does not depend on the record length, so one window gives it. The
phase of each FM noise is white or flicker PM summed once or twice, and
the flicker rows are the mean of the two flicker models, which lie within
2.4 % of it in the first flicker row and 1.2 % in the second.

The values depart most from the published ones at the smallest ``m``. On
white PM the modified total factor is 0.938 at ``m = 3``, where the
published value is 0.995, and on white FM it is 0.814 at ``m = 2`` against
0.771 ([Bias correction](bias.md) has the table).

## EDF

| | Model | Largest difference from the exact EDF |
|:---|:---|:---|
| Total, white PM | ``\dfrac{72 N^2 B^2}{140N + 225m + 64m^2}`` | 0.8 % |
| Total, flicker PM | ``\dfrac{B^2}{g(m)/N + 0.93\,(m/N)^2}`` | 1.4 % up to ``τ = T/3``, 2.5 % at ``T/2`` |
| Hadamard total | ``B^2\,\dfrac{N}{g}\;\dfrac{1 - a\,e^{-(N-3m)/(bN)}}{1 + c\,\min(s, s_0)}`` | see below |

Here ``g(m) = m/(0.789 + 0.276 \ln m + 0.204 \ln^2 m)`` is ``N/ν`` of the
overlapping Allan variance on sampled flicker PM for ``N ≫ m``, fitted
within 0.9 %.

In the white PM form of the total variance the terms have a meaning. The
interior of the record contributes ``140N``: second differences of white
noise have covariances 6, −4 and 1 at lags 0, ``m`` and ``2m``, and
``2(6^2 + 2·4^2 + 2·1^2) = 140``. The end samples shared through the
reflections contribute ``64m^2``, which is why the EDF falls so fast at
long ``τ``. The coefficient 225 is fitted.

For the Hadamard total, ``g`` is 2.45 on white PM and ``g(m)`` on flicker
PM, and ``(c, s_0, a, b)`` is (2.0, 0.125, 1.0, 0.03) and
(2.5, 0.16, 0.6, 0.04). These constants were set by hand against the exact
EDF for records of 512 and 1024 points. The model differs from the exact
EDF by:

| | ``m < 16`` | ``m ≥ 16``, ``τ ≤ T/4`` | ``T/4`` to ``0.32\,T`` | beyond |
|:---|:---:|:---:|:---:|:---:|
| White PM | −22 to +4 % | −6 to +2 % | −8 to −3 % | −27 to −7 % |
| Flicker PM | −4 to +5 % | −1 to +2 % | −2 to +7 % | −5 to +10 % |

A negative difference means an EDF that is too small and an interval that
is too wide. The white PM model is low for odd ``m`` below 16 and at the
longest averaging times, where the last subsequences run out.

The models replace the EDF of the parent, which the two estimators would
otherwise take. For a record of 1024 points:

| | ``m = 16`` | 64 | 128 | 256 |
|:---|---:|---:|---:|---:|
| Total, white PM, exact | 481 | 210 | 84 | 30 |
| Allan variance model | 514 | 478 | 432 | 354 |
| Total, flicker PM, exact | 196 | 74 | 38 | 18 |
| Allan variance model | 195 | 78 | 45 | 23 |
| Hadamard total, white PM, exact | 756 | 708 | 641 | 677 |
| Hadamard variance model | 427 | 379 | 318 | 256 |
| Hadamard total, flicker PM, exact | 274 | 116 | 66 | 40 |
| Hadamard variance model | 164 | 63 | 34 | 17 |

On white PM the total variance has far fewer degrees of freedom than the
Allan variance at long ``τ``, because a few shared end samples dominate
it, and the parent's EDF would give an interval three times too narrow at
a quarter of the record.
The Hadamard total variance has about twice the degrees of freedom of the
Hadamard variance on both PM noises, and the parent's EDF would give an
interval too wide.

## In a default call

The factors are exact for the noise type they are given. A default call
identifies the noise type first, and with few averaged points the label is
uncertain ([Noise identification](noise_identification.md)). In simulation
(1024 points, 2000 records per noise type), the mean of the variance over
that of the parent:

| | ``m = 16`` | 64 | 128 | 256 |
|:---|---:|---:|---:|---:|
| `totdev`, white PM, uncorrected | 1.02 | 1.08 | 1.17 | 1.33 |
| default | 1.00 | 1.00 | 1.00 | 1.00 |
| `totdev`, flicker PM, uncorrected | 1.01 | 1.05 | 1.10 | 1.20 |
| default | 1.00 | 1.00 | 1.00 | 1.19 |
| `totdev`, white FM, uncorrected | 1.00 | 1.00 | 1.01 | 1.01 |
| default | 1.00 | 1.01 | 1.02 | 1.06 |
| `htotdev`, white PM, uncorrected | 1.37 | 1.38 | 1.38 | 1.38 |
| default | 1.00 | 1.00 | 1.00 | 1.00 |
| `htotdev`, flicker PM, uncorrected | 1.22 | 1.26 | 1.27 | 1.28 |
| default | 1.01 | 1.04 | 1.04 | 1.25 |
| `htotdev`, white FM, uncorrected | 1.00 | 1.00 | 1.00 | 1.00 |
| default | 1.00 | 1.00 | 1.01 | 1.00 |

On white PM the default removes the bias at every averaging time of the
table. On flicker PM it removes all but 4 % of it down to 8 averaged
points (``m = 128``). With 4 points flicker PM is labeled white FM in most
records, and its bias stays. White FM pays up to 6 % in the total variance
at the longest averaging times of the table, for the records labeled
flicker or random-walk FM. Stating the noise type removes both effects:

```julia
htotdev(p; noise = :WHPM)
```

## Other conventions

Stable32 and the published tables apply no factor on PM noise.
`correct_bias=:published` and `correct_bias=:stable32` do the same: on the
rows labeled white or flicker PM they return the uncorrected deviation.
The default is ``\sqrt{B}`` times smaller: 10 to 15 % for the Hadamard
total on white PM, and 7 % for the total on white PM at ``τ = T/8``.

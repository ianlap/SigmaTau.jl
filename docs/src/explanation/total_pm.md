# Total estimators on PM noise

The published bias and EDF tables of TOTDEV and HTOTDEV stop at white FM:
TOTDEV's cover white, flicker and random-walk FM
[greenhall-1999-totvar](@cite), [riley-2008-sp1065](@cite), HTOTDEV's white
FM to random-run FM [howe-2001-total-estimators](@cite). On white PM and flicker
PM both estimators read high, and the EDF of their unextended parents is the
wrong size for them. SigmaTau corrects both from the exact moments of the
two statistics, computed by
[`tools/total_pm_models.jl`](https://github.com/ianlap/SigmaTau.jl/blob/main/tools/total_pm_models.jl)
with no random numbers, in under a minute.

## Method

Both statistics are quadratic forms x′Qx in the phase record. For Gaussian
noise with phase covariance Σ = FF′,

```math
E = \operatorname{tr}(F'QF), \qquad
\mathrm{edf} = \frac{\operatorname{tr}(F'QF)^2}{\lVert F'QF\rVert^2} ,
```

the EDF of Greenhall and Riley [greenhall-2003-edf-stability](@cite). The
bias factor B is E over the expectation of the parent, AVAR or HVAR, on the
same noise. The tool builds Q from the definitions; the test suite extracts
it from the package kernels and holds the shipped models to the same moments.

White PM is Σ = I. Flicker PM has no unique discrete form, and two are used:

- the **sampled** process, a continuous 1/f phase noise read every τ₀, whose
  structure function D(k) is twice the second central difference of t² ln t
  at t = k. This is the model behind the Greenhall–Riley EDF of every other
  estimator in the package.
- the **filtered** process, white noise through the half-order filter of
  [`noise_gen`](@ref).

The two agree on the bias to 2 % from m = 16. On the EDF they do not: the
filtered process has less power near the Nyquist frequency and up to 22 %
more degrees of freedom. The shipped EDF follows the sampled process, as the
rest of the package does, so on `noise_gen` flicker PM it reads up to 17 %
low and the interval is up to 9 % wider than it need be.

## Bias

With s = τ/T, T = (N − 1)τ₀:

| | White PM | Flicker PM |
|:---|:---|:---|
| TOTVAR | 1 + (8/3)·⌊m/2⌋/(N − 2), exact | 1 + (0.593 + 0.068 ln m)·s |
| TOTHVAR, even m | 1.3838 − 0.2759/m | (5.148 + 1.440 ln m)/(4.75 + ln m) − 0.114/m |
| TOTHVAR, odd m | 1.3697 − 0.8198/m | (5.148 + 1.440 ln m)/(4.75 + ln m) − 0.396/m |

TOTVAR reads high in proportion to τ/T, as it reads low on flicker and
random-walk FM: 33 % in variance on white PM and 24 % on flicker PM at
τ = T/4. The odd reflection about each end puts the end sample, doubled,
into every extended term.

TOTHVAR's bias does not depend on the record length, as on FM. It is large
at every τ: the variance reads 25 % high on white PM at m = 2, 37 % at
m = 16 and 38 % beyond; on flicker PM 4 to 11 % at m = 2, 21 % at m = 16
and 28 % at m = 256. On PM the half-average slope removed from each 3m-point
subsequence is set by three or four phase samples, and its error is as large
as the noise being measured. Odd and even m follow different curves because
the two half-averages share the middle phase sample only when 3m is even.

The white PM forms are within 0.1 % of the exact value at every m. The
flicker PM forms are fitted to the mean of the two flicker models and are
within 1.0 % of either for TOTVAR through τ = T/4 (1.6 % at T/2), and from
m = 16 for TOTHVAR. Below that the two models part, and at m = 2 the TOTHVAR
form is 3.5 % from each.

## EDF

| | Model | Largest residual |
|:---|:---|---:|
| TOTVAR, white PM | 72N²B² / (140N + 225m + 64m²) | 0.8 % |
| TOTVAR, flicker PM | B² / (g(m)/N + 0.93·(m/N)²) | 2.0 % to τ = T/3, 3.7 % at T/2 |
| TOTHVAR | B²·(N/g)·(1 − a·e^(−(N − 3m)/(bN))) / (1 + c·min(s, s₀)) | below |

Here g(m) = m / (0.789 + 0.276 ln m + 0.204 ln² m) is N/edf of overlapping
AVAR on sampled flicker PM for N ≫ m. In the white PM form, 140N is the
interior of the record (second differences with covariances 6, −4, 1) and
64m² the samples shared through the reflected ends; 225m is fitted.

For TOTHVAR, g is 2.45 on white PM and g(m) on flicker PM, and
(c, s₀, a, b) = (2.0, 0.125, 1.0, 0.03) and (2.5, 0.16, 0.6, 0.04). The
constants were set by hand against the exact EDF at N = 512, 1024 and 2048.
Residuals of the shipped EDF, low to high:

| | m < 16 | m ≥ 16, τ ≤ T/4 | to 0.32 T | last few m |
|:---|---:|---:|---:|---:|
| white PM | −22 to +4 % | −7 to +2 % | −8 to −3 % | −27 to −4 % |
| flicker PM (sampled) | −4 to +5 % | −1 to +3 % | −2 to +7 % | −5 to +11 % |

The white PM form reads low, the safe side, most of all for odd m below 16
and at the reach limit τ = T/3, where one subsequence is left.

What the correction replaces is the EDF of the parent. For a 1024-point
record:

| | m = 16 | 64 | 128 | 256 |
|:---|---:|---:|---:|---:|
| TOTVAR white PM, exact | 481 | 210 | 84 | 30 |
| AVAR model | 514 | 478 | 432 | 354 |
| TOTVAR flicker PM, exact | 196 | 74 | 38 | 18 |
| AVAR model | 195 | 78 | 45 | 23 |
| TOTHVAR white PM, exact | 756 | 708 | 641 | 677 |
| HVAR model | 427 | 379 | 318 | 256 |
| TOTHVAR flicker PM, exact | 274 | 116 | 66 | 40 |
| HVAR model | 164 | 63 | 34 | 17 |

On white PM the reflected ends of TOTVAR are a few shared samples that
dominate the variance at long τ, so its interval was far too narrow; TOTHVAR
has about twice the degrees of freedom of HVAR, so its interval was too wide.

## When the label is wrong

The factor follows the identified noise type, and with fewer than about 32
averaged points white FM is labelled as a phase noise in 20 to 50 % of
records. In simulation (N = 1024, 2000 records per noise, default
identification), the mean of σ over its true value:

| | m = 16 | 64 | 128 | 256 |
|:---|---:|---:|---:|---:|
| HTOTDEV, white PM, without the PM factor | 1.17 | 1.18 | 1.18 | 1.19 |
| with it | 1.00 | 1.03 | 1.07 | 1.08 |
| HTOTDEV, flicker PM, without | 1.10 | 1.12 | 1.13 | 1.14 |
| with it | 1.01 | 1.03 | 1.05 | 1.06 |
| HTOTDEV, white FM, without | 1.00 | 1.00 | 0.99 | 0.97 |
| with it | 0.99 | 0.96 | 0.93 | 0.92 |
| TOTDEV, white PM, without | 1.01 | 1.04 | 1.08 | 1.16 |
| with it | 1.00 | 1.01 | 1.03 | 1.06 |
| TOTDEV, white FM, without | 1.00 | 0.99 | 0.99 | 0.98 |
| with it | 1.00 | 0.99 | 0.97 | 0.95 |

The PM factor removes most of a bias of 10 to 19 % on phase noise at every τ
of HTOTDEV, and what is left at long τ is the records labelled FM. It costs
up to 5 % on white FM and 3 % on flicker FM at the longest τ, where the
scatter of a single record is 20 to 30 %.

## Stable32

Stable32 applies no bias factor on PM, so on rows labelled white or flicker
PM the default TOTDEV and HTOTDEV are lower than its values, HTOTDEV by 10
to 15 %. `correct_bias=false` returns the raw statistic, which matches
Stable32's on those rows.

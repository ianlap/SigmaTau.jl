# Theory: Noise Identification

Each averaging-factor `m` of an estimator may correspond to a different
dominant noise type. The EDF formula and the χ² CI both depend on the
α exponent of the local power-law process, so SigmaTau classifies the
noise *per τ* before assigning a confidence interval.

## Why classify?

The Greenhall–Riley EDF formulas (see [Confidence](confidence.md))
parameterize on α explicitly. Pretending the noise is white when it is
flicker yields an overconfident CI; pretending it is random-walk when
it is white yields an underconfident one. Classification matters
quantitatively, not just descriptively.

The five power-law families SigmaTau recognises are white phase
modulation (WHPM, α=2), flicker phase modulation (FLPM, α=1), white
frequency modulation (WHFM, α=0), flicker frequency modulation
(FLFM, α=−1), and random-walk frequency modulation (RWFM, α=−2). The
α exponent feeds straight into the EDF and bias-correction tables.

## Lag-1 ACF on differenced phase

Riley & Greenhall 2004 [@cite Riley2004] show that the lag-1
autocorrelation `r₁` of an appropriately differenced phase series maps
monotonically onto the spectral exponent α. Their procedure:

1. Detrend (quadratic) and difference the input phase `m` times.
2. Compute `r₁` on the result.
3. Map `r₁` to α via the published thresholds.

For `r₁` near 0.5 the noise is white; near zero, flicker; near −0.5,
random-walk. The α-to-`r₁` mapping is implemented in
`lib/SigmaTauStability/src/noise/lag1.jl`.

```julia
identify_noise(PhaseData(x, 1.0); m=8)
```

The result is a `Symbol`: one of `:WHPM`, `:FLPM`, `:WHFM`, `:FLFM`,
`:RWFM`.

Internally the kernel iterates the difference operator: at each step it
recomputes `r₁`, converts to ρ = r₁ / (1 + r₁), and stops once ρ falls
below 0.25 (or `dmax` differences have been taken). The terminal
difference order `d` plus ρ recovers α via α = −2(ρ + d) + 2. The
iterative criterion is what makes the method robust at boundaries
between adjacent power-law families.

## B1/R(n) fallback

For short records or borderline cases, SigmaTau falls back to the
Allan-variance B1 ratio with the R(n) factor for WPM/FPM
disambiguation, per SP1065 §6 [@cite RileyHowe2008]. This path is
triggered automatically when the lag-1 method's confidence is low.

B1 compares the classical (standard) variance of the averaged
frequency to the Allan variance and discriminates RWFM/FLFM/WHFM/FPM
boundaries. When B1 lands in the WPM/FPM region, an additional R(n)
test — the ratio of MDEV² to ADEV² — separates the two phase-modulated
families, since MDEV's extra phase-averaging window collapses WHPM
faster than FLPM.

## NEFF_RELIABLE = 30

The lag-1 method needs a minimum effective sample count for `r₁` to
be a reliable estimator. SigmaTau uses `NEFF_RELIABLE = 30` as the
threshold; below this, classification falls back to B1/R(n) regardless
of `r₁`. The threshold tracks the Riley R 2020 recommendation
[@cite Riley_R_2020] and the historical SigmaTau policy mandate.

## The `noise_type` field

`adev`, `mdev`, etc. populate `StabilityResult.noise_type` with one
classification per τ. Empty when classification was not requested or
not possible. When a particular m fails (e.g. degenerate variance,
too few points after differencing), the kernel carries forward the
last reliable classification rather than emitting `:unknown` mid-run.

## Implementation notes

- 5σ outlier rejection plus linear detrend run before classification.
- `_noise_id_lag1acf` and `_noise_id_b1rn` are the two internal kernels;
  the public `identify_noise` chooses between them based on `N_eff`.
- Per-m quadratic detrending is on by default to remove frequency
  offset and drift from each decimated subseries; pass `detrend=false`
  to match Stable32's no-detrend convention when comparing fixtures
  point-for-point.

## See also

- [Theory: Confidence](confidence.md).
- [API: SigmaTauStability](../reference/stability.md).

## References

- Riley & Greenhall, *Power-law noise identification using the lag-1
  autocorrelation*, PTTI 2004 [@cite Riley2004].
- SP1065 §6 [@cite RileyHowe2008].
- Riley, *Frequency Stability Analysis Using R*, 2020
  [@cite Riley_R_2020].

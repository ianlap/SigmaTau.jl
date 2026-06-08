# Theory: PDEV (PVAR) confidence intervals

`pdev` computes the parabolic deviation σ_PDEV(τ), the square root of the
parabolic variance PVAR — a wavelet variance built from a least-squares
straight-line (equivalently, parabolic-weight) fit of the phase over each
analysis window [Vernotte–Lenczner–Bourgeois–Rubiola 2016](@cite vernotte-2016-pvar).
Like every other SigmaTau deviation, when `ci=true` it reports a per-τ noise
type, an equivalent number of degrees of freedom (EDF) `ν`, and a χ²-based
confidence interval. This page documents the EDF model behind that interval.

Unlike MHTOTDEV (see ["MHTOTDEV bias and EDF"](mhtotdev_bias_edf.md)), PVAR's EDF
is **published**, so SigmaTau reproduces the literature model rather than
measuring its own.

## χ² confidence interval

A variance estimator that is `k·χ²ν`-distributed satisfies

```math
\nu = \frac{2\,E^2[\hat\sigma^2(\tau)]}{V[\hat\sigma^2(\tau)]}
```

([Vernotte–Chen–Rubiola 2020](@cite vernotte-2020-pvar-noninteger) Eq. 17; the
same estimator SigmaTau uses for every family). Given `ν`, the two-sided interval
at confidence `p` on the *variance* is `[ν/χ²_{hi}, ν/χ²_{lo}]·\hat\sigma^2`, and
the deviation bounds are the square roots — exactly the `confidence_intervals`
path shared with ADEV/HDEV/etc.

## The EDF model

For white PM noise the EDF has a closed form
[(Vernotte et al. 2016](@cite vernotte-2016-pvar), Eq. 24). The 2020 paper
generalizes it to all five power-law noises as

```math
\nu \approx \frac{35}{A(\alpha)\,(m/M) - B\,(m/M)^2},
\qquad M = N - 2m, \qquad B = 12,
```

([Vernotte–Chen–Rubiola 2020](@cite vernotte-2020-pvar-noninteger), Eqs. 16, 18),
where `m = τ/τ0` is the averaging factor, `N` the number of phase samples, and
`M = N − 2m` the number of PVAR estimates (the same `M` used in the kernel). The
amplitude coefficient is the cubic

```math
A(\alpha) = 27 + \tfrac{1}{4}\alpha + \tfrac{5}{14}\alpha^2 - \tfrac{3}{4}\alpha^3 ,
```

determined from massive Monte Carlo with `A(+2) = 23` fixed by the white-PM
closed form. Evaluated at the five integer noise types it gives

| α (noise) | +2 (WHPM) | +1 (FLPM) | 0 (WHFM) | −1 (FLFM) | −2 (RWFM) |
|:---------:|:---------:|:---------:|:--------:|:---------:|:---------:|
| `A(α)`    | 23        | 27        | 27       | 28        | 34        |

SigmaTau stores these integer-α values (`_pvar_A`), since noise identification
only yields integer α.

### Endpoints

Two regimes need separate handling:

- **`m = 1`.** PVAR(τ0) ≡ AVAR(τ0): the kernel `_pdev_core` returns overlapping
  ADEV at `m = 1`, so the EDF uses ADEV's Greenhall–Riley value
  (`_calc_edf_core` with `d = 2`). This is the exact result for the identity and
  sidesteps the paper's note that the closed-form approximation is least accurate
  at `m = 1, 2`. (At `m = 2` no exact fallback exists, so the formula is used as
  published; it reads somewhat high there.)
- **Large τ (`m ≳ N/4`).** The Eq. (16) approximation degrades beyond `m ≈ N/4`.
  The paper fills `N/4 < m ≤ N/2` with a semi-logarithmic interpolation
  ([Vernotte–Chen–Rubiola 2020](@cite vernotte-2020-pvar-noninteger), Eqs. 22–23)

  ```math
  \nu(m) = a\,\ln m + b, \qquad
  a = \frac{\nu(m_1) - 1}{\ln m_1 - \ln m_2},
  ```

  anchored at the Eq.-16 value at `m_1 = \mathrm{round}(1.11\,N/4)` and at `ν = 1`
  (a single estimate) at `m_2 = \mathrm{round}(0.901\,N/2) \approx N/2`. SigmaTau
  holds `ν = 1` for any `m ≥ m_2` rather than extrapolating below one.

## Accuracy and scope

The model reproduces the paper's Monte Carlo EDF table (Table III, N = 2048) to
within the stated ≈ ±10 % across the noise types — e.g. white PM at m = 64/128/256
gives 46.5/22.1/10.0 against the tabulated 46.9/22.0/10.0. The in-suite tripwire
`test/stab/pdev_edf_mc.jl` re-measures `ν = 2E[V]²/Var[V]` from a small fresh
Monte Carlo and checks it against `_pvar_edf` at a few cells.

PVAR is an **unbiased** wavelet variance, so — unlike the total family — no bias
correction is applied to `pdev`.

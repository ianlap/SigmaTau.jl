# Theory: Allan Family

The Allan family covers the six estimators built from second- and
third-difference operators on phase data: ADEV, MDEV, TDEV, HDEV, MHDEV,
and HTDEV. They share the same boundary handling (no extension) and
differ in difference order, phase-averaging, and final scaling.

## ADEV — overlapping Allan deviation

The original two-sample variance, in its overlapping estimator form
(SP1065 Eq. 14):

```math
\sigma^2_y(\tau) \;=\; \frac{1}{2(N - 2m)\,(m\tau_0)^2}
\sum_{i=1}^{N-2m} \bigl(x_{i+2m} - 2\,x_{i+m} + x_i\bigr)^2
```

with `τ = m·τ₀` and `N` phase samples. ADEV slope versus α:
WPM/FPM both give μ_dev = −1 (degeneracy), WFM gives μ_dev = −1/2, FFM
gives μ_dev = 0, RWFM gives μ_dev = +1/2.

In SigmaTau:

```julia
adev(PhaseData(x, τ₀), [1, 2, 4, 8, 16])
```

(Cite [@cite Sullivan_NBS_TN_1337] for origin; SP1065 §5
[@cite RileyHowe2008].)

## MDEV — modified Allan deviation

A phase-averaged second difference. SP1065 Eq. 16:

```math
\mathrm{MVAR}(\tau) \;=\; \frac{1}{2 m^4 \tau_0^2 \, N_e}
\sum_{j=1}^{N_e} \biggl[\sum_{k=0}^{m-1}
\bigl(x_{j+k+2m} - 2\,x_{j+k+m} + x_{j+k}\bigr)\biggr]^2,
\qquad N_e = N - 3m + 1.
```

The inner phase-averaging step decouples WPM (μ_dev = −3/2) from FPM
(μ_dev = −1). In `lib/SigmaTauStability/src/core/allan.jl` this is
implemented in third-difference form via prefix sums (Greenhall 1997
[@cite Greenhall1997]) — algebraically identical to the SP1065 form
above; see `legdocs/equations/allan.md` for the equivalence proof.

```julia
mdev(PhaseData(x, τ₀), τs)
```

## TDEV — time deviation

A scaled MDEV (SP1065 Eq. 17):

```math
\mathrm{TVAR}(\tau) \;=\; \frac{\tau^2}{3}\,\mathrm{MVAR}(\tau).
```

`tdev(...)` wraps `mdev(...)` and applies the scaling.

## HDEV — overlapping Hadamard deviation

A third-difference variant that is insensitive to linear frequency
drift. SP1065 §5.4:

```math
\sigma^2_{H,y}(\tau) \;=\; \frac{1}{6\,(N - 3m)\,(m\tau_0)^2}
\sum_{i=1}^{N-3m} \bigl(x_{i+3m} - 3\,x_{i+2m} + 3\,x_{i+m} - x_i\bigr)^2
```

The `1/6` prefactor arises from the third-difference variance
normalization. Drift insensitivity follows because a linear-in-`t` term
is annihilated by the third difference.

```julia
hdev(PhaseData(x, τ₀), τs)
```

## MHDEV — modified Hadamard deviation

Phase-averaged third difference; same relation to HDEV that MDEV has
to ADEV.

```julia
mhdev(PhaseData(x, τ₀), τs)
```

## HTDEV — Hadamard time deviation

HTDEV is to MHDEV what TDEV is to MDEV: a τ-scaled time-domain
deviation built on the third-difference (Hadamard) kernel instead of
the second-difference (Allan) kernel. SigmaTau implements it as
`mhdev` followed by the scaling `τ / √(10/3)`:

```math
\sigma_{HT,y}(\tau) \;=\; \frac{\tau}{\sqrt{10/3}} \, \sigma_{MH,y}(\tau).
```

```julia
htdev(PhaseData(x, τ₀), τs)
```

The `√(10/3)` factor mirrors the `√3` factor in TDEV and follows from
the third-difference Hadamard kernel variance.

**Provenance.** The construction is original to this package; the
standard time-and-frequency references — SP1065 [@cite RileyHowe2008],
IEEE 1139-2022 [@cite IEEE1139_2022], NBS-TN-1337
[@cite Sullivan_NBS_TN_1337] — do not define it. The earlier name
`ldev` is retained as a deprecated alias for one release.

## Slope vs noise table

For each estimator, the deviation slope `μ_dev` versus the spectral
exponent `α` of `S_y(f)`:

| α  | ADEV | MDEV | HDEV | MHDEV | TDEV | HTDEV |
|----|------|------|------|-------|------|------|
| +2 (WPM)  | −1     | −3/2 | −1     | −3/2 | −1/2  | −1/2 |
| +1 (FPM)  | −1*    | −1   | −1*    | −1   | 0     | 0    |
|  0 (WFM)  | −1/2   | −1/2 | −1/2   | −1/2 | +1/2  | +1/2 |
| −1 (FFM)  |  0     |  0   |  0     |  0   | +1    | +1   |
| −2 (RWFM) | +1/2   | +1/2 | +1/2   | +1/2 | +3/2  | +3/2 |

(*) ADEV/HDEV at FPM include a `log(2π·m)` factor (see SP1065 §5).

## Demonstration

```@example allan
using SigmaTau, Random
Random.seed!(0)

# Pure WPM noise
N = 8192
x = randn(N)
τs = [1, 4, 16, 64, 256]

a  = adev(PhaseData(x, 1.0), τs).dev
m  = mdev(PhaseData(x, 1.0), τs).dev

# Slopes (deviation log-log)
slope(σ) = (log10(σ[end]) - log10(σ[1])) / (log10(τs[end]) - log10(τs[1]))
round.([slope(a), slope(m)]; sigdigits=3)
```

ADEV slope should be near −1; MDEV slope should be near −3/2 — the
characteristic split that makes MDEV able to disambiguate WPM from FPM.

## Implementation notes

- All cores in `lib/SigmaTauStability/src/core/{allan,hadamard}.jl` take
  `Vector{Float64}` and return raw arrays; the public API in
  `lib/SigmaTauStability/src/api/` wraps them and returns
  `StabilityResult`.
- `MDEV/MHDEV` use a prefix-sum form algebraically equivalent to the
  textbook `1/m⁴` form; see `legdocs/equations/allan.md` for the proof.
- `tdev` and `htdev` are scaling wrappers; they do no extra kernel work.

## See also

- [Theory: Total family](total_family.md) — boundary-extended estimators.
- [Theory: Confidence](confidence.md) — EDF and CI for these estimators.
- [API: SigmaTauStability](../reference/stability.md) — function signatures.

## References

- SP1065 §5 [@cite RileyHowe2008].
- NBS Technical Note 1337 [@cite Sullivan_NBS_TN_1337].
- Greenhall, *Third-difference approach to MVAR*, IEEE T-IM 1997
  [@cite Greenhall1997].
- IEEE 1139-2022 [@cite IEEE1139_2022] for canonical ADEV / MDEV / HDEV
  / TDEV / MHDEV definitions; HTDEV is **not** defined there.

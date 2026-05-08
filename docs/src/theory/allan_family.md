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

Time deviation. TDEV measures time-error stability, so its units are
seconds (it is a `\sigma_x` quantity, not `\sigma_y`). SP1065 Eq. 17
defines it as a scaled MDEV:

```math
\sigma_x^2(\tau) \;=\; \mathrm{TVAR}(\tau) \;=\; \frac{\tau^2}{3}\,\mathrm{MVAR}(\tau).
```

`tdev(...)` wraps `mdev(...)` and applies the `\tau / \sqrt{3}` scaling.
The factor of `\sqrt{3}` comes from the integral of the second-difference
sampling function squared against the fractional-frequency-to-time-error
kernel.

**Use case.** TDEV is the standard tool for characterizing the time
error of a clock or a time-distribution system — telecommunications
synchronization networks (PTP / SyncE / Recommendation G.810), GNSS
time transfer, atomic-clock comparison links, and any infrastructure
where the operative question is "how badly will a downstream
consumer's clock drift away from the reference over τ seconds?"
SP1065 §5 [@cite RileyHowe2008] frames it as the way to characterize
"the time error of a time source (clock) or distribution system." For
records dominated by white PM noise, TVAR reduces to the standard
variance of the time deviations themselves; for the other power-law
noises it remains a convergent estimator and inherits MDEV's WPM/FPM
disambiguation.

## HDEV — overlapping Hadamard deviation

A third-difference variant with two advantages over ADEV: linear
frequency drift is filtered out automatically, and the variance
integral converges over a wider range of low-frequency noise. SP1065
§5.4:

```math
\sigma^2_{H,y}(\tau) \;=\; \frac{1}{6\,(N - 3m)\,(m\tau_0)^2}
\sum_{i=1}^{N-3m} \bigl(x_{i+3m} - 3\,x_{i+2m} + 3\,x_{i+m} - x_i\bigr)^2 .
```

The `1/6` prefactor arises from the third-difference variance
normalization.

**Drift insensitivity.** A linear-in-`t` term in `y(t)` (constant
frequency drift) is annihilated by the third difference, so HDEV is
not contaminated by drift the way ADEV is. SP1065 demonstrates this
on a simulated rubidium record: ADEV picks up a `+τ` slope at long τ
without prior detrending, while HDEV gives essentially the same
answer as drift-removed ADEV [@cite RileyHowe2008].

**Noise-type convergence.** ADEV's variance integral diverges for
`α ≤ −3` (flicker walk FM, random run FM) at the low-frequency end.
HDEV remains finite down to `α = −4`, because the third difference
adds an extra factor of `f²` to the kernel that suppresses the
`f → 0` singularity [@cite Greenhall1997]. In practice this matters
only for records with very-low-frequency power-law content; most
laboratory clocks are well-described in `−2 ≤ α ≤ +2` and ADEV
suffices.

HDEV is the preferred Allan-family estimator for clocks with known
drift (Cs, Rb, H-maser) and for any analysis that needs to
characterize the longest power-law tails.

```julia
hdev(PhaseData(x, τ₀), τs)
```

## MHDEV — modified Hadamard deviation

A phase-averaged third difference. Combining the SP1065 form with the
prefix-sum / third-difference equivalence from Greenhall 1997
[@cite Greenhall1997]:

```math
\mathrm{MHVAR}(\tau) \;=\; \frac{1}{6\,m^4 \tau_0^2 \, N_e}
\sum_{j=1}^{N_e} \biggl[\sum_{k=0}^{m-1}
\bigl(x_{j+k+3m} - 3\,x_{j+k+2m} + 3\,x_{j+k+m} - x_{j+k}\bigr)\biggr]^2,
\qquad N_e = N - 4m + 1.
```

The relationship between MHDEV and HDEV mirrors the MDEV / ADEV
relationship: phase-averaging the inner third-difference window splits
the WPM/FPM degeneracy (slope `μ_dev` is `−3/2` under WPM, `−1` under
FPM), while the third-difference kernel preserves HDEV's drift
insensitivity. MHDEV is the right choice when a record contains
linear frequency drift *and* phase noise that ADEV / MDEV cannot
disambiguate. SigmaTau's kernel uses the prefix-sum form for
performance; the equivalence to the textbook expression above is in
`legdocs/equations/hadamard.md`.

```julia
mhdev(PhaseData(x, τ₀), τs)
```

## HTDEV — Hadamard time deviation

HTDEV is to MHDEV what TDEV is to MDEV: a τ-scaled time-domain
deviation built on the third-difference (Hadamard) kernel instead of
the second-difference (Allan) kernel. Like TDEV, HTDEV has units of
seconds (it is a `\sigma_x` quantity). SigmaTau implements it as
`mhdev` followed by the scaling `τ / √(10/3)`:

```math
\sigma_{x,\mathrm{HT}}(\tau) \;=\; \frac{\tau}{\sqrt{10/3}} \, \sigma_{y,\mathrm{MH}}(\tau).
```

```julia
htdev(PhaseData(x, τ₀), τs)
```

The `\sqrt{10/3}` factor mirrors the `\sqrt{3}` factor in TDEV: each
arises from the integral of the kernel's sampling function squared
against the frequency-to-time-error map. The second-difference (Allan)
kernel produces `1/3`; the third-difference (Hadamard) kernel produces
`3/10`.

**Use case.** HTDEV inherits TDEV's purpose — characterizing the time
error of a clock or time-distribution system — and adds two
properties TDEV cannot offer:

- **Linear frequency-drift insensitivity.** The third-difference
  kernel annihilates terms linear in `t`, so a record from a Cesium,
  Rubidium, or H-maser clock with significant drift can be analyzed
  without first removing the drift. SP1065 §5 makes the equivalent
  point for HDEV vs ADEV: "the Hadamard deviation may be used to
  reject linear frequency drift when a stability analysis is performed"
  [@cite RileyHowe2008]. HTDEV carries that benefit into the
  time-domain.
- **Wider noise-type convergence.** The Hadamard family converges
  over `α ∈ {−4, −3}` (frequency walk-walk and random-run FM) where
  the Allan family diverges [@cite Greenhall1997]. HTDEV is the
  time-domain extension of that range.

For records dominated by white-PM noise without drift, TDEV is fine
and slightly more efficient (tighter CI per τ on shorter records).
For records with drift, divergent low-frequency noise, or both, HTDEV
gives a numerically valid time-stability summary where TDEV would be
contaminated.

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

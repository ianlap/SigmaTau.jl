# The Hadamard time deviation

[`htdev`](@ref) is defined in this package:
[SP1065](@cite riley-2008-sp1065) and
[IEEE Std 1139-2022](@cite ieee1139-2022-definitions) give a time deviation
for the modified Allan deviation only. This page derives the factor
``\sqrt{10/3}`` of its definition from the rule that gives TDEV its
``\sqrt{3}``, and shows what the two time deviations measure as errors of
prediction.

## Differences of averages

Let ``\bar x(t)`` be the phase averaged over ``τ = mτ_0``, the mean of
``m`` consecutive samples starting at ``t``, with second and third
differences at spacing ``τ``

```math
\begin{aligned}
Δ^2\bar x(t) &= \bar x(t+2τ) - 2\bar x(t+τ) + \bar x(t), \\
Δ^3\bar x(t) &= \bar x(t+3τ) - 3\bar x(t+2τ) + 3\bar x(t+τ) - \bar x(t).
\end{aligned}
```

The modified Allan variance sums the second difference of single samples
over ``m`` adjacent starts
([SP1065 §5.2.5, Eq. (14)](@cite riley-2008-sp1065)), and that sum divided
by ``m`` is ``Δ^2\bar x``. The modified Hadamard variance does the same
with the third difference
([SP1065 §5.2.10, Eq. (23)](@cite riley-2008-sp1065)). With
``\langle\cdot\rangle`` the mean over every start the record allows, the
estimators of [`mdev`](@ref) and [`mhdev`](@ref) are therefore

```math
\mathrm{Mod}\,σ_y^2(τ) = \frac{\langle (Δ^2\bar x)^2 \rangle}{2τ^2},
\qquad
\mathrm{Mod}\,σ_H^2(τ) = \frac{\langle (Δ^3\bar x)^2 \rangle}{6τ^2} .
```

## The factor 10/3

The time variance is ``σ_x^2(τ) = τ^2\,\mathrm{Mod}\,σ_y^2(τ)/3``
([SP1065 §5.2.6, Eq. (15)](@cite riley-2008-sp1065)), which is
``\langle (Δ^2\bar x)^2 \rangle / 6``. The 3 comes from white PM. With
variance ``σ^2`` per sample, averages that do not overlap are independent
with variance ``σ^2/m``, and the coefficients ``(1, -2, 1)`` of the second
difference have squares that sum to 6:

```math
\langle (Δ^2\bar x)^2 \rangle = \frac{6σ^2}{m}, \qquad
\mathrm{Mod}\,σ_y^2(τ) = \frac{3σ^2}{mτ^2}, \qquad
σ_x^2(τ) = \frac{σ^2}{m} .
```

The 3 is what the variance's own 2 leaves of the 6. Dividing it out makes
TDEV on white PM the standard deviation of the ``τ``-average of the phase
([SP1065 §5.2.6](@cite riley-2008-sp1065)).

The third difference has coefficients ``(-1, 3, -3, 1)``, whose squares sum
to 20:

```math
\langle (Δ^3\bar x)^2 \rangle = \frac{20σ^2}{m}, \qquad
\mathrm{Mod}\,σ_H^2(τ) = \frac{10}{3}\,\frac{σ^2}{mτ^2} .
```

So 10/3 takes the place of TDEV's 3, and the Hadamard time variance is

```math
σ_{x,H}^2(τ) = \frac{τ^2\,\mathrm{Mod}\,σ_H^2(τ)}{10/3}
             = \frac{\langle (Δ^3\bar x)^2 \rangle}{20} .
```

HTDEV is then ``σ/\sqrt{m}`` on white PM, as TDEV is, and the two can be
read on one plot. On a drifting clock they part: a linear frequency drift
is a quadratic in the phase, which the third difference removes and the
second does not.

The same rule holds at any difference order ``d``. The variance divides by
``d!`` ([SP1065 §5.2.10, Eq. (24)](@cite riley-2008-sp1065)), the squared
binomial coefficients sum to ``\binom{2d}{d}``, and the time form divides
by ``\binom{2d}{d}/d!``, which is 3 at ``d = 2`` and 10/3 at ``d = 3``.

## Two prediction errors

Each difference is a value less a prediction of it:

```math
\begin{aligned}
Δ^2\bar x(t) &= \bar x(t+2τ) - \bigl[\,2\bar x(t+τ) - \bar x(t)\,\bigr], \\
Δ^3\bar x(t) &= \bar x(t+3τ) - \bigl[\,3\bar x(t+2τ) - 3\bar x(t+τ) + \bar x(t)\,\bigr].
\end{aligned}
```

The first bracket is the line through two successive averages, carried one
step on. The second is the parabola through three. With the two time
variances above,

```math
\mathrm{rms}\,Δ^2\bar x = \sqrt{6}\;σ_x(τ), \qquad
\mathrm{rms}\,Δ^3\bar x = \sqrt{20}\;σ_{x,H}(τ) = \sqrt{6}\;τ\,\mathrm{Mod}\,σ_H(τ) .
```

``\sqrt{6}`` TDEV is the RMS error of predicting the next ``τ``-average of
the phase with the frequency learned from the two before it, and
``\sqrt{20}`` HTDEV the error with the frequency and the drift learned from
the three before it. The second needs no factor beyond the modified
Hadamard variance's own 6. [Holdover and the time deviation](@ref)
compares both with the measured holdover test, which predicts single
readings from a least-squares fit.

## In numbers

On white PM of 1 ns per sample, both time deviations follow
``σ/\sqrt{m}``, here in nanoseconds:

```@example htdev
using SigmaTau, Random

p = PhaseData(1e-9 .* randn(Xoshiro(1), 100_000), 1.0)
m = [1, 4, 16, 64]
round.(1e9 .* [1e-9 ./ sqrt.(m) tdev(p, m; ci = false).dev htdev(p, m; ci = false).dev];
       digits = 4)
```

The third difference of the averages, over ``\sqrt{20}``, is HTDEV on any
record:

```@example htdev
x̄ = [sum(p.x[j:j+15]) / 16 for j in 1:length(p.x)-15]
Δ³ = x̄[49:end] .- 3 .* x̄[33:end-16] .+ 3 .* x̄[17:end-32] .- x̄[1:end-48]
sqrt(sum(abs2, Δ³) / 20length(Δ³)) ≈ htdev(p, [16]; ci = false).dev[1]
```

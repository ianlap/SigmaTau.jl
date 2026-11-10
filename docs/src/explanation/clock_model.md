# The clock model

The Kalman filters and timescales describe a clock by a small state driven
by white noise [zucca-2005-clock-model](@cite). With phase x, frequency y and
drift d,

```math
\dot x = y + \xi_1, \qquad \dot y = d + \xi_2, \qquad \dot d = \xi_3,
```

where ξ₁, ξ₂, ξ₃ are independent white noises with two-sided spectral
densities q₁, q₂, q₃. The reading is x plus white phase noise of variance R
per sample. Without drift, d and ξ₃ are dropped.

## From the h_α

Each drive gives one power law in the one-sided spectrum
S_y(f) = Σ h_α f^α. White FM is ξ₁ itself, so h₀ = 2q₁. Random-walk FM is
ξ₂ integrated once, S_y = 2q₂/(2πf)², so h₋₂ = q₂/2π². Random-run FM is ξ₃
integrated twice, so h₋₄ = q₃/8π⁴. White PM sampled every τ₀ has
S_x = 2τ₀R up to f_h = 1/2τ₀, so h₂ = 8π²τ₀R. [`clock_model`](@ref) inverts
these:

| noise | α | model |
|:------|:-:|:------|
| white PM | 2 | R = h₂ / 8π²τ₀ |
| white FM | 0 | q₁ = h₀ / 2 |
| random-walk FM | −2 | q₂ = 2π² h₋₂ |
| random-run FM | −4 | q₃ = 8π⁴ h₋₄ |

Flicker noises (α = 1, −1, −3) have a spectrum no finite state produces.
[`clock_model`](@ref) leaves flicker FM and flicker walk out;
`powerlaw_sigma(model, τ)` next to `powerlaw_sigma(fit, τ)` shows how far
that moves σ(τ).

## Flicker PM

Flicker PM is reading noise that is not white. A filter whose R is the
white PM alone takes its readings for quieter than they are: on simulated
clocks whose PM is nine tenths flicker, its innovations have 7 to 20 times
the variance it expects and its phase estimate is 2 to 3 times noisier at
τ₀. So [`clock_model`](@ref) folds the flicker PM into R. A filter averages
the PM that is faster than its own response and follows the rest as if it
were the clock, and R stands for the part it averages: the measurement
variance of [van-dierendonck-1984-allan-kalman](@cite), the phase spectrum
above the tracker's corner frequency. With σ²_PM(τ) the Allan variance of
the white and flicker PM together, and m the last lag at which it exceeds
that of the FM terms (in the Hadamard variance when the model has flicker
walk or random-run FM, and at most 65 536 samples),

```math
R = \frac{(m\tau_0)^2\,\sigma^2_{\mathrm{PM}}(m\tau_0)}{3}.
```

For white PM alone σ² = 3R/τ² at every lag, so this is h₂/8π²τ₀ again. For
flicker PM it grows with the logarithm of m: it is close to the variance of
the flicker PM between f = 1/(9mτ₀) and the Nyquist frequency.

No white R makes the filter right about flicker PM. Its slow part is still
followed as clock phase, the innovations stay correlated, and the phase
uncertainty the filter reports is too small. Read
[`innovation_check`](@ref) with that in mind.

The fit has to carry the flicker PM term (α = 1) for this to work. Fitted
with white PM, white FM and random-walk FM only, flicker PM passes for white
FM, and q₁ comes out several times too large. The exact response of a
sampled record separates white from flicker PM in the Allan deviation, from
m = 1; how the two share the level at τ₀ depends on the flicker's spectrum
near the Nyquist frequency, their sum does not.

## Over one step

Over a step of δ seconds the state moves by
X(t + δ) = Φ(δ) X(t) + w, with

```math
\Phi(\delta) = \begin{pmatrix} 1 & \delta & \delta^2/2 \\ 0 & 1 & \delta \\ 0 & 0 & 1 \end{pmatrix}
```

and w the drives integrated over the step, ∫₀^δ Φ(δ − s) ξ(s) ds. Its
covariance Q(δ) = ∫₀^δ Φ(s) diag(q₁, q₂, q₃) Φ(s)ᵀ ds is
[greenhall-2003-jones-tryon](@cite)

```math
Q(\delta) = \begin{pmatrix}
q_1\delta + q_2\delta^3/3 + q_3\delta^5/20 & q_2\delta^2/2 + q_3\delta^4/8 & q_3\delta^3/6 \\
q_2\delta^2/2 + q_3\delta^4/8 & q_2\delta + q_3\delta^3/3 & q_3\delta^2/2 \\
q_3\delta^3/6 & q_3\delta^2/2 & q_3\delta
\end{pmatrix}.
```

Because the drives are white, two steps compose into one:
Q(a + b) = Φ(b) Q(a) Φ(b)ᵀ + Q(b). The filters use this to step over gaps
and irregular sampling with one call to [`process_noise`](@ref).

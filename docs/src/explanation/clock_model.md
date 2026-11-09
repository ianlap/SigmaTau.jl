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
[`clock_model`](@ref) leaves them out; `powerlaw_sigma(model, τ)` next to
`powerlaw_sigma(fit, τ)` shows how far that moves σ(τ).

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

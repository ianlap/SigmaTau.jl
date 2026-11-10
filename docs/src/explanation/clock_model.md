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
By default [`clock_model`](@ref) leaves flicker FM and flicker walk out and
folds flicker PM into R; `powerlaw_sigma(model, τ)` next to
`powerlaw_sigma(fit, τ)` shows how far that moves σ(τ). With
`flicker = true` it keeps flicker FM and flicker PM as the Markov states
described [below](@ref "Flicker states").

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
[`innovation_check`](@ref) with that in mind, or give the flicker PM
[states](@ref "Flicker states"): R is then the white PM alone.

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

## Flicker states

A finite state cannot make a 1/f spectrum, but a few single-pole processes
can follow one over the span that matters. Davis, Greenhall and Stacey
[davis-greenhall-stacey-2005-flicker](@cite) model flicker FM as a sum of
independent Markov frequency processes mⱼ that add to the frequency,

```math
\dot x = y + \sum_j m_j + \xi_1, \qquad \dot m_j = -m_j / T_j + \eta_j,
```

each stationary with the same variance σ², and with time constants Tⱼ that
rise geometrically (their eqs 2.1 and 2.10). One process has the one-sided
spectrum 4σ²T/(1 + (2πfT)²) and the Allan variance σ² g(x)/x², with
x = τ/T and g(x) = 2x − 3 + 4e⁻ˣ − e⁻²ˣ (their eq. 2.5), which peaks at
τ = 1.89 T. Over a step, with a = e^(−δ/T), the state decays by a and moves
the phase by T(1 − a) times itself (eq. 2.7), and it adds

```math
Q_{11} \mathrel{+}= \sigma^2 T^2 g(\delta/T), \qquad Q_{1m} = \sigma^2 T (1 - a)^2, \qquad Q_{mm} = \sigma^2 (1 - a^2)
```

to the process noise (eqs 2.8–2.9, written there with σ_M² = 2σ²/T). The
filter starts each state at zero with variance σ².

The paper spaces the time constants by 8 or 16 and fits the levels to the
measured Hadamard deviation. Here they are a decade apart and the level
comes from h₋₁: a ladder with ratio ρ has the spectrum σ²/(f ln ρ) inside
its span, so S_y = h₋₁/f needs

```math
\sigma^2 = h_{-1} \ln 10 ,
```

with a ripple of ±2/cosh(π²/ln ρ), 5.5 % in the spectrum at a decade. In
σ(τ) the sum stays within about 3 % of flicker FM across the span. The same
relation gives the level and the ripple of the paper's figure 4.

Flicker PM gets the same construction one integration up, which is this
package's own: no published clock filter carries it. Markov phase processes
pⱼ of variance σ² = h₁ ln 10 / 4π² add to each reading and not to the phase,

```math
z = x + \sum_j p_j + v ,
```

so flicker PM stays measurement noise, as white PM is, but with memory. The
filter estimates the pⱼ and keeps them out of the phase, and in a timescale
they leave the weights to the clocks' own noise.

The states sit only where a flicker matters: across the τ, from τ₀ to
`tau_max`, at which it holds at least a quarter of the model's Allan
variance, a decade apart from the first such τ to at least the last. A
flicker that never reaches a quarter gets no states. What the states leave
out at τ₀ becomes white noise: white FM that brings the model's Allan
variance at τ₀ to the 2h₋₁ of sampled flicker FM, and white PM that does
the same for flicker PM. `fm_band` and `pm_band` place the states by hand.

In simulation, on a record with a flicker floor between white and
random-walk FM, the white-only fit predicts a third worse at 1000 s and a
half worse or more at 10 000 s than the model with states, and claims up to
four times the error it makes. Narrowing the states from the whole record to the
quarter band leaves the predictions where they were; the innovations keep a
small correlation from the flicker below the band. With flicker PM states
the innovations have the variance the filter claims, where the fold into R
claims nearly twice it, and the next reading is predicted 10 to 20 % better.
In a scale of three equal clocks with flicker PM on two channels, the
states leave the weights equal and give a scale 10 to 15 % quieter at 16 to
64 s than the fold, 15 % noisier at 1 s, and the same from 1000 s. The levels are
only as good as the fit: at the longest τ of a short record, flicker FM and
random-walk FM are hard to tell apart.

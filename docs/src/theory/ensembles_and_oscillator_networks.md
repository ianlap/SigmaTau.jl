# Theory: Time-Scale Algorithms and Oscillator Networks

A clock ensemble combines multiple physical oscillators into a single
realised timescale that is more uniform than any of its members. Two
paradigms appear in the literature:

1. **Globally coupled** schemes — the time-scale equation
   ([Stein 2003](@cite stein-2003-timescales)) and Kalman-filter
   ensembles ([Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation),
   [Breakiron 2001](@cite breakiron-2001-kalman-timescales)) — in
   which every clock's correction depends on a global quantity (the
   ensemble mean, or the covariance of the full state vector). These
   require pairwise difference measurements among all clocks at each
   epoch.
2. **Locally coupled** schemes — Kuramoto coupling and
   nearest-neighbour averaging
   ([Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto)) —
   in which each clock applies a phase correction that depends only
   on its directly connected neighbours, and every clock is itself a
   physical realisation of the ensemble timescale. These are designed
   for SWaP-constrained PNT timekeeping where global state-tracking
   across thousands of clocks is prohibitive.

This page covers the basic time-scale equation and its closure
constraints, the three- and N-cornered hat for separating individual-
clock variances, dynamic clock weights for picking one solution out of
the ambiguous family, the Kuramoto scheme and its EWFA baseline,
nearest-neighbour coupling, constellation fragmentation /
defragmentation, and three machine-learning extensions: ANN-learned
ensemble weights, ML-based clock-bias forecasting, and
telemetry-based stability estimation.

## Basic time-scale equation

A time-scale algorithm estimates the time error of each clock in an
ensemble; the corrected ensemble time is more uniform than the time
of any individual clock
[Stein 2003](@cite stein-2003-timescales). Because only pairwise
difference measurements are available — only `z_{ij}(t_k) =
x_i(t_k) − x_j(t_k) + v(t_k)` is observable — individual clock
corrections are formally unobservable, and the algorithm picks one
solution from an infinite ambiguous family by imposing closure
constraints
[Stein 2003](@cite stein-2003-timescales). The standard closure is a
weighted-sum-zero on the estimated phase, frequency, and
frequency-aging shocks across the ensemble:

```math
\sum_{i=1}^{N} a_i(t_k)\,\hat\varepsilon_i(t_k) \;=\; 0, \qquad
\sum_{i=1}^{N} b_i(t_k)\,\hat\eta_i(t_k) \;=\; 0, \qquad
\sum_{i=1}^{N} c_i(t_k)\,\hat\alpha_i(t_k) \;=\; 0,
```

with weights $a_i, b_i, c_i$ chosen to reflect each clock's noise
level
[Stein 2003](@cite stein-2003-timescales). The closure constraints
satisfy a limit theorem: the weighted sum of random shocks converges
to zero in the infinite-ensemble limit
[Stein 2003](@cite stein-2003-timescales).

The phase noise of a clock is the sum of direct phase shocks plus
integrated frequency and frequency-aging shocks, introducing
additional ambiguities beyond simple voltage averaging
[Stein 2003](@cite stein-2003-timescales). The unexplained behaviour
of historical time-scale algorithms — NIST AT1, USNO maser ensembles,
PTB time-scale, BIPM Algos — originated from implicit, ad-hoc choices
of one solution among the infinite family, not from algorithmic
idiosyncrasies
[Stein 2003](@cite stein-2003-timescales).

The same observability obstruction motivates Tryon and Jones's NBS
Kalman ensemble, which constrains drifts to sum to zero across the
seven-clock ensemble: because clock readings are differential, drifts
must be constrained to sum to zero in the ensemble fit, leaving
common-mode drift unobservable
[Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation).
USNO maser-ensemble Kalman timescales follow the same pattern,
realised through inverse-variance weights with an upper bound on any
single clock's contribution
[Breakiron 2001](@cite breakiron-2001-kalman-timescales). The NIST
AT1 family more generally couples a per-clock weight choice to a
target stability metric
[Sullivan, Allan, Howe & Walls 1990](@cite sullivan-1990-tn1337).

## Three-cornered hat

The three-cornered hat is the classical technique for separating
individual-clock variances from a triple of pairwise frequency-
stability measurements when no fourth, more-stable reference is
available
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability). It is
the `M = 3` instance of a more general M-cornered closure that
recovers each clock's variance from all `\binom{M}{2}` pairwise
variances under the assumption that the clocks are statistically
independent (see [N-cornered hat](#n-cornered-hat) below).

Frequency stability cannot be assessed with a lone oscillator using
classical phase techniques; verification requires pairwise comparisons
or three-cornered-hat comparisons among `≥ 3` clocks
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability). The NIST
SP1065 analysis workflow places the three-cornered hat as a late
stage: data precision check, gap/outlier/jump preprocessing, drift
analysis, variance analysis, spectral analysis, outlier recognition,
plotting, variance selection, then the three-cornered hat for
separating the oscillator under test from the reference
[Riley & Howe 2008](@cite riley-2008-sp1065).

## N-cornered hat

The N-cornered hat generalises the three-cornered hat to `M > 3`
ensemble clocks: each clock's variance is recovered from the full set
of pairwise variances by inverting an `M`-clock closure that assumes
statistical independence of the clocks
[Nandita et al. 2020](@cite nandita-2020-annensemble):

```math
\sigma_i^{2} \;=\; \frac{1}{M-2}\Bigl[\sum_{j=1}^{M}\sigma_{ij}^{2} - B\Bigr],
\qquad B \;=\; \frac{1}{2(M-1)}\sum_{k=1}^{M}\sum_{j=1}^{M}\sigma_{kj}^{2}.
```

The closure subtracts a normalised total `B` that captures the
cross-coupling and divides by `M − 2`, recovering `σ_i²` given the
full pairwise input set
[Nandita et al. 2020](@cite nandita-2020-annensemble).

The four deviation kernels typically used as pairwise inputs in this
pipeline are overlapping Allan deviation, modified Allan deviation,
Hadamard deviation, and overlapping Hadamard deviation
[Nandita et al. 2020](@cite nandita-2020-annensemble). Standard
statistical deviation is unsuitable for clock data because the
underlying noise processes are non-stationary; M-sample deviations
such as ADEV and HDEV are required instead
[Nandita et al. 2020](@cite nandita-2020-annensemble).

## Dynamic clock weights

Dynamic clock weights `a_i(t_k), b_i(t_k), c_i(t_k)` are the
time-varying coefficients applied in the closure constraints of a
Kalman-based time-scale algorithm; their choice picks one solution out
of the infinite ambiguous family that the under-determined ensemble
problem admits
[Stein 2003](@cite stein-2003-timescales). Selecting a single weight
family that matches one noise type produces a time-scale that tracks
either the short-term-best clock or the long-term-best clock
[Stein 2003](@cite stein-2003-timescales).

When both phase and frequency weights match the random-walk-phase-
noise levels, the two-clock time-scale tracks the better short-term
clock and improves performance by about 5 % short-term and about
10 % long-term. When both weight families instead match the
random-walk-frequency-noise levels, the time-scale tracks the better
long-term clock; improvement is about 18 % short-term and about 6 %
long-term
[Stein 2003](@cite stein-2003-timescales). **Splitting** the weight
choice — phase weights from random-walk-phase-noise levels and
frequency weights from random-walk-frequency-noise levels — yields a
time-scale 3 % better than the best clock at short `τ` and 8 %
better at long `τ`
[Stein 2003](@cite stein-2003-timescales).

In the USNO maser ensemble the same principle motivates an upper limit
on any single clock's weight to keep real-time operation robust
[Breakiron 2001](@cite breakiron-2001-kalman-timescales). Inverse-
variance weighting in a Kalman ensemble mean systematically
underestimates clock variance because each clock contributes to the
ensemble it is compared to; the bias is corrected by
`σ_u² = σ² / (1 − w)` or by referencing each clock to a mean of the
rest
[Breakiron 2001](@cite breakiron-2001-kalman-timescales).

## Kuramoto clock synchronisation

Kuramoto clock synchronisation steers each clock in an ensemble by a
phase correction proportional to the sine of its phase difference with
neighbouring clocks, applied as an additional term in the underlying
two-state Tavella stochastic clock model
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto). The
generic Kuramoto correction applied by clock `q` over its `n`
neighbours is

```math
\Delta\phi_{q,\text{corr}} \;=\; \sum_{j=1}^{n} a_{q,j}\,\sin(\phi_j - \phi_q), \qquad q \neq j,\; n \le N,
```

with coupling strengths `a_{q,j}`. In the small-phase-difference
linearisation with equal weights `a_{q,j} = 1/n` this reduces to the
arithmetic mean-deviation form

```math
\Delta\phi_{q,\text{corr}} \;=\; \frac{1}{n}\sum_{j=1}^{n}(\phi_j - \phi_q).
```

The two-state SDE with the Kuramoto correction term added to clock
`i` is

```math
\begin{aligned}
dX_{1,i}(t) &= X_{2,i}(t)\,dt + \sigma_{1,i}\,dW_{1,i}(t) + \sum_{j=1}^{n} a_{i,j}\,\sin(\phi_j - \phi_i)\,dt, \\
dX_{2,i}(t) &= \sigma_{2,i}\,dW_{2,i}(t).
\end{aligned}
```

Unlike conventional ensemble algorithms (Kalman variants, NIST
AT1/AT2, PTB time-scale, Leader-Follower), it does **not** require
global information flow: each clock only needs the present-epoch phase
of its nearest neighbours, and every member of the ensemble is itself
a physical realisation of the ensemble timescale
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto). The
motivating application is SWaP-constrained PNT timekeeping for large
cislunar / Martian / pLEO satellite constellations where global
state-tracking across thousands of clocks is prohibitive
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

After roughly 100 s of averaging, every clock in a 5-clock Kuramoto
chain achieves an MDEV equivalent to Equal-Weights Frequency Averaging
(EWFA, below), with the expected `1/√N` improvement over a single
free-running PRS-10
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).
Monte-Carlo experiments show the 1σ time-error buildup of the
Kuramoto algorithm matches EWFA and is approximately `2×` better
than a free-running PRS-10
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto). The
simulated PRS-10 Rb noise model uses `S_y(f) = h₀ f⁰ + h₋₂ f⁻²` with
`h₀ = 8 × 10⁻¹²` (white FM) and `h₋₂ = 4 × 10⁻¹⁵` (random-walk FM)
measured in the authors' laboratory
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

The Kuramoto–EWFA equivalence is established empirically at the MDEV
level after roughly 100 s of averaging rather than via an analytic
stability proof
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

!!! note "Planned implementation"
    [`KuramotoOscillator`](@ref) is currently exported as a stub in
    SigmaTauEnsemble. The mathematical form above (the per-link
    correction, the linearisation, the SDE augmentation) is
    documented; no time-stepped chain simulator is wired up yet.

## Nearest-neighbour coupling

Nearest-neighbour coupling is the topology specialisation of the
Kuramoto scheme in which each clock applies its phase correction
using only the present-epoch phases of its directly connected
neighbours, rather than the global ensemble
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto). On
a linear chain with two nearest neighbours, the Kuramoto correction
reduces to the symmetric two-neighbour average about each interior
node:

```math
\Delta\phi_{q,\text{corr}} \;=\; \tfrac{1}{2}\bigl(\phi_{q+1} - \phi_q + \phi_{q-1} - \phi_q\bigr).
```

Other connection topologies — ring, barbell, complete graph, and
application-specific shapes — are admissible, with the caveat that
some topologies have multiple stable equilibria. The locality of the
data flow is what makes the scheme robust against changes in topology
mid-mission: clocks are mutually steered, so each ensemble member is
itself a physical realisation of the ensemble timescale even when the
topology changes
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

## EWFA baseline

Equal-Weights Frequency Averaging (EWFA) is the global-information-
flow comparison baseline used to benchmark the Kuramoto clock-
synchronisation scheme: each clock `q` applies a frequency correction
equal to the negative deviation of its frequency from the
ensemble-mean frequency:

```math
\Delta f_{q,\text{corr}} \;=\; -\,\Bigl(f_q - \tfrac{1}{N}\sum_{j=1}^{N} f_j\Bigr).
```

EWFA sums over all ensemble members, requiring every clock to know
the present frequencies of every other member, in contrast to
Kuramoto's nearest-neighbour sum
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).
Empirically the Kuramoto-coupled ensemble matches EWFA at the MDEV
level after roughly 100 s of averaging in a 5-clock chain, with the
expected `1/√N` improvement over a single free-running clock
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto). In
its simplest form the Kuramoto ensemble's frequency stability is
equivalent to EWFA but requires only nearest-neighbour data flow,
providing scalability advantages for large constellations
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

## Constellation fragmentation and defragmentation

Constellation fragmentation is the operational scenario in which
crosslinks between satellite clocks fail and the ensemble breaks into
disconnected subsets, each running as an independent Kuramoto fragment
until the crosslinks are re-established
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto). The
Kuramoto scheme handles fragmentation gracefully because each clock
only needs its current neighbours: a fragment connected to a
higher-stability lunar-ground reference clock starts to follow the
reference at `τ ≈ 300 s`, while an isolated fragment performs `√2`
better than a single member but cannot inherit the ground-reference
stability
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).
After defragmentation — the re-establishment of crosslinks — the
full Kuramoto-coupled constellation re-synchronises to the
ground-station clock within a few hundred seconds of averaging
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

By contrast, classical globally coupled algorithms (Kalman variants,
NIST AT1/AT2, PTB timescale, Leader-Follower) require statistically
independent member clocks and global information flow, so every
satellite must track every other member's correction history —
prohibitive for 1000+ satellite constellations and brittle under
fragmentation
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

## Kuramoto coupling strength and loop attack time

The per-link weight `a_{i,j}` in the Kuramoto correction term controls
the trade-off between synchronisation speed and steady-state stability
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

Define the Kuramoto **loop-attack time** `τ_K` as the smallest `τ` at
which the Kuramoto deviation `σ_{Kuramoto}(τ)` falls below the EWFA
TIE: `\mathrm{TIE}_{\text{EWFA}}(\tau) ≥ σ_{\text{Kuramoto}}(τ_K)`.
Simulations show that the loop-attack time decreases monotonically
with the coupling strength up to a regime where over-correction drives
oscillatory instability; coupling strengths roughly between `0.01`
and `1` produce stable synchronisation in the Ristoff–Kettering–
Camparo simulations, with strengths above ~`1` driving instability
through over-correction
[Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto).

## ANN ensemble timescale

An ANN ensemble timescale uses a small feed-forward perceptron per
averaging time to learn positive, sum-to-one clock weights from
observed deviations rather than relying on inverse-deviation or
Kalman recursion
[Nandita et al. 2020](@cite nandita-2020-annensemble). The training
inputs are the four deviation estimators listed under
[N-cornered hat](#n-cornered-hat) — overlapping Allan, modified
Allan, Hadamard, and overlapping Hadamard — each using their textbook
formulas
[Nandita et al. 2020](@cite nandita-2020-annensemble). The
constraints encode the requirement that a weighted-mean ensemble must
not introduce additive bias and that weights must remain non-negative
[Nandita et al. 2020](@cite nandita-2020-annensemble).

The N-cornered-hat closure supplies per-clock variances as inputs:

```math
\sigma_i^{2} \;=\; \frac{1}{M-2}\Bigl[\sum_{j=1}^{M}\sigma_{ij}^{2} - B\Bigr],
\qquad B \;=\; \frac{1}{2(M-1)}\sum_{k=1}^{M}\sum_{j=1}^{M}\sigma_{kj}^{2}.
```

The proposed IRNWT ensemble enforces positive weights summing to one
and zero bias, since a weighted-mean ensemble cannot have an additive
offset
[Nandita et al. 2020](@cite nandita-2020-annensemble). On the
demonstration dataset the best three-clock ensemble (using OHDEV
inputs) reaches `σ_y(1 s) = 2.99 × 10⁻¹³` and
`σ_y(1 day) = 6.99 × 10⁻¹⁶` with learned weights `A = 0.003,
C = 0.871, D = 0.125`; the reported long-term stability `~10⁻¹⁶` at
one day is roughly an order of magnitude better than reported Kalman
and METS implementations on the same dataset
[Nandita et al. 2020](@cite nandita-2020-annensemble). Each per-`τ`
perceptron trains in about 20 s over 1000 epochs (~20 ms/epoch)
[Nandita et al. 2020](@cite nandita-2020-annensemble). Beyond
`τ ≈ 10⁶ s` the ensemble curve diverges from its trend because
individual-clock frequency drift dominates and would need explicit
removal
[Nandita et al. 2020](@cite nandita-2020-annensemble).

## ML-based clock-bias forecasting

ML-based clock-bias forecasting predicts future GPS satellite clock
bias by training neural-network models on the residuals of a
parametric clock model rather than fitting a polynomial directly
[Song et al. 2025](@cite song-2025-mlclockbias). The bare quadratic-
polynomial clock-bias model is

```math
x \;=\; a_0 + a_1\, t + a_2\, t^{2} + \varepsilon,
```

with offset, drift, and drift-rate coefficients fit by least squares.
The QPMwPT pipeline removes outliers, fits the quadratic polynomial
trend, FFT-extracts dominant residual frequencies, and adds a periodic
correction term to give the stronger baseline

```math
x \;=\; a_0 + a_1\, t + a_2\, t^{2} + \sum_{i=1}^{p}\bigl(A_i \sin(2\pi f_i t) + B_i \cos(2\pi f_i t)\bigr) + \varepsilon,
```

with frequencies `f_i` extracted from the QP-residual sequence
[Song et al. 2025](@cite song-2025-mlclockbias). QPMwPT improves on
the bare QP by about 13.5 % in average 1-day prediction accuracy
[Song et al. 2025](@cite song-2025-mlclockbias). Average 1-day
prediction-accuracy improvements over QP across six test satellites
are 39.45 % (BPNN), 57.57 % (WNN), 27.28 % (LSTM), 29.14 % (GRU); the
wavelet neural network achieves the lowest MAE of the four
architectures
[Song et al. 2025](@cite song-2025-mlclockbias). All four
architectures share identical hyperparameters: 15-sample lookback
window, 10 hidden units, learning rate `5 × 10⁻⁴`, 1000 epochs,
single feature
[Song et al. 2025](@cite song-2025-mlclockbias).

Training data span 9 days (1–9 January 2024) of IGS precision
clock-bias data sampled every 5 minutes; the held-out test day is
10 January 2024
[Song et al. 2025](@cite song-2025-mlclockbias). The six evaluation
satellites (G1, G2, G4, G5, G10, G24) cover all five active GPS
atomic-clock types: Block IIR Rb, IIR-M Rb, III-A Rb, IIF Rb, IIF Cs
[Song et al. 2025](@cite song-2025-mlclockbias).

## Telemetry-based stability estimation

Telemetry-based stability estimation is a supervised-learning method
that estimates the frequency stability of a single oscillator using
only its internal telemetry — cavity temperature, supply voltages,
oven currents, and similar engineering channels — bypassing the
classical requirement of a pairwise phase comparison against a second
comparable clock
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability). The
regression target is the dynamic Allan deviation (DADEV) of a sliding
window of phase residuals, and the input features are the Allan
deviations of each telemetry channel computed over the same window
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability).

Pairwise phase comparison between two comparably stable oscillators
is the classical prerequisite for any frequency-stability assessment,
and a lone oscillator in isolation cannot be characterised this way
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability). Dynamic
Allan deviation `σ_y(t, τ)` extends the static ADEV by a time axis,
repeatedly computing static ADEV over sliding windows of length `T`
at advancement interval `Δt`
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability). After
removing linear drift, hydrogen-maser frequency residuals are
dominated by random-walk frequency noise associated with quasi-
independent processes inside the maser physics package
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability).

Each maser sample is arranged as input shape `(N, 21 channels, 24 τ)`
and target shape `(N, 24 ADEV values)`, with the test set held out
from a non-overlapping time period
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability). On the
oscillator-specific data set, a CNN achieves `R² = 0.85` and
`MAE = 5.61 × 10⁻¹⁶`, outperforming elastic net (`R² = 0.77`),
random forest (`R² = 0.78`), and the constant-mean baseline
(`R² = 0.73`)
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability). On the
oscillator-agnostic data set (trained on Masers A, B, C, D, X; tested
on Maser Y), the CNN achieves `R² = 0.75` and
`MAE = 6.53 × 10⁻¹⁶`, against a baseline `R²` of 0.28
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability).

Gradient-based CNN saliency identifies cavity-temperature-control
telemetry as dominant for `τ` from `~100` to `~10 000 s`; receiver
electronics (IF amplitude, cavity register, VCO) dominate at
`τ > 10 000 s`. Hydrogen masers exhibit a frequency-temperature
sensitivity of order `1 × 10⁻¹⁵` per °C, motivating the prominence of
thermal-control telemetry in the saliency maps
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability). The
technique provides only coarse stability estimates at present but is
potentially valuable for optical clocks where running multiple units
in parallel for pairwise comparison is impractical
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability).
Telemetry sampled at 1 sample per minute yields a minimum `τ` of
60 s for stability estimation; faster telemetry would extend
predictions to shorter `τ`
[McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability).

## See also

- [Theory: Clock State-Space Models](ensemble_overview.md) — per-clock
  SDE underlying both the time-scale equation and the Kuramoto
  augmentation.
- [Theory: Kalman Filter and Variants](kalman.md) — supplies the
  per-clock state estimates the closure constraints combine.
- [Theory: Clock Steering with PID Controllers](steering.md) — a
  steered ensemble realisation can use the same controller against
  the time-scale output.
- [Allan family](allan_family.md) — the deviation kernels used as
  pairwise inputs to the cornered-hat closures.
- [API: `SigmaTau.Est`](../reference/est.md) —
  [`KuramotoOscillator`](@ref).

## References

- [Sullivan et al. 1990](@cite sullivan-1990-tn1337) — NIST
  Technical Note 1337, primary-frequency-standard timescale context.
- [Stein 2003](@cite stein-2003-timescales) — basic time-scale
  equation and the closure-ambiguity unification.
- [Tryon & Jones 1983](@cite tryon-1983-cesium-parameter-estimation) —
  NBS Kalman-ensemble drift-sum-zero constraint.
- [Breakiron 2001](@cite breakiron-2001-kalman-timescales) — USNO
  maser-ensemble inverse-variance weighting and bias correction.
- [Riley & Howe 2008](@cite riley-2008-sp1065) — NIST SP1065
  three-cornered-hat workflow.
- [Nandita et al. 2020](@cite nandita-2020-annensemble) — IRNWT
  ANN-learned ensemble weights and N-cornered-hat closure.
- [Song et al. 2025](@cite song-2025-mlclockbias) — ML-based GPS
  clock-bias forecasting.
- [McKelvy et al. 2025](@cite mckelvy-2025-telemetrystability) —
  CNN-based telemetry-only frequency-stability estimation.
- [Ristoff, Kettering & Camparo 2026](@cite ristoff-2026-kuramoto) —
  Kuramoto clock synchronisation, EWFA baseline, fragmentation /
  defragmentation, loop-attack time.

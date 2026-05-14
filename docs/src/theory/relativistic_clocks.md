# Theory: Relativistic Clock Models

A relativistic clock model corrects an idealised atomic clock for
general-relativistic and special-relativistic rate effects so that
its ticks can be related to a chosen coordinate timescale
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time)
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt). For
**terrestrial** timekeeping the IAU TCB / TCG transformations and
the Sagnac correction already handle the operationally significant
effects, validated in orbit by ACES at the $10^{-17}$ level
[Cacciapuoti & Salomon 2009](@cite cacciapuoti-2009-aces-space-clocks).
For **cislunar PNT** — a lunar-surface clock, an Earth-Moon-Lagrangian
relay, or a frozen-eccentricity lunar orbiter — the same machinery has
to be carried into a Moon-centred frame, which introduces TCL and
the TL fixed-rate constant alongside the familiar TCB–TCG pair
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales)
[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).

This page is the **hub** for the relativistic-PNT documentation. It
states the model in broad strokes, anchors the operational headline
numbers, and points to three companion pages where the mathematical
machinery is set out in full:

- [Theory: Relativistic Frames and Time Scales](relativistic_frames_and_timescales.md)
  — BCRS / GCRS / LCRS reference systems and the six time scales
  TT / TCG / TCB / TDB / TCL / TL with their defining constants and
  conversions.
- [Theory: Relativistic Clock Corrections](relativistic_corrections.md)
  — 1PN proper-time mapping, gravitational redshift on the Moon,
  $L_{Gm}$ Earth–Moon rate constant and Lagrange-point offsets,
  cislunar orbit drift regimes (vLLO / LLO / ELFO / L1 / NRHO), and
  the Shapiro / Sagnac light-time corrections.
- [Theory: Lunar PNT Systems](lunar_pnt_systems.md) — two-way time
  transfer (synchronous and asynchronous) in the ESA Moonlight LCNS
  architecture, and the relativistic positioning system (emission
  coordinates / Autonomous Basis of Coordinates) tradition.

## Operational headline numbers

The few numbers worth memorising, all from the same multi-source
ingest:

- A clock on the Moon's equator ticks **faster** than one on Earth's
  geoid by `56.0199(12) − 0.108 434 17(89) cos f` µs/day, with $f$
  the true anomaly of the Moon's orbit
  [Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).
- A static observer on the lunar surface ticks **slower** than TCG
  by `−2.7119 µs/day`
  [Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).
- An Earth–Moon L1 clock leads an Earth-geoid clock by
  `58.612 420(12) − 0.107 361 06(12) cos f` µs/day; the dominant
  60.2 µs/day term is just the depth of L1 in Earth's potential
  [Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).
- An ELFO satellite at $e \approx 0.69$ drifts 58.1152 µs/day vs TT,
  with periodic terms $\{0.115, 0.040, 0.018\}$ µs at the orbital
  frequency and its first two harmonics
  [Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).
- An NRHO (Gateway-like, $e = 0.9088$, period 7.49 d) drifts
  58.5431 µs/day vs TT, with a TT-mapping geometry term up to 0.81 µs
  near apoapsis
  [Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

## Operational design pattern

Across the lunar PNT literature the recommended stack is the same
[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic):

1. **Internal bookkeeping in TDB.** Internal ephemeris dynamics use
   BCRS / TDB so the dynamics core is decoupled from later
   TL-rate-constant choices.
2. **Deterministic relativistic correction first.** Apply the 1PN
   proper-time correction (and Shapiro / Sagnac when relevant)
   *before* the stochastic clock-state estimator runs. Failing to do
   this absorbs the relativistic signature into the clock estimate
   as if it were stochastic noise, degrading both estimation
   accuracy and time-dissemination fidelity.
3. **TL mapping at the API layer.** The choice of $L_L$, $TL_0$,
   and any TL standardisation lives at the broadcast / mapping
   interface, not in the dynamics or estimation code. SigmaTau's
   future implementation should parameterise these constants
   accordingly.

## Three look-alike $\vec v \cdot \vec R / c^{2}$ terms

Three terms in the cislunar relativistic toolbox look superficially
identical but refer to different velocities and reference frames, and
must not be conflated
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt):

1. The $\vec v_L \cdot \vec r / c^{2}$ term in the **TCB–TCL**
   integrand (and its analogue $\vec v_e \cdot \vec R / c^{2}$ in
   TCB–TCG).
2. The satellite proper-time periodic term $2\,\vec R \cdot \vec V / c^{2}$
   in the Keplerian-orbit proper-time formula.
3. The Sagnac path correction $2\,\vec\Omega \cdot \vec A / c^{2}$,
   which is a path-area integral in a rotating frame rather than a
   point-velocity scalar.

Concept and implementation notes — and the eventual unit-test harness
for [`RelativisticClock`](@ref) — must keep these three terms
independently named.

## SigmaTau implementation

The relativistic clock corrections are surfaced in the ensemble
subpackage as [`RelativisticClock`](@ref), an `AbstractClockModel`
intended to wrap an underlying SDE clock model and apply the
proper-time differential at the boundary between the SDE state and
the chosen coordinate timescale.

!!! note "Planned implementation"
    [`RelativisticClock`](@ref) is currently exported as a stub in
    `SigmaTau.Est`. The mathematical form of every correction —
    1PN proper-time differential, TCB / TCG / TCL transformations,
    surface redshift, Sagnac, Shapiro, ELFO orbit formulas — is
    documented across this page and the three companion pages above;
    no numerical implementation of the corrections is wired up yet,
    and no DE440 / INPOP21a ephemeris reader is present. The
    near-term recommended implementation path is the closed-form
    Keplerian-Orbital form described in
    [Theory: Relativistic Clock Corrections](relativistic_corrections.md)
    (per [Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt) §6),
    which avoids any ephemeris dependency. A DE440 path may
    eventually justify a separate SigmaTauLunar subpackage; the
    closed-form path does not. The three velocity-dot-radius warning
    above should be carried into the eventual implementation as a
    unit-test harness rather than left as documentation alone.

## See also

- [Theory: Clock State-Space Models](ensemble_overview.md) — the
  underlying stochastic clock model that the relativistic correction
  post-processes.
- [Theory: Kalman Filter and Variants](kalman.md) — Kalman recursion
  for the noise-driven phase residuals downstream of the proper-time
  transformation.
- [Theory: Single-Clock Steering](steering.md) — the PID controller
  that closes the loop on a disciplined clock.
- [Allan family](allan_family.md) — verification of the
  post-correction stability with $\sigma_y(\tau)$ slope identification.
- [API: `SigmaTau.Est`](../reference/est.md) —
  [`RelativisticClock`](@ref).

## Primary sources

The relativistic-PNT pages draw on a multi-source ingest of recent
lunar-PNT and space-clock literature:

- [Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time) —
  *A Relativistic Framework to Establish Coordinate Time on the Moon
  and Beyond.* The proof-of-concept reference; introduces $L_{Gm}$ and
  the headline 56.02 µs/day Earth–Moon rate.
- [Turyshev 2025](@cite turyshev-2025-cislunar-time-scales) —
  *High-Precision Relativistic Time Scales for Cislunar Navigation.*
  The most comprehensive worked-out realisation of IAU Resolution II
  (2024); orbit-by-orbit drift catalogue.
- [Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic)
  — *Relativistic Framework for Disseminating Lunar Time for the
  Intuitive Machines LCRNS PNT Constellation.* The deployment-side
  companion: TDB internal, TL broadcast, ELFO operational.
- [Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt) —
  *Relativistic Time Modeling for Lunar PNT* (master's thesis); the
  closed-form Keplerian-Orbital path and the three-look-alike-terms
  warning.
- [Iess et al. 2025](@cite iess-2025-cislunar-od-time-sync) and
  [Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer)
  — ESA Moonlight LCNS architecture (NAVI 2025 long form, PTTI 2026
  synchronous-vs-asynchronous follow-up).
- [Cacciapuoti & Salomon 2009](@cite cacciapuoti-2009-aces-space-clocks)
  — *Space Clocks and Fundamental Tests: The ACES Experiment.* The
  Earth-orbit precursor for in-orbit relativistic-clock validation.
- [Reinhardt, Hartwig & Heinzel 2024](@cite reinhardt-2024-lisa-clock-sync)
  — heliocentric BCRS clock synchronisation for LISA; the
  cross-cutting reference that demonstrates the same machinery
  outside the lunar context.
- [Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning) —
  the relativistic-positioning-system tradition; emission coordinates
  and Autonomous Basis of Coordinates as the conceptual foreground
  to the lunar work.

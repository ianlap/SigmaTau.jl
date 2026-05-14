# Theory: Relativistic Clock Corrections

The corrections in this page ride on top of the frames and time scales
collected in
[Theory: Relativistic Frames and Time Scales](relativistic_frames_and_timescales.md).
The leading 1PN proper-time mapping is the universal foundation
([Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time),
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales),
[Reinhardt, Hartwig & Heinzel 2024](@cite reinhardt-2024-lisa-clock-sync),
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt)); the
gravitational redshift on the Moon, the orbit-averaged Earth–Moon
rate constant $L_{Gm}$ with Earth–Moon Lagrange-point offsets, the
cislunar orbit drift regimes, and the Shapiro / Sagnac light-time
corrections specialise that framework to the sub-µs / 0.1 ps targets
of cislunar PNT.

## 1PN proper-time mapping

The first post-Newtonian proper-time mapping is the order-$1/c^{2}$
relation between the proper time $\tau$ recorded by a clock and the
coordinate time $T$ of a chosen post-Newtonian frame, parameterised
by the clock's gravitational potential $\Phi$ (or $U$) and squared
velocity $V^{2}$:

```math
d\tau \;=\; \Bigl(1 + \frac{\Phi}{c^{2}} - \frac{V^{2}}{2 c^{2}}\Bigr)\,dT + \mathcal{O}(c^{-4})
```

[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time)
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt). The two
contributions are physically distinct: the gravitational-redshift term
$\Phi/c^{2}$ causes deeper-potential clocks to tick slower; the
kinematic term $-V^{2}/(2c^{2})$ causes faster-moving clocks to tick
slower
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time)
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

The Schwarzschild-plus-kinematic specialisation for a heliocentric
spacecraft (Reinhardt LISA form) is

```math
\delta\dot\tau_i \;=\; -\,\frac{1}{2}\,\frac{|\vec v_i|^{2}}{c^{2}} - \frac{r_s}{2\,|\vec x_i|},
```

with $r_s$ the solar Schwarzschild radius
[Reinhardt, Hartwig & Heinzel 2024](@cite reinhardt-2024-lisa-clock-sync).
The same equation appears in the Ashby, Turyshev, and Seyffert
references up to notation; numerical $L_L$ estimates differ by 0.4–1.4
ns/d depending on the chosen lunar reference radius rather than from
any derivation conflict
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales)
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

2PN corrections (squared potential and gravitomagnetic cross terms)
are negligible for current cislunar clock accuracies and may be
dropped from operational implementations
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt). At
$5 \times 10^{-18}$ stability targets, however, Turyshev retains
explicit 2PN bookkeeping in the LCRS metric truncation
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

ACES on the ISS is the Earth-orbit operational validator for the 1PN
mapping at the $10^{-17}$ level, with a target gravitational
frequency-shift sensitivity of 2 ppm
[Cacciapuoti & Salomon 2009](@cite cacciapuoti-2009-aces-space-clocks).

!!! warning "Apply deterministic corrections before estimation"
    Failing to remove the deterministic 1PN signature before the
    stochastic clock-state estimator runs absorbs the relativistic
    signature into the clock estimate as if it were stochastic,
    degrading both estimation accuracy and dissemination fidelity
    [Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).

## Gravitational redshift on the Moon

For a static observer on the lunar surface, the redshift relative to
TCG, the corresponding fractional-frequency shift, and the
chronometric-levelling height resolution follow from the proper-time
differential applied at fixed `v = 0`:

```math
\frac{\Delta t}{t} \;=\; -\,\frac{V_{\mathrm{surface}}}{c^{2}}, \qquad
\frac{\Delta\nu}{\nu} \;=\; \frac{\Delta U}{c^{2}}, \qquad
\Delta h_b \;=\; \frac{\Delta U}{g_b}
```

[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt). The
Moon's TCG-companion rescaling constant in the Ashby form is

```math
L_m \;=\; -\,\frac{\Phi_{0m}}{c^{2}} \;=\; 3.138\,81(15) \times 10^{-11} \;\approx\; 2.71 \;\mu\mathrm{s/d},
```

referenced to the lunar selenoid
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time). This
is consistent with Turyshev's selenoid-derived
$L_L = 3.139\,05 \times 10^{-11}$ ≈ 2.7121 µs/d to within the
0.84–2.17 ns/d offsets between selenoid- and South-Pole-anchored
evaluations
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

Numerical operating points worth retaining:

- A static observer on the lunar surface ticks slower than TCG by
  `−2.7119 µs/day`
  [Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).
- Topographic potential variations across the lunar surface produce
  `±15 ns/day` of clock-rate variation, with a maximum 28.7 ns/day
  from low to high terrain
  [Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).
- Because lunar gravity is roughly 16.5 % of Earth's, the same
  fractional-frequency clock sensitivity yields about `6×` poorer
  chronometric-levelling height resolution on the Moon than on Earth
  [Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

ACES PHARAO + SHM target instabilities are
$1 \times 10^{-13}\,\tau^{-1/2}$ short-term combined with
$1.5 \times 10^{-15}$ at $10^{4}$ s long-term, supporting a 2-ppm
gravitational frequency-shift test in Earth orbit; the MWL link
targets a time deviation of 0.3 ps at 300 s and 23 ps at 10 days,
1–2 orders of magnitude better than TWSTFT or GPS time transfer
[Cacciapuoti & Salomon 2009](@cite cacciapuoti-2009-aces-space-clocks).

## $L_{Gm}$ — the Earth–Moon rate constant and Lagrange-point offsets

$L_{Gm}$ is the orbit-averaged 1PN rate constant between the lunar
selenoid and Earth's geoid, derived in the Ashby–Patla post-Newtonian
framework as the headline numerical signature of any cislunar
timekeeping system
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time). The
result is independent of the choice among barycentric,
centre-of-mass-rotating, and Earth-origin-rotating coordinate
systems; centrifugal-potential contributions cancel any apparent
frame dependence
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).

The fractional rate difference between a clock on the Moon's equator
and a clock on Earth's geoid (Ashby–Patla) is

```math
\frac{d\tau_m - d\tau_e}{d\tau_e} \;=\; \frac{G M_m - G M_e}{c^{2} D} + \frac{\Phi_{0m} - \Phi_0}{c^{2}}
       - (1 - 2\mu)\,\frac{G M_T \,(1 + 2 e \cos f + e^{2})}{2 a c^{2} (1 - e^{2})},
```

with $D$ the Earth–Moon distance, $f$ the true anomaly, and
$\mu = M_m / M_T$
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).
Numerically, this is $6.483\,78(15) \times 10^{-10} - 1.255\,025\,18(89) \times 10^{-12}\cos f$,
i.e. $56.0199(12) - 0.108\,434\,17(89)\cos f$ µs/day
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).

The $L_{Gm}$ rate constant (TCL relative to TCG) is

```math
L_{Gm} \;=\; \frac{G M_m - G M_e}{c^{2} D} - \frac{V_m^{2} - V_e^{2}}{2 c^{2}},
```

varying as $-1.49373 - 0.10967 \cos f$ µs/day; the Keplerian model
accumulates up to about 75 ns of error per lunar orbit relative to a
DE440 numerical integration
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).

Earth–Moon Lagrange-point clocks (L1, L2, L4 / L5) carry the same
redshift-plus-kinematic structure but evaluated at the equilibrium
points; their rate differences relative to Earth's surface are
dominated by being higher in Earth's potential (~60.2 µs/day from
$-\Phi_0/c^{2}$) plus orbit-specific corrections
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time):

| Point | Fractional rate difference vs Earth geoid (µs/day) |
|---|---|
| L1 | $58.612\,420(12) - 0.107\,361\,06(12)\cos f$ |
| L2 | $58.619\,639(12) - 0.124\,455\,90(12)\cos f$ |

L4 / L5 are equidistant from Earth and Moon; their rate difference
reduces to a sum of geopotential terms at $1/r_e$ and $1/r_m$ plus
second-order Doppler
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).

Three classes of position-dependent corrections were omitted from
Ashby's orbit average: an Earth-rotation term $\sim$ 0.0002 µs/day,
a Moon-orientation term $\sim$ 0.0045 µs/day, and an Earth-orbit /
spin cross term $\sim$ 0.0055 µs/day
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time). Tidal
effects (perturbed potential, perturbed position, perturbed
second-order Doppler) are not included in Ashby's analytical
estimates; residuals computed against DE440 reach a few parts in
$10^{13}$ over a lunar orbit
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).
Turyshev's L1-specific orbit drift (58.6182 µs/d) agrees with Ashby's
L1 result to within the operational-precision bookkeeping
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

## Cislunar proper-time drift regimes

[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales) catalogues
five canonical cislunar regimes — vLLO, LLO, ELFO, Earth–Moon L1,
NRHO — each with a characteristic secular rate $L_{CL}$ vs TT and a
periodic-term spectrum dominated by the orbital frequency and its
harmonics. The generic τ-vs-TT mapping is

```math
\tau - \mathrm{TT} \;=\; -\,L_{CL}\,(\tau - \tau_0) - P_{CL}(\tau) + (\text{LCRS geometry terms}).
```

| Regime | $L_{CL}$ | Drift vs TT | Periodic structure |
|---|---|---|---|
| vLLO (10 km) | $4.6818 \times 10^{-11}$ | 54.6926 µs/day | Lunar gravity needs $\ell_{\max} \gtrsim 300$ |
| LLO (100 km) | $4.4521 \times 10^{-11}$ | 54.8912 µs/day | $J_{2M}$ at 2.28 ps; $C_{22}$ at 0.46–0.50 ps |
| ELFO ($e = 0.6917$) | $7.2372 \times 10^{-12}$ | 58.1152 µs/day | K+M $\{0.115, 0.040, 0.018\}$ µs at $\{\omega, 2\omega, 3\omega\}$ |
| Earth–Moon L1 | $1.3827 \times 10^{-12}$ | 58.6182 µs/day | Monthly K+M 25.3 ns; geometry $\lesssim 36$ ns |
| NRHO ($e = 0.9088$, 7.49 d) | $2.2537 \times 10^{-12}$ | 58.5431 µs/day | $\{0.137, 0.062, 0.038\}$ µs harmonics; geometry up to 0.81 µs near apoapsis |

[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). Lunar
harmonics through $\ell = 9$ and external tides through $\ell = 8$
suffice for $5 \times 10^{-18}$ stability at the deep cislunar
regimes (L1, NRHO); near-surface and very low orbits require
$\ell_{\max} \gtrsim 300$
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

For the IM LDN ELFO ($e \approx 0.69$), the secular drift is
consistent with Turyshev's $L_{CL} = 7.24 \times 10^{-12}$, and
periodic terms at orbital frequency dominate at sub-µs amplitude
[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).

The Seyffert closed-form Keplerian periodic correction for a
satellite orbit is

```math
\Delta\tau \;=\; \Bigl(1 - \frac{3 G M}{2 c^{2} a}\Bigr)\,\Delta t \;-\; \frac{2}{c^{2}}\sqrt{G M\,a}\,e\,\sin E,
```

with $a$ the semi-major axis, $e$ the eccentricity, and $E$ the
eccentric anomaly
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt). In
simulation, the Cartesian-Orbital and Keplerian-Orbital proper-time
formulas agree at machine precision; the simpler Lander-Like form
differs by $\sim 10^{-7}$ s/year ($\sim 10^{-3}$ of the main signal)
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt). Most
circular lunar orbits are dynamically unstable (mascons / large
$C_{20}, C_{22}$ for low orbits, third-body Earth perturbations for
high orbits); ELFOs with $e \approx 0.6$ are the practical stable
design for navigation constellations
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

The choice between a closed-form Keplerian path (no ephemeris
dependency) and a DE440 numerical-integration path is an
architectural commitment for [`RelativisticClock`](@ref); the
closed-form path is the recommended near-term implementation in
`SigmaTau.Est` (per Seyffert §6) and avoids any ephemeris
dependency, while a DE440 path may eventually justify a separate
`SigmaTauLunar` subpackage
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

## Light-time corrections — Shapiro and Sagnac

Total Earth–Moon coordinate-time light propagation decomposes as

```math
\Delta t_{1 \to 2} \;=\; \frac{R_{12}}{c} + \sum_{B \in \{S, E, M\}}\Delta_B^{\mathrm{Sh}} + \Delta_{(1)}^{\mathrm{Sag}} + \mathcal{O}(c^{-4})
```

[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales), with the
Shapiro term per gravitating body $B$ given by the standard
Schwarzschild log form

```math
\Delta_B^{\mathrm{Sh}} \;=\; \frac{2 G M_B}{c^{3}} \ln\!\left(\frac{r_{1B} + r_{2B} + R_{12}}{r_{1B} + r_{2B} - R_{12}}\right)
```

[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). At the
0.1 ps target of an LCRS-class architecture, Shapiro contributions
from the Sun (~20–30 ns), Earth (~0.1–0.2 ns), and Moon (1–3 ps) are
all individually relevant
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

The first-order Sagnac correction in the GCRS is

```math
\Delta_{(1)}^{\mathrm{Sag}} \;=\; -\,\frac{\vec\Omega_{\oplus}}{c^{2}} \cdot (\vec r_2 \times \vec r_1)_{\mathrm{GCRS}}
```

[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). Seyffert
gives the equivalent closed-path form for a one-way signal:

```math
\Delta t_{\mathrm{Sagnac}} \;\approx\; \frac{2\,\vec\Omega \cdot \vec A}{c^{2}} \;=\; \frac{2 \omega A_z}{c^{2}} \;=\; \frac{2 \omega}{c^{2}}\bigl(x_r y_s - y_r x_s\bigr).
```

Equatorial Sagnac for a signal circumnavigating Earth is 207.4 ns;
a GPS satellite-to-equator path can reach 130 ns — about 26 m of
range error if uncorrected
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

The Sagnac term carries opposite signs for a signal travelling east
vs west around the rotating frame; closed-path topology determines
the sign by the area-vector convention. The closed-form Seyffert
Sagnac and the Turyshev GCRS Sagnac differ in sign convention by the
choice of area-vector orientation, so an implementation must commit
to a consistent orientation across the code base
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt)
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

## Three look-alike $\vec v \cdot \vec R / c^{2}$ terms

Three terms in the cislunar relativistic toolbox look superficially
identical but refer to different velocities and reference frames, and
must not be conflated
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt):

1. The $\vec v_L \cdot \vec r / c^{2}$ term in the **TCB–TCL**
   integrand (and its analogue $\vec v_e \cdot \vec R / c^{2}$ in
   TCB–TCG), where the velocity is the orbital velocity of the
   central body's barycentric motion and $\vec r$ (or $\vec R$) is
   the position of the clock in the local frame.
2. The satellite proper-time periodic term $2\,\vec R \cdot \vec V / c^{2}$
   in the Keplerian-orbit proper-time formula, where the velocity is
   the satellite's velocity relative to the central body and $\vec R$
   is the satellite's position vector in the same frame.
3. The Sagnac path correction $2\,\vec\Omega \cdot \vec A / c^{2}$,
   which is a path-area integral in a rotating frame rather than a
   point-velocity scalar.

Concept and implementation notes must not conflate these three
terms; they are best carried in software as separately named
variables and exercised in a unit-test harness that verifies their
independent meanings
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

## See also

- [Theory: Relativistic Frames and Time Scales](relativistic_frames_and_timescales.md)
  — BCRS / GCRS / LCRS and the time scales these corrections live in.
- [Theory: Relativistic Clock Models](relativistic_clocks.md) — hub
  page for [`RelativisticClock`](@ref) and the orbit-averaged
  Earth–Moon clock-rate machinery in operational form.
- [Theory: Lunar PNT Systems](lunar_pnt_systems.md) — two-way time
  transfer and relativistic positioning systems built on these
  corrections.
- [Allan family](allan_family.md) — verification of post-correction
  stability with $\sigma_y(\tau)$ slope identification.

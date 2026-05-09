# Theory: Relativistic Frames and Time Scales

A cislunar PNT stack works with three nested post-Newtonian reference
systems and six coordinate / realised time scales. The Earth side
(BCRS / GCRS, with TCB / TT / TCG / TDB) is the IAU 2000–2006
framework already used by the IERS and modern GNSS; the Moon side
(LCRS, with TCL and TL) was added to the IAU programme by
Resolution II (2024). This page collects definitions, the defining
constants, and the explicit transformation forms — drawing on a
multi-source ingest from the lunar PNT literature
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time),
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales),
[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic),
[Reinhardt, Hartwig & Heinzel 2024](@cite reinhardt-2024-lisa-clock-sync),
[Cacciapuoti & Salomon 2009](@cite cacciapuoti-2009-aces-space-clocks),
and [Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

The companion page
[Theory: Relativistic Clock Corrections](relativistic_corrections.md)
covers the 1PN proper-time mapping, gravitational redshift on the
Moon, the orbit-averaged $L_{Gm}$ rate constant, the cislunar drift
regimes, and the Shapiro / Sagnac light-time corrections.
[Theory: Lunar PNT Systems](lunar_pnt_systems.md) covers two-way time
transfer and relativistic positioning systems.

## Reference systems

### BCRS — Barycentric Celestial Reference System

The BCRS is the IAU global post-Newtonian coordinate system whose
spatial origin is the Solar-System barycentre and whose coordinate
time is TCB
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). Its metric
is constructed by post-Newtonian techniques, with a 1PN proper-time
rate equation governing how clocks anywhere in the Solar System map
to TCB
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time)
[Reinhardt, Hartwig & Heinzel 2024](@cite reinhardt-2024-lisa-clock-sync).
In cislunar PNT operations the BCRS is the natural frame for ephemeris
bookkeeping, light-time solutions, and any computation that crosses
the Earth–Moon system; the Earth-centred (GCRS) and Moon-centred
(LCRS) systems are constructed as nested local frames inside it
[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).

### GCRS — Geocentric Celestial Reference System

The GCRS is the IAU local post-Newtonian frame centred on Earth's
centre of mass, with coordinate time TCG
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales)
[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).
It is the natural frame for clocks that orbit Earth or sit on Earth's
surface; the realised ground-clock scale TT is obtained from TCG by a
single fixed-rate rescaling
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time). In a
cislunar PNT stack, the GCRS is the bridging frame between BCRS
ephemerides and ground-segment clocks; the Earth-side endpoint
correction $(\vec v_E \cdot \vec r_E)/c^{2}$ enters every TT–TDB
conversion
[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).

### LCRS — Lunicentric Celestial Reference System

The LCRS is the IAU-Resolution-II (2024) local post-Newtonian frame
centred on the Moon's centre of mass, constructed analogously to the
GCRS by the IAU B1.5 prescription, with coordinate time TCL
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). Turyshev
truncates its metric at the level needed for $5 \times 10^{-18}$
fractional-frequency / 0.1 ps timing accuracy — retaining lunar
gravity field harmonics through degree $\ell = 9$ and external (Earth,
Sun, planetary) tides through degree $\ell = 8$ in deep cislunar
regimes, with much higher cutoffs near the lunar surface
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). An
equivalent first-order construction in Fermi-normal coordinates
centred on the Moon, with Earth and Sun as external potentials, gives
identical fractional rates — confirming that the proper-time
difference between Moon-fixed and Earth-fixed clocks is independent of
the choice among barycentric, centre-of-mass-rotating, and
Earth-origin-rotating coordinate systems
[Ashby & Patla 2024](@cite ashby-2024-lunar-coordinate-time).

The dominant frame-conversion endpoint correction at the lunar-vicinity
clock's location is

```math
\frac{\vec v_M \cdot \vec r_M}{c^{2}},
```

which reaches $\pm 0.58$ µs annually with $\pm 21$ ps at the lunar
sidereal period of 27.32166 d
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). The
endpoint $-(\vec v_{EM} \cdot \vec X)/c^{2}$ in the TT mapping is a
coordinate artefact, not a physical clock effect, and must be
distinguished in operational notes from the LCRS periodic terms and
from Sagnac
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

## Earth-system time scales

### TT, TCG, TCB, TDB and the defining constants

TT (Terrestrial Time) is a uniformly rescaled TCG using the IAU 2000
defining constant $L_G$:

```math
\mathrm{TCG} - \mathrm{TT} \;=\; \frac{L_G}{1 - L_G}\,(\mathrm{TT} - T_0), \qquad
L_G \;=\; 6.969\,290\,134 \times 10^{-10}.
```

The constant corresponds to a TT-vs-TCG drift of about 60.2147 µs per
day and is fixed by IAU 2000 Resolution B1.9
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

TCB versus TCG is governed by the long-time average of Earth's total
orbital energy:

```math
L_C \;=\; \left\langle\,\frac{1}{c^{2}}\Bigl[\,\frac{v_E^{2}}{2} + \sum_{B \neq E}\frac{G M_B}{r_{BE}}\,\Bigr]\right\rangle
\;=\; 1.480\,826\,854\,55 \times 10^{-8}
\;\approx\; 1.279\,434\,4 \;\mathrm{ms/d}
```

[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

TDB (Barycentric Dynamical Time) is a fixed-rate rescaling of TCB
defined by IAU 2006 Resolution B3:

```math
L_B \;=\; L_G + L_C - L_G\,L_C \;=\; 1.550\,519\,768 \times 10^{-8}
\;\approx\; 1.339\,65 \;\mathrm{ms/d}
```

[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

TT and TDB have no secular rate difference; only periodic terms remain
at the $2 \times 10^{-19}$ / sub-ns level over years
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

### TT vs TDB — ephemeris-plus-endpoint structure

The TT–TDB transformation decomposes into a precomputable time
ephemeris plus an endpoint correction at the GCRS clock's
instantaneous location:

```math
t_{\mathrm{TT}} - t_{\mathrm{TDB}} \;=\; \mathrm{TimeEph}(t) + \frac{1}{c^{2}}\,\vec v_E \cdot \vec r_E + \mathcal{O}(c^{-4})
```

[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).
This is the IERS-2010 / Petit & Luzum operational form; the time
ephemeris is precomputed (table or polynomial expansion), the endpoint
is a per-clock correction.

Internal navigation-system bookkeeping is done in TDB so that the
dynamics / estimation core is decoupled from later time-standard
choices (e.g., how TL is anchored on the Moon)
[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).

## Moon-system time scales

### TCL — Lunicentric Coordinate Time

TCL is the coordinate time of the LCRS, constructed by the same
B1.5-style local-frame prescription that defines TCG inside the GCRS
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). Its
long-time secular rate vs TCB is governed by the Moon's BCRS
orbital-energy average:

```math
L_H \;=\; 1.482\,536\,24 \times 10^{-8} \;\approx\; 1.280\,913\,2 \;\mathrm{ms/d}
```

[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

The Moon-centred TCB–TCL conversion mirrors the IAU TCB–TCG
transformation one-for-one, replacing Earth's barycentric quantities
with the Moon's:

```math
\mathrm{TCB} - \mathrm{TCL} \;=\; \frac{1}{c^{2}}\int_{t_0}^{t}\!\Bigl[\,\frac{v_L^{2}}{2} + \sum_{A \neq L}\frac{G M_A}{r_{LA}}\,\Bigr]\,dt + \frac{\vec v_L \cdot \vec r}{c^{2}} + \mathcal{O}(c^{-4})
```

[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

The fixed-rate TCL–TCG ratio is

```math
\frac{d\mathrm{TCL}}{d\mathrm{TCG}} \;=\; 1 - \frac{\alpha_{LE}}{c^{2}} + \mathcal{O}(c^{-4}),\qquad
\alpha_{LE} \;=\; \frac{v_{LE}^{2}}{2} + \sum_{A \neq L}\frac{G M_A}{r_{LA}} - \sum_{A \neq E}\frac{G M_A}{r_{EA}}.
```

Numerically, TCL leads TCB by 1.2808 ms/day and TCG by 1.4769 µs/day
over a 10-year DE440 ephemeris
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt). TCL–TCG
agreement between the Kopeikin (DE440) and Fienga (INPOP21a)
ephemerides is at the sub-µs level over the 2020–2022 window, with
the dominant periodic line at the anomalistic month
$M \approx 27.55\,\mathrm{d}$, amplitude 0.4778 µs
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt).

### TL — Lunar Time and the $L_L$ defining constant

TL is a surface-realisable lunar time scale defined as a fixed-rate
rescaling of TCL by a constant $L_L$, in direct analogy with the
TT-vs-TCG construction on Earth:

```math
t_{\mathrm{TL}} \;=\; t_{\mathrm{TCL}} - L_L\,(t_{\mathrm{TCL}} - TL_0)
```

[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).

Turyshev recommends promoting $L_L$ to a defining constant for the
LCRS, fixed at the selenoid value
$L_L^{(\mathrm{def})} = 3.139\,05 \times 10^{-11}$, so that operational
realisations at South-Pole or other surface sites are documented as
small offsets rather than as competing definitions
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). The
combined TL-vs-TCB rate constant is

```math
L_M \;=\; L_L + L_H - L_L\,L_H \;=\; 1.485\,675\,294 \times 10^{-8}
\;\approx\; 1.283\,62 \;\mathrm{ms/d}
```

— mirroring $L_B = L_G + L_C - L_G L_C$ on Earth
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

[Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic)
documents three operational choices for $L_L$: $L_L = 0$ (TL = TCL);
$L_L$ chosen to match a selenoid clock's mean rate; or $L_L$ chosen
to minimise the long-term TT–TL drift. Numerical $L_L$ values differ
by 0.84–2.17 ns/d depending on whether the selenoid (RMQ = 1738.0 km)
or the South Pole is the reference site; the IAU-recommended defining
value selects the selenoid
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

A consequence of the $L_L$ definition is that a static observer on
the lunar surface ticks slower than TCG by `−2.7119 µs/day`
[Seyffert 2025](@cite seyffert-2025-relativistic-lunar-pnt); see
[Theory: Relativistic Clock Corrections](relativistic_corrections.md)
for the surface-redshift derivation.

## Time-scale chain

The chained relationship from TL all the way back to TT is

```text
TL ──L_L──► TCL ──L_H──► TCB ──L_C──► TCG ──L_G──► TT
```

with each $L$ a fixed-rate rescaling and the coordinate-time ↔
coordinate-time legs (TCL ↔ TCB, TCB ↔ TCG) carrying their own
endpoint and periodic terms. The chain TL → TCL → TCB → TDB → TT
yields a TL–TT relation with a secular rate constant
$L_{EM} \approx 1.7093906 \times 10^{-11}$ and a series of periodic
corrections at sub-µs amplitudes; the dominant monthly term has
amplitude 0.473 µs
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales). The
geometry term $-(\vec v_{EM} \cdot \vec X)/c^{2}$ contributes one-way
amplitudes scaling with orbit altitude — a few tens of ns at low
orbits up to ~0.81 µs near NRHO apoapsis
[Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

## Practitioner notes

- **Bookkeeping in TDB.** Lunar-PNT internal ephemeris bookkeeping is
  performed in BCRS / TDB so that the dynamics core is decoupled from
  later TL-standard choices; TCL and any TL-mapped product appear at
  the API / mapping layer
  [Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).
- **Parameterise, don't hardcode.** Whether the IAU adopts
  $L_L^{(\mathrm{def})} = 3.139\,05 \times 10^{-11}$ or a differently
  anchored value is still an open community decision; SigmaTauEnsemble's
  [`RelativisticClock`](@ref) should parameterise $L_L$ rather than
  hardcode it
  [Turyshev 2025](@cite turyshev-2025-cislunar-time-scales)
  [Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).
- **Don't bake the rate into the estimator.** Failing to remove
  deterministic relativistic effects before the stochastic clock-state
  estimator runs absorbs the relativistic signature into the clock
  estimate as if it were stochastic, degrading both estimation
  accuracy and dissemination fidelity
  [Leonard, Stewart & Gaylor 2026](@cite leonard-2026-lcrns-relativistic).

## Open questions

- Two upstream primary sources for the LCRS construction
  (Kopeikin & Kaplan 2024, and Turyshev / Williams / Boggs / Park 2025)
  are not yet captured by the SigmaTau bibliography; numerical results
  here cite them transitively via
  [Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).
- An alternative lunar-reference-timescale formulation
  (Bourgoin, Defraigne & Meynadier 2025) is cited by Turyshev and by
  Seyffert; ingest is optional but advisable if the alternative
  becomes operationally relevant
  [Turyshev 2025](@cite turyshev-2025-cislunar-time-scales).

## See also

- [Theory: Relativistic Clock Corrections](relativistic_corrections.md)
  — the corrections (1PN, gravitational redshift, $L_{Gm}$, cislunar
  drift, Shapiro, Sagnac) that ride on top of these frames.
- [Theory: Relativistic Clock Models](relativistic_clocks.md) —
  hub page for the [`RelativisticClock`](@ref) Julia stub and for the
  Earth–Moon clock-rate machinery in operational form.
- [Theory: Lunar PNT Systems](lunar_pnt_systems.md) — two-way time
  transfer and relativistic positioning systems built on the framework
  above.

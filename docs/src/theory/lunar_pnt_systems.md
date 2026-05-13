# Theory: Lunar PNT Systems

This page collects the operational PNT-system patterns that consume
the relativistic clock framework on the preceding pages
([Theory: Relativistic Frames and Time Scales](relativistic_frames_and_timescales.md),
[Theory: Relativistic Clock Corrections](relativistic_corrections.md)):
two-way time transfer in synchronous and asynchronous variants
([Iess et al. 2025](@cite iess-2025-cislunar-od-time-sync),
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer))
and relativistic positioning systems with emission coordinates
([Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning)).

## Two-way time transfer (TWSTFT)

Two-way time transfer compares a remote clock to a reference clock by
exchanging coded pseudorange signals in both directions over the same
path, so that the differential desynchronisation observable cancels
the largest portion of the path delay
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).
In ESA's Moonlight Lunar Communication and Navigation Service (LCNS)
architecture, two TWSTFT variants are formalised: an asynchronous
mode using two independent one-way pseudorange measurements (eight
world-line events $t_1$–$t_8$) and a novel synchronous mode in which
the onboard transponder timestamps the uplink-code epoch with its
local clock and coherently retransmits a downlink code synchronised
in chip rate and code epochs ($t_4 = t_5$ enforced)
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).

The asynchronous TWSTFT desynchronisation is

```math
\Delta \;=\; \tfrac{1}{2}\!\left[\,(t_4 + t_5) - (t_1 + t_8) + \delta_{SC-GS} + \delta_{Media\downarrow\uparrow} + \delta_{up-dn} + \delta_{SC\downarrow\uparrow}\,\right].
```

The synchronous (coherent retransmission) form collapses $t_4, t_5$
to a single epoch:

```math
\Delta \;=\; \tfrac{1}{2}\!\left[\,2\,t_4 - (t_1 + t_8) + \delta_{SC-GS} + \delta_{Media\downarrow\uparrow} + \delta_{up-dn} + \delta_{SC\downarrow\uparrow}\,\right],\qquad t_4 = t_5.
```

[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).
Both methods reach about 1 ns absolute accuracy when internal delays
are calibrated, and the synchronous method is preferred for Moonlight
because it preserves coherent radiometric tracking — the asynchronous
method requires a switch to a non-coherent transponder mode that
interrupts the radiometric observables used for orbit determination
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).

### ATLAS K-band signal structure

The ATLAS K-band SS code is a truncated maximum-length sequence with
about 262 000 chips at about 24 Mcps, giving a 10.9 ms repetition
period for unambiguous epoch identification
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).
OD accuracy of < 10 cm LOS (0.3 ns light-time equivalent) supports
about 1 ns ground-to-space clock synchronisation when internal delays
are calibrated to sub-ns
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).

For an internal-delay calibration error of about 1 ns, the resulting
desynchronisation error scales as $\sim 10^{-5}\,\delta(t_3 - t_4) \approx 10$ fs
— negligible at the 1 ns target
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).
Differential troposphere over the round-trip light time scales as
$\mathrm{RTLT} \times \Delta_T \times dE/dt$, negligible above 20°
elevation; ionospheric contribution at K-band is < 0.2 ns even at
100 TECU
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).

Asynchronous suppression is effective only when $t_4$ and $t_5$ are
made as close as possible (the "L-configuration" of the ACES
experiment, [Cacciapuoti & Salomon 2009](@cite cacciapuoti-2009-aces-space-clocks));
the synchronous scheme is always in L-configuration by construction
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer).

### MSPA ground architecture

The Moonlight architecture comprises three ground stations separated
by about 120° longitude, each with a ~30 cm K-band / ~90 cm X-band
antenna, simultaneously tracking the entire constellation via
coherent two-way microwave links — the Multiple-Spacecraft-per-Aperture
(MSPA) pattern
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer)
[Iess et al. 2025](@cite iess-2025-cislunar-od-time-sync). The ATLAS
architecture range / range-rate performance is about 30 cm and
0.01 mm/s at 10 / 60 s integration with two-way K-band SS signals
[Iess et al. 2025](@cite iess-2025-cislunar-od-time-sync).

The Moonlight payload baseline is a miniRAFS oscillator with
$\sigma_y(\tau = 1000\,\mathrm{s}) = 3 \times 10^{-13}$, sufficient
short-term but requiring frequent ground sync to meet the SISE
specification of 25 m RMS in initial deployment, 10 m RMS fully
deployed; time synchronisation 1–5 ns between satellite clocks (RAFS
or miniRAFS) and the constellation reference
[Iess, Boscagli & Di Benedetto 2026](@cite iess-2026-moonlight-time-transfer)
[Iess et al. 2025](@cite iess-2025-cislunar-od-time-sync).

## Relativistic positioning system

A Relativistic Positioning System (RPS) replaces conventional
Newtonian-time GNSS triangulation with positioning by emission
coordinates $(\tau_1, \tau_2, \tau_3, \tau_4)$ — the proper times of
four emitting satellites at the moments their signals reach the
receiver — making the navigation framework covariant and
observer-independent
[Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning). In
the Autonomous Basis of Coordinates (ABC) construction, two satellites
can determine each other's constants of motion via inter-satellite
proper-time exchange; additional satellites improve accuracy without
requiring ground tracking
[Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning).

An RPS treats relativity as primary in system definition rather than
as a corrective layer added on top of a Newtonian-time framework —
the design pattern explored by the European PECS "Relativistic Global
Navigation System" project (2011–2014)
[Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning).
Two ESA Ariadna projects (2010, 2011) demonstrated RPS feasibility
in the idealised Schwarzschild metric; the PECS project extended to
a perturbed metric with Earth multipoles, tides, rotation, and
Moon / Sun / planet gravity
[Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning).

The linearised Schwarzschild perturbation theory underlying the RPS
work writes

```math
g_{\mu\nu} \;=\; g_{\mu\nu}^{(0)} + h_{\mu\nu},
```

with the background $g_{\mu\nu}^{(0)}$ Schwarzschild and
$h_{\mu\nu} \ll g_{\mu\nu}^{(0)}$ admitting linear perturbation
theory
[Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning).
The vacuum Einstein equation at first order yields a
Regge–Wheeler decomposition of $h_{\mu\nu}$ into odd-parity
($h_0, h_1$) and even-parity ($H_0, H_1, K, H_2$) tensor
spherical-harmonic components, with $H_0 = H_2 = H$ in vacuum
[Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning).

### Perturbation hierarchy at GNSS altitudes

At GNSS altitudes (~20 000 km), the gravitational-perturbation
hierarchy from largest to smallest is
[Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning):

1. Earth multipoles
2. Moon and Sun gravity
3. Solar radiation pressure / Earth albedo (not modelled in RPS)
4. Earth tides
5. Relativistic non-Schwarzschild effects, Jupiter / Venus gravity
6. Earth-rotation gravitomagnetic effects

This is the Earth-orbit analogue of Ashby & Patla's lunar-orbit
"$L_{Gm}$ + 75 ns/orbit residuals" perturbation analysis
([Theory: Relativistic Clock Corrections](relativistic_corrections.md));
both authors are NIST / ESA-adjacent (Ashby NIST; Defraigne SYRTE).

### Status

No current cislunar-PNT architecture in the ingested literature plans
an RPS-style emission-coordinate realisation; Moonlight, LCRNS, and
LDN all use ground-tracked TWSTFT-style architectures with
relativistic corrections layered on
[Gomboc et al. 2013](@cite gomboc-2013-relativistic-positioning).
The RPS concept remains a research-direction reference rather than
an operational blueprint.

## See also

- [Theory: Relativistic Frames and Time Scales](relativistic_frames_and_timescales.md)
  — BCRS / GCRS / LCRS and TCB / TCG / TT / TDB / TCL / TL.
- [Theory: Relativistic Clock Corrections](relativistic_corrections.md)
  — 1PN proper-time mapping, gravitational redshift, $L_{Gm}$,
  cislunar drift regimes, Shapiro and Sagnac.
- [Theory: Relativistic Clock Models](relativistic_clocks.md) — hub
  page for the [`RelativisticClock`](@ref) Julia stub.
- [Theory: Single-Clock Steering](steering.md) — the PID controller
  on which any operational time-transfer / disciplining loop in
  `SigmaTau.Est` would close.

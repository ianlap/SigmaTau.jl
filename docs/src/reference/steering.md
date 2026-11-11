# Steering

A steered clock is kept on its reference by small changes of its frequency,
each computed from an estimate of how far off the clock is. Tight steering
holds the phase close and carries the noise of the reference and of the
readings into the clock; loose steering leaves the clock its own
short-term stability and lets it wander. [`steering_gains`](@ref) sets the
gains from the time constant that the loop should have and
[`lqg_gains`](@ref) from a cost, [`steering_poles`](@ref) shows what any
gains do, and [`steer`](@ref) replays a measured record through the loop,
or predicts the steady state of the loop from a clock model.

```@example steering
using SigmaTau, Random

# the clock against its reference: white FM of 1e-11 at 1 s and random-walk
# FM, read with 0.1 ns of white PM
m = ClockModel(; tau0 = 1.0, R = 1e-20, q1 = 1e-22, q2 = 1e-26)
p, _ = simulate_clock(m, 65_536; rng = Xoshiro(1))
# a loop that settles in 300 s, with a steer each second
g = steering_gains(300.0, 1.0)
```

## The loop

The state is the phase ``p`` and the frequency ``f`` of the clock against
its reference, estimated by the [Kalman filter](kalman.md). After each
reading the frequency of the clock is changed by ``u = -G\hat x``, with
``G = (g_P, g_D)`` the gains, and the change stays in the clock, so over
the interval ``τ`` to the next steer the state moves as

```math
x_{k+1} = Φ\,x_k + B\,u_k + w_k, \qquad
Φ = \begin{pmatrix} 1 & τ \\ 0 & 1 \end{pmatrix}, \qquad
B = \begin{pmatrix} τ \\ 1 \end{pmatrix},
```

and the filter, which knows the steer, predicts
``\hat x^-_{k+1} = (Φ - BG)\,\hat x_k``
([Koppang and Leland (1999), Eqs. (1)–(3)](@cite koppang-1999-lqg-masers);
[Matsakis (2019), Eqs. (A.11)–(A.18)](@cite matsakis-2019-steering)).
``g_P`` is in 1/s: the frequency steer per second of phase error. ``g_D``
is a pure number: the share of the frequency error that one steer removes.

The steer is a change of frequency that persists, the second of the three
actuators that
[Farina et al. (2010), Eqs. (28)–(38)](@cite farina-2010-steering-control)
set side by side. Their Eqs. (42)–(45) turn such steers into those of
another actuator, a phase microstepper for one. For a loop whose control
sets the frequency offset itself, as the tuning voltage of an oscillator
does, the frequency gain is ``g_D - 1``
([Novick et al. (2023), §I](@cite novick-2023-gpsdc-control)).

## Gains from a time constant

The eigenvalues of ``Φ - BG`` are the poles of the loop, the roots of
``z^2 + (τ g_P + g_D - 2)\,z + (1 - g_D)``
([Matsakis (2019), Eqs. (3.8)–(3.9)](@cite matsakis-2019-steering)). A
pole at ``z`` decays as ``e^{-t/T}`` with ``|z| = e^{-τ/T}``. Pole
placement chooses the poles and solves for the gains. With both poles at
``r = e^{-τ/T}`` the loop is critically damped, the fastest response
without ringing, and the gains are those of Eqs. (3.10)–(3.11) of that
paper,

```math
g_P = (1 - r)^2/τ, \qquad g_D = 1 - r^2 .
```

```@example steering
s = steering_poles(g)
s.stable, round.(s.time_constants)
```

Poles anywhere else come from the other method, and gains found any other
way go straight into [`SteeringGains`](@ref):

```@example steering
# a pair of poles, with a steer each minute
steering_gains([0.9 + 0.1im, 0.9 - 0.1im], 60.0)
```

A loop of two gains is stable for ``g_P > 0``, ``g_D > 0`` and
``τ g_P + 2 g_D < 4``
([Matsakis (2019), Eq. (3.17) and the text after Eq. (3.18)](@cite matsakis-2019-steering)).

### Gains from a cost

Pole placement fixes how the loop answers a disturbance and leaves the
variances to follow. The linear quadratic Gaussian (LQG) design
([Koppang and Leland (1999), §III](@cite koppang-1999-lqg-masers);
[Farina et al. (2010), §III-C-2](@cite farina-2010-steering-control))
does the reverse: [`lqg_gains`](@ref) returns the gains that minimize the
sum of

```math
W_p\,p^2 + W_f\,f^2 + W_u\,u^2
```

over the steers, and the response is what it turns out to be.

```@example steering
gc = lqg_gains(1.0; phase = 1.0, steer = 1e10)
```

```@example steering
round.(steering_poles(gc).time_constants)
```

Only the ratios ``W_p τ^2/W_u`` and ``W_f/W_u`` matter, and the gains do
not depend on the noise of the clock or of the readings
([Matsakis (2019), §1](@cite matsakis-2019-steering)): the weights say
what is tolerated, and [`steer`](@ref) on the clock model then gives the
variances that result. A heavier steer weight is a slower loop.
[`lqg_gains`](@ref) reproduces the three designs printed in
[Koppang and Leland (1999), §IV](@cite koppang-1999-lqg-masers), each gain
within 3 parts in ``10^4`` of the printed figure.

### The integral gain

A third gain ``g_I`` acts on the running sum of the phase estimates, the
current one included
([Matsakis (2020), Eq. (3.1)](@cite matsakis-2020-pid)). The poles are
then the roots of

```math
(z - 1)^3 + τ g_I\,z^2 + τ g_P\,z\,(z - 1) + g_D\,(z - 1)^2 = 0 ,
```

Eqs. (3.10) and (4.15) of that paper. At ``g_I = 0`` this is ``z - 1``
times the quadratic above. Three poles at ``r`` give its
Eqs. (3.15)–(3.17),

```math
g_I = (1 - r)^3/τ, \qquad g_P = (1 - 3r^2 + 2r^3)/τ, \qquad g_D = 1 - r^3 ,
```

from `steering_gains(T, τ; integral = true)`. The integral gain is there
for what the model leaves out, a frequency drift above all, and it costs
variance: with no such error the RMS of phase, frequency and steers is
least at ``g_I = 0``
([Matsakis (2020), §§II–III](@cite matsakis-2020-pid)). It drives the
mean of the estimated phase to zero. It does not correct what the filter
itself gets wrong: under a drift that the clock model does not hold, the
filter lags the clock and the steered phase keeps that lag.

## Replaying a record

`steer(data, model, g)` takes the phase record of the free-running clock
against the reference and returns the record as the loop would have left
it.

```@example steering
r = steer(p, m, g)
```

`r.steered` is the steered record, a [`PhaseData`](@ref) for any
deviation, `r.steers` the frequency steer made at each sample, and
`r.correction` the frequency correction in force.

The filter runs once over the record. Because it knows every steer
exactly, its estimate of the steered clock is its estimate of the
free-running clock plus the accumulated correction, so the loop itself is
a short recursion on top of the filter, and trying other gains costs
little. The filter runs at every sample. The steers can come less often,
at `g.interval`, with the gains computed for that interval
([Gödel et al. (2019), §3.1](@cite goedel-2019-steering-comparison)).

```@example steering
using Plots

ms = 2 .^ (0:14)
plot(adev(p, ms); label = "free-running", color = :gray)
plot!(adev(r.steered, ms); label = "steered, T = 300 s")
plot!(ms, steer(m, g; tau = ms, readings = true).adev; label = "predicted",
      linestyle = :dash, color = :black, legend = :bottomleft)
```

The Allan deviation of the steered record follows the clock at averaging
times well below ``T`` and falls as ``1/τ`` well above it, where the clock
is held to the reference.

`sigmatau steer` runs the replay from the [command line](cli.md), with the
model fitted to the first stretch of the record alone and the results
taken from the end of that stretch on.

## Predicted steady state

On a [`ClockModel`](@ref) in place of a record, `steer` returns the steady
state: the RMS of the time error, of the frequency error and of the
steers, and the Allan deviation of the steered clock against its
reference.

```@example steering
steer(m, g)
```

It answers what a time constant costs before any record is replayed:

```@example steering
T = [30, 300, 3000]
q = [steer(m, steering_gains(t, 1.0)) for t in T]
# T (s), the RMS time error (s) and the RMS steer
[T  [x.phase for x in q]  [x.steers for x in q]]
```

Steering is deterministic, so the Kalman gain ``K`` and the covariance
``P`` of the estimate are those of the free-running clock. The estimate
then obeys ``\hat x_{k+1} = (Φ - BG)\,\hat x_k + K ν_{k+1}`` with white
innovations ``ν``, its covariance solves a Lyapunov equation, and the
covariance of the clock itself is that plus ``P``
([Matsakis (2019), §2 and Appendix B](@cite matsakis-2019-steering);
[Matsakis (2020), §II](@cite matsakis-2020-pid)). The Allan deviation
comes from the covariance of the phase at two lags. In simulation, on 40
records of 100 000 samples of the clock above steered at ``T`` = 30 s with
two gains and with three, the RMS of phase, frequency and steers and the
Allan deviation at 1, 16 and 256 s were within 0.4 % of the prediction in
the mean and within 6 % in every record.

The prediction is for the clock itself, without the noise of the readings,
and for a steer at every sample. For a longer interval, replay a record
from [`simulate_clock`](@ref).

A replayed record is of readings, ``z_k = h^\mathsf{T} x_k + v_k``, and
`readings = true` gives `phase` and `adev` for those, to set beside it:

```@example steering
τ = [1, 16, 256, 4096]
predicted = steer(m, g; tau = τ, readings = true).adev
[τ  predicted  adev(r.steered, τ; ci = false).dev]    # τ (s), predicted, replayed
```

The noise ``v`` of a reading, of variance ``R``, is not in the clock at
that instant, so the variance of the reading is that of
``h^\mathsf{T} x`` plus ``R``. It is in the clock afterwards, through the
steers: ``j`` steps on, the joint state of estimate and error has the
covariance ``M^j b_v R`` with it, where ``M`` is the transition of that
state and ``b_v`` what one reading adds to it. With ``C_j`` the covariance
of ``h^\mathsf{T} x`` at lag ``j``, ``c`` the row that picks
``h^\mathsf{T} x`` out of the joint state, and
``ρ_j = c^\mathsf{T} M^j b_v R``,

```math
σ_y^2(jτ) = \frac{3C_0 - 4C_j + C_{2j} + 3R - 4ρ_j + ρ_{2j}}{(jτ)^2} .
```

On 40 simulated records each of five clocks and loops (white PM with white
and random-walk FM under two gains, three gains and gains from a cost,
white PM far above the clock, and flicker PM and FM states), the median
Allan deviation of the steered record at 1 to 256 samples stayed within
1 % of this. Without the terms of the readings the prediction at one
sample was low by a factor of 3 to 17.

## Flicker and drift states

The published loops act on a phase and a frequency state. Two extensions
here are this package's own.

With flicker FM states in the model, the frequency gain acts on the whole
frequency estimate, the frequency state plus the Markov states. In
simulation, on a clock with white FM of ``10^{-11}`` at 1 s, a flicker
floor of ``3 × 10^{-13}`` and 0.1 ns of white PM, steered every second
(five records of 300 000 samples), the phase RMS was 1.10, 1.65 and 2.39
times smaller at ``T`` = 30, 300 and 3000 s than with the gain on the
frequency state alone, for steers 1.08, 1.6 and 2.7 times larger.

A drift state ``d`` is steered out ahead of the feedback: the steer
carries ``-τ\hat d``, and the frequency that the gain acts on is
``\hat f - τ\hat d/2``. Phase and that frequency then obey the two-state
loop above exactly, so a steady drift leaves no standing phase offset,
with or without an integral gain.

```@docs
SteeringGains
steering_gains
lqg_gains
steering_poles
steer
SteeringResult
SteeringPrediction
```

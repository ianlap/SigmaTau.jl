# Steering

A steered clock is kept on its reference by small changes of its frequency,
each computed from an estimate of how far off the clock is. Tight steering
holds the phase close and carries the noise of the reference and of the
readings into the clock; loose steering leaves the clock its own short-term
stability and lets it wander. [`steering_gains`](@ref) sets the gains from
the time constant the loop should have and [`lqg_gains`](@ref) from a cost,
[`steering_poles`](@ref) shows what any gains do, and [`steer`](@ref)
replays a measured record through the loop, or predicts the steady state of
the loop from a clock model.

```julia
m = clock_model(powerlaw_fit(p; alphas=[2, 0, -2]))   # p: clock against reference
g = steering_gains(3600.0, p.tau0)                    # settle in an hour
r = steer(p, m, g)
adev(r.steered)                                       # the clock as steered
steer(m, g)                                           # what the model predicts
```

## The loop

The state is the phase ``p`` and frequency ``f`` of the clock against its
reference, estimated by the Kalman filter of
[`kalman_filter`](@ref). After each reading the frequency of the clock is
changed by ``u = -G\hat x``, with ``G = (g_P, g_D)`` the gains, and the
change stays in the clock, so over the interval ``τ`` to the next steer the
state moves as

```math
x_{k+1} = Φ\,x_k + B\,u_k + w_k, \qquad
Φ = \begin{pmatrix} 1 & τ \\ 0 & 1 \end{pmatrix}, \qquad
B = \begin{pmatrix} τ \\ 1 \end{pmatrix},
```

and the filter, which knows the steer, predicts
``\hat x^-_{k+1} = (Φ - BG)\,\hat x_k``
[koppang-1999-lqg-masers, matsakis-2019-steering](@cite). ``g_P`` is in 1/s: the frequency steer per
second of phase error. ``g_D`` is a pure number: the share of the frequency
error that one steer removes. The steer here is a change of frequency that
persists, the second of the three actuators that
[farina-2010-steering-control](@cite) sets side by side; its Eqs. (42)–(45)
turn such steers into those of another actuator, a phase microstepper for
one. For a loop whose control sets the frequency offset itself, as the
tuning voltage of an oscillator does, the frequency gain is ``g_D - 1``
[novick-2023-gpsdc-control](@cite).

## Gains from a time constant

The eigenvalues of ``Φ - BG`` are the poles of the loop, the roots of
``z^2 + (τ g_P + g_D - 2)\,z + (1 - g_D)``. A pole at ``z`` decays as
``e^{-t/T}`` with ``|z| = e^{-τ/T}``. Pole placement chooses the poles and
solves for the gains. With both poles at ``r = e^{-τ/T}`` the loop is
critically damped, the fastest response without ringing, and

```math
g_P = (1 - r)^2/τ, \qquad g_D = 1 - r^2 .
```

```julia
g = steering_gains(3600.0, 60.0)        # T = 1 h, a steer each minute
steering_poles(g)                       # poles, time constants, ringing, stability
steering_gains([0.9 + 0.1im, 0.9 - 0.1im], 60.0)    # any pair of poles
SteeringGains(1e-5, 0.05, 0, 60.0)      # or gains from elsewhere
```

A loop of two gains is stable for ``g_P > 0``, ``g_D > 0`` and
``τ g_P + 2 g_D < 4``. Gains found any other way go straight into
[`SteeringGains`](@ref).

### Gains from a cost

Pole placement fixes how the loop answers a disturbance and leaves the
variances to follow. The linear quadratic Gaussian (LQG) design
[koppang-1999-lqg-masers, farina-2010-steering-control](@cite) does the
reverse: [`lqg_gains`](@ref) returns the gains that minimize the sum of

```math
W_p\,p^2 + W_f\,f^2 + W_u\,u^2
```

over the steers, and the response is what it turns out to be; read it with
[`steering_poles`](@ref).

```julia
g = lqg_gains(3600.0; phase=1e-3, frequency=1e-3, steer=1e12)   # a gentle loop, steered hourly
steering_poles(g).time_constants
steer(m, g)                              # what it does to the clock
```

Only the ratios ``W_p τ^2/W_u`` and ``W_f/W_u`` matter, and the gains do not
depend on the noise of the clock or of the readings
[matsakis-2019-steering](@cite): the weights say what is tolerated, and
[`steer`](@ref) on the clock model then gives the variances that result. A
heavier steer weight is a slower loop. [`lqg_gains`](@ref) reproduces the
three designs printed in [koppang-1999-lqg-masers](@cite), §IV, each gain
within 3 parts in 10⁴ of the printed figure.

### The integral gain

A third gain ``g_I`` acts on the running sum of the phase estimates, the
current one included [matsakis-2020-pid](@cite). The poles are then the
roots of

```math
(z - 1)^3 + τ g_I\,z^2 + τ g_P\,z\,(z - 1) + g_D\,(z - 1)^2 = 0 ,
```

Eqs. (3.10) and (4.15) of that paper. At ``g_I = 0`` this is ``z - 1`` times the
quadratic above. Three poles at ``r`` give its Eqs. (3.15)–(3.17),

```math
g_I = (1 - r)^3/τ, \qquad g_P = (1 - 3r^2 + 2r^3)/τ, \qquad g_D = 1 - r^3 ,
```

from `steering_gains(T, τ; integral=true)`. The integral gain is there for
what the model leaves out, a frequency drift above all, and it costs
variance: with no such error the RMS of phase, frequency and steers is
least at ``g_I = 0`` [matsakis-2020-pid](@cite). It drives the mean of the
estimated phase to zero. It does not correct what the filter itself gets
wrong: under a drift that the clock model does not hold, the filter lags the
clock and the steered phase keeps that lag.

## Replaying a record

`steer(data, model, g)` takes the phase record of the free-running clock
against the reference and returns the record as the loop would have left
it. The filter runs once over the record. Because the filter knows every
steer exactly, its estimate of the steered clock is its estimate of the
free-running clock plus the accumulated correction, so the loop itself is a
short recursion on top of the filter, and trying other gains costs little.

```julia
r = steer(p, m, g)
r.steered          # PhaseData: pass it to adev, mdev, mtie, ...
r.steers           # the frequency steer made at each sample
r.correction       # the frequency correction in force
```

The Allan deviation of the steered record follows the clock at averaging
times well below ``T`` and falls as ``1/τ`` well above it, where the clock
is held to the reference. The filter runs at every sample; the steers can
come less often, at `g.interval`, with the gains computed for that interval
[goedel-2019-steering-comparison](@cite).

`sigmatau steer` runs the replay from the [command line](../guide/cli.md),
with the model fitted to the first stretch of the record alone and the
results taken from the end of that stretch on.

## Predicted steady state

On a [`ClockModel`](@ref) in place of a record, `steer` returns the steady
state: the RMS of the time error, of the frequency error and of the steers,
and the Allan deviation of the steered clock against its reference.

```julia
q = steer(m, g; tau=[1, 10, 100, 1000])
q.phase, q.frequency, q.steers
q.adev
```

Steering is deterministic, so the Kalman gain ``K`` and the covariance
``P`` of the estimate are those of the free-running clock. The estimate
then obeys ``\hat x_{k+1} = (Φ - BG)\,\hat x_k + K ν_{k+1}`` with white
innovations ``ν``, its covariance solves a Lyapunov equation, and the
covariance of the clock itself is that plus ``P``
[matsakis-2019-steering, matsakis-2020-pid](@cite). The Allan deviation comes from the
covariance of the phase at two lags. In simulation, on 40 records of
100 000 samples of a clock with white and random-walk FM steered every
sample at ``T = 30`` samples, the measured RMS of phase, frequency and
steers and the Allan deviation at 1, 16 and 256 samples all stayed within
6 % of the prediction.

The prediction is for the clock itself, without the noise of the readings,
and for a steer at every sample. For a longer interval, replay a record
from [`simulate_clock`](@ref).

## Flicker and drift states

The published loops act on a phase and a frequency state. Two extensions
here are this package's own.

With flicker FM states in the model, the frequency gain acts on the whole
frequency estimate, the frequency state plus the Markov states. In
simulation, on a clock with white FM `q1 = 1e-22` s, a flicker floor of
3 × 10⁻¹³ and 0.1 ns of white PM, steered every second (five records of
300 000 samples), the phase RMS was 1.10, 1.65 and 2.39 times smaller at
``T`` = 30, 300 and 3000 s than with the gain on the frequency state alone,
for steers 1.08, 1.6 and 2.7 times larger.

A drift state ``d`` is steered out ahead of the feedback: the steer carries
``-τ\hat d``, and the frequency the gain acts on is ``\hat f - τ\hat d/2``.
Phase and that frequency then obey the two-state loop above exactly, so a
steady drift leaves no standing phase offset, with or without an integral
gain.

```@docs
SteeringGains
steering_gains
lqg_gains
steering_poles
steer
SteeringResult
SteeringPrediction
```

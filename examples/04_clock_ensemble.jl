# # From one clock to an ensemble: noise model, Kalman filter, timescale
#
# A clock's stability plot is also the input to its state model. This
# tutorial fits a noise model to a simulated clock, turns it into a
# [`ClockModel`](@ref), tracks the clock with a Kalman filter, asks how long it
# can hold over, and then combines four clocks into a timescale that is
# steadier than any of them. The data are simulated, so every estimate can
# be checked against the truth.

using SigmaTau
using Random
using Statistics

# ## One clock: from ADEV to a state model
#
# White FM at short τ and random-walk FM at long τ, the noise of a typical
# atomic frequency standard, plus a little white PM from the measurement.

τ₀ = 1.0
m_true = ClockModel(; tau0 = τ₀, R = 1e-22, q1 = 1e-24, q2 = 1e-30)
p, truth = simulate_clock(m_true, 100_000; rng = Xoshiro(1));
r = mhdev(p)

# `noise_fit` builds the power-law model one noise type at a time, each
# fitted on the τ range where it dominates, with the slope fixed. The
# indices pick points of `r.tau`.

f = noise_fit(r)
add_component!(f, :WHPM, 1:2)
add_component!(f, :WHFM, 4:9)
add_component!(f, :RWFM, length(r.tau)-4:length(r.tau))
f

# `clock_model` turns the levels into the filter's process noise. Each
# should land within the scatter of a single realization of the truth.

m = clock_model(f)
[m.R m_true.R; m.q1 m_true.q1; m.q2 m_true.q2]

# ## The Kalman filter
#
# `kalman_filter` estimates phase and frequency at every sample;
# `innovation_check` tests whether the model explains the record (its
# normalized innovations should be white with unit variance).

kf = kalman_filter(p, m)
innovation_check(kf).ok

# The frequency estimate tracks the true random-walk frequency to within its
# own 1σ error bar:

err = kf.x[:, 2] .- truth[:, 2]
std(err[1000:end]), sqrt(mean(kf.P[2, 2, 1000:end]))

# `holdover_error` answers the operational question: once the reference is
# lost, how far does the clock's time wander, 1σ, after τ seconds?

holdover_error(kf, [60.0, 3600.0, 86400.0])

# ## Four clocks: a timescale
#
# Now four clocks of different quality, each read against an ideal
# reference (`reference = 0`), so the scale's own stability can be seen
# directly.

ms = [ClockModel(; tau0 = τ₀, q1 = 1e-24 * 4.0^k, q2 = 1e-30 * 4.0^k) for k in 0:3]
cc, _ = simulate_ensemble(ms, 50_000; reference = 0, rng = Xoshiro(2));

# EWFA averages the clocks' predicted phases with equal weights, so the
# noisiest clock pulls the scale down. AT1 learns each clock's
# prediction-error variance and weights it by the inverse.

ewfa = timescale(cc; algorithm = :ewfa)
at1  = timescale(cc; algorithm = :at1)
vec(mean(at1.weights[end÷2:end, :]; dims = 1))

# AT1's weights fall from the best clock to the worst. Against the ideal
# reference, the AT1 scale is steadier than both the EWFA scale and the best
# single clock:

τ = [1, 10, 100, 1000]
[adev(cc[1], τ; ci = false).dev adev(timescale_phase(ewfa), τ; ci = false).dev adev(timescale_phase(at1), τ; ci = false).dev]

# In practice the reference is one of the clocks or a laboratory standard,
# and `timescale_phase` is the scale against it. The same scale comes from
# the command line with `sigmatau timescale FILE --value-cols 2,3,4,5`, and
# from the explorer's **Ensemble** page.

using Plots

p1 = plot(adev(cc[1]); label = "best clock", ylabel = "σ_y(τ)", lw = 1.5)
plot!(p1, adev(cc[4]); label = "worst clock", lw = 1.5)
plot!(p1, adev(timescale_phase(ewfa)); label = "EWFA", ls = :dash)
plot!(p1, adev(timescale_phase(at1)); label = "AT1", lw = 2.5)

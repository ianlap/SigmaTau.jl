# # From one clock to a time scale
#
# A stability curve describes the noise of a clock. The same noise, written
# as a state model, lets a Kalman filter follow the clock, predict it once
# its reference is lost, and weigh it against other clocks. This tutorial
# takes one simulated clock from its record to a noise fit, a
# [`ClockModel`](@ref), a filter and a holdover test, and then forms a time
# scale from four clocks. The data are simulated, so that each estimate can
# be set beside the truth.

using SigmaTau
using Random

# ## A day of one clock
#
# The clock has white FM of ``10^{-12}`` at 1 s, a flicker floor of
# ``10^{-13}`` and random-walk FM that takes over near 1000 s, and it is
# read every second against a reference with 1 ps of white PM.
# [`simulate_clock`](@ref) draws the record and also returns the states
# behind it, the true phase of the clock among them:

τ₀ = 1.0
N  = 86_400
truth = ClockModel(; tau0 = τ₀, R = 1e-24, q1 = 1e-24, q2 = 3e-29,
                   flicker_fm = 1e-26 / 2log(2), tau_max = N)
p, x = simulate_clock(truth, N; rng = Xoshiro(1))
nothing #hide

# [`powerlaw_fit`](@ref) fits the levels ``h_α`` of the five power-law noise
# types to the modified Hadamard deviation of the record:

fit = powerlaw_fit(p)

#-

using Plots

r = mhdev(p)
plot(r; label = "measured")
plot!(r.tau, powerlaw_sigma(fit, r.tau; deviation = :mhdev); label = "fit",
      ls = :dash, color = :black, legend = :bottomleft)

# ## From the fit to a state model
#
# [`clock_model`](@ref) turns the levels into the noise of a state model:
# white PM becomes the variance `R` of a reading, and white and random-walk
# FM become the levels `q1` and `q2` that drive the phase and the frequency
# ([Noise levels](@ref)). Flicker noise has no finite state. With
# `flicker = true` the flicker FM is kept as a few Markov processes that
# add to the frequency ([Flicker states](@ref)); by default it is left out.

model = clock_model(fit; flicker = true)

# The fit also finds a little flicker PM, which this clock does not have.
# It never holds a quarter of the Allan variance, so the model gives it no
# states. The columns are the levels of the fitted model and those that the
# record was made with:

[model.R           truth.R
 model.q1          truth.q1
 model.flicker_fm  truth.flicker_fm
 model.q2          truth.q2]

# The four levels are within 13 % of the truth.

# ## The filter
#
# [`kalman_filter`](@ref) estimates the phase and the frequency of the
# clock at every sample:

kf = kalman_filter(p, model)
nothing #hide

# A reading is the phase of the clock plus the noise of the measurement.
# The simulation kept the phase itself, so the error of the readings and
# that of the estimate can both be measured. The columns are the RMS error
# of the readings, that of the estimate, and the standard deviation that
# the filter gives its estimate, in seconds:

rms(v) = sqrt(sum(abs2, v) / length(v))
k = 1001:N
[rms(p.x[k] .- x[k, 1])  rms(kf.x[k, 1] .- x[k, 1])  sqrt(kf.P[1, 1, end])]

# The estimate is closer to the clock than the readings are, and its error
# is the one that the filter states.
#
# ## Holdover
#
# A measured clock has no truth to compare with, and the record itself has
# to test the model. [`holdover`](@ref) does it where it matters: from
# every sample past the first quarter of the record it carries the estimate
# ahead with no further readings and compares with the reading there
# ([Holdover](../reference/holdover.md)). With `limit` it also answers how
# long the clock keeps time to 1 ns once its reference is lost:

tau = [10, 100, 1000, 10_000]
h = holdover(p, model; tau, limit = 1e-9)

# `rms` is the error made and `claimed` the standard deviation that the
# filter itself gives it. The columns are the holdover time and the ratio
# of `rms` to `claimed`, for this model and for the default one, which
# leaves the flicker FM out:

h0 = holdover(p, clock_model(fit); tau)
[tau  h.rms ./ h.claimed  h0.rms ./ h0.claimed]

# With the states the errors are those claimed, within 4 %. Without them
# the errors are nearly the same, and up to 1.5 times what the filter
# claims.
#
# [`innovation_check`](@ref) tests the one-step predictions of the filter
# instead. It is the stricter test: on a record of this length it fails
# both models, the first for innovations 3 % larger in variance than the
# model expects ([Testing the model](@ref)).
#
# ## Four clocks
#
# Four clocks now, each twice as noisy as the one before it, with white FM
# and random-walk FM. [`simulate_ensemble`](@ref) reads
# them against one reference. With `reference = 0` the reference is
# perfect, so that the stability of each clock and of the scale shows
# directly:

clocks = [ClockModel(; tau0 = τ₀, q1 = 1e-24 * 4.0^i, q2 = 1e-30 * 4.0^i)
          for i in 0:3]
cc, _ = simulate_ensemble(clocks, N; reference = 0, rng = Xoshiro(2))
cc

# A time scale is a weighted mean of its clocks. Equal weights (`:ewfa`)
# let the noisiest clock set the level. AT1 learns the prediction error of
# each clock from the record and weighs by its inverse
# ([Weiss and Weissert (1991), Eqs. (1)–(10)](@cite weiss-weissert-1991-at2)):

ewfa = timescale(cc; algorithm = :ewfa)
at1  = timescale(cc; algorithm = :at1)

# The Kalman scale `:greenhall` takes its weights from the clock models
# ([Greenhall (2003), §4.2](@cite greenhall-2003-jones-tryon)). A measured
# ensemble has to fit them, which [`clock_models`](@ref) does from the
# comparisons. The values are the fitted random-walk FM levels over the
# true ones:

fitted = clock_models(cc)
round.([c.q2 for c in fitted] ./ [c.q2 for c in clocks]; sigdigits = 2)

#-

kalman = timescale(cc; algorithm = :greenhall, models = fitted)
nothing #hide

# [`timescale_phase`](@ref) returns the reference against the scale, which
# against a perfect reference is the error of the scale. The columns are
# ``τ`` and the Allan deviation of the best clock, of the equal-weight
# scale, of AT1 and of the Kalman scale:

m = [1, 10, 100, 1000, 10_000]
σ(scale) = adev(timescale_phase(scale), m; ci = false).dev
[m  adev(cc[1], m; ci = false).dev  σ(ewfa)  σ(at1)  σ(kalman)]

# AT1 is more stable than its best clock at every ``τ`` here. The Kalman
# scale equals it up to 100 s and falls behind the best clock from 1000 s
# on. A day shows too little of the random-walk FM of these clocks to fit
# its level, as the ratios above say, and the scale is as good as its
# models ([Fitted models](@ref)).

m = 2 .^ (0:14)
plot(adev(cc[1], m); label = "best clock", color = :gray)
plot!(adev(timescale_phase(ewfa), m); label = "equal weights")
plot!(adev(timescale_phase(at1), m); label = "AT1")
plot!(adev(timescale_phase(kalman), m); label = "Kalman, fitted models",
      legend = :bottomleft)

# In a laboratory the reference is one of the clocks or a standard, and
# `timescale_phase` is then that reference against the scale.
#
# ## Where next
#
# | For | See |
# |:---|:---|
# | the states and noise levels of the model | [The clock model](../explanation/clock_model.md) |
# | the filter, the smoother and tuning | [Kalman filter](../reference/kalman.md) |
# | holdover with a least-squares fit, and predicted from a noise model | [Holdover](../reference/holdover.md) |
# | the six time scales compared, steps and missing readings | [Clock ensembles](../reference/ensemble.md) |
# | a clock kept on its reference | [Steering](../reference/steering.md) |
#
# `sigmatau holdover --kalman` and `sigmatau timescale` run the same tests
# from the [command line](../reference/cli.md), and the explorer has them
# on its **Holdover** and **Ensemble** pages.

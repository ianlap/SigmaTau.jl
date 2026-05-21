# # Phase data and frequency data
#
# `SigmaTau` operates on two timing-data containers: [`PhaseData`](@ref)
# and [`FrequencyData`](@ref). Both wrap a regularly-sampled vector and
# the sample interval ``\tau_0``. Every public deviation function (`adev`,
# `mdev`, `hdev`, the total family, …) accepts either container — the
# package converts internally as needed.
#
# This tutorial walks through:
#
# 1. Constructing a [`PhaseData`](@ref) from synthetic samples.
# 2. The same operation on a [`FrequencyData`](@ref) record.
# 3. The phase ↔ frequency identity that ties the two together.
#
# The companion tutorial [Computing the Allan deviation](02_compute_adev.md)
# picks up where this one ends, running `adev` on the same fixture.

using SigmaTau
using Random

# ## Synthetic phase residuals
#
# Phase residuals ``x[k]`` (in seconds) sampled every ``\tau_0 = 1`` s.
# We synthesise a realistic mix: white phase modulation (WPM) plus a
# random-walk frequency component (RWFM, shows up as a curving drift).

Random.seed!(20260509)

τ₀ = 1.0
N  = 4096
wpm  = randn(N) .* 1e-9
rwfm = cumsum(cumsum(randn(N) .* 1e-12))
x    = wpm .+ rwfm

pd = PhaseData(x, τ₀)

# A `PhaseData` exposes the sample vector and the cadence:

length(pd.x), pd.tau0

# ## The same data as fractional frequency
#
# Frequency residuals ``y[k] = (x[k+1] - x[k]) / \tau_0`` are produced
# by differencing the phase. Construct a [`FrequencyData`](@ref) directly
# when you already have ``y[k]`` (e.g. from a counter that reports
# fractional-frequency offsets):

y  = randn(N) .* 1e-10
fd = FrequencyData(y, τ₀)

length(fd.y), fd.tau0

# ## Phase ↔ frequency equivalence
#
# `adev(pd)` and `adev(FrequencyData(diff(pd.x) ./ τ₀, τ₀))` produce
# bit-identical numbers — `SigmaTau` runs the same kernel on both
# representations. This is verified in the test suite under
# `test/stab/runtests.jl` (`FrequencyData ↔ PhaseData equivalence`).
# In day-to-day use, pick whichever container matches your input.

# ## Plotting the raw phase residual
#
# A simple time-series plot of ``x[k]`` against ``t = k\tau_0`` shows
# the WPM jitter on top of the RWFM drift. We use `Plots.jl` here;
# the docs build renders this through PGFPlotsX, so the figure below
# is publication-quality vector output.

using Plots
t = (0:N-1) .* τ₀
plot(t, pd.x;
     xlabel = raw"Time \(t\) (s)",
     ylabel = raw"Phase residual \(x(t)\) (s)",
     title  = "Synthetic clock-phase residual (WPM + RWFM)",
     legend = false,
     lw     = 0.8)

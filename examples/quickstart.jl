# quickstart.jl — minimal SigmaTau.jl walkthrough
#
# Run with:
#   julia --project=. examples/quickstart.jl

using SigmaTau
using Random

# ── 1. Synthetic phase data ──────────────────────────────────────────────────
Random.seed!(2026)

tau0 = 1.0                                       # base sample interval, seconds
N    = 4096
# Mix of WPM (white phase) and RWFM (random walk frequency).
x    = randn(N) .* 1e-9 .+ cumsum(cumsum(randn(N) .* 1e-12))

data = PhaseData(x, tau0)

# ── 2. Deviations across an octave-spaced τ-grid ─────────────────────────────
m_grid = [1, 2, 4, 8, 16, 32, 64]

a = adev(data, m_grid)
m = mdev(data, m_grid)
t = tdev(data, m_grid)
h = hdev(data, m_grid)
total = totdev(data, m_grid)

println("τ  ", join(round.(a.tau; digits=2), " "))
println("ADEV    ", a.dev)
println("Noise:  ", a.noise_type)
println("EDF:    ", round.(a.edf; digits=1))

# CI bounds are χ² where EDF is finite, Gaussian fallback otherwise.
@assert all(a.ci_lower .<= a.dev .<= a.ci_upper)

# ── 3. FrequencyData entry point ─────────────────────────────────────────────
# adev(FrequencyData(y, τ₀), …) is identical to adev on the equivalent phase.
y  = randn(N) .* 1e-10
fd = FrequencyData(y, tau0)
a_freq = adev(fd, m_grid)
println("\nADEV from FrequencyData: ", a_freq.dev)

# ── 4. Plot recipes (when Plots is loaded) ──────────────────────────────────-
# Uncomment if Plots.jl is installed:
#
#   using Plots
#   plt = plot(adev(data, m_grid); label="ADEV")
#   plot!(plt, mdev(data, m_grid); label="MDEV")
#   savefig(plt, "stability.png")

# ── 5. Kalman filter (SigmaTauEnsemble) ──────────────────────────────────────
using LinearAlgebra: I

clock = ThreeStateClock(tau=tau0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
est   = StandardKalmanFilter([x[1], 0.0, 0.0], 1e-12 * Matrix(I(3)))

for k in 1:N
    predict!(est, clock, tau0)
    update!(est, clock, x[k])
end

println("\nFinal Kalman state estimate (phase, freq, drift): ", est.x)

# clock_steering.jl — Kalman filter + PID steering on a drifting clock
#
# Demonstrates the optional steering loop on a `StandardKalmanFilter`. The
# clock has a constant fractional-frequency offset; the PID controller pulls
# the residual phase back toward zero by injecting a steering correction
# into the state predict step.
#
# Run with:
#   julia --project=. examples/clock_steering.jl

using SigmaTau
using LinearAlgebra: I
using Random

Random.seed!(7)

N        = 500          # samples
tau      = 1.0          # s
f_offset = 1e-9         # fractional-frequency offset

# Phase data: linear drift from a constant frequency offset, plus noise.
data = [k * tau * f_offset + 1e-12 * randn() for k in 0:N-1]

model = ThreeStateClock(tau=tau, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
est   = StandardKalmanFilter([data[1], 0.0, 0.0], 1e-12 * Matrix(I(3)))
pid   = PIDController(g_p=0.5, g_i=0.05, g_d=0.1)

phase_history = zeros(N)
steer_history = zeros(N)

for k in 1:N
    corr = steer_to_correction(pid.last_steer, 3, tau)
    predict!(est, model, tau; steering=corr)
    update!(est, model, data[k])
    step!(pid, est.x)

    phase_history[k] = est.x[1]
    steer_history[k] = pid.last_steer
end

unsteered_endpoint = N * tau * f_offset

println("Unsteered final phase (no PID): $unsteered_endpoint s")
println("Steered final phase:             $(phase_history[end]) s")
println("Final accumulated steer:         $(steer_history[end])")
println("Reduction factor:                $(round(abs(phase_history[end]) / unsteered_endpoint; digits=4))")

# With Plots.jl loaded:
#   using Plots
#   plot((1:N) .* tau, phase_history; label="phase residual", xlabel="t (s)")

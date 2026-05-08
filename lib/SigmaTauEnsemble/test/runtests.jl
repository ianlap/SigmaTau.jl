using Test
using LinearAlgebra
using Statistics
using StaticArrays
using SigmaTauEnsemble

# Wrap the legacy reference code in its own module so its `step!` and
# `PIDController` definitions don't collide with ours when present locally.
const LEGACY_DIR = joinpath(@__DIR__, "..", "..", "..", "legacy", "julia", "src")
module LegacyKF
    using LinearAlgebra
    using Statistics
    const LEGACY_DIR = joinpath(@__DIR__, "..", "..", "..", "legacy", "julia", "src")
    include(joinpath(LEGACY_DIR, "clock_model.jl"))
    include(joinpath(LEGACY_DIR, "filter.jl"))
end

using .LegacyKF: ClockNoiseParams, ClockModel2, ClockModel3,
                 build_phi, build_Q,
                 kalman_filter, PhaseOnlyMeasurement

@testset "SigmaTauEnsemble.jl" begin

    # ── Φ / Q matrix parity ──────────────────────────────────────────────────
    @testset "Phi / Q matrix parity" begin
        tau = 1.0
        q0 = 1e-22; q1 = 1e-23; q2 = 1e-33; q3 = 1e-43

        legacy_noise = ClockNoiseParams(q_wpm=q0, q_wfm=q1, q_rwfm=q2, q_irwfm=q3)
        legacy3 = ClockModel3(noise=legacy_noise, tau=tau)
        new3    = ThreeStateClock(tau=tau, q0=q0, q1=q1, q2=q2, q3=q3)

        @test Matrix(state_transition(new3)) ≈ build_phi(legacy3) atol=1e-30
        @test Matrix(process_noise(new3))    ≈ build_Q(legacy3)   atol=1e-30

        legacy2 = ClockModel2(noise=legacy_noise, tau=tau)
        new2    = TwoStateClock(tau=tau, q0=q0, q1=q1, q2=q2)

        @test Matrix(state_transition(new2)) ≈ build_phi(legacy2) atol=1e-30
        @test Matrix(process_noise(new2))    ≈ build_Q(legacy2)   atol=1e-30
    end

    # ── Full Kalman filter parity (legacy_compat=true) ───────────────────────
    @testset "StandardKalmanFilter Parity (legacy_compat)" begin
        # Use a fixed seed for deterministic data
        import Random; Random.seed!(42)

        N   = 100
        tau = 1.0
        data = cumsum(randn(N) * 1e-10) # random walk phase

        # Realistic noise parameters (tiny — triggers safe_sqrt clamping)
        q0 = 1e-22 # WPM
        q1 = 1e-23 # WFM
        q2 = 1e-33 # RWFM
        q3 = 1e-43 # IRWFM

        x0_init = [data[1], 0.0, 0.0]
        P0_init = Matrix(1e-12 * I(3))

        # ── Legacy run (steering disabled) ───────────────────────────────────
        legacy_noise = ClockNoiseParams(q_wpm=q0, q_wfm=q1, q_rwfm=q2, q_irwfm=q3)
        legacy_model = ClockModel3(noise=legacy_noise, tau=tau)

        legacy_res = kalman_filter(data, legacy_model, PhaseOnlyMeasurement();
                                   x0=copy(x0_init), P0=copy(P0_init),
                                   g_p=0.0, g_i=0.0, g_d=0.0)

        # ── New filter run (legacy_compat=true) ──────────────────────────────
        new_model = ThreeStateClock(tau=tau, q0=q0, q1=q1, q2=q2, q3=q3)
        est = StandardKalmanFilter(x0_init, P0_init; legacy_compat=true)

        new_phase_est = zeros(N)
        new_freq_est  = zeros(N)
        new_drift_est = zeros(N)
        new_P_history = zeros(3, 3, N)

        for k in 1:N
            predict!(est, new_model, tau)
            update!(est, new_model, data[k])

            new_phase_est[k]       = est.x[1]
            new_freq_est[k]        = est.x[2]
            new_drift_est[k]       = est.x[3]
            new_P_history[:, :, k] .= est.P
        end

        @test legacy_res.phase_est ≈ new_phase_est atol=1e-25 rtol=1e-12
        @test legacy_res.freq_est  ≈ new_freq_est  atol=1e-25 rtol=1e-12
        @test legacy_res.drift_est ≈ new_drift_est atol=1e-25 rtol=1e-12
        @test legacy_res.P_history ≈ new_P_history atol=1e-25 rtol=1e-12
    end

    # ── AD-clean path (legacy_compat=false, no clamping) ─────────────────────
    @testset "StandardKalmanFilter AD-clean (no clamping)" begin
        import Random; Random.seed!(42)

        N   = 100
        tau = 1.0
        data = cumsum(randn(N) * 1e-10)

        # Larger noise so P stays well above 1e-10 threshold naturally
        q0 = 1e-2; q1 = 1e-3; q2 = 1e-4; q3 = 1e-5

        x0_init = [data[1], 0.0, 0.0]
        P0_init = Matrix(1.0 * I(3))

        # Legacy (steering off, large noise → safe_sqrt is a no-op)
        legacy_noise = ClockNoiseParams(q_wpm=q0, q_wfm=q1, q_rwfm=q2, q_irwfm=q3)
        legacy_model = ClockModel3(noise=legacy_noise, tau=tau)
        legacy_res   = kalman_filter(data, legacy_model, PhaseOnlyMeasurement();
                                     x0=copy(x0_init), P0=copy(P0_init),
                                     g_p=0.0, g_i=0.0, g_d=0.0)

        # New (default legacy_compat=false)
        new_model = ThreeStateClock(tau=tau, q0=q0, q1=q1, q2=q2, q3=q3)
        est       = StandardKalmanFilter(x0_init, P0_init)

        new_phase_est = zeros(N)
        new_freq_est  = zeros(N)
        new_drift_est = zeros(N)
        new_P_history = zeros(3, 3, N)

        for k in 1:N
            predict!(est, new_model, tau)
            update!(est, new_model, data[k])

            new_phase_est[k]       = est.x[1]
            new_freq_est[k]        = est.x[2]
            new_drift_est[k]       = est.x[3]
            new_P_history[:, :, k] .= est.P
        end

        @test legacy_res.phase_est ≈ new_phase_est atol=1e-25 rtol=1e-12
        @test legacy_res.freq_est  ≈ new_freq_est  atol=1e-25 rtol=1e-12
        @test legacy_res.drift_est ≈ new_drift_est atol=1e-25 rtol=1e-12
        @test legacy_res.P_history ≈ new_P_history atol=1e-25 rtol=1e-12
    end

    # ── PID steering controller ──────────────────────────────────────────────
    @testset "PIDController.step!" begin
        # Zero state → zero steer.
        pid = PIDController()
        @test step!(pid, [0.0, 0.0, 0.0]) == 0.0

        # Positive phase → negative steer (drives toward zero).
        pid = PIDController(g_p=0.1, g_i=0.0, g_d=0.0)
        s1 = step!(pid, [1.0, 0.0, 0.0])
        @test s1 == -0.1
        @test pid.last_steer == -0.1

        # Integral term accumulates.
        pid = PIDController(g_p=0.0, g_i=0.5, g_d=0.0)
        step!(pid, [1.0, 0.0])    # sumx = 1.0 → -0.5
        step!(pid, [1.0, 0.0])    # sumx = 2.0 → -1.0
        @test pid.last_steer == -1.0

        # Derivative term picks up frequency state.
        pid = PIDController(g_p=0.0, g_i=0.0, g_d=0.2)
        s = step!(pid, [0.0, 1.5, 0.0])
        @test s ≈ -0.3
    end

    @testset "Predict with steering correction" begin
        import Random; Random.seed!(7)

        N   = 200
        tau = 1.0
        # Constant-frequency-offset clock: phase grows linearly with τ·offset.
        f_offset = 1e-9
        data = [k * tau * f_offset + 1e-12 * randn() for k in 0:N-1]

        model = ThreeStateClock(tau=tau, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        est   = StandardKalmanFilter([data[1], 0.0, 0.0], 1e-12 * Matrix(I(3)))
        pid   = PIDController(g_p=0.5, g_i=0.05, g_d=0.1)

        for k in 1:N
            corr = steer_to_correction(pid.last_steer, 3, tau)
            predict!(est, model, tau; steering=corr)
            update!(est, model, data[k])
            step!(pid, est.x)
        end

        # PID should have driven the residual phase down to a small fraction
        # of the unsteered drift after N=200 steps.
        unsteered_endpoint = N * tau * f_offset
        @test abs(est.x[1]) < 1e-2 * unsteered_endpoint
    end

    # ── TwoStateClock smoke test ─────────────────────────────────────────────
    @testset "TwoStateClock basic run" begin
        import Random; Random.seed!(99)

        N   = 50
        tau = 1.0
        data = cumsum(randn(N) * 1e-10)

        model = TwoStateClock(tau=tau, q0=1e-2, q1=1e-3, q2=1e-4)
        est   = StandardKalmanFilter([data[1], 0.0], Matrix(1.0 * I(2)))

        for k in 1:N
            predict!(est, model, tau)
            update!(est, model, data[k])
        end

        @test length(est.x) == 2
        @test size(est.P) == (2, 2)
        @test est.k == N
    end
end

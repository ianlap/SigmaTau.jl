using Test
using LinearAlgebra
using Statistics
using StaticArrays
using SigmaTau

# Legacy reference code lives under `legacy/julia/src/` and is gitignored,
# so it is only present on Ian's working tree. CI checkouts skip the
# legacy-parity testsets entirely; everything else still runs.
const LEGACY_DIR     = joinpath(@__DIR__, "..", "..", "legacy", "julia", "src")
const LEGACY_PRESENT = isfile(joinpath(LEGACY_DIR, "clock_model.jl")) &&
                       isfile(joinpath(LEGACY_DIR, "filter.jl"))

if LEGACY_PRESENT
    # Wrap the legacy reference in its own module so its `step!` and
    # `PIDController` definitions don't collide with ours.
    @eval module LegacyKF
        using LinearAlgebra
        using Statistics
        const LEGACY_DIR = joinpath(@__DIR__, "..", "..", "legacy", "julia", "src")
        include(joinpath(LEGACY_DIR, "clock_model.jl"))
        include(joinpath(LEGACY_DIR, "filter.jl"))
    end
    using .LegacyKF: ClockNoiseParams, ClockModel2, ClockModel3,
                     build_phi, build_Q,
                     kalman_filter, PhaseOnlyMeasurement
else
    @info "legacy/julia/src not present, skipping legacy-KF parity testsets"
end

@testset "SigmaTau.Est" begin

    # ── Φ / Q matrix parity ──────────────────────────────────────────────────
    if LEGACY_PRESENT
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
    end

    # ── Full Kalman filter parity (legacy_compat=true) ───────────────────────
    if LEGACY_PRESENT
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
    end  # if LEGACY_PRESENT

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

    # ── prop!: covariance-only propagation ───────────────────────────────────
    @testset "state_transition / process_noise dt-overload parity" begin
        # The single-arg methods must equal the two-arg form at dt = model.tau,
        # bit-exact — locks that the dt refactor introduced no drift.
        m2 = TwoStateClock(tau=2.5, q0=1e-22, q1=1e-23, q2=1e-33)
        m3 = ThreeStateClock(tau=2.5, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)

        @test state_transition(m2)             == state_transition(m2, m2.tau)
        @test state_transition(m3)             == state_transition(m3, m3.tau)
        @test Matrix(process_noise(m2))        == Matrix(process_noise(m2, m2.tau))
        @test Matrix(process_noise(m3))        == Matrix(process_noise(m3, m3.tau))

        # Φ scales linearly in dt for the polynomial integrator.
        @test state_transition(m3, 1.0)[1, 2]  == 1.0
        @test state_transition(m3, 5.0)[1, 2]  == 5.0
        @test state_transition(m3, 5.0)[1, 3]  == 5.0^2 / 2.0
    end

    @testset "prop! Q-integration parity" begin
        # Hand-derive Q(dt) for ThreeStateClock and check the SMatrix output.
        q1 = 1e-23; q2 = 1e-33; q3 = 1e-43
        dt = 0.7
        Q11 = q1*dt + q2*dt^3/3 + q3*dt^5/20
        Q12 = q2*dt^2/2 + q3*dt^4/8
        Q13 = q3*dt^3/6
        Q22 = q2*dt    + q3*dt^3/3
        Q23 = q3*dt^2/2
        Q33 = q3*dt
        Q_expected = [Q11 Q12 Q13; Q12 Q22 Q23; Q13 Q23 Q33]

        m = ThreeStateClock(tau=1.0, q0=1e-22, q1=q1, q2=q2, q3=q3)
        @test Matrix(process_noise(m, dt)) ≈ Q_expected atol=0.0 rtol=1e-14
    end

    @testset "prop! single-step propagation" begin
        # Φ(dt) x and Φ(dt) P Φ(dt)' + Q(dt) — exact match against manual math.
        m = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        x0 = [3.0, 1e-10, 1e-15]
        P0 = Matrix(1e-18 * I(3))

        est = StandardKalmanFilter(copy(x0), copy(P0))
        prop!(est, m, 1.5)

        Phi = Matrix(state_transition(m, 1.5))
        Q   = Matrix(process_noise(m, 1.5))
        @test Vector(est.x) ≈ Phi * x0       atol=0.0 rtol=1e-14
        @test Matrix(est.P) ≈ Phi * P0 * Phi' + Q  atol=0.0 rtol=1e-14
        # prop! must NOT increment k.
        @test est.k == 0
    end

    @testset "prop! group / additivity composition" begin
        # Two prop!s of dt₁ then dt₂ must equal one prop! of dt₁+dt₂ exactly,
        # because Φ has the group property and Q is additive under it:
        #   Q(dt₁+dt₂) = Φ(dt₂) Q(dt₁) Φ(dt₂)' + Q(dt₂)
        m = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        x0 = [1.0, 1e-10, 1e-16]
        P0 = Matrix(1e-18 * I(3))

        a = StandardKalmanFilter(copy(x0), copy(P0))
        prop!(a, m, 0.4)
        prop!(a, m, 0.6)

        b = StandardKalmanFilter(copy(x0), copy(P0))
        prop!(b, m, 1.0)

        @test Vector(a.x)  ≈ Vector(b.x)  atol=0.0 rtol=1e-14
        @test Matrix(a.P)  ≈ Matrix(b.P)  atol=0.0 rtol=1e-14
    end

    @testset "prop! parity with predict! after k>0" begin
        # Once est.k > 0, predict! and prop! must produce identical state/covariance
        # when called with dt = model.tau (since the gate is the only difference).
        import Random; Random.seed!(123)

        m = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        x0 = [0.0, 0.0, 0.0]
        P0 = Matrix(1e-12 * I(3))

        a = StandardKalmanFilter(copy(x0), copy(P0))
        b = StandardKalmanFilter(copy(x0), copy(P0))

        # Drive both filters past the k>0 gate with one update.
        update!(a, m, 1e-9)
        update!(b, m, 1e-9)
        @test a.k == 1 && b.k == 1

        # Step forward identically with predict! vs prop!.
        for _ in 1:10
            predict!(a, m, m.tau)
            prop!(b, m, m.tau)
        end

        @test Vector(a.x) ≈ Vector(b.x) atol=0.0 rtol=1e-14
        @test Matrix(a.P) ≈ Matrix(b.P) atol=0.0 rtol=1e-14
        # Crucially, prop! did not bump k.
        @test a.k == 1 && b.k == 1
    end

    @testset "prop! steering correction" begin
        # Steering vector adds to the predicted state mean exactly as in predict!.
        m = TwoStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33)
        x0 = [0.0, 0.0]
        P0 = Matrix(1e-18 * I(2))

        u = -3.0e-10              # PID-style frequency correction
        dt = 0.5
        steer = steer_to_correction(u, 2, dt)

        est = StandardKalmanFilter(copy(x0), copy(P0))
        prop!(est, m, dt; steering=steer)

        # Φ(dt)·x₀ = 0; steering adds [u·dt, u].
        @test est.x[1] ≈ u * dt atol=0.0 rtol=1e-14
        @test est.x[2] ≈ u      atol=0.0 rtol=1e-14
    end

    @testset "prop! covariance band over horizons (holdover example pattern)" begin
        # Operationally: build a 1σ covariance band around a deterministic
        # forward projection by prop!ing a side-channel estimator from a fixed
        # starting P0 over each horizon h·τ. Check monotonic growth of σ_x(τ),
        # which is the property the holdover plot relies on.
        m = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=0.0)
        x0 = [0.0, 0.0, 0.0]
        P0 = Matrix(1e-24 * I(3))

        horizons = [1.0, 2.0, 5.0, 10.0, 50.0, 100.0]
        sigmas = Float64[]
        for h in horizons
            est = StandardKalmanFilter(copy(x0), copy(P0))
            prop!(est, m, h)
            push!(sigmas, sqrt(est.P[1, 1]))
        end

        @test all(diff(sigmas) .> 0.0)
    end

    @testset "Stub type ArgumentError" begin
        rc = RelativisticClock()
        @test_throws ArgumentError nstates(rc)
        @test_throws ArgumentError state_transition(rc, 1.0)
        @test_throws ArgumentError process_noise(rc, 1.0)
        @test_throws ArgumentError measurement_matrix(rc)
        @test_throws ArgumentError measurement_noise(rc)

        m = TwoStateClock(tau=1.0, q0=1e-22)
        @test_throws ArgumentError predict!(UDFactorizedFilter(), m, 1.0)
        @test_throws ArgumentError update!(UDFactorizedFilter(), m, 1e-9)
        @test_throws ArgumentError predict!(KuramotoOscillator(), m, 1.0)
        @test_throws ArgumentError update!(KuramotoOscillator(), m, 1e-9)
    end

    # ── ClockEnsemble: stacked time-scale Kalman model ───────────────────────
    @testset "ClockEnsemble construction & operator algebra" begin
        cA = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        cB = ThreeStateClock(tau=1.0, q0=1e-22, q1=2e-23, q2=2e-33, q3=2e-43)
        cC = ThreeStateClock(tau=1.0, q0=1e-22, q1=4e-23, q2=4e-33, q3=4e-43)

        # Binary + → ClockEnsemble
        e2 = cA + cB
        @test e2 isa ClockEnsemble{2, ThreeStateClock}
        @test e2.clocks === (cA, cB)
        @test e2.ref == 1
        @test e2.tau == 1.0

        # Left-associative chaining
        e3 = (cA + cB) + cC
        @test e3 isa ClockEnsemble{3, ThreeStateClock}
        @test e3.clocks === (cA, cB, cC)

        # Right-associative chaining preserves the original reference clock
        # (semantic, not syntactic, associativity).
        e3r = cA + (cB + cC)
        @test e3r isa ClockEnsemble{3, ThreeStateClock}
        @test e3r.clocks === (cA, cB, cC)

        # Varargs explicit constructor
        e3v = ClockEnsemble(cA, cB, cC)
        @test e3v isa ClockEnsemble{3, ThreeStateClock}
        @test e3v.clocks === (cA, cB, cC)

        # Heterogeneous mix throws
        t2 = TwoStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33)
        @test_throws ArgumentError t2 + cA
        @test_throws ArgumentError cA + t2

        # Mismatched tau throws
        cD = ThreeStateClock(tau=2.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        @test_throws ArgumentError ClockEnsemble((cA, cD))

        # Ref-index bounds
        @test_throws ArgumentError ClockEnsemble((cA, cB); ref=0)
        @test_throws ArgumentError ClockEnsemble((cA, cB); ref=3)

        # N ≥ 2 enforcement (single-clock "ensemble" is meaningless)
        @test_throws ArgumentError ClockEnsemble((cA,))
    end

    @testset "ClockEnsemble matrix shapes" begin
        cA = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        cB = ThreeStateClock(tau=1.0, q0=2e-22, q1=2e-23, q2=2e-33, q3=2e-43)
        cC = ThreeStateClock(tau=1.0, q0=3e-22, q1=3e-23, q2=3e-33, q3=3e-43)
        e3 = cA + cB + cC

        @test nstates(e3) == 9

        Φ = state_transition(e3, 1.0)
        Q = process_noise(e3, 1.0)
        @test size(Φ) == (9, 9)
        @test size(Q) == (9, 9)

        # Block-diagonal structure of Φ
        for k in 1:3
            off = (k - 1) * 3
            @test Matrix(Φ[off+1:off+3, off+1:off+3]) ≈ Matrix(state_transition(e3.clocks[k], 1.0))
        end
        # Off-block entries identically zero
        @test all(Φ[i, j] == 0.0 for i in 1:3, j in 4:9)
        @test all(Φ[i, j] == 0.0 for i in 4:6, j in [1:3..., 7:9...])

        # Block-diagonal structure of Q
        for k in 1:3
            off = (k - 1) * 3
            @test Matrix(Q[off+1:off+3, off+1:off+3]) ≈ Matrix(process_noise(e3.clocks[k], 1.0))
        end

        # H selects phase differences vs ref=1
        H = measurement_matrix(e3)
        @test size(H) == (2, 9)
        @test H[1, 1] == -1.0 && H[1, 4] == 1.0   # x_B − x_A
        @test H[2, 1] == -1.0 && H[2, 7] == 1.0   # x_C − x_A
        @test sum(abs, H) == 4.0   # exactly four ±1 entries

        # R diagonal: q0_ref + q0_i per row
        R = measurement_noise(e3)
        @test size(R) == (2, 2)
        @test R[1, 1] ≈ cA.q0 + cB.q0
        @test R[2, 2] ≈ cA.q0 + cC.q0
        @test R[1, 2] == 0.0 && R[2, 1] == 0.0

        # Non-default ref: x_A − x_B, x_C − x_B
        e3_ref2 = ClockEnsemble((cA, cB, cC); ref=2)
        H2 = measurement_matrix(e3_ref2)
        @test H2[1, 4] == -1.0 && H2[1, 1] == 1.0
        @test H2[2, 4] == -1.0 && H2[2, 7] == 1.0
    end

    @testset "ClockEnsemble auto-weights (Stein §VI–VII)" begin
        # Identical clocks → equal weights
        cI = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        eI = cI + cI
        @test eI.weights.a ≈ [0.5, 0.5]
        @test eI.weights.b ≈ [0.5, 0.5]
        @test eI.weights.c ≈ [0.5, 0.5]

        # 10-dB lower q1 on clock B → a_B is 10× a_A (inverse-variance weight)
        cA = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=1e-43)
        cB = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-24, q2=1e-33, q3=1e-43)
        e = cA + cB
        @test e.weights.a[2] / e.weights.a[1] ≈ 10.0 rtol=1e-12
        @test sum(e.weights.a) ≈ 1.0
        @test sum(e.weights.b) ≈ 1.0
        @test sum(e.weights.c) ≈ 1.0

        # All q3 == 0 → c-vector is zeros (TwoStateClock-like sentinel)
        cN1 = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=0.0)
        cN2 = ThreeStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33, q3=0.0)
        eN = cN1 + cN2
        @test all(eN.weights.c .== 0.0)

        # TwoStateClock ensemble: c always zero
        t1 = TwoStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33)
        t2 = TwoStateClock(tau=1.0, q0=1e-22, q1=1e-23, q2=1e-33)
        et = t1 + t2
        @test all(et.weights.c .== 0.0)
        @test et.weights.a ≈ [0.5, 0.5]

        # q1 == 0 → ArgumentError (auto-weights undefined)
        c0 = ThreeStateClock(tau=1.0, q0=1e-22, q1=0.0, q2=1e-33, q3=1e-43)
        @test_throws ArgumentError c0 + cA

        # Explicit weights override
        w = EnsembleWeights{2}(SVector(0.3, 0.7), SVector(0.4, 0.6), SVector(0.0, 0.0))
        eW = ClockEnsemble((cA, cB); weights=w)
        @test eW.weights === w
    end

    @testset "ClockEnsemble runs through StandardKalmanFilter" begin
        import Random; Random.seed!(2025)

        # Two ThreeStateClocks, mild noise — full Kalman over a synthesised
        # phase-difference record. Tests that the generic predict!/update!
        # path consumes the ensemble model unmodified.
        τ  = 1.0
        cA = ThreeStateClock(tau=τ, q0=1e-20, q1=1e-22, q2=1e-30, q3=0.0)
        cB = ThreeStateClock(tau=τ, q0=1e-20, q1=2e-22, q2=2e-30, q3=0.0)
        e  = cA + cB

        # Synthesise independent free-running phase records via cumulative
        # noise — good enough for filter-stability checking.
        N      = 256
        xA     = cumsum(randn(N) * sqrt(cA.q1 * τ))
        xB     = cumsum(randn(N) * sqrt(cB.q1 * τ))
        z      = xB .- xA .+ sqrt(cA.q0 + cB.q0) .* randn(N)

        ns = nstates(e)
        x0 = zeros(ns); x0[1] = xA[1]; x0[4] = xB[1]
        P0 = Matrix(1e-14 * I(ns))
        kf = StandardKalmanFilter(x0, P0)

        for k in 1:N
            predict!(kf, e, τ)
            update!(kf, e, z[k])
        end

        # Numerical health
        @test all(isfinite, kf.x)
        @test all(isfinite, Matrix(kf.P))
        # Covariance stays symmetric
        Pm = Matrix(kf.P)
        @test Pm ≈ Pm' atol=1e-30 rtol=1e-12
        # PSD — all eigenvalues ≥ 0 (within fp noise)
        @test minimum(eigvals(Symmetric(Pm))) > -1e-20

        # Innovation drift: at the end the predicted difference should track
        # the last measurement within a few σ of the observation noise.
        H        = measurement_matrix(e)
        z_hat    = (H * kf.x)[1]
        z_sigma  = sqrt(cA.q0 + cB.q0)
        @test abs(z_hat - z[end]) < 10 * z_sigma   # loose — KF lag at the boundary

        # Counter advanced once per update
        @test kf.k == N
    end

    @testset "ClockEnsemble observable vs non-observable modes" begin
        # With ref=1, the differential mode (x_B − x_A) is observable
        # and its posterior variance stays bounded, while the common
        # mode (x_A absolute phase) is non-observable and its variance
        # grows monotonically — Stein §II / Galleani–Tavella's free-clock
        # non-observability statement, made concrete on the filter output.
        import Random; Random.seed!(7)

        τ   = 1.0
        c   = ThreeStateClock(tau=τ, q0=1e-20, q1=1e-22, q2=1e-30, q3=0.0)
        e   = c + c
        ns  = nstates(e)

        # Drive z = 0 to isolate the structural growth, not noise.
        Nshort = 16
        Nlong  = 512
        kf_short = StandardKalmanFilter(zeros(ns), Matrix(1e-14 * I(ns)))
        for _ in 1:Nshort
            predict!(kf_short, e, τ); update!(kf_short, e, 0.0)
        end
        kf_long = StandardKalmanFilter(zeros(ns), Matrix(1e-14 * I(ns)))
        for _ in 1:Nlong
            predict!(kf_long, e, τ); update!(kf_long, e, 0.0)
        end

        H = measurement_matrix(e)
        Pdiff_short = (H * SMatrix(kf_short.P) * H')[1, 1]
        Pdiff_long  = (H * SMatrix(kf_long.P)  * H')[1, 1]
        Pref_short  = kf_short.P[1, 1]
        Pref_long   = kf_long.P[1, 1]

        # Reference-clock absolute-phase variance grows substantially
        # over the run.
        @test Pref_long > 10 * Pref_short
        # Differential-mode variance grows far less (close to its
        # steady-state bound — Riccati fixed point dominated by R).
        @test Pdiff_long < 2 * Pdiff_short
        # And on absolute scale, the non-observable mode dominates the
        # observable mode after sufficient steps.
        @test Pref_long > Pdiff_long
    end
end

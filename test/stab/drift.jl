# drift.jl — frequency drift fits. scipy reference fits and the published
# slopes of the 1000-point test suite are in test/validation; these check
# the models on exact data and the estimators against their published forms.

using Test
using SigmaTau

@testset "drift_fit" begin
    tau0 = 60.0
    t = collect(0:1999) .* tau0
    T = t[end]

    @testset "$model recovers an exact record" for (model, c, yfun, rate) in (
            (:linear,    [1e-9, 3e-16],            t -> 1e-9 + 3e-16t,              3e-16),
            (:quadratic, [1e-9, 3e-16, 1e-22],     t -> 1e-9 + 3e-16t + 1e-22t^2,   3e-16 + 2e-22T),
            (:log,       [2e-10, 1e-4, 5e-10],     t -> 2e-10log1p(1e-4t) + 5e-10,  2e-10 * 1e-4 / (1 + 1e-4T)),
            (:diffusion, [1e-10, 3e-13, 5000.0],   t -> 1e-10 + 3e-13sqrt(t + 5000), 3e-13 / 2sqrt(T + 5000)))
        f = drift_fit(FrequencyData(yfun.(t), tau0); model=model)
        @test f.coef ≈ c rtol=1e-6
        @test f.rate ≈ rate * 86400 rtol=1e-6
        @test f.rms < 1e-6 * maximum(abs, f.y)
    end

    @testset "gaps are skipped and phase is differenced" begin
        y = 1e-9 .+ 3e-16 .* t .+ 1e-12 .* sin.(t ./ 7000)
        yg = copy(y); yg[300:400] .= NaN
        f = drift_fit(FrequencyData(yg, tau0))
        keep = .!isnan.(yg)
        X = hcat(ones(count(keep)), t[keep])
        b = (X' * X) \ (X' * y[keep])
        @test f.coef ≈ b rtol=1e-10
        @test isnan(f.y[350]) && isfinite(f.fit[350])
        x = [0.0; cumsum(y)] .* tau0
        @test drift_fit(PhaseData(x, tau0)).coef ≈ drift_fit(FrequencyData(y, tau0)).coef rtol=1e-9
    end

    @testset "bad input" begin
        @test_throws ArgumentError drift_fit(FrequencyData(ones(10), 1.0); model=:cubic)
        @test_throws ArgumentError drift_fit(FrequencyData([1.0, 2.0, NaN], 1.0); model=:linear)
        @test_throws ArgumentError drift_fit(FrequencyData(ones(10), 1.0); method=:five_point)
        @test_throws ArgumentError drift_fit(FrequencyData(ones(10), 1.0); model=:log, method=:bisection)
    end

    # A drifting record with a deterministic wobble on it, as phase and as frequency.
    N = 601
    k = 0:N-1
    x = @. 1e-6 + 2e-9 * k * tau0 + 1.5e-16 * (k * tau0)^2 + 1e-10 * sinpi(k / 37) + 3e-11 * cospi(k / 5.3)
    p, fd = PhaseData(x, tau0), FrequencyData(diff(x) ./ tau0, tau0)
    rate(data, method) = drift_fit(data; method).coef[2]

    @testset "$method recovers an exact drift" for method in SigmaTau.DRIFT_METHODS
        f = drift_fit(FrequencyData(1e-9 .+ 3e-16 .* t, tau0); method)
        @test f.method === method
        @test f.coef ≈ [1e-9, 3e-16] rtol=1e-9
        @test f.rate ≈ 3e-16 * 86400 rtol=1e-9
        @test rate(p, method) ≈ rate(fd, method) rtol=1e-6
    end

    @testset "estimators as published" begin
        T = (N - 1) * tau0
        # Three-point estimator: Weiss & Hackman (1992), Eq. (1); SP1065 §5.21.
        @test rate(p, :bisection) ≈ 4 * (x[end] - 2x[(N+1)÷2] + x[1]) / T^2 rtol=1e-9
        # Mean second difference and the standard deviation of its mean: Barnes (1983), §IV.
        d2 = diff(diff(x)) ./ tau0^2
        f = drift_fit(p; method=:second_difference)
        @test f.coef[2] ≈ sum(d2) / length(d2) rtol=1e-9
        @test f.stderr[2] ≈ sqrt(sum(abs2, d2 .- f.coef[2]) / (length(d2) - 1) / length(d2)) rtol=1e-6
        @test f.rate_stderr == f.stderr[2] * 86400
        # Four-point estimator: Greenhall (1997), Eq. (3), with n₁ = N/10 to the nearest integer.
        w = [0.0; cumsum(x)]
        n1, r = 60, 60 / N
        @test rate(p, :four_point) ≈
              6 / (N^3 * tau0^2 * r * (1 - r)) * (w[N+1] - w[1] - (w[N-n1+1] - w[n1+1]) / (1 - 2r)) rtol=1e-6
        # Quadratic fit to phase: drift 2c, its standard error from Barnes (1983), App. A.
        f = drift_fit(p; method=:phase_fit)
        r2 = detrend(p; method=:quadratic).x
        s2 = sum(abs2, r2) / (N - 3)
        @test f.stderr[2] ≈ 2sqrt(s2 * 180 / (N * (N - 1) * (N - 2) * (N + 1) * (N + 2))) / tau0^2 rtol=1e-6
        @test all(isnan, drift_fit(p; method=:bisection).stderr)
        @test isnan(drift_fit(p; method=:four_point).rate_stderr)
    end

    @testset "detrend by a fit" begin
        f = drift_fit(p; method=:phase_fit)
        r = detrend(p, f).x
        @test r .- r[1] ≈ (q = detrend(p; method=:quadratic).x; q .- q[1]) atol=1e-15
        @test detrend(fd, f).y ≈ f.y .- f.fit
        # The second difference that the three-point estimator takes is zero afterwards.
        r = detrend(p, drift_fit(p; method=:bisection)).x
        @test abs(r[end] - 2r[(N+1)÷2] + r[1]) < 1e-18
        @test_throws ArgumentError detrend(PhaseData(x[1:100], tau0), f)
        @test_throws ArgumentError detrend(PhaseData(x, 2tau0), f)
    end

    @testset "gaps and the estimators" begin
        yg = 1e-9 .+ 3e-16 .* t
        yg[[5, 700, 701, 1500]] .= NaN
        @test rate(FrequencyData(yg, tau0), :bisection) ≈ 3e-16 rtol=1e-9
        @test rate(FrequencyData(yg, tau0), :second_difference) ≈ 3e-16 rtol=1e-9
        @test_throws ArgumentError drift_fit(FrequencyData(yg, tau0); method=:phase_fit)
        @test_throws ArgumentError drift_fit(FrequencyData(yg, tau0); method=:four_point)
        xg = copy(x); xg[[50, 300]] .= NaN
        keep = isfinite.(xg)
        X = hcat(ones(count(keep)), k[keep] .* tau0, (k[keep] .* tau0) .^ 2)
        @test rate(PhaseData(xg, tau0), :phase_fit) ≈ 2 * ((X' * X) \ (X' * x[keep]))[3] rtol=1e-6
        @test_throws ArgumentError drift_fit(PhaseData(xg, tau0); method=:four_point)
        @test_throws ArgumentError drift_fit(FrequencyData([NaN, NaN, NaN, 1.0, 2.0, 3.0, 4.0], 1.0); method=:bisection)
        @test_throws ArgumentError drift_fit(FrequencyData([1.0, NaN, 2.0, NaN, 3.0, NaN, 4.0], 1.0); method=:second_difference)
    end
end

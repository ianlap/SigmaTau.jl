# drift.jl — frequency drift fits. scipy reference fits are in
# test/validation (drift_rows); these check the models on exact data.

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
    end
end

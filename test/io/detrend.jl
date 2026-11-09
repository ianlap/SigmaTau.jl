@testset "detrend" begin
    Random.seed!(20260510)

    @testset ":linear removes a planted ramp" begin
        N = 1024
        noise = randn(N) .* 1e-9
        ramp  = 1e-3 .* (0:N-1)
        x     = noise .+ ramp
        pd    = PhaseData(x, 1.0)

        pd_d = detrend(pd; method=:linear)

        # After OLS detrend, both the residual mean and the residual best-fit
        # slope must be zero to machine precision.
        ns      = collect(0:N-1)
        sx      = sum(ns); sxx = sum(abs2, ns)
        sy      = sum(pd_d.x); sxy = sum(ns .* pd_d.x)
        slope   = (N * sxy - sx * sy) / (N * sxx - sx * sx)
        @test abs(sum(pd_d.x) / N) < 1e-12
        @test abs(slope) < 1e-15

        # tau0 preserved, original untouched
        @test pd_d.tau0 == pd.tau0
        @test pd.x == x
    end

    @testset ":endpoint" begin
        N = 256
        x = collect(0.0:N-1) .* 0.5 .+ 7.0
        pd = PhaseData(x, 1.0)
        pd_d = detrend(pd; method=:endpoint)
        @test maximum(abs, pd_d.x) < 1e-12        # exact line removal
    end

    @testset ":mean zeroes the mean (FrequencyData)" begin
        y = randn(512) .+ 42.0
        fd = FrequencyData(y, 1.0)
        fd_d = detrend(fd; method=:mean)
        @test abs(sum(fd_d.y) / length(fd_d.y)) < 1e-12
        @test fd_d.tau0 == fd.tau0
    end

    @testset ":none returns a copy" begin
        x = randn(64)
        pd = PhaseData(copy(x), 1.0)
        pd_d = detrend(pd; method=:none)
        @test pd_d.x == x
        @test pd_d.x !== pd.x                     # not aliased
    end

    @testset "gaps stay NaN and do not enter the fit" begin
        N = 512
        t = collect(0.0:N-1)
        x = 1e-9 .* randn(N) .+ 2e-6 .* t .+ 3e-9 .* t .^ 2
        x[100:140] .= NaN
        ok = isfinite.(x)
        for (method, X) in ((:linear, hcat(ones(N), t)), (:quadratic, hcat(ones(N), t, t .^ 2)))
            d = detrend(PhaseData(x, 1.0); method).x
            b = (X[ok, :]' * X[ok, :]) \ (X[ok, :]' * x[ok])
            @test isequal(isnan.(d), .!ok)
            @test d[ok] ≈ x[ok] .- X[ok, :] * b atol=1e-12
        end
        @test sum(filter(isfinite, detrend(PhaseData(x, 1.0); method=:mean).x)) ≈ 0 atol=1e-9
        e = detrend(PhaseData([NaN; x[2:end]], 1.0); method=:endpoint).x
        @test isnan(e[1]) && abs(e[2]) < 1e-12 && abs(e[end]) < 1e-12
    end

    @testset "unknown method throws" begin
        pd = PhaseData(randn(32), 1.0)
        @test_throws ArgumentError detrend(pd; method=:cubic)
    end
end

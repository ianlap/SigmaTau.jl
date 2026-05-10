@testset "fillgaps (Howe)" begin

    @testset "no-NaN passthrough" begin
        Random.seed!(20260510)
        x = cumsum(randn(256))
        pd = PhaseData(copy(x), 1.0)
        pd_f = fillgaps(pd)
        @test pd_f.x == x
        @test pd_f.tau0 == pd.tau0
    end

    @testset "single-sample gap → mean of neighbours" begin
        x = collect(1.0:64.0)
        x[32] = NaN
        pd = PhaseData(x, 1.0)
        pd_f = fillgaps(pd)
        @test pd_f.x[32] ≈ 0.5 * (x[31] + x[33])
        @test count(isnan, pd_f.x) == 0
    end

    @testset "multi-sample gap is filled and AVAR is preserved" begin
        # Plant a fairly long gap in a power-law-noisy random walk and
        # check that the Allan curve of the filled record stays in family
        # with the gap-free original. The Howe paper claims sub-decade
        # agreement; we use 2x as a robust ceiling.
        Random.seed!(20260510)
        N    = 4096
        tau0 = 1.0
        x_clean = cumsum(randn(N))
        x_gap   = copy(x_clean)
        x_gap[1500:1600] .= NaN

        pd_filled = fillgaps(PhaseData(x_gap, tau0))
        @test count(isnan, pd_filled.x) == 0
        @test length(pd_filled.x) == N

        ms = [1, 2, 4, 8, 16, 32, 64]
        r_clean  = adev(PhaseData(x_clean,         tau0), ms; calc_ci=false)
        r_filled = adev(pd_filled,                       ms; calc_ci=false)

        # Howe's claim: filled AVAR is "close" to clean AVAR. Allow a 2x
        # ratio either way at every τ in the test grid.
        ratios = r_filled.dev ./ r_clean.dev
        @test all(0.5 .<= ratios .<= 2.0)
    end

    @testset "FrequencyData entry point" begin
        Random.seed!(20260510)
        y = randn(512) .* 1e-9
        y[100:120] .= NaN
        fd = FrequencyData(y, 1.0)
        fd_f = fillgaps(fd)
        @test count(isnan, fd_f.y) == 0
        @test fd_f.tau0 == fd.tau0
    end

    @testset "_make_equispaced snap to min(diff)" begin
        # Irregular sampling, dt=1 nominal; one missing slot at index 5.
        t = [0.0, 1.0, 2.0, 3.0, 5.0, 6.0]
        x = [0.0, 1.0, 2.0, 3.0, 5.0, 6.0]
        tfilled, xfilled = SigmaTau._make_equispaced(t, x)
        @test tfilled == collect(0.0:1.0:6.0)
        @test isnan(xfilled[5])
        @test xfilled[[1,2,3,4,6,7]] == [0.0,1.0,2.0,3.0,5.0,6.0]
    end
end

# test/stab/samplestats.jl — sample_statistics, sample_histogram and acf
using Test, Random, Statistics, SigmaTau

@testset "Sample statistics, histogram and ACF" begin
    rng = MersenneTwister(5)
    y = randn(rng, 4096) .* 1e-12
    fd = FrequencyData(y, 1.0)
    pd = PhaseData([0.0; cumsum(y)], 1.0)

    @testset "sample_statistics" begin
        s = sample_statistics(fd)
        @test (s.n, s.gaps, s.min, s.max) == (4096, 0, minimum(y), maximum(y))
        @test (s.mean, s.median, s.std) == (mean(y), median(y), std(y))
        @test s.std_y == s.std
        @test s.adev ≈ adev(fd, [1]; ci = false).dev[1]
        @test s.noise === :WHFM
        # The phase of the same record: its own samples, the same σ_y(τ₀).
        p = sample_statistics(pd)
        @test p.n == 4097 && p.max == maximum(pd.x)
        @test p.adev ≈ s.adev
        @test p.std_y ≈ s.std
        # A gap: left out of the statistics and of every ADEV term reading it.
        x = [0.0, 1.0, 3.0, NaN, 4.0, 4.0, 6.0]
        g = sample_statistics(PhaseData(x, 2.0))
        @test (g.n, g.gaps, g.median) == (6, 1, 3.5)
        @test g.std_y == std([1.0, 2.0, 0.0, 2.0] ./ 2)
        @test g.adev ≈ sqrt(((1 / 2)^2 + (2 / 2)^2) / 4)
        @test_throws ArgumentError sample_statistics(PhaseData([1.0, NaN, 2.0], 1.0))
    end

    @testset "sample_histogram" begin
        h = sample_histogram(fd)
        s = std(y)
        @test length(h.counts) == ceil(Int, (maximum(y) - minimum(y)) / (3.49s * 4096^(-1 / 3)))
        @test (h.edges[1], h.edges[end]) == extrema(y)
        @test sum(h.counts) == h.n == 4096
        @test (h.mean, h.std) == (mean(y), s)
        h3 = sample_histogram(PhaseData([0.0, 1.0, NaN, 3.0, 2.9], 1.0); bins = 3)
        @test (h3.edges, h3.counts) == ([0.0, 1.0, 2.0, 3.0], [1, 1, 2])
        @test length(sample_histogram(fd; bins = 7).counts) == 7
        # One sample far out on a long record: at most 1000 bins.
        @test length(sample_histogram(FrequencyData([zeros(30000); 1.0], 1.0)).counts) == 1000
        @test sample_histogram(FrequencyData(fill(2.0, 5), 1.0)).counts == [5]
        # Samples on the edges, which the division alone puts one bin low.
        @test sample_histogram(FrequencyData((0:10) ./ 10, 1.0); bins = 10).counts == [1, 1, 1, 1, 1, 1, 1, 1, 1, 2]
        @test_throws ArgumentError sample_histogram(fd; bins = 0)
    end

    @testset "acf" begin
        a = acf(fd)
        c(z, k) = sum((z[1:end-k] .- mean(z)) .* (z[1+k:end] .- mean(z)))
        @test a.lag == 0:2048 && a.n == 4096
        @test a.r[1:20] ≈ [c(y, k) / c(y, 0) for k in 0:19]
        @test abs(a.alpha) < 0.2 && a.noise === :WHFM
        # m: frequency averaged without overlap, phase taken every m samples.
        a8 = acf(fd; m = 8, lags = 5)
        z = vec(mean(reshape(y, 8, :); dims = 1))
        @test (a8.tau, a8.n, a8.lag) == (8.0, 512, 0:5)
        @test a8.r ≈ [c(z, k) / c(z, 0) for k in 0:5]
        @test acf(pd; m = 8, lags = 5).r ≈ [c(pd.x[1:8:end], k) / c(pd.x[1:8:end], 0) for k in 0:5]
        @test a8.noise === identify_noise(fd, [8])[1]
        # White PM read as frequency: r₁ near −1/2, and α near 2.
        w = acf(FrequencyData(diff(randn(rng, 4097)), 1.0))
        @test w.r[2] ≈ -0.5 atol = 0.05
        @test w.alpha ≈ 2 atol = 0.2
        @test isnan(acf(fd; m = 256).alpha)                 # 16 points
        # α as identify_noise counts: 239 phase samples hold 29 of m = 8.
        @test isnan(acf(PhaseData(cumsum(y[1:239]), 1.0); m = 8).alpha)
        @test_throws ArgumentError acf(PhaseData([0.0, 1.0, NaN, 2.0, 3.0], 1.0))
        @test_throws ArgumentError acf(fd; m = 0)
        @test_throws ArgumentError acf(fd; lags = 4096)
        @test_throws ArgumentError acf(fd; lags = -5)
    end
end

# Barnes's bias functions B1, B2, B3 and the dead-time correction of ADEV.

using SigmaTau: _dt_g

# Expected N-sample variance of averages whose phase has the structure function
# |t|^(μ+2)/(4 − 2^(μ+2)) (TN 1318 Appendix), summed pair by pair from the
# phase sample times: an average of M measurements over 1 that start every r,
# or with `distributed = false` one measurement over M that starts every M·r.
function structure_variance(N, M, r, μ; distributed=true)
    D(t) = t == 0 ? 0.0 : abs(t)^(μ + 2) / (4 - 2.0^(μ + 2))
    ends(i) = distributed ? [((i * M + j) * r + e, (2e - 1) / M) for j in 0:M-1 for e in (1, 0)] :
                            [(i * M * r + e * M, (2e - 1) / M) for e in (1, 0)]
    s = 0.0
    for i in 0:N-2, j in i+1:N-1
        c = vcat(ends(i), [(t, -w) for (t, w) in ends(j)])
        s -= sum(w1 * w2 * D(t1 - t2) for (t1, w1) in c, (t2, w2) in c) / 2
    end
    return s / (N * (N - 1))
end

@testset "dead-time bias functions" begin
    sig4(x) = round(x; sigdigits=4)

    @testset "published tables" begin
        # NBS TN 375, B1(N, r, μ) for r = 0.001 (p. 11) and B2(r, μ) (pp. 35–37).
        @test sig4(bias_b1(16, 0.001, 0.4)) == 42.80
        @test sig4(bias_b1(1024, 0.001, -1.6)) == 9.946
        @test sig4(bias_b1(Inf, 0.001, -0.6)) == 1.592e4
        @test sig4(bias_b2(0.2, 0.2)) == 0.07996
        @test sig4(bias_b2(0.8, 1.2)) == 0.6910
        # NIST TN 1318: B1 for r = 1.01 and 2 (pp. 21–22), B2 (pp. 29–30) and
        # B3(2, M, r, μ) for r = 0.03, 1.01, 2 and 32 (pp. 31–36).
        @test sig4(bias_b1(64, 1.01, 0.6)) == 10.89
        @test sig4(bias_b1(256, 2, -0.2)) == 1.934
        @test sig4(bias_b1(Inf, 2, -1.4)) == 0.9532
        @test sig4(bias_b1(1024, 2, 0)) == 3.515
        @test sig4(bias_b2(0.0003, -0.8)) == 6.956e-5
        @test sig4(bias_b2(1.1, -1.6)) == 0.8708
        @test sig4(bias_b2(4096, 0)) == 7.082
        @test sig4(bias_b2(64, 0.2)) == 6.834
        @test sig4(bias_b2(4096, 1.8)) == 3.406e6
        @test sig4(bias_b3(2, 64, 0.03, -1.8)) == 488.0
        @test sig4(bias_b3(2, 1024, 1.01, -1.6)) == 9.057
        @test sig4(bias_b3(2, 8, 2, 0)) == 0.6684
        @test sig4(bias_b3(2, 1024, 2, -0.6)) == 0.5757
        @test sig4(bias_b3(2, 16, 32, 0.2)) == 0.3821
        @test sig4(bias_b3(2, 512, 32, 1.4)) == 0.8049
        # SP1065 Fig. 24: √B2 at r = 36 for flicker FM (printed 1.92, here
        # 1.915), random-walk FM and white PM.
        @test all(isapprox.(sqrt.(bias_b2.(36, [:FLFM, :RWFM, :WHPM])), [1.92, 7.31, 0.82]; atol=0.006))
    end

    @testset "identities of TN 1318 Table 1" begin
        @test bias_b1(2, 3.7, 0.3) ≈ 1
        @test bias_b1(9, 2.5, 2) ≈ 9 * 10 / 6
        @test bias_b1(9, 1, 1) ≈ 9 / 2
        @test bias_b1(9, 1, 0) ≈ 9 * log(9) / (16 * log(2))
        @test bias_b1(9, 1, -0.7) ≈ 9 * (1 - 9^-0.7) / (16 * (1 - 2^-0.7))
        @test bias_b1(9, 5, -1) ≈ 1
        @test bias_b1(9, 5, -2) == 1
        @test bias_b2(0, 0.5) == 0
        @test bias_b2(1, -1.3) == 1
        @test bias_b2(7, 2) ≈ 49
        @test bias_b2(7, 1) ≈ (3 * 7 - 1) / 2
        @test bias_b2(0.3, -1) ≈ 0.3
        @test bias_b2(7, -1) ≈ 1
        @test bias_b2(1, -2) == 1 && bias_b2(7, -2) ≈ 2 / 3 && bias_b2(0.3, -2) ≈ 2 / 3
        @test bias_b3(2, 6, 1, 0.4) ≈ 1
        @test bias_b3(2, 6, 3, -2) == 6
        @test bias_b3(2, 1, 3, 0.4) == 1
        @test bias_b3(2, 6, 3, 2) ≈ 1
        @test bias_b3(2, 6, 3, -1) ≈ 1
        # A noise type stands for its μ.
        @test bias_b2(7, :FLFM) == bias_b2(7, 0) && bias_b3(2, 6, 3, :FLPM) == 6
        @test bias_b1(9, 1, :WHFM) ≈ 1
    end

    @testset "accuracy" begin
        # μ = 0, where the published forms are 0/0, is the limit from either side.
        @test bias_b2(4, 1e-9) ≈ bias_b2(4, 0) ≈ bias_b2(4, -1e-9)
        @test bias_b3(2, 16, 4, 1e-9) ≈ bias_b3(2, 16, 4, 0)
        # The series taken from A = 32 joins the three-term form there.
        @test _dt_g(prevfloat(32.0), 0.7) ≈ _dt_g(32.0, 0.7) rtol = 1e-12
        # TN 1318 Eq. (12) in 256-bit arithmetic, at a ratio where its terms
        # cancel to one part in 10¹⁴.
        F(A, μ) = 2A^(μ + 2) - (A + 1)^(μ + 2) - abs(A - 1)^(μ + 2)
        @test bias_b2(1e7, 0.5) ≈ Float64((1 + F(big(1e7), big(0.5)) / 2) / (2 * (1 - sqrt(big(2.0))))) rtol = 1e-12
    end

    @testset "against the structure function" begin
        b3(N, M, r, μ) = structure_variance(N, M, r, μ) / structure_variance(N, M, r, μ; distributed=false)
        @test bias_b3(2, 5, 1.7, -0.6) ≈ b3(2, 5, 1.7, -0.6)
        # More than two samples, which the reference defines and does not tabulate.
        @test bias_b3(5, 3, 1.7, -0.6) ≈ b3(5, 3, 1.7, -0.6)
        @test bias_b3(7, 6, 36, 1) ≈ b3(7, 6, 36, 1)
        @test bias_b1(5, 0.3, -1.5) ≈ structure_variance(5, 1, 0.3, -1.5) / structure_variance(2, 1, 0.3, -1.5)
    end

    @testset "white PM averages that share an end point" begin
        # At r = 1/k the average k places on starts where the first ends. The
        # tables print the value for any other r (1 for B1, M for B3).
        @test bias_b1(3, 1 / 2, -2) ≈ 7 / 6
        @test bias_b3(2, 128, 0.01, -2) ≈ 150
    end

    @testset "arguments" begin
        @test_throws ArgumentError bias_b2(2, 2.1)
        @test_throws ArgumentError bias_b2(-1, 0)
        @test_throws ArgumentError bias_b1(4, 0, 0)
        @test_throws ArgumentError bias_b1(1, 2, 0)
        @test_throws ArgumentError bias_b1(2.5, 2, 0)
        @test_throws ArgumentError bias_b3(1, 4, 2, 0)
        @test_throws ArgumentError bias_b3(2, 0, 2, 0)
        @test_throws ArgumentError bias_b2(2, :FWFM)
        @test bias_b1(Inf, 2, 0.5) == Inf
    end
end

@testset "ADEV with dead time" begin
    y = noise_gen(FrequencyData, 2^17, 1.0; sigma1=Dict(:RWFM => 1.0), rng=MersenneTwister(11)).y
    # 16 s averages started every 64 s: three quarters of the record is dead time.
    fd = FrequencyData([sum(@view y[i:i+15]) / 16 for i in 1:64:length(y)-15], 16.0)
    ms = [1, 2, 4]
    raw = adev(fd, ms; noise=:RWFM)
    res = adev(fd, ms; dead_time=48, noise=:RWFM)
    bias = bias_b2(4, 1) .* bias_b3.(2, ms, 4, 1)

    @test res.dev ≈ raw.dev ./ sqrt.(bias)
    @test [c.hi for c in res.ci] ≈ [c.hi for c in raw.ci] ./ sqrt.(bias)
    @test (res.tau, res.edf, res.noise_type) == (raw.tau, raw.edf, raw.noise_type)
    @test adev(fd, ms; dead_time=48, noise=:RWFM, ci=false).dev == res.dev
    @test deviation(:adev, fd, ms; dead_time=48, noise=:RWFM).dev == res.dev
    @test adev(fd, ms; dead_time=0).dev == adev(fd, ms).dev

    # Uncorrected, the deviation reads 2.3 to 2.0 times high. The corrected one
    # scatters by 1 to 2 % about the truth over seeds at this length.
    truth = adev(FrequencyData(y, 1.0), 16 .* ms; ci=false).dev
    @test isapprox(res.dev, truth; rtol=0.1)
    @test all(raw.dev .> 1.8 .* truth)

    @test_throws ArgumentError adev(fd, ms; dead_time=48)                       # noise not stated
    @test_throws ArgumentError adev(fd, ms; dead_time=-16, noise=:RWFM)         # r = 0
    @test_throws ArgumentError adev(PhaseData(cumsum(y)), ms; dead_time=48, noise=:RWFM)
    @test_throws ArgumentError adev(FrequencyData([fd.y[1:9]; NaN; fd.y[11:40]], 16.0), ms; dead_time=48, noise=:RWFM)
    @test_throws ArgumentError deviation(:mdev, fd, ms; dead_time=48, noise=:RWFM)
end

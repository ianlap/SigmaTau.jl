# nch.jl — N-cornered-hat noise separation tests.
#
# The primary contract is that `nch` reproduces the manual Gray–Allan
# three-cornered-hat solution of examples/06_three_cornered_hat.jl on
# the tutorial's own fixture (same seeds, same grids, same NaN policy),
# plus exact-algebra recovery checks for the general N-clock formula
#   σ̂²ᵢ = (Rᵢ − S/(N−1)) / (N−2),  Rᵢ = Σ_{j≠i} s²ᵢⱼ,  S = Σ_{j<k} s²ⱼₖ.

using Test
using Random
using FFTW                                  # AbstractFFTs backend for noise synth
using SigmaTau
using SigmaTau: _gen_powerlaw_phase, _CIBound

# Synthetic StabilityResult with prescribed deviations on a fixed grid.
_nch_mk(devs; taus = Float64.(eachindex(devs)), kind = :adev,
        neff = fill(100, length(devs))) =
    StabilityResult(kind, collect(Float64, taus), collect(Float64, devs),
                    Symbol[], _CIBound[], Float64[], neff)

@testset "N-cornered hat (nch)" begin

    @testset "reproduces tutorial 06's manual three-cornered hat" begin
        # Same fixture as examples/06_three_cornered_hat.jl: three
        # independent white-FM clocks, pairwise ADEV, manual TCH inversion.
        N    = 8192
        tau0 = 1.0
        m_values = unique(round.(Int, exp10.(range(0, log10(N ÷ 128); length = 12))))

        function clock_record(alpha, scale, seed)
            Random.seed!(seed)
            return scale .* _gen_powerlaw_phase(alpha, N; tau0 = tau0)
        end
        x1 = clock_record(0.0, 1.00e-12, 51)
        x2 = clock_record(0.0, 1.05e-12, 152)
        x3 = clock_record(0.0, 1.10e-12, 253)

        a12 = adev(PhaseData(x1 .- x2, tau0), m_values; ci = false)
        a13 = adev(PhaseData(x1 .- x3, tau0), m_values; ci = false)
        a23 = adev(PhaseData(x2 .- x3, tau0), m_values; ci = false)

        # Manual solution, transcribed from the tutorial's tch_solve.
        v12 = a12.dev .^ 2; v13 = a13.dev .^ 2; v23 = a23.dev .^ 2
        positive_sqrt(v) = v > 0 ? sqrt(v) : NaN
        s1 = positive_sqrt.((v12 .+ v13 .- v23) ./ 2)
        s2 = positive_sqrt.((v12 .+ v23 .- v13) ./ 2)
        s3 = positive_sqrt.((v13 .+ v23 .- v12) ./ 2)

        # 3-arg convenience form: (A−B, B−C, C−A); x13 stands in for the
        # C−A record since the deviation is orientation-invariant.
        clocks = nch(a12, a23, a13)
        @test length(clocks) == 3
        for (manual, got) in zip((s1, s2, s3), clocks)
            @test isnan.(manual) == isnan.(got.dev)
            @test isapprox(got.dev, manual; rtol = 1e-12, nans = true)
        end

        # Output conventions on the realistic fixture.
        for r in clocks
            @test r.deviation_type == :adev          # input type is kept
            @test r.tau == a12.tau
            @test isempty(r.ci) && isempty(r.edf) && isempty(r.noise_type)
            @test r.neff == min.(a12.neff, a13.neff, a23.neff)
        end

        # Matrix form (upper triangle) agrees with the 3-arg form exactly.
        P = Matrix{StabilityResult}(undef, 3, 3)
        P[1, 2] = a12; P[1, 3] = a13; P[2, 3] = a23
        clocks_m = nch(P)
        for (a, b) in zip(clocks, clocks_m)
            @test isapprox(a.dev, b.dev; rtol = 0, nans = true)
        end

        # Lower-triangle fill is accepted too (symmetric access).
        Q = Matrix{StabilityResult}(undef, 3, 3)
        Q[2, 1] = a12; Q[3, 1] = a13; Q[3, 2] = a23
        clocks_q = nch(Q)
        @test isapprox(clocks_q[1].dev, clocks[1].dev; rtol = 0, nans = true)

        # Recovered deviations track the synthesised ground truth within
        # finite-N scatter at the well-populated short τ end.
        truth1 = adev(PhaseData(x1, tau0), m_values; ci = false).dev
        @test isapprox(clocks[1].dev[1], truth1[1]; rtol = 0.1)
    end

    @testset "exact algebra: N = 3, 4, 5 recovery" begin
        # Noise-free pairwise variances s²ᵢⱼ = σ²ᵢ + σ²ⱼ must invert to the
        # exact per-clock variances at every N (validates the general
        # formula and its N=3 reduction to var_A = ½(s²AB + s²AC − s²BC)).
        for n_clocks in (3, 4, 5)
            vars = collect(1.0:n_clocks) .* 0.7
            taus = [1.0, 2.0, 4.0]
            P = Matrix{StabilityResult}(undef, n_clocks, n_clocks)
            for i in 1:(n_clocks - 1), j in (i + 1):n_clocks
                P[i, j] = _nch_mk(fill(sqrt(vars[i] + vars[j]), 3); taus)
            end
            out = nch(P)
            @test length(out) == n_clocks
            for i in 1:n_clocks
                @test all(isapprox.(out[i].dev .^ 2, vars[i]; rtol = 1e-12))
            end
        end

        # Explicit N=3 hand check of the classic formula.
        s2ab, s2bc, s2ca = 5.0, 13.0, 10.0
        out = nch(_nch_mk([sqrt(s2ab)]; taus = [1.0]),
                  _nch_mk([sqrt(s2bc)]; taus = [1.0]),
                  _nch_mk([sqrt(s2ca)]; taus = [1.0]))
        @test out[1].dev[1]^2 ≈ (s2ab + s2ca - s2bc) / 2   # var_A = 1
        @test out[2].dev[1]^2 ≈ (s2ab + s2bc - s2ca) / 2   # var_B = 4
        @test out[3].dev[1]^2 ≈ (s2bc + s2ca - s2ab) / 2   # var_C = 9
    end

    @testset "negative variance estimates become NaN" begin
        # v12 = v13 = 1, v23 = 4 → var₁ = (1 + 1 − 4)/2 = −1 → NaN;
        # var₂ = var₃ = 2 stay finite. NaN inputs also propagate to NaN.
        out = nch(_nch_mk([1.0]; taus = [1.0]),     # A−B
                  _nch_mk([2.0]; taus = [1.0]),     # B−C (dev 2 → var 4)
                  _nch_mk([1.0]; taus = [1.0]))     # C−A
        @test isnan(out[1].dev[1])
        @test out[2].dev[1]^2 ≈ 2.0
        @test out[3].dev[1]^2 ≈ 2.0

        out_nan = nch(_nch_mk([NaN]; taus = [1.0]),
                      _nch_mk([1.0]; taus = [1.0]),
                      _nch_mk([1.0]; taus = [1.0]))
        @test all(isnan(r.dev[1]) for r in out_nan)
    end

    @testset "argument validation" begin
        good = _nch_mk([1.0, 2.0])
        # Non-square / too-small matrices.
        @test_throws ArgumentError nch(Matrix{StabilityResult}(undef, 3, 4))
        P2 = Matrix{StabilityResult}(undef, 2, 2)
        P2[1, 2] = good
        @test_throws ArgumentError nch(P2)
        # Missing pair.
        P = Matrix{StabilityResult}(undef, 3, 3)
        P[1, 2] = good; P[1, 3] = good            # (2,3) unassigned
        @test_throws ArgumentError nch(P)
        # Mismatched deviation type.
        P[2, 3] = _nch_mk([1.0, 2.0]; kind = :hdev)
        @test_throws ArgumentError nch(P)
        # Mismatched τ grid.
        P[2, 3] = _nch_mk([1.0, 2.0]; taus = [1.0, 3.0])
        @test_throws ArgumentError nch(P)
        # neff propagates the elementwise minimum across all pairs.
        P[2, 3] = _nch_mk([1.0, 2.0]; neff = [7, 200])
        out = nch(P)
        @test all(r.neff == [7, 100] for r in out)
    end
end

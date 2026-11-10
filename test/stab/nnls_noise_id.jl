# nnls_noise_id.jl — the weighted-NNLS power-law decomposition (`noise_id_nnls`).

@testset "_nnls_exact" begin
    # Feasible unconstrained solution is returned as is.
    A = [1.0 0.0; 0.0 1.0; 1.0 1.0]; b = [1.0, 2.0, 3.0]
    a, rss = SigmaTau._nnls_exact(A, b)
    @test a ≈ A \ b
    @test rss ≈ sum(abs2, A * (A \ b) .- b) atol = 1e-14
    # Infeasible one ([2, −1]) drops the negative column; all-negative gives zero.
    I2 = [1.0 0.0; 0.0 1.0]
    @test SigmaTau._nnls_exact(I2, [2.0, -1.0]) == ([2.0, 0.0], 1.0)
    @test SigmaTau._nnls_exact(I2, [-1.0, -2.0]) == ([0.0, 0.0], 5.0)
    # Collinear columns: an exact nonnegative fit.
    a, rss = SigmaTau._nnls_exact([1.0 2.0; 1.0 2.0; 1.0 2.0], [3.0, 3.0, 3.0])
    @test rss ≈ 0.0 atol = 1e-20
    @test all(>=(0.0), a) && a[1] + 2a[2] ≈ 3.0
    # Three columns: the residual is the minimum over every feasible support.
    rng = Xoshiro(7)
    A = rand(rng, 12, 3) .+ 0.1
    b = A * [1.0, 0.0, 2.0] .+ 0.01 .* randn(rng, 12)
    a, rss = SigmaTau._nnls_exact(A, b)
    best = minimum(begin s = A[:, S] \ b; all(>=(0.0), s) ? sum(abs2, A[:, S] * s .- b) : Inf end
                   for S in ([1], [2], [3], [1, 2], [1, 3], [2, 3], [1, 2, 3]))
    @test all(>=(0.0), a)
    @test rss ≈ sum(abs2, A * a .- b) atol = 1e-12
    @test rss <= best + 1e-12
    @test_throws ArgumentError SigmaTau._nnls_exact(rand(4, 21), rand(4))
    @test_throws DimensionMismatch SigmaTau._nnls_exact(rand(4, 2), rand(3))
end

@testset "noise_id_nnls" begin
    N = 8192
    pure = (α, σ; seed = 20260807) -> noise_gen(PhaseData, N, 1.0; sigma1 = Dict(α => σ), rng = Xoshiro(seed + α))
    cases = ((2, :WHPM, 1e-11), (0, :WHFM, 1e-12), (-2, :RWFM, 1e-14))

    @testset "basis bookkeeping" begin
        # μ(α): modified family μ = −α−1; unmodified collapses both PM noises onto −2.
        @test all(SigmaTau._nnls_mu(:modified, α) == -α - 1 for α in 2:-1:-4)
        @test [SigmaTau._nnls_mu(:unmodified, α) for α in (2, 1, 0, -2, -4)] == [-2, -2, -1, 1, 3]
        # The Allan family stops at α = −2, the Hadamard family reaches −4 and
        # an unmodified basis carries one PM column.
        @test SigmaTau._nnls_default_alphas(:mdev) == [2, 1, 0, -1, -2]
        @test SigmaTau._nnls_default_alphas(:mhdev) == SigmaTau._nnls_default_alphas(:mhtotdev) == [2, 1, 0, -1, -2, -3, -4]
        @test SigmaTau._nnls_default_alphas(:hdev) == [2, 0, -1, -2, -3, -4]
    end

    @testset "pure noise is identified at every τ" begin
        # The analytic basis is asymptotic, so the assertion is on the dominant
        # type per τ, not on the individual amplitudes.
        for deviation in (:mdev, :mhdev, :mhtotdev, :hdev), (α, sym, σ) in cases
            r = noise_id_nnls(pure(α, σ); deviation)
            @test r.deviation === deviation && all(>=(0.0), r.amplitudes)
            @test r.fit ≈ vec(sum(r.contributions; dims = 2))
            @test all(==(sym), r.dominant)
        end
        # α = −4 needs a third-difference estimator; the Allan family cannot carry it.
        p4 = noise_gen(PhaseData, N, 1.0; h = Dict(-4 => 1e-40), rng = Xoshiro(20260808))
        for deviation in (:mhdev, :hdev, :mhtotdev)
            @test all(==(:RRFM), noise_id_nnls(p4; deviation).dominant)
        end
        @test_throws ArgumentError noise_id_nnls(p4; deviation = :mdev, alphas = [0, -4])
    end

    @testset "amplitudes are on the right scale" begin
        # Against the naive per-τ estimate v(τ)·(τ/τ₀)^−μ, unbiased for a pure
        # power law; catches any uniform-scaling regression.
        for deviation in (:mdev, :hdev), (α, _, σ) in cases
            r = noise_id_nnls(pure(α, σ; seed = 20260810); deviation)
            μ = SigmaTau._nnls_mu(SigmaTau._nnls_family(deviation), α)
            naive = sort(r.var ./ r.tau .^ μ)
            @test isapprox(r.amplitudes[findfirst(==(α), r.alphas)], naive[(end + 1) ÷ 2]; rtol = 0.25)
        end
        p4 = noise_gen(PhaseData, N, 1.0; h = Dict(-4 => 1e-40), rng = Xoshiro(20260811))
        r = noise_id_nnls(p4; deviation = :mhdev)
        naive = sort(r.var ./ r.tau .^ 3)
        @test isapprox(r.amplitudes[findfirst(==(-4), r.alphas)], naive[(end + 1) ÷ 2]; rtol = 0.25)
    end

    @testset "two-component record" begin
        # WFM below the crossover, RWFM above it; the dominant α falls
        # monotonically with τ (the octave at the crossover may read FLFM).
        p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14), rng = Xoshiro(20260809))
        r = noise_id_nnls(p; deviation = :mhdev)
        @test r.dominant[1] === :WHFM && r.dominant[end] === :RWFM
        α_of = Dict(zip(r.noise_types, r.alphas))
        @test issorted([α_of[s] for s in r.dominant]; rev = true)
    end

    @testset "grid, weights and FrequencyData entry" begin
        p = pure(0, 1e-12; seed = 5)
        r = noise_id_nnls(p; deviation = :mhdev)
        @test r.tau == Float64.(filter(>=(4), SigmaTau._default_m_values(N, :mhdev)))
        @test length(r.weights) == length(r.tau) && all(>(0.0), r.weights)
        @test noise_id_nnls(p, 8:8:64; deviation = :mhdev).tau == 8.0:8.0:64.0
        for w in (:relative, fill(1.0, length(r.tau)))
            @test all(>=(0.0), noise_id_nnls(p; deviation = :mhdev, weights = w).amplitudes)
        end
        rf = noise_id_nnls(FrequencyData(diff(p.x) ./ p.tau0, p.tau0); deviation = :mhdev)
        @test rf.deviation === :mhdev && all(>=(0.0), rf.amplitudes)
    end

    @testset "argument validation" begin
        p = noise_gen(PhaseData, 2048, 1.0; sigma1 = Dict(0 => 1e-12), rng = Xoshiro(2))
        @test_throws ArgumentError noise_id_nnls(p; deviation = :adev)
        @test_throws ArgumentError noise_id_nnls(p; alphas = Int[])
        @test_throws ArgumentError noise_id_nnls(p; alphas = [3, 0])
        @test_throws ArgumentError noise_id_nnls(p; deviation = :mdev, alphas = [0, -3])   # below the Allan floor
        @test_throws ArgumentError noise_id_nnls(p; deviation = :hdev, alphas = [2, 1, 0]) # two PM columns
        @test_throws ArgumentError noise_id_nnls(p; weights = :nope)
        @test_throws DimensionMismatch noise_id_nnls(p; weights = [1.0, 2.0])
        @test_throws ArgumentError noise_id_nnls(p, [64, 128]; deviation = :mhdev)         # more columns than τ
    end
end

# nnls_noise_id.jl — α = -4 (RRFM) synthesis and the weighted-NNLS power-law
# decomposition (`noise_id_nnls`).

@testset "RRFM (α = -4) synthesis" begin
    @testset "h₋₄ normalization: driver variance q = 8π⁴τ₀³h₋₄" begin
        # y = √q · cumsum(cumsum(w)) with w ~ N(0, 1), so the exact second
        # difference of the generated y recovers the white driver:
        # Δ²y[k] = √q · w[k]. Its sample variance is the sharpest possible check
        # on the PSD normalization S_y(f) → h₋₄ f⁻⁴ (the whole derivation is the
        # map h₋₄ ↦ q; everything else is an exact discrete identity).
        for (tau0, h) in ((1.0, 1e-40), (0.5, 3e-38), (10.0, 1e-45))
            q = 8 * π^4 * tau0^3 * h
            f = noise_gen(FrequencyData, 200_000, tau0; h = Dict(-4 => h),
                          rng = Xoshiro(20260805))
            d2 = diff(diff(f.y))
            μ  = sum(d2) / length(d2)
            @test sum(abs2, d2 .- μ) / (length(d2) - 1) ≈ q rtol = 0.02
        end
    end

    @testset "MHDEV log-log slope is +3/2" begin
        # σ_MH(τ) ∝ τ^(-α-1)/2 = τ^(+3/2) for α = -4. The Allan family has no
        # such law — its integral diverges — which is why the slope check runs
        # on the Hadamard-family estimator.
        p = noise_gen(PhaseData, 16_384, 1.0; h = Dict(-4 => 1e-40),
                      rng = Xoshiro(20260806))
        ms = [8, 16, 32, 64, 128, 256, 512]
        r  = mhdev(p, ms; ci = false)
        lx = log.(r.tau); ly = log.(r.dev); n = length(lx)
        slope = (n * sum(lx .* ly) - sum(lx) * sum(ly)) / (n * sum(lx .^ 2) - sum(lx)^2)
        @test isapprox(slope, 1.5; atol = 0.15)
    end

    @testset "mixes with the spectral-shaper components" begin
        p = noise_gen(PhaseData, 4096, 1.0; h = Dict(0 => 2e-24, -4 => 1e-42),
                      rng = Xoshiro(11))
        @test length(p.x) == 4096
        @test all(isfinite, p.x)
        # Zero amplitude contributes nothing.
        a = noise_gen(FrequencyData, 512, 1.0; h = Dict(0 => 2e-24), rng = Xoshiro(3))
        b = noise_gen(FrequencyData, 512, 1.0; h = Dict(0 => 2e-24, -4 => 0.0),
                      rng = Xoshiro(3))
        @test a.y == b.y
    end

    @testset "sigma1 at α ≤ -3 raises, h does not" begin
        # The Allan-variance σ_y(τ₀) convention has no finite value for α ≤ -3,
        # so `sigma1` must refuse it — with a message that points at `h`.
        err = try
            noise_gen(PhaseData, 512, 1.0; sigma1 = Dict(-4 => 1e-18))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("diverges", err.msg)
        @test occursin("h = Dict(-4", err.msg)
        @test_throws ArgumentError noise_gen(PhaseData, 512, 1.0; sigma1 = Dict(-3 => 1e-18))
        # α = -3 has no generator at all, in either mode.
        @test_throws ArgumentError noise_gen(PhaseData, 512, 1.0; h = Dict(-3 => 1e-30))
        @test_throws ArgumentError noise_gen(PhaseData, 512, 1.0; h = Dict(-5 => 1e-30))
        @test_throws ArgumentError noise_gen(PhaseData, 512, 1.0; h = Dict(-4 => -1e-40))
        # …but the h path at α = -4 works.
        @test noise_gen(PhaseData, 512, 1.0; h = Dict(-4 => 1e-40)) isa PhaseData
    end
end

@testset "_nnls_exact — exact nonnegative least squares" begin
    @testset "matches the unconstrained solution when it is feasible" begin
        A = [1.0 0.0; 0.0 1.0; 1.0 1.0]
        b = [1.0, 2.0, 3.0]
        a, rss = SigmaTau._nnls_exact(A, b)
        @test a ≈ A \ b
        @test rss ≈ sum(abs2, A * (A \ b) .- b) atol = 1e-14
    end

    @testset "clips an infeasible unconstrained solution" begin
        # Unconstrained: a = [2, -1]; the NNLS optimum drops the second column.
        A = [1.0 0.0; 0.0 1.0]
        b = [2.0, -1.0]
        a, rss = SigmaTau._nnls_exact(A, b)
        @test a ≈ [2.0, 0.0]
        @test rss ≈ 1.0
        # All-negative target → the zero solution.
        a0, rss0 = SigmaTau._nnls_exact(A, [-1.0, -2.0])
        @test a0 == [0.0, 0.0]
        @test rss0 ≈ 5.0
    end

    @testset "collinear columns: picks the better single column" begin
        A = [1.0 2.0; 1.0 2.0; 1.0 2.0]
        b = [3.0, 3.0, 3.0]
        a, rss = SigmaTau._nnls_exact(A, b)
        @test rss ≈ 0.0 atol = 1e-20
        @test all(>=(0.0), a)
        @test A * a ≈ b
    end

    @testset "hand-checkable three-column problem" begin
        # Exhaustive enumeration is the reference: the returned residual must be
        # the minimum over every feasible support, checked by brute force here
        # against an independent projected-gradient descent.
        rng = Xoshiro(7)
        A = rand(rng, 12, 3) .+ 0.1
        b = A * [1.0, 0.0, 2.0] .+ 0.01 .* randn(rng, 12)
        a, rss = SigmaTau._nnls_exact(A, b)
        @test all(>=(0.0), a)
        x = fill(0.1, 3)
        L = 2 * maximum(abs, A' * A)
        for _ in 1:200_000
            x = max.(x .- (2 / L) .* (A' * (A * x .- b)), 0.0)
        end
        @test sum(abs2, A * a .- b) <= sum(abs2, A * x .- b) + 1e-10
        @test rss ≈ sum(abs2, A * a .- b) atol = 1e-12
    end

    @testset "argument validation" begin
        @test_throws ArgumentError SigmaTau._nnls_exact(rand(4, 21), rand(4))
        @test_throws DimensionMismatch SigmaTau._nnls_exact(rand(4, 2), rand(3))
    end
end

@testset "noise_id_nnls" begin
    N = 8192

    @testset "basis bookkeeping" begin
        # μ(α) per family. Modified: μ = -α-1 everywhere. Unmodified: the two PM
        # noises collapse onto μ = -2.
        for α in (2, 1, 0, -1, -2, -3, -4)
            @test SigmaTau._nnls_mu(:modified, α) == -α - 1
        end
        @test SigmaTau._nnls_mu(:unmodified, 2) == -2
        @test SigmaTau._nnls_mu(:unmodified, 1) == -2
        @test SigmaTau._nnls_mu(:unmodified, 0) == -1
        @test SigmaTau._nnls_mu(:unmodified, -2) == 1
        @test SigmaTau._nnls_mu(:unmodified, -4) == 3

        # Defaults: the Allan family stops at α = -2 (its second-difference
        # kernel diverges below that); the Hadamard family reaches α = -4; an
        # unmodified basis carries only one PM column.
        @test SigmaTau._nnls_default_alphas(:mdev)     == [2, 1, 0, -1, -2]
        @test SigmaTau._nnls_default_alphas(:mhdev)    == [2, 1, 0, -1, -2, -3, -4]
        # MHTOTDEV's per-window detrend is quadratic — matched to the third
        # difference's annihilation order — so its measured α = -4 response is a
        # genuine τ³ and it carries the full modified-family basis like MHDEV.
        @test SigmaTau._nnls_default_alphas(:mhtotdev) == [2, 1, 0, -1, -2, -3, -4]
        @test SigmaTau._nnls_default_alphas(:hdev)     == [2, 0, -1, -2, -3, -4]
        @test 1 ∉ SigmaTau._nnls_default_alphas(:hdev)
    end

    @testset "pure-noise identification, modified family" begin
        # A pure α record must be classified as that α at every fitted τ. The
        # analytic basis is asymptotic, so the assertion is on the dominant type
        # per τ (robust) rather than on the individual amplitudes (adjacent
        # power-law columns trade amplitude under sampling noise).
        for deviation in (:mdev, :mhdev, :mhtotdev)
            for (α, sym, σ) in ((2, :WHPM, 1e-11), (0, :WHFM, 1e-12), (-2, :RWFM, 1e-14))
                p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(α => σ),
                              rng = Xoshiro(20260807 + α))
                r = noise_id_nnls(p; deviation = deviation)
                @test r.deviation === deviation
                @test all(>=(0.0), r.amplitudes)
                @test r.fit ≈ vec(sum(r.contributions; dims = 2))
                @test count(==(sym), r.dominant) == length(r.dominant)
            end
        end
    end

    @testset "pure-noise identification, unmodified family" begin
        for (α, sym, σ) in ((2, :WHPM, 1e-11), (0, :WHFM, 1e-12), (-2, :RWFM, 1e-14))
            p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(α => σ),
                          rng = Xoshiro(20260807 + α))
            r = noise_id_nnls(p; deviation = :hdev)
            @test count(==(sym), r.dominant) == length(r.dominant)
        end
    end

    @testset "RRFM is identified by the third-difference estimators" begin
        # α = -4 lives outside the Allan family's convergence range, so `:mdev`
        # cannot even carry the column — asking for it is an error. The
        # third-difference estimators pass the τ³ response through and
        # identify it — including MHTOTDEV, whose per-window quadratic detrend
        # preserves the Hadamard range through the boundary extension.
        for deviation in (:mhdev, :hdev, :mhtotdev)
            p = noise_gen(PhaseData, N, 1.0; h = Dict(-4 => 1e-40),
                          rng = Xoshiro(20260808))
            r = noise_id_nnls(p; deviation = deviation)
            @test :RRFM in r.noise_types
            @test count(==(:RRFM), r.dominant) == length(r.dominant)
        end
        p2048 = noise_gen(PhaseData, 2048, 1.0; h = Dict(-4 => 1e-40), rng = Xoshiro(1))
        @test_throws ArgumentError noise_id_nnls(p2048; deviation = :mdev, alphas = [0, -4])

    end

    @testset "recovered amplitude is on the right scale" begin
        # Dominance says "the right column won"; it says nothing about the
        # magnitude. Cross-check the NNLS amplitude against the naive per-τ
        # estimate v(τ)·(τ/τ₀)^-μ, which for a pure power law is an unbiased
        # amplitude estimate at every point and involves no least squares at all.
        # Catches any uniform-scaling regression (weighting, column
        # equilibration, the unscale after the solve).
        for (deviation, cases) in ((:mdev,     ((2, 1e-11), (0, 1e-12), (-2, 1e-14))),
                                   (:mhtotdev, ((2, 1e-11), (0, 1e-12), (-2, 1e-14))),
                                   (:mhdev,    ((2, 1e-11), (0, 1e-12), (-2, 1e-14))),
                                   (:hdev,     ((2, 1e-11), (0, 1e-12), (-2, 1e-14))))
            family = SigmaTau._nnls_family(deviation)
            for (α, σ) in cases
                p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(α => σ),
                              rng = Xoshiro(20260810 + α))
                r = noise_id_nnls(p; deviation = deviation)
                i = findfirst(==(α), r.alphas)
                μ = SigmaTau._nnls_mu(family, α)
                naive = sort([r.var[k] / r.tau[k]^μ for k in eachindex(r.tau)])
                ref = naive[(length(naive) + 1) ÷ 2]        # median, no Statistics
                @test isapprox(r.amplitudes[i], ref; rtol = 0.25)
            end
        end
        # Same check on the α = -4 column, which only the third-difference
        # estimators carry.
        for deviation in (:mhdev, :hdev)
            p = noise_gen(PhaseData, N, 1.0; h = Dict(-4 => 1e-40),
                          rng = Xoshiro(20260811))
            r = noise_id_nnls(p; deviation = deviation)
            i = findfirst(==(-4), r.alphas)
            naive = sort([r.var[k] / r.tau[k]^3 for k in eachindex(r.tau)])
            ref = naive[(length(naive) + 1) ÷ 2]
            @test isapprox(r.amplitudes[i], ref; rtol = 0.25)
        end
    end

    @testset "two-component record: both regimes appear" begin
        # WFM below the crossover, RWFM above it (σ_WFM/σ_RWFM = 100 s).
        p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14),
                      rng = Xoshiro(20260809))
        r = noise_id_nnls(p; deviation = :mhdev)
        @test r.dominant[1] === :WHFM
        @test r.dominant[end] === :RWFM
        # The dominant α must fall monotonically with τ: a mixture of power laws
        # crosses over once, steepest-rising component last. (The octave nearest
        # the crossover can legitimately report the intermediate FLFM column —
        # a τ⁻¹ and a τ¹ variance summed do look locally flicker-like there —
        # so the assertion is monotonicity, not a hard two-valued split.)
        # (`_alpha_from_noise` stops at α = -2, so map through the basis instead.)
        α_of = Dict(zip(r.noise_types, r.alphas))
        αdom = [α_of[s] for s in r.dominant]
        @test issorted(αdom; rev = true)
    end

    @testset "default grid, weights and FrequencyData entry" begin
        p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(0 => 1e-12), rng = Xoshiro(5))
        r = noise_id_nnls(p; deviation = :mhdev)
        # Default grid = the estimator's octave grid from m = 4 up.
        @test r.tau == Float64.(filter(>=(4), SigmaTau._default_m_values(N, :mhdev)))
        @test length(r.weights) == length(r.tau)
        @test all(>(0.0), r.weights)
        # An explicit grid is honoured verbatim.
        r_ex = noise_id_nnls(p, [2, 4, 8, 16, 32, 64, 128]; deviation = :mhdev)
        @test r_ex.tau == [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0]
        # Alternative weightings run and stay nonnegative.
        for w in (:relative, fill(1.0, length(r.tau)))
            rw = noise_id_nnls(p; deviation = :mhdev, weights = w)
            @test all(>=(0.0), rw.amplitudes)
        end
        # FrequencyData dispatches through _freq_to_phase.
        f = FrequencyData(diff(p.x) ./ p.tau0, p.tau0)
        rf = noise_id_nnls(f; deviation = :mhdev)
        @test rf.deviation === :mhdev
        @test all(>=(0.0), rf.amplitudes)
    end

    @testset "argument validation" begin
        p = noise_gen(PhaseData, 2048, 1.0; sigma1 = Dict(0 => 1e-12), rng = Xoshiro(2))
        @test_throws ArgumentError noise_id_nnls(p; deviation = :adev)
        @test_throws ArgumentError noise_id_nnls(p; alphas = Int[])
        @test_throws ArgumentError noise_id_nnls(p; alphas = [3, 0])
        # Allan family below its convergence floor.
        @test_throws ArgumentError noise_id_nnls(p; deviation = :mdev, alphas = [0, -3])
        # Unmodified family with both PM columns — degenerate basis.
        @test_throws ArgumentError noise_id_nnls(p; deviation = :hdev, alphas = [2, 1, 0])
        @test_throws ArgumentError noise_id_nnls(p; weights = :nope)
        @test_throws DimensionMismatch noise_id_nnls(p; weights = [1.0, 2.0])
        # More columns than usable τ points.
        @test_throws ArgumentError noise_id_nnls(p, [64, 128]; deviation = :mhdev)
    end
end

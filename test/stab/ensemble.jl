# ensemble.jl — N-cornered hat and Groslambert covariance: the inversion
# σ̂²ᵢ = (Rᵢ − S/(N−1)) / (N−2), Rᵢ = Σ_{j≠i} s²ᵢⱼ, S = Σ_{j<k} s²ⱼₖ, and the
# KLTG uncertainties of both, checked against their scatter over many
# realizations.

using Test
using Random
using Statistics
using SigmaTau
using SigmaTau: _CIBound

# Synthetic StabilityResult with prescribed deviations on a fixed grid.
_nch_mk(devs; taus = Float64.(eachindex(devs)), kind = :adev,
        neff = fill(100, length(devs))) =
    StabilityResult(kind, collect(Float64, taus), collect(Float64, devs),
                    Symbol[], _CIBound[], Float64[], neff)

# Observed scatter over R realizations against the predicted spread
# Var(σ̂²) = 2σ̂⁴/edf, for three clocks of one power-law noise type.
function ensemble_scatter(α, ms; N=1024, R=200, rng=Xoshiro(20261108))
    σ = [1.0, 1.5, 0.8]
    est = [Float64[] for _ in 1:2, _ in ms]; sd = [Float64[] for _ in 1:2, _ in ms]
    gl1 = zeros(length(ms)); truth = zeros(length(ms))
    for _ in 1:R
        x = [σ[i] .* noise_gen(PhaseData, N, 1.0; h=Dict(α => 1.0), rng=rng).x for i in 1:3]
        d(i, j) = PhaseData(x[i] .- x[j], 1.0)
        h = ncornered_hat(adev(d(1, 2), ms), adev(d(2, 3), ms), adev(d(3, 1), ms))[2]
        g = groslambert(d(1, 2), d(2, 3), d(3, 1), ms)[2]
        for (e, r) in enumerate((h, g)), k in eachindex(ms)
            push!(est[e, k], r.dev[k]^2); push!(sd[e, k], r.dev[k]^2 * sqrt(2 / r.edf[k]))
        end
        gl1 .+= g.dev .^ 2 ./ R
        truth .+= adev(PhaseData(x[2], 1.0), ms; ci=false).dev .^ 2 ./ R
    end
    return [std(est[c]) / mean(sd[c]) for c in eachindex(est)], gl1 ./ truth
end

@testset "ncornered_hat" begin
    @testset "three white-FM clocks against the classic formula" begin
        N = 8192
        ms = unique(round.(Int, exp10.(range(0, log10(N ÷ 128); length = 12))))
        x = [s .* noise_gen(PhaseData, N, 1.0; h=Dict(0 => 1.0), rng=Xoshiro(k)).x
             for (s, k) in ((1.0, 51), (1.05, 152), (1.1, 253))]
        d(i, j) = adev(PhaseData(x[i] .- x[j], 1.0), ms; ci = false)
        a12, a13, a23 = d(1, 2), d(1, 3), d(2, 3)
        v12, v13, v23 = a12.dev .^ 2, a13.dev .^ 2, a23.dev .^ 2
        root(v) = v > 0 ? sqrt(v) : NaN
        manual = (root.((v12 .+ v13 .- v23) ./ 2), root.((v12 .+ v23 .- v13) ./ 2),
                  root.((v13 .+ v23 .- v12) ./ 2))

        # 3-arg form (A−B, B−C, C−A); the deviation is orientation-invariant.
        clocks = ncornered_hat(a12, a23, a13)
        @test length(clocks) == 3
        @test all(isapprox(got.dev, m; rtol = 1e-12, nans = true) for (m, got) in zip(manual, clocks))
        for r in clocks
            @test r.deviation_type == :adev && r.tau == a12.tau
            @test isempty(r.ci) && isempty(r.edf) && isempty(r.noise_type)   # inputs without edf
            @test r.neff == min.(a12.neff, a13.neff, a23.neff)
        end
        # Matrix form, upper or lower triangle, agrees exactly.
        P = Matrix{StabilityResult}(undef, 3, 3); P[1, 2] = a12; P[1, 3] = a13; P[2, 3] = a23
        Q = Matrix{StabilityResult}(undef, 3, 3); Q[2, 1] = a12; Q[3, 1] = a13; Q[3, 2] = a23
        @test all(isequal(a.dev, b.dev) for (a, b) in zip(clocks, ncornered_hat(P)))
        @test isequal(ncornered_hat(Q)[1].dev, clocks[1].dev)
        # Recovered deviation tracks the truth at the well-populated short τ end.
        @test isapprox(clocks[1].dev[1], adev(PhaseData(x[1], 1.0), ms; ci = false).dev[1]; rtol = 0.1)
    end

    @testset "exact algebra: N = 3, 4, 5 recovery" begin
        for n in (3, 4, 5)
            vars = collect(1.0:n) .* 0.7
            P = Matrix{StabilityResult}(undef, n, n)
            for i in 1:n-1, j in i+1:n
                P[i, j] = _nch_mk(fill(sqrt(vars[i] + vars[j]), 3); taus = [1.0, 2.0, 4.0])
            end
            out = ncornered_hat(P)
            @test length(out) == n
            @test all(isapprox(out[i].dev .^ 2, fill(vars[i], 3); rtol = 1e-12) for i in 1:n)
        end
        # Hand check: s²AB = 5, s²BC = 13, s²CA = 10 gives variances 1, 4, 9.
        out = ncornered_hat(_nch_mk([sqrt(5.0)]), _nch_mk([sqrt(13.0)]), _nch_mk([sqrt(10.0)]))
        @test [r.dev[1]^2 for r in out] ≈ [1.0, 4.0, 9.0]
    end

    @testset "negative variances and NaN inputs become NaN" begin
        # v12 = v13 = 1, v23 = 4: var₁ = −1 → NaN, var₂ = var₃ = 2.
        out = ncornered_hat(_nch_mk([1.0]), _nch_mk([2.0]), _nch_mk([1.0]))
        @test isnan(out[1].dev[1]) && out[2].dev[1]^2 ≈ 2.0 && out[3].dev[1]^2 ≈ 2.0
        @test all(isnan(r.dev[1]) for r in ncornered_hat(_nch_mk([NaN]), _nch_mk([1.0]), _nch_mk([1.0])))
    end

    @testset "argument validation" begin
        good = _nch_mk([1.0, 2.0])
        @test_throws ArgumentError ncornered_hat(Matrix{StabilityResult}(undef, 3, 4))
        P2 = Matrix{StabilityResult}(undef, 2, 2); P2[1, 2] = good
        @test_throws ArgumentError ncornered_hat(P2)
        P = Matrix{StabilityResult}(undef, 3, 3); P[1, 2] = good; P[1, 3] = good
        @test_throws ArgumentError ncornered_hat(P)                        # (2,3) missing
        P[2, 3] = _nch_mk([1.0, 2.0]; kind = :hdev)
        @test_throws ArgumentError ncornered_hat(P)                        # mixed type
        P[2, 3] = _nch_mk([1.0, 2.0]; taus = [1.0, 3.0])
        @test_throws ArgumentError ncornered_hat(P)                        # mixed grid
        P[2, 3] = _nch_mk([1.0, 2.0]; neff = [7, 200])
        @test all(r.neff == [7, 100] for r in ncornered_hat(P))           # elementwise minimum
    end
end

@testset "ensemble uncertainties" begin
    # White PM at short τ, random-walk FM out to τ = N/16, where the EDF of
    # each record's noise type carries the long-lag correlation.
    for (α, ms) in ((2, [1, 8]), (-2, [16, 64]))
        ratios, bias = ensemble_scatter(α, ms)
        @test all(0.8 .< ratios .< 1.25)
        @test all(abs.(bias .- 1) .< 0.1)
    end
end

@testset "groslambert" begin
    rng = Xoshiro(7)
    x = [s .* randn(rng, 2048) for s in (1.0, 2.0, 0.5, 1.5)]
    d(i, j) = PhaseData(x[i] .- x[j], 1.0)
    g3 = groslambert(d(1, 2), d(2, 3), d(3, 1), [1, 2, 4])
    P = Matrix{PhaseData}(undef, 3, 3)
    P[1, 2] = d(1, 2); P[3, 2] = d(3, 2); P[1, 3] = d(1, 3)   # mixed orientation
    @test all(groslambert(P, [1, 2, 4])[i].dev ≈ g3[i].dev for i in 1:3)
    P[1, 1] = PhaseData(zeros(16), 1.0)                         # ignored, also for the grid
    @test groslambert(P)[1].tau == groslambert(d(1, 2), d(2, 3), d(3, 1))[1].tau
    g32 = groslambert(d(1, 2), PhaseData(Float32.(x[2] .- x[3]), 1.0), d(3, 1), [1, 2, 4])
    @test isapprox(g32[1].dev, g3[1].dev; rtol=1e-5)
    # CI without noise types: shown and saved with an empty noise column.
    @test occursin("ci_lo", sprint(show, MIME"text/plain"(), g3[1]))
    path = save(tempname() * ".tsv", g3[1])
    @test load(path).dev == g3[1].dev && isempty(load(path).noise_type)
    # Four clocks: each clock's estimate averages the three couples of the others.
    Q = Matrix{PhaseData}(undef, 4, 4)
    for i in 1:3, j in i+1:4
        Q[i, j] = d(i, j)
    end
    g4 = groslambert(Q, [1, 4])
    @test all(isapprox(g4[i].dev[1], adev(PhaseData(x[i], 1.0), [1]; ci=false).dev[1]; rtol=0.1) for i in 1:4)
    @test g4[1].deviation_type === :adev && length(g4[1].ci) == 2 && g4[1].confidence == DEFAULT_CONFIDENCE
    P2 = Matrix{PhaseData}(undef, 2, 2); P2[1, 2] = d(1, 2)
    @test_throws ArgumentError groslambert(P2, [1])
    @test_throws ArgumentError groslambert(d(1, 2), d(2, 3), PhaseData(x[3][1:100], 1.0), [1])
    # KLTS: the χ² interval for the dominant clock; a lower bound of 0 for the
    # quiet clock once it is not resolved (m = 512), where χ² stays above zero.
    k3 = groslambert(d(1, 2), d(2, 3), d(3, 1), [1, 64, 512]; interval=:klts)
    c3 = groslambert(d(1, 2), d(2, 3), d(3, 1), [1, 64, 512])
    @test all(isapprox(collect(k3[2].ci[q]), collect(c3[2].ci[q]); rtol=0.01) for q in 1:3)
    @test k3[3].ci[3].lo == 0 && c3[3].ci[3].lo > 0 && k3[3].dev == c3[3].dev
    @test_throws ArgumentError groslambert(Q, [1]; interval=:klts)
    @test_throws ArgumentError groslambert(P, [1]; interval=:kltg)
    # The χ² interval is withheld below 1 EDF.
    @test all(isnan, values(SigmaTau._gl_chi2_ci([1.0], [0.5], 0.95)[1]))
end

@testset "KLTS against SigmaTheta" begin
    # SigmaTheta's GCUncert inputs at τ = 100 s and 500 s of a synthetic
    # three-clock run (Groslambert codeviations in units of 1e-13, a negative
    # one for a negative covariance; the closure-noise deviation; the EDF) and
    # its KLTS 97.7 % bounds, with the lower bound of the noisiest clock.
    function st(g, noise, ν)
        a, b, c = (sign(v) * v^2 for v in g); vN = noise^2
        S = [a+b+vN -b -a; -b b+c+vN -c; -a -c c+a+vN]
        q, zero_lower = SigmaTau._klts(S, vN, ν, (0.023, 0.977))
        return sqrt.(q), zero_lower
    end
    q, z = st((1.7971408, 2.2816487, 7.4142426), 0.230, 44.028683)
    @test isapprox(q[:, 2], [3.391053, 3.542078, 9.690462]; rtol=0.02)
    @test isapprox(q[3, 1], 6.234503; rtol=0.02) && z == [true, true, false]
    q, z = st((2.7856007, -2.2631499, 8.7271228), 0.0993, 8.5863312)
    @test isapprox(q[:, 2], [2.857472, 2.565684, 16.37095]; rtol=0.02)
    @test isapprox(q[3, 1], 5.840700; rtol=0.02) && z == [true, true, false]
end

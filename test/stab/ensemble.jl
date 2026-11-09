# ensemble.jl — uncertainties of the N-cornered hat and the Groslambert
# covariance, checked against their scatter over many realizations.

using Test
using Random
using Statistics
using SigmaTau

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
end

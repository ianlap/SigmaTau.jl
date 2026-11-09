# ensemble.jl — uncertainties of the N-cornered hat and the Groslambert
# covariance, checked against their scatter over many realizations.

using Test
using Random
using Statistics
using SigmaTau

@testset "ensemble uncertainties" begin
    rng = Xoshiro(20261108)
    N, R, ms = 1024, 200, [1, 8]
    σ = [1.0, 1.5, 0.8]                        # white PM levels of clocks A, B, C
    hat = [Float64[] for _ in 1:3, _ in ms]; hat_sd = similar(hat)
    gl = [Float64[] for _ in 1:3, _ in ms]; gl_sd = similar(gl)
    for c in eachindex(hat)
        hat_sd[c] = Float64[]; gl_sd[c] = Float64[]
    end
    truth = zeros(3, length(ms))
    for _ in 1:R
        x = [σ[i] .* randn(rng, N) for i in 1:3]
        d(i, j) = PhaseData(x[i] .- x[j], 1.0)
        h = ncornered_hat(adev(d(1, 2), ms), adev(d(2, 3), ms), adev(d(3, 1), ms))
        g = groslambert(d(1, 2), d(2, 3), d(3, 1), ms)
        for i in 1:3, k in eachindex(ms)
            # Predicted 1σ spread of the variance from the equivalent edf.
            push!(hat[i, k], h[i].dev[k]^2); push!(hat_sd[i, k], h[i].dev[k]^2 * sqrt(2 / h[i].edf[k]))
            push!(gl[i, k], g[i].dev[k]^2); push!(gl_sd[i, k], g[i].dev[k]^2 * sqrt(2 / g[i].edf[k]))
            truth[i, k] += adev(PhaseData(x[i], 1.0), ms[k:k]; ci=false).dev[1]^2 / R
        end
    end
    ratio(v, s) = std(v) / mean(s)
    @test all(0.8 < ratio(hat[c], hat_sd[c]) < 1.25 for c in eachindex(hat))
    @test all(0.8 < ratio(gl[c], gl_sd[c]) < 1.25 for c in eachindex(gl))
    @test all(abs(mean(gl[i, k]) / truth[i, k] - 1) < 0.05 for i in 1:3, k in eachindex(ms))
end

@testset "groslambert" begin
    rng = Xoshiro(7)
    x = [s .* randn(rng, 2048) for s in (1.0, 2.0, 0.5, 1.5)]
    d(i, j) = PhaseData(x[i] .- x[j], 1.0)
    g3 = groslambert(d(1, 2), d(2, 3), d(3, 1), [1, 2, 4])
    P = Matrix{PhaseData}(undef, 3, 3)
    P[1, 2] = d(1, 2); P[3, 2] = d(3, 2); P[1, 3] = d(1, 3)   # mixed orientation
    @test all(groslambert(P, [1, 2, 4])[i].dev ≈ g3[i].dev for i in 1:3)
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

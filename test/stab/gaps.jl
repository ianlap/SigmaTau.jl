# gaps.jl — estimators on records with NaN gaps, against direct sums over
# the terms whose samples were all measured. The allantools gradev check is
# in test/validation (gap_rows).

using Test
using Random
using SigmaTau

# Term at start i of each estimator, reading x[i .. i+S].
_bf_span(est, m) = est === :adev ? 2m : est === :hdev ? 3m : est === :mdev ? 3m - 1 :
                   est === :mhdev ? 4m - 1 : m
function _bf_term(est, x, i, m)
    est === :adev  && return x[i+2m] - 2x[i+m] + x[i]
    est === :hdev  && return x[i+3m] - 3x[i+2m] + 3x[i+m] - x[i]
    est === :tierms && return x[i+m] - x[i]
    w = est === :mdev ? (1, -2, 1) : (-1, 3, -3, 1)
    return sum(sum(w[k+1] * x[i+j+k*m] for k in 0:length(w)-1) for j in 0:m-1)
end
_bf_pts(est, m) = est === :adev ? (0, m, 2m) : est === :hdev ? (0, m, 2m, 3m) : (0, m)

# Direct gap-tolerant deviation: `x` phase with NaN at missing samples, or
# `y` frequency (then `x` = its running sum with 0 at gaps and `ybad` marks
# the missing steps).
function _bf_dev(est, x, m, tau0; ybad=nothing)
    s, n = 0.0, 0
    S = _bf_span(est, m)
    for i in 1:length(x)-S
        ok = if ybad !== nothing
            !any(ybad[i:i+S-1])
        elseif est in (:mdev, :mhdev)
            all(isfinite, x[i:i+S])
        else
            all(o -> isfinite(x[i+o]), _bf_pts(est, m))
        end
        ok || continue
        s += _bf_term(est, x, i, m)^2
        n += 1
    end
    c = est === :adev ? 2m^2 : est === :hdev ? 6m^2 : est === :mdev ? 2m^4 : est === :mhdev ? 6m^4 : 1
    return sqrt(s / (c * n * (est === :tierms ? 1 : tau0^2))), n
end

# Pairs at most m apart inside a run that holds a full window of m + 1 samples.
function _bf_mtie(x, m; ybad=nothing)
    best = 0.0
    run = ybad === nothing ? zeros(Int, length(x)) : [0; cumsum(ybad)]
    for a in eachindex(x), b in a+1:min(a + m, length(x))
        isfinite(x[a]) && isfinite(x[b]) && run[a] == run[b] || continue
        count(==(run[a]), run) > m || continue
        best = max(best, abs(x[b] - x[a]))
    end
    return best
end

@testset "gap-tolerant estimators" begin
    rng = MersenneTwister(20261108)
    N, tau0 = 1500, 0.5
    x = cumsum(randn(rng, N)) .* 1e-9 .+ randn(rng, N) .* 1e-9 .+ 1e-12 .* (1:N) .^ 2
    xg = copy(x); xg[90:130] .= NaN; xg[400] = NaN; xg[1001:1006] .= NaN
    y = diff(x) ./ tau0
    yg = copy(y); yg[50:70] .= NaN; yg[700] = NaN
    ybad = isnan.(yg)
    xy = [0.0; cumsum(replace(yg, NaN => 0.0))] .* tau0
    ms = [1, 2, 3, 8, 30, 100]
    pg, fg = PhaseData(xg, tau0), FrequencyData(yg, tau0)

    @testset "$est matches the direct sum" for (est, f) in
            ((:adev, adev), (:mdev, mdev), (:hdev, hdev), (:mhdev, mhdev), (:tierms, tierms))
        rp, rf = f(pg, ms; ci=false), f(fg, ms; ci=false)
        bp = [_bf_dev(est, xg, m, tau0) for m in ms]
        bf = [_bf_dev(est, xy, m, tau0; ybad=ybad) for m in ms]
        @test rp.dev ≈ first.(bp) rtol=1e-9
        @test rp.neff == last.(bp)
        @test rf.dev ≈ first.(bf) rtol=1e-9
        @test rf.neff == last.(bf)
    end

    @testset "mtie over the measured samples" begin
        @test mtie(pg, ms).dev ≈ [_bf_mtie(xg, m) for m in ms]
        @test mtie(fg, ms).dev ≈ [_bf_mtie(xy, m; ybad=ybad) for m in ms]
        # A run shorter than the window adds nothing, however large its swing.
        ys = copy(yg); ys[300:303] .= NaN; ys[305:307] .= NaN; ys[304] = 1.0
        sbad = isnan.(ys)
        xs = [0.0; cumsum(replace(ys, NaN => 0.0))] .* tau0
        @test mtie(FrequencyData(ys, tau0), ms).dev ≈ [_bf_mtie(xs, m; ybad=sbad) for m in ms]
    end

    @testset "tdev and htdev scale the gap-tolerant mdev and mhdev" begin
        for d in (pg, fg)
            @test tdev(d, ms).dev ≈ mdev(d, ms).dev .* ms .* tau0 ./ sqrt(3)
            @test htdev(d, ms).dev ≈ mhdev(d, ms).dev .* ms .* tau0 ./ sqrt(10 / 3)
        end
    end

    @testset "confidence from the terms used" begin
        full, gapped = adev(PhaseData(x, tau0), ms), adev(pg, ms)
        @test all(isfinite, gapped.edf)
        @test all(gapped.edf .< full.edf)
        @test all(c.lo < d < c.hi for (c, d) in zip(gapped.ci, gapped.dev))
        @test gapped.noise_type == full.noise_type
    end

    @testset "a gap-free record takes the same path either way" begin
        g = SigmaTau._gap_record(PhaseData(x, tau0))
        for est in (:adev, :mhdev, :mtie)
            a = SigmaTau._gap_result(est, g, ms)
            b = getfield(SigmaTau, est)(PhaseData(x, tau0), ms)
            @test a.dev ≈ b.dev rtol=1e-12
            @test a.neff == b.neff
            @test a.edf ≈ b.edf
        end
    end

    @testset "too few terms give NaN" begin
        short = fill(NaN, 40); short[1:5] .= 1.0:5.0; short[30:32] .= 0.0
        r = adev(PhaseData(short, 1.0), [1, 2, 4])
        @test r.neff == [4, 0, 0]
        @test isnan(r.dev[2]) && isnan(r.dev[3])
    end
end

using Test
using Random
using SigmaTau
using Tables

@testset "Tables.jl extension" begin
    pd = PhaseData(cumsum(randn(MersenneTwister(23), 256)) .* 1e-9, 1.0)
    r  = adev(pd; ci=false)
    rc = adev(pd; ci=true)
    @test Tables.istable(typeof(r)) && Tables.rowaccess(typeof(r))
    @test Tables.schema(r).names ==
          (:deviation_type, :tau, :dev, :noise_type, :ci_lower, :ci_upper, :edf, :neff)

    # Without CI the CI columns are missing; neff is always there.
    c = Tables.columntable(r)
    @test (c.tau, c.dev, c.neff) == (r.tau, r.dev, r.neff)
    @test all(c.deviation_type .=== :adev) && eltype(c.neff) == Int
    @test all(ismissing, c.noise_type) && all(ismissing, c.ci_lower) &&
          all(ismissing, c.ci_upper) && all(ismissing, c.edf)
    cc = Tables.columntable(rc)
    @test (cc.noise_type, cc.ci_lower, cc.ci_upper, cc.edf) ==
          (rc.noise_type, ci_lower(rc), ci_upper(rc), rc.edf)

    # A suite stacks its results in order.
    s = stability(pd; devs=(:adev, :mdev), m=Octave, ci=false)
    @test Tables.istable(typeof(s)) && Tables.rowaccess(typeof(s))
    sc = Tables.columntable(s)
    @test sc.deviation_type == [fill(:adev, length(s[:adev].tau)); fill(:mdev, length(s[:mdev].tau))]
    @test (sc.tau, sc.neff) == ([s[:adev].tau; s[:mdev].tau], [s[:adev].neff; s[:mdev].neff])
    @test length(Tables.rows(stability(pd; devs=(), ci=false))) == 0

    # A hat stacks its clocks, each row named.
    h = ncornered_hat(ClockComparisons(randn(MersenneTwister(4), 256, 2), 1.0; names=["a", "b"]), [1, 2, 4])
    hc = Tables.columntable(h)
    @test Tables.schema(h).names[1] === :clock && hc.clock == repeat(["a", "b", "reference"]; inner=3)
    @test (hc.dev, hc.edf) == (reduce(vcat, r.dev for r in h), reduce(vcat, r.edf for r in h))
end

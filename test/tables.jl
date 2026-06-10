using Test
using Random
using SigmaTau
using Tables

@testset "Tables.jl extension" begin
    Random.seed!(23)
    pd = PhaseData(cumsum(randn(256)) .* 1e-9, 1.0)

    result = adev(pd; ci=false)
    @test Tables.istable(typeof(result))
    @test Tables.rowaccess(typeof(result))
    @test Tables.schema(result).names == (
        :deviation_type,
        :tau,
        :dev,
        :noise_type,
        :ci_lower,
        :ci_upper,
        :edf,
        :neff,
    )

    rows = collect(Tables.rows(result))
    @test length(rows) == length(result.tau)
    @test rows[1].deviation_type === :adev
    @test rows[1].tau == result.tau[1]
    @test rows[1].dev == result.dev[1]
    @test ismissing(rows[1].noise_type)
    @test ismissing(rows[1].ci_lower)
    @test ismissing(rows[1].ci_upper)
    @test ismissing(rows[1].edf)
    # neff is populated even on a ci=false result.
    @test rows[1].neff == result.neff[1]
    @test rows[1].neff isa Int

    result_ci = adev(pd; ci=true)
    rows_ci = collect(Tables.rows(result_ci))
    @test rows_ci[1].noise_type == result_ci.noise_type[1]
    @test rows_ci[1].ci_lower == result_ci.ci[1].lo
    @test rows_ci[1].ci_upper == result_ci.ci[1].hi
    @test rows_ci[1].edf == result_ci.edf[1]
    @test rows_ci[1].neff == result_ci.neff[1]
    @test [row.neff for row in rows_ci] == result_ci.neff

    suite = stability(pd; devs=(:adev, :mdev), taus=Octave, ci=false)
    @test Tables.istable(typeof(suite))
    @test Tables.rowaccess(typeof(suite))

    suite_rows = collect(Tables.rows(suite))
    @test length(suite_rows) == length(suite[:adev].tau) + length(suite[:mdev].tau)
    @test suite_rows[1].deviation_type === :adev
    @test suite_rows[length(suite[:adev].tau) + 1].deviation_type === :mdev

    cols = Tables.columntable(suite)
    @test cols.deviation_type == [row.deviation_type for row in suite_rows]
    @test cols.tau == [row.tau for row in suite_rows]
    @test cols.dev == [row.dev for row in suite_rows]
    @test cols.neff == vcat(suite[:adev].neff, suite[:mdev].neff)

    # A bare StabilityResult also round-trips through columntable (rowaccess +
    # schema fallback on ResultRows).
    result_cols = Tables.columntable(result)
    @test result_cols.tau == result.tau
    @test result_cols.dev == result.dev

    # An empty suite (devs=()) reports zero rows and must not throw — the
    # SuiteRows length reduction seeds init=0 for the empty case.
    empty_suite = stability(pd; devs=(), ci=false)
    @test length(Tables.rows(empty_suite)) == 0
    @test length(Tables.columntable(empty_suite).tau) == 0
end

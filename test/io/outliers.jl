# test/io/outliers.jl — find_outliers (Stable32-style MAD check) and
# remove_outliers (flag → NaN → Howe fill). Fixtures are deterministic:
# a sinusoidal base record with planted glitches.

@testset "find_outliers / remove_outliers" begin
    # Deterministic base records (no RNG): a slow sinusoid at the 1e-12 level.
    N      = 128
    y_base = sinpi.((1:N) ./ 7.3) .* 1e-12
    x_base = cumsum(y_base)               # tau0 = 1 ⇒ phase is the running sum

    @testset "FrequencyData: planted glitch is flagged" begin
        y     = copy(y_base)
        y[50] = 1e-9                       # counter glitch, ~1000x the base level
        fd    = FrequencyData(y, 1.0)
        @test find_outliers(fd) == [50]
        # Clean record flags nothing.
        @test isempty(find_outliers(FrequencyData(y_base, 1.0)))
        # Tighter threshold flags at least as much.
        @test 50 in find_outliers(fd; nsigma=3)
        @test_throws ArgumentError find_outliers(fd; nsigma=0)
    end

    @testset "PhaseData: glitch flags the two samples around the steps" begin
        # An isolated phase glitch at x[40] produces a step up in Δ₃₉ and a
        # step down in Δ₄₀; the documented convention flags sample i+1 per
        # flagged difference Δᵢ, i.e. samples 40 and 41.
        x     = copy(x_base)
        x[40] += 5e-9
        pd    = PhaseData(x, 1.0)
        @test find_outliers(pd) == [40, 41]

        # A persistent phase step (all samples shifted from index 70 on) is a
        # single frequency outlier: only sample 70 — where the step lands —
        # is flagged.
        xs = copy(x_base)
        xs[70:end] .+= 5e-9
        @test find_outliers(PhaseData(xs, 1.0)) == [70]

        # Clean record flags nothing.
        @test isempty(find_outliers(PhaseData(x_base, 1.0)))
    end

    @testset "NaN samples are excluded and never flagged" begin
        y     = copy(y_base)
        y[10] = NaN
        y[80] = 2e-9
        idx   = find_outliers(FrequencyData(y, 1.0))
        @test idx == [80]                  # the gap is not an outlier
    end

    @testset "remove_outliers round-trip (frequency)" begin
        y     = copy(y_base)
        y[50] = 1e-9
        fd    = FrequencyData(y, 2.0; source="glitchy.csv")
        clean = remove_outliers(fd)
        @test clean isa FrequencyData
        @test all(isfinite, clean.y)
        @test isempty(find_outliers(clean))
        @test clean.y[49] == y[49]                 # neighbours untouched
        @test abs(clean.y[50]) < 1e-10             # glitch replaced
        @test clean.tau0 == 2.0                    # metadata carries over
        @test clean.source == "glitchy.csv"
        @test fd.y[50] == 1e-9                     # input untouched
        # Nothing flagged ⇒ the input comes back as-is.
        fd_ok = FrequencyData(y_base, 1.0)
        @test remove_outliers(fd_ok) === fd_ok
    end

    @testset "remove_outliers round-trip (phase)" begin
        x     = copy(x_base)
        x[40] += 5e-9
        pd    = PhaseData(x, 1.0; source="glitchy.dat")
        clean = remove_outliers(pd)
        @test clean isa PhaseData
        @test all(isfinite, clean.x)
        @test isempty(find_outliers(clean))
        @test abs(clean.x[40] - x_base[40]) < 1e-9  # glitch gone
        @test clean.source == "glitchy.dat"
        @test pd.x[40] == x[40]                     # input untouched
    end
end

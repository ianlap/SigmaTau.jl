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
        @test clean.y[49] == y[49]                 # neighbors untouched
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

    @testset "remove_outliers removes a phase step (frequency space)" begin
        # A persistent phase step is one outlying frequency sample. Cleaning in
        # frequency space and reintegrating removes the step from every later
        # sample; refilling the phase sample it lands on would leave it there.
        step = 5e-9
        xs   = copy(x_base)
        xs[70:end] .+= step
        clean = remove_outliers(PhaseData(xs, 1.0; source="step.dat"))
        @test clean isa PhaseData
        @test all(isfinite, clean.x)
        @test clean.x[1:69] ≈ xs[1:69] rtol=1e-12      # before the step: unchanged
        @test maximum(abs.(clean.x .- x_base)) < 1e-2 * step
        @test isempty(find_outliers(clean))
        @test clean.source == "step.dat"

        # ADEV recovers to the clean record's level instead of carrying the step.
        ms     = [1, 2, 4, 8]
        a_base = adev(PhaseData(x_base, 1.0), ms; ci=false).dev
        a_step = adev(PhaseData(xs, 1.0), ms; ci=false).dev
        a_fix  = adev(clean, ms; ci=false).dev
        @test all(a_step .> 100 .* a_base)
        @test all(a_fix .< 2 .* a_base)
    end

    @testset "remove_outliers: phase result matches the frequency path" begin
        # A frequency spike integrated to phase is a phase step. Cleaning the
        # phase record must equal cleaning the frequency record and integrating.
        tau0  = 2.0
        y     = copy(y_base)
        y[50] = 1e-9
        x     = vcat(0.0, cumsum(y) .* tau0)
        via_phase = remove_outliers(PhaseData(x, tau0))
        fd_clean  = remove_outliers(FrequencyData(y, tau0))
        via_freq  = vcat(0.0, cumsum(fd_clean.y) .* tau0)
        @test via_phase.x ≈ via_freq rtol=1e-12 atol=1e-24
        @test via_phase.tau0 == tau0
    end

    @testset "a quantized record has no scale to test against" begin
        steps = repeat([0.0, 0.0, 0.0, 1e-9, -1e-9], 40)      # MAD of the steps is 0
        @test isempty(find_outliers(FrequencyData(steps, 1.0)))
        @test isempty(find_outliers(PhaseData(cumsum(steps), 1.0)))
    end

    @testset "remove_outliers: samples already missing stay missing" begin
        xs = copy(x_base)
        xs[90:end] .+= 5e-9
        xs[20] = NaN
        clean = remove_outliers(PhaseData(xs, 1.0))
        @test findall(isnan, clean.x) == [20]
        @test abs(clean.x[end] - x_base[end]) < 5e-11
        y = copy(y_base)
        y[20] = NaN
        y[50] = 1e-9
        clean = remove_outliers(FrequencyData(y, 1.0))
        @test findall(isnan, clean.y) == [20]
        @test abs(clean.y[50]) < 1e-10
        # Nothing flagged ⇒ the input comes back as-is.
        pd_ok = PhaseData(x_base, 1.0)
        @test remove_outliers(pd_ok) === pd_ok
    end
end

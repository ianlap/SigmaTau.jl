@testset "fillgaps (Howe and Schlossberger)" begin
    N, gap = 4096, 1401:2600                      # a 30 % gap
    holed(v) = (w = copy(v); w[gap] .= NaN; w)
    # Over 300 seeds the ratios bounded below run from 0.95 to 1.05 and the
    # flicker FM contrast from 1.09 to 1.37; the bounds stand clear of both.
    ratio(f, o, ms) = adev(f, ms; ci=false).dev ./ adev(o, ms; ci=false).dev

    @testset "nothing to fill, and a gap of one sample" begin
        x = cumsum(randn(Xoshiro(1), 256))
        @test fillgaps(PhaseData(x, 2.0)).x == x && fillgaps(PhaseData(x, 2.0)).tau0 == 2.0
        x[32] = NaN
        @test fillgaps(PhaseData(x, 1.0)).x[32] ≈ (x[31] + x[33]) / 2
        @test_throws ArgumentError fillgaps(PhaseData(x, 1.0); form=:spline)
    end

    @testset "white PM on a frequency offset: filled as phase, around its trend" begin
        x = randn(Xoshiro(1), N) .+ 100.0 .* (0:N-1)
        d = PhaseData(holed(x), 1.0)
        f = fillgaps(d)
        @test f.x == fillgaps(d; form=:phase).x
        @test f.x[.!isnan.(d.x)] == x[.!isnan.(d.x)]              # measured samples as they came
        @test all(r -> 0.9 < r < 1.1, ratio(f, PhaseData(x, 1.0), [1, 16, 256]))
    end

    @testset "random-walk FM phase: filled on its first differences" begin
        x = cumsum(cumsum(randn(Xoshiro(1), N)))
        d = PhaseData(holed(x), 1.0)
        f = fillgaps(d)
        @test f.x == fillgaps(d; form=:frequency).x
        @test f.x[.!isnan.(d.x)] == x[.!isnan.(d.x)]              # the differences sum to the measured change
        @test 0.9 < only(ratio(f, PhaseData(x, 1.0), [1])) < 1.1
    end

    @testset "flicker FM past many short gaps: its deviation does not fall, filled as frequency" begin
        x = noise_gen(PhaseData, N, 1.0; sigma1=Dict(-1 => 1.0), rng=Xoshiro(1)).x
        d = PhaseData(copy(x), 1.0)
        d.x[[k + j for k in 20:40:N-40 for j in 0:11]] .= NaN     # 12 of every 40 lost
        f = fillgaps(d)
        @test f.x == fillgaps(d; form=:frequency).x
        @test 0.9 < only(ratio(f, PhaseData(x, 1.0), [1])) < 1.1
        @test only(ratio(fillgaps(d; form=:phase), PhaseData(x, 1.0), [1])) > 1.05   # a step at every join
    end

    @testset "white PM frequency record: filled as phase, the step across the gap free" begin
        y = diff(randn(Xoshiro(1), N))
        d = FrequencyData(holed(y), 1.0)
        f = fillgaps(d)
        @test f.y == fillgaps(d; form=:phase).y
        @test all(r -> 0.9 < r < 1.1, ratio(f, FrequencyData(y, 1.0), [1, 16, 256]))
    end

    @testset "gaps at the record ends join the data they touch" begin
        x = cumsum(randn(Xoshiro(8), 1000))
        x[1:30] .= NaN; x[500:520] .= NaN; x[961:1000] .= NaN
        for form in (:phase, :frequency)
            f = fillgaps(PhaseData(x, 1.0); form=form).x
            @test all(isfinite, f)
            @test abs(f[31] - f[30]) < 4 && abs(f[961] - f[960]) < 4
        end
        y = copy(x)
        @test all(isfinite, fillgaps(FrequencyData(y, 1.0); form=:phase).y)
    end

    # A fill is the samples beside it reversed and inverted, plus the line
    # that joins it and the parabola set aside as the trend: added to the
    # samples in their order, it leaves no third difference.
    bend(fill, source) = maximum(abs, diff(diff(diff(fill .+ reverse(source)))))

    @testset "outward from the longest run: the gap after it from the left, the gap before it from the right" begin
        x = randn(Xoshiro(2), 512)
        d = copy(x); d[61:100] .= NaN; d[401:440] .= NaN          # runs of 60, 300 and 72
        f = fillgaps(PhaseData(d, 1.0); form=:phase).x
        @test bend(f[401:440], x[361:400]) < 1e-9
        @test bend(f[61:100], x[101:140]) < 1e-9
    end

    @testset "the longest run may be the first, and a fill serves the next gap" begin
        d = randn(Xoshiro(2), 305)
        d[101:160] .= NaN; d[211:265] .= NaN                      # runs of 100, 50 and 40
        f = fillgaps(PhaseData(d, 1.0); form=:phase).x
        @test bend(f[101:160], d[41:100]) < 1e-9
        @test bend(f[211:265], f[156:210]) < 1e-9                 # 5 filled samples and the run of 50
    end

    @testset "a gap longer than the runs beside it is returned to" begin
        d = fill(NaN, 2000)
        for run in (1:150, 601:640, 1101:1400, 1801:2000)         # gaps of 450, 460 and 400
            d[run] = randn(Xoshiro(first(run)), length(run))
        end
        @test all(isfinite, fillgaps(PhaseData(d, 1.0); form=:phase).x)
    end

    @testset "an unfillable layout throws instead of hanging" begin
        # Neither gap has enough data on either side to mirror.
        x = fill(NaN, 40); x[1:5] .= 1.0:5.0; x[30:32] .= 0.0
        @test_throws ErrorException fillgaps(PhaseData(x, 1.0))
    end
end

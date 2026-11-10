@testset "fillgaps (Howe and Schlossberger)" begin
    N, gap = 4096, 1401:2600                      # a 30 % gap
    holed(v) = (w = copy(v); w[gap] .= NaN; w)
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
        @test all(r -> 0.95 < r < 1.05, ratio(f, PhaseData(x, 1.0), [1, 16, 256]))
    end

    @testset "random-walk FM phase: filled on its first differences" begin
        x = cumsum(cumsum(randn(Xoshiro(1), N)))
        d = PhaseData(holed(x), 1.0)
        f = fillgaps(d)
        @test f.x == fillgaps(d; form=:frequency).x
        @test f.x[.!isnan.(d.x)] == x[.!isnan.(d.x)]              # the differences sum to the measured change
        @test 0.95 < only(ratio(f, PhaseData(x, 1.0), [1])) < 1.05
        @test only(ratio(fillgaps(d; form=:phase), PhaseData(x, 1.0), [1])) > 1.2   # the joins show
    end

    @testset "flicker FM past many short gaps: its deviation does not fall, filled as frequency" begin
        x = noise_gen(PhaseData, N, 1.0; sigma1=Dict(-1 => 1.0), rng=Xoshiro(1)).x
        d = PhaseData(copy(x), 1.0)
        d.x[[k + j for k in 20:40:N-40 for j in 0:11]] .= NaN     # 12 of every 40 lost
        f = fillgaps(d)
        @test f.x == fillgaps(d; form=:frequency).x
        @test 0.95 < only(ratio(f, PhaseData(x, 1.0), [1])) < 1.05
        @test only(ratio(fillgaps(d; form=:phase), PhaseData(x, 1.0), [1])) > 1.1
    end

    @testset "white PM frequency record: filled as phase, the step across the gap free" begin
        y = diff(randn(Xoshiro(1), N))
        d = FrequencyData(holed(y), 1.0)
        f = fillgaps(d)
        @test f.y == fillgaps(d; form=:phase).y
        @test all(r -> 0.95 < r < 1.05, ratio(f, FrequencyData(y, 1.0), [1, 16, 256]))
        @test ratio(fillgaps(d; form=:frequency), FrequencyData(y, 1.0), [256])[1] > 2   # a level error, integrated
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

    @testset "an unfillable layout throws instead of hanging" begin
        # Neither gap has enough data on either side to mirror.
        x = fill(NaN, 40); x[1:5] .= 1.0:5.0; x[30:32] .= 0.0
        @test_throws ErrorException fillgaps(PhaseData(x, 1.0))
    end
end

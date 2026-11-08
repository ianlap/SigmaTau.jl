const FMT = joinpath(@__DIR__, "..", "fixtures", "formats")
fmt(f) = joinpath(FMT, f)

@testset "instrument formats" begin
    @testset "detect_format" begin
        @test detect_format(fmt("timelab.tim")) === :timelab
        @test detect_format(fmt("stable32.phd")) === :stable32
        @test detect_format(fmt("stable32.frd")) === :stable32
        @test detect_format(fmt("counter_53230a_freq.csv")) === :counter
        @test detect_format(fmt("counter_53230a_ti.txt")) === :counter
        @test detect_format(fmt("gpsdo.csv")) === :table
        @test detect_format(joinpath(@__DIR__, "..", "fixtures", "validation", "stable32gen.DAT")) === :columns
    end

    @testset "TimeLab .tim" begin
        p = read_data(fmt("timelab.tim"))
        raw = [parse(Float64, l) for l in eachline(fmt("timelab.tim")) if occursin(r"^\s+[-\d]", l)]
        @test p isa PhaseData
        @test p.tau0 == 0.1
        @test p.x == raw && length(p) == 64
        @test read_data(fmt("timelab.tim"); tau0=0.2).tau0 == 0.2
    end

    @testset "Stable32 .phd / .frd" begin
        p = read_data(fmt("stable32.phd"))
        @test p isa PhaseData && length(p) == 32
        @test p.tau0 == 1.0                        # header Tau:
        f = read_data(fmt("stable32.frd"); tau0=10.0)
        @test f isa FrequencyData && f.tau0 == 10.0
        @test_throws ArgumentError read_data(fmt("stable32.frd"))   # no τ₀ anywhere
    end

    @testset "counter logs" begin
        y = read_data(fmt("counter_53230a_freq.csv"); tau0=1.0)
        f = [parse(Float64, l) for l in eachline(fmt("counter_53230a_freq.csv"))]
        @test y isa FrequencyData
        @test isnan(y.y[18])                       # 9.91E+37 = no reading
        @test y.y[1] ≈ (f[1] - 10e6) / 10e6 rtol=1e-12
        @test read_data(fmt("counter_53230a_freq.csv"); tau0=1.0, nominal=1e7).y[1] == y.y[1]
        x = read_data(fmt("counter_53230a_ti.txt"); tau0=1.0)
        @test x isa PhaseData
        @test x.x[1] ≈ 1.2345e-7 rtol=1e-3
        @test_throws ArgumentError read_data(fmt("counter_53230a_ti.txt"))  # gate time needed
    end

    @testset "named columns (GPSDO log)" begin
        p = read_data(fmt("gpsdo.csv"))
        rows = [split(l, ',') for l in Iterators.drop(eachline(fmt("gpsdo.csv")), 1)]
        @test p isa PhaseData
        @test p.tau0 == 1.0                        # ISO timestamps
        @test p.x ≈ parse.(Float64, getindex.(rows, 2)) .* 1e-9
        g = read_data(fmt("gpsdo.csv"); fillgaps=true)
        @test length(g) == length(p) + 1           # the missing second, filled
        o = read_data(fmt("gpsdo.csv"); column="OSC offset (ppb)")
        @test o isa FrequencyData && o.y ≈ parse.(Float64, getindex.(rows, 3)) .* 1e-9
        @test read_data(fmt("gpsdo.csv"); column=4, kind=:phase) isa PhaseData
        @test_throws ArgumentError read_data(fmt("gpsdo.csv"); column="nope")
    end
end

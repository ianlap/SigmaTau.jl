using Dates

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
        @test detect_format(fmt("bipm_utc_lab.txt")) === :bipm_utc
        @test detect_format(fmt("rinex_clk_synth.clk")) === :rinex_clk
        @test detect_format(joinpath(@__DIR__, "..", "fixtures", "validation", "stable32gen.DAT")) === :columns
    end

    @testset "BIPM UTC - UTC(k)" begin
        # CRLF, a "+/-" uncertainty column in the later rows, a DOS end mark, two rows missing.
        p = read_data(fmt("bipm_utc_lab.txt"))
        @test p isa PhaseData && p.tau0 == 5 * 86400.0
        @test length(p) == 60 && count(isnan, p.x) == 2          # the two missing rows stay empty slots
        rows = [split(l) for l in eachline(fmt("bipm_utc_lab.txt")) if occursin(r"^\d", l)]
        @test p.x[1:3] == [parse(Float64, r[2]) * 1e-9 for r in rows[1:3]]     # ns -> s
        q = read_data(fmt("bipm_utc_lab.txt"); fillgaps=true)
        @test !any(isnan, q.x) && length(q) == 60
        @test start_time(p) == DateTime(1998, 1, 1)               # MJD 50814
    end

    @testset "RINEX clock" begin
        # Three clocks (G01 out of view for two epochs) on a 30 s grid.
        g = read_data(fmt("rinex_clk_synth.clk"))                       # the first clock
        @test g isa PhaseData && g.tau0 == 30 && length(g) == 40 && count(isnan, g.x) == 2
        @test g.x[1] == 1.5e-4
        @test start_time(g) == DateTime(2025, 8, 17, 12)
        e = read_data(fmt("rinex_clk_synth.clk"); column="e02")          # by name, any case
        @test e.tau0 == 30 && !any(isnan, e.x) && e.x[1] == 8.0e-5
        @test read_data(fmt("rinex_clk_synth.clk"); column=2).x[1] == -3.2e-7   # by position: the station
        @test !any(isnan, read_data(fmt("rinex_clk_synth.clk"); fillgaps=true).x)
        @test_throws "among the 3 clocks" read_data(fmt("rinex_clk_synth.clk"); column="G99")
        @test_throws "exceeds the 3 clocks" read_data(fmt("rinex_clk_synth.clk"); column=9999)
        # Old IGS files write the exponent with D; the record runs on past midnight.
        old = tempname() * ".clk"
        write(old, """
                 2.00           C                   G                   RINEX VERSION / TYPE
                                                                        END OF HEADER
            AS G01  2025 08 17 23 59 30.000000  2   -.123000000000D-03  .250000000000D-10
            AS G01  2025 08 18 00 00  0.000000  2   -.123100000000D-03  .250000000000D-10
            AS G01  2025 08 18 00 00 30.000000  2   -.123200000000D-03  .250000000000D-10
            """)
        o = read_data(old)
        @test o.tau0 == 30 && o.x == [-1.230e-4, -1.231e-4, -1.232e-4]
        @test start_time(o) == DateTime(2025, 8, 17, 23, 59, 30)
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
        @test start_time(p) == DateTime(2026, 11, 8)    # MJD time tags
        f = read_data(fmt("stable32.frd"); tau0=10.0)
        @test f isa FrequencyData && f.tau0 == 10.0
        @test_throws ArgumentError read_data(fmt("stable32.frd"))   # no τ₀ anywhere
    end

    @testset "counter logs" begin
        y = read_data(fmt("counter_53230a_freq.csv"); tau0=1.0)
        f = [parse(Float64, l) for l in eachline(fmt("counter_53230a_freq.csv"))]
        @test y isa FrequencyData
        @test isnan(y.y[18])                       # 9.91E+37 = no reading
        @test all(isfinite, read_data(fmt("counter_53230a_freq.csv"); tau0=1.0, fillgaps=true).y)
        @test y.y[1] ≈ (f[1] - 10e6) / 10e6 rtol=1e-12
        @test read_data(fmt("counter_53230a_freq.csv"); tau0=1.0, nominal=1e7).y[1] == y.y[1]
        x = read_data(fmt("counter_53230a_ti.txt"); tau0=1.0)
        @test x isa PhaseData
        @test x.x[1] ≈ 1.2345e-7 rtol=1e-3
        @test_throws ArgumentError read_data(fmt("counter_53230a_ti.txt"))  # gate time needed
        # Sub-hertz readings look like time intervals unless the caller says.
        slow = tempname() * ".txt"
        write(slow, "+5.00000001E-01\n+4.99999998E-01\n+5.00000002E-01\n")
        @test read_data(slow; tau0=1.0) isa PhaseData
        @test read_data(slow; tau0=1.0, kind=:frequency).y ≈ [2e-9, -4e-9, 4e-9] rtol=1e-6
    end

    @testset "named columns (GPSDO log)" begin
        p = read_data(fmt("gpsdo.csv"))
        rows = [split(l, ',') for l in Iterators.drop(eachline(fmt("gpsdo.csv")), 1)]
        @test p isa PhaseData
        @test p.tau0 == 1.0                        # ISO timestamps
        @test start_time(p) == DateTime(2026, 11, 8, 12)
        @test isnan(p.x[24]) && length(p) == length(rows) + 1    # the missing second
        @test filter(!isnan, p.x) ≈ parse.(Float64, getindex.(rows, 2)) .* 1e-9
        g = read_data(fmt("gpsdo.csv"); fillgaps=true)
        @test length(g) == length(p) && all(isfinite, g.x)
        o = read_data(fmt("gpsdo.csv"); column="OSC offset (ppb)")
        @test o isa FrequencyData && filter(!isnan, o.y) ≈ parse.(Float64, getindex.(rows, 3)) .* 1e-9
        @test read_data(fmt("gpsdo.csv"); column=4, kind=:phase) isa PhaseData
        @test_throws ArgumentError read_data(fmt("gpsdo.csv"); column="OSC offset (ppb)", kind=:phase)
        @test_throws ArgumentError read_data(fmt("gpsdo.csv"); column="PPS offset (ns)", kind=:frequency)
        @test_throws ArgumentError read_data(fmt("gpsdo.csv"); column="nope")
        @test SigmaTau._read_data(fmt("gpsdo.csv"), :auto; fillgaps=true).filled == (1:41 .== 24)
        semi = tempname() * ".txt"
        write(semi, "time;phase (ns)\n0;1.0\n1;2.0\n2;4.0\n")
        @test read_data(semi; delim=';').x == [1e-9, 2e-9, 4e-9]      # detected with the delimiter
        write(semi, "0;1\n1;2\n2;4\n")
        @test read_data(semi; delim=';').x == [1, 2, 4]               # plain columns keep it too
        write(semi, "0,1\n1,2\n2,4\n")
        @test read_data(semi).x == [1, 2, 4]                          # a comma found in a .txt
        # The unit in a name fixes the record type over its words, and a time
        # unit scales the time column.
        units = tempname() * ".csv"
        write(units, "time (ms),drift ppm,osc (ms),freq (ns/s)\n0,1,1,1\n100,2,2,2\n200,4,4,4\n")
        d = read_data(units; column="drift ppm")
        @test d isa FrequencyData && d.y ≈ [1e-6, 2e-6, 4e-6] && d.tau0 == 0.1
        @test read_data(units; column="osc (ms)") isa PhaseData
        @test read_data(units; column="freq (ns/s)").y ≈ [1e-9, 2e-9, 4e-9]
        @test_throws ArgumentError read_data(units; column="drift ppm", kind=:phase)
        # A time column found by its ISO timestamps rather than its name.
        iso = tempname() * ".csv"
        write(iso, "measurement_time,phase (ns)\n2026-11-08T12:00:00Z,1\n2026-11-08T12:00:01Z,2\n2026-11-08T12:00:03Z,4\n")
        @test length(read_data(iso; fillgaps=true)) == 4
    end

    @testset "time grid" begin
        # Jittered timestamps keep every reading in order; a missing one
        # leaves a NaN slot; a record with none reads as written.
        log = tempname() * ".csv"
        write(log, "time,phase (s)\n0.02,1\n0.98,2\n2.01,3\n4.03,5\n4.97,6\n")
        p = read_data(log)
        @test p.tau0 == 1.0                        # the fitted step, to a round value
        @test isequal(p.x, [1, 2, 3, NaN, 5, 6])
        @test start_time(p) === nothing
        @test read_data(log; fillgaps=true).x[[1, 2, 3, 5, 6]] == [1, 2, 3, 5, 6]
        @test SigmaTau._time_grid([0.0, 1.0, 3.1, 4.0]) == ([1, 2, 4, 5], 1.0)
        # Exact timestamps keep their own step on a short record, not a rounder one.
        @test SigmaTau._time_grid(86400.0 .* (0:9))[2] == 86400.0     # not 86000
        @test SigmaTau._time_grid(604800.0 .* (0:99))[2] == 604800.0  # not 605000
        # MJD to six decimals (0.0864 s): 1 s steps read as 0.9504 or 1.0368 s.
        mjd = round.(59646.9583 .+ [0:999; 1500:2999] ./ 86400; digits=6)
        tbl = tempname() * ".csv"
        write(tbl, "MJD,ratio\n" * join(string.(mjd, ",", 1:2500), "\n"))
        q = read_data(tbl; column=2, kind=:frequency)
        @test q.tau0 == 1.0 && length(q) == 3000 && count(isnan, q.y) == 500
        # MJD by its range under a column headed `t`, as optical clock
        # comparisons write it, both as plain columns and as a table; a
        # name that says seconds keeps seconds.
        ratio = string.(mjd, " ", 1e-15 .* (1:2500), " 1")
        dat = tempname() * ".dat"
        write(dat, "# A-B comparison\n# t  ratio  flag\n" * join(ratio, "\n"))
        c = read_data(dat)
        @test c isa FrequencyData && c.tau0 == 1.0 && length(c) == 3000
        @test start_time(c) == DateTime(2022, 3, 8, 22, 59, 57, 120)
        write(tbl, "t,ratio\n" * join(string.(mjd, ",", 1e-15 .* (1:2500)), "\n"))
        @test read_data(tbl) isa FrequencyData && read_data(tbl).tau0 == 1.0
        # Seconds of day at 10 Hz lie in the MJD range too, and stay seconds.
        sod = string.(43200 .+ (0:99) ./ 10, " ", 1e-9 .* (1:100))
        write(dat, join(sod, "\n"))
        @test read_data(dat).tau0 == 0.1 && start_time(read_data(dat)) === nothing
        for name in ("time", "t", "UTC")
            write(tbl, "$name,phase\n" * join(replace.(sod, " " => ","), "\n"))
            @test read_data(tbl).tau0 == 0.1
        end
        write(tbl, "t (s),phase\n" * join(string.(50000 .+ (0:9) ./ 2, ",", 1:10), "\n"))
        @test read_data(tbl).tau0 == 0.5
        # A column sampled every third row keeps its own step; a lone column
        # beside the time column is the clock.
        write(tbl, "MJD,a (ns),b (ns)\n" * join(["$(60000 + i / 86400),$(i % 3 == 0 ? i : NaN),$i" for i in 0:29], "\n"))
        r = read_data(tbl; column="a (ns)")
        @test r.tau0 ≈ 3.0 && r.x ≈ (0:3:27) .* 1e-9
        write(tbl, "MJD,UTC-UTC(k) (ns)\n60000,1\n60005,2\n60010,4\n")
        @test read_data(tbl).x ≈ [1e-9, 2e-9, 4e-9]
        @test isempty(SigmaTau._time_grid([0.0, 1.0, 0.5])[1])        # not a time axis
        @test isempty(SigmaTau._time_grid([0.0, 1.0, 2.5, 3.5])[1])   # fits no grid
        # A counter log with ISO timestamps.
        ctr = tempname() * ".txt"
        write(ctr, "2026-11-08T12:00:00Z +1.0000000001E+07\n2026-11-08T12:00:01Z +1.0000000002E+07\n" *
                   "2026-11-08T12:00:03Z +1.0000000001E+07\n")
        y = read_data(ctr)
        @test y isa FrequencyData && y.tau0 == 1.0 && isnan(y.y[3])
        @test start_time(y) == DateTime(2026, 11, 8, 12)
    end

    @testset "ISO 8601 offsets" begin
        utc = SigmaTau._iso_seconds("2026-11-08T12:00:00Z")
        @test SigmaTau._iso_seconds("2026-11-08T14:00:00+02:00") == utc
        @test SigmaTau._iso_seconds("2026-11-08 06:30:00-0530") == utc
        @test SigmaTau._iso_seconds("2026-11-08T12:00:00") == utc
    end
end

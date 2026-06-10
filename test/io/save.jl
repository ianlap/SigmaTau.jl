# test/io/save.jl — save(data) two-column text format, comment-header
# round-trip through read_phase/read_frequency, tau0 pickup, source
# propagation, and the save() delegation to the result/suite writers.

@testset "save / data round-trip" begin

    @testset "PhaseData round-trip with no read keywords" begin
        Random.seed!(20260610)
        N  = 256
        x  = cumsum(randn(N)) .* 1e-10
        pd = PhaseData(x, 0.5)
        mktempdir() do dir
            path = joinpath(dir, "phase.txt")
            @test save(path, pd) == path

            lines = readlines(path)
            @test lines[1] == "# SigmaTau phase data"
            @test lines[2] == "# source: user"
            @test lines[3] == "# tau0: 0.5"
            @test length(lines) == 3 + N
            # Two columns: sample time (i-1)·tau0, then the value.
            c1 = split(lines[4], '\t')
            @test parse(Float64, c1[1]) == 0.0
            @test parse(Float64, c1[2]) == x[1]
            c2 = split(lines[5], '\t')
            @test parse(Float64, c2[1]) == 0.5

            # Round-trip: comments auto-skipped, tau0 picked up from header.
            pd2 = read_phase(path)
            @test pd2 isa PhaseData
            @test pd2.x ≈ x
            @test pd2.tau0 == 0.5
            @test pd2.source == path
        end
    end

    @testset "FrequencyData round-trip, csv delimiter" begin
        Random.seed!(20260610)
        N  = 128
        y  = randn(N) .* 1e-12
        fd = FrequencyData(y, 10.0)
        mktempdir() do dir
            path = joinpath(dir, "freq.csv")
            save(path, fd)
            lines = readlines(path)
            @test lines[1] == "# SigmaTau frequency data"
            @test occursin(',', lines[4])     # csv extension → comma delimiter

            fd2 = read_frequency(path)
            @test fd2 isa FrequencyData
            @test fd2.y ≈ y
            @test fd2.tau0 == 10.0
            @test fd2.source == path
        end
    end

    @testset "explicit tau0 kwarg overrides the header value" begin
        pd = PhaseData([1.0, 2.0, 3.0, 4.0] .* 1e-9, 2.0)
        mktempdir() do dir
            path = joinpath(dir, "p.txt")
            save(path, pd)
            @test read_phase(path).tau0 == 2.0
            @test read_phase(path; tau0=7.0).tau0 == 7.0
        end
    end

    @testset "source written to header and re-read records new path" begin
        mktempdir() do dir
            first_path  = joinpath(dir, "a.txt")
            second_path = joinpath(dir, "b.txt")
            pd  = PhaseData(randn(16), 1.0)
            save(first_path, pd)
            pd1 = read_phase(first_path)
            @test pd1.source == first_path
            # Saving a read record writes its provenance into the header …
            save(second_path, pd1)
            @test any(==("# source: $first_path"), readlines(second_path))
            # … and reading the new file records the new path.
            @test read_phase(second_path).source == second_path
        end
    end

    @testset "source propagation through detrend / fillgaps" begin
        # Constructor default
        pd = PhaseData(collect(1.0:32.0), 1.0)
        @test pd.source == "user"
        @test detrend(pd).source == "user"

        mktempdir() do dir
            path = joinpath(dir, "p.txt")
            save(path, PhaseData(cumsum(randn(64)), 1.0))
            pdr = read_phase(path)
            @test pdr.source == path
            @test detrend(pdr; method=:quadratic).source == path
            xg = copy(pdr.x); xg[20:22] .= NaN
            @test fillgaps(PhaseData(xg, 1.0; source=pdr.source)).source == path

            fdr = read_frequency(path)
            @test fdr.source == path
            @test detrend(fdr).source == path
            @test fillgaps(fdr).source == path
        end
    end

    @testset "leading comments do not break plain files or .DAT headers" begin
        mktempdir() do dir
            # A file with hand-written leading comments but no tau0 entry:
            # comments are skipped, tau0 inferred from the time column.
            path = joinpath(dir, "commented.txt")
            open(path, "w") do io
                println(io, "# my counter log")
                println(io, "# exported 2026-06-10")
                for k in 0:31
                    println(io, Float64(k), '\t', k * 1e-10)
                end
            end
            pd = read_phase(path)
            @test length(pd.x) == 32
            @test pd.tau0 ≈ 1.0

            # A non-leading '#' line is NOT treated as a comment header:
            # the Stable32 .DAT convention ("# Header End" inside a plain-text
            # header) still needs the explicit header kwarg.
            datpath = joinpath(pkgdir(SigmaTau), "test", "fixtures",
                               "validation", "stable32gen.DAT")
            pd_dat = read_phase(datpath; header=10, time_col=0, value_col=1,
                                tau0=1.0)
            @test length(pd_dat.x) == 8192
            @test pd_dat.source == datpath
        end
    end

    @testset "save delegates to result / suite writers" begin
        Random.seed!(20260610)
        pd = PhaseData(randn(128) .* 1e-9, 1.0)
        mktempdir() do dir
            r     = adev(pd, [1, 2, 4]; ci=true)
            rpath = joinpath(dir, "r.tsv")
            @test save(rpath, r) == rpath
            r2 = load_result(rpath)
            @test r2.dev ≈ r.dev
            @test ci_upper(r2) ≈ ci_upper(r)

            s     = stability(pd; devs=(:adev, :hdev))
            spath = joinpath(dir, "s.tsv")
            @test save(spath, s; source_file="p.txt") == spath
            s2 = load_suite(spath)
            @test s2[:hdev].dev ≈ s[:hdev].dev
            @test any(l -> l == "# source_file=p.txt", readlines(spath))
        end
    end
end

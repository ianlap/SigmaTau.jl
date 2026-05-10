@testset "read_phase / read_frequency" begin

    @testset "basic 2-column TSV round-trip" begin
        Random.seed!(20260510)
        N = 256
        t = collect(0.0:N-1)
        x = randn(N)
        mktempdir() do dir
            path = joinpath(dir, "phase.tsv")
            open(path, "w") do io
                for i in 1:N
                    println(io, t[i], '\t', x[i])
                end
            end

            pd = read_phase(path)
            @test length(pd.x) == N
            @test pd.x ≈ x
            @test pd.tau0 ≈ 1.0          # inferred from time column
        end
    end

    @testset "scaling and explicit tau0 (no time column)" begin
        N = 64
        x = randn(N)
        mktempdir() do dir
            path = joinpath(dir, "phase_only.txt")
            open(path, "w") do io
                for i in 1:N
                    println(io, x[i])
                end
            end

            pd = read_phase(path; time_col=0, value_col=1, tau0=0.1, scaling=1e-9)
            @test pd.tau0 == 0.1
            @test pd.x ≈ x .* 1e-9
        end
    end

    @testset "detrend kwarg removes a planted ramp" begin
        N = 1024
        ramp = 1e-3 .* (0:N-1)
        x = randn(N) .+ ramp
        mktempdir() do dir
            path = joinpath(dir, "trended.tsv")
            open(path, "w") do io
                for i in 1:N
                    println(io, i - 1, '\t', x[i])
                end
            end

            pd = read_phase(path; detrend=:linear)
            @test abs(sum(pd.x) / N) < 1e-12              # zero mean
            ns    = collect(0:N-1)
            sx    = sum(ns); sxx = sum(abs2, ns)
            sy    = sum(pd.x); sxy = sum(ns .* pd.x)
            slope = (N * sxy - sx * sy) / (N * sxx - sx * sx)
            @test abs(slope) < 1e-15                      # zero residual slope
        end
    end

    @testset "fillgaps kwarg requires a time column" begin
        mktempdir() do dir
            path = joinpath(dir, "p.txt")
            open(path, "w") do io
                for v in 1:32
                    println(io, Float64(v))
                end
            end
            @test_throws ArgumentError read_phase(path; time_col=0, tau0=1.0, fillgaps=true)
        end
    end

    @testset "invalid value_col errors" begin
        mktempdir() do dir
            path = joinpath(dir, "p.tsv")
            open(path, "w") do io
                println(io, "0\t1.0")
                println(io, "1\t2.0")
            end
            @test_throws ArgumentError read_phase(path; value_col=5)
        end
    end

    @testset "read_frequency returns FrequencyData" begin
        Random.seed!(20260510)
        N = 128
        y = randn(N) .* 1e-9
        mktempdir() do dir
            path = joinpath(dir, "freq.csv")
            open(path, "w") do io
                for i in 1:N
                    println(io, i - 1, ',', y[i])
                end
            end
            fd = read_frequency(path)
            @test fd isa FrequencyData
            @test fd.y ≈ y
            @test fd.tau0 ≈ 1.0
        end
    end
end

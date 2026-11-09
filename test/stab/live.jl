# live.jl — live sources and LiveStability: streamed σ(τ) equals the batch
# estimators on the readings so far, for each reading quantity.

using Test
using Random
using Sockets
using SigmaTau

@testset "Live measurement" begin
    tau0 = 0.5
    ms = [1, 2, 4, 16]
    devs = [:adev, :mdev, :hdev, :mhdev]
    batch = Dict(:adev => adev, :mdev => mdev, :hdev => hdev, :mhdev => mhdev)
    same(r, ref) = isequal(isnan.(r.dev), isnan.(ref.dev)) &&
                   all(isapprox(p, q; rtol=1e-10) for (p, q) in zip(r.dev, ref.dev) if !isnan(q))

    @testset "simulated phase matches batch" begin
        src = SimulatedSource(white_pm=1e-9, white_fm=1e-11, rw_fm=1e-13, tau0=tau0,
                              rng=Xoshiro(1))
        live = LiveStability(devs, tau0, ms)
        x = [read_reading(src) for _ in 1:300]
        foreach(v -> push!(live, v), x)
        @test nsamples(live) == 300
        for (d, r) in zip(devs, snapshot(live))
            @test same(r, batch[d](PhaseData(x, tau0), ms; ci=false))
        end
    end

    @testset "frequency and hertz readings integrate like FrequencyData" begin
        y = 1e-11 .* randn(Xoshiro(2), 200)
        live = LiveStability([:adev], tau0, ms; quantity=:frequency)
        foreach(v -> push!(live, v), y)
        ref = adev(FrequencyData(y, tau0), ms; ci=false)
        @test same(only(snapshot(live)), ref)

        hz = LiveStability([:adev], tau0, ms; quantity=:hz)
        foreach(v -> push!(hz, 10e6 * (1 + v)), y)
        @test hz.nominal == 10e6
        @test only(snapshot(hz)).dev ≈ ref.dev rtol = 1e-5
    end

    @testset "white FM level" begin
        src = SimulatedSource(white_fm=2e-12, rng=Xoshiro(3))
        live = LiveStability([:adev], 1.0, [1])
        for _ in 1:20_000
            push!(live, read_reading(src))
        end
        @test only(snapshot(live)).dev[1] ≈ 2e-12 rtol = 0.03
    end

    @testset "instrument lines" begin
        io = IOBuffer("+1.000000000012E+07\n\n10000000.5,VOLT\r\n  9.9999999E6 Hz\n")
        src = InstrumentSource(io)
        @test [read_reading(src) for _ in 1:3] == [1.000000000012e7, 10000000.5, 9.9999999e6]
        @test isnan(read_reading(src))
        @test_throws ArgumentError read_reading(InstrumentSource(IOBuffer("OVLD\n")))

        # An LXI-style SCPI socket: the query goes out before every reading.
        server = listen(Sockets.localhost, 0)
        port = Int(getsockname(server)[2])
        replies = @async begin
            c = accept(server)
            n = 0
            while readline(c) == "READ?"
                n += 1
                println(c, n / 4)
            end
            close(c)
            n
        end
        sock = Sockets.connect(port)
        scpi = InstrumentSource(sock; query="READ?")
        @test [read_reading(scpi) for _ in 1:3] == [0.25, 0.5, 0.75]
        close(sock); close(server)
        @test fetch(replies) == 3
    end

    @test_throws ArgumentError LiveStability([:adev], 1.0, [1]; quantity=:volts)
    @test_throws ArgumentError push!(LiveStability([:adev], 1.0, [1]), NaN)
end

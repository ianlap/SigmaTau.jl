# streaming.jl — streaming accumulator equivalence tests.
#
# The correctness contract of `StreamingStability` is exact batch parity:
# after pushing x[1:k], `snapshot(acc).dev` must reproduce the batch
# `adev/mdev/hdev/mhdev(PhaseData(x[1:k], tau0), m_values; ci=false)` at
# every supported m (NaN where the batch gives NaN, identical neff). The
# scheme is Dobrogowski & Kasznia 2007 (IEEE FCS, pp. 877–882): running sums
# of squared second differences (eqs. 6–9) and the inner-sum/overall-sum
# form for the modified family (eqs. 10–14), generalized here to the
# third-difference family for :hdev/:mhdev.

using Test
using Random
using SigmaTau

@testset "Streaming accumulators (Dobrogowski–Kasznia)" begin
    # Fixed-seed mixed power-law fixture (WPM + RWFM), same recipe as the
    # legacy-parity testset. tau0 ≠ 1 so normalization errors can't hide.
    Random.seed!(20260610)
    N    = 600
    tau0 = 0.5
    x    = randn(N) .* 1e-9 .+ cumsum(cumsum(randn(N) .* 1e-12))

    # Non-octave m's included; checkpoints straddle each kernel's activation
    # boundary (e.g. mhdev at m=33 needs t = 4m+1 = 133 for its 2nd window;
    # at m=64 it needs t = 257).
    ms          = [1, 2, 5, 8, 33, 64]
    checkpoints = (2, 3, 7, 50, 132, 133, 257, 400, N)
    batch_fns   = Dict(:adev => adev, :mdev => mdev, :hdev => hdev, :mhdev => mhdev)

    @testset "sample-by-sample batch equivalence: $dev" for dev in
            (:adev, :mdev, :hdev, :mhdev)
        acc = StreamingStability(dev, tau0, ms)
        for k in 1:N
            push!(acc, x[k])
            k in checkpoints || continue
            @test nsamples(acc) == k
            r  = snapshot(acc)
            rb = batch_fns[dev](PhaseData(x[1:k], tau0), ms; ci=false)
            @test r.deviation_type == dev
            @test r.tau == rb.tau
            @test r.neff == rb.neff
            for j in eachindex(ms)
                if isnan(rb.dev[j])
                    @test isnan(r.dev[j])
                else
                    @test isapprox(r.dev[j], rb.dev[j]; rtol=1e-10)
                end
            end
            # snapshot() carries no CI machinery (on-demand EDF/CI is future
            # work) — those fields are empty, exactly like ci=false batch.
            @test isempty(r.ci) && isempty(r.edf) && isempty(r.noise_type)
        end
    end

    @testset "append! chunked feeding ≡ sample-by-sample: $dev" for dev in
            (:adev, :mdev, :hdev, :mhdev)
        acc_chunked = StreamingStability(dev, tau0, ms)
        # Irregular chunk sizes, including a 1-sample chunk and a big tail.
        i = 1
        for len in (1, 7, 64, 100, 13, 250, N)   # last chunk clipped to N
            hi = min(i + len - 1, N)
            hi < i && break
            append!(acc_chunked, x[i:hi])
            i = hi + 1
        end
        @test nsamples(acc_chunked) == N

        acc_single = StreamingStability(dev, tau0, ms)
        for v in x
            push!(acc_single, v)
        end

        r_chunked = snapshot(acc_chunked)
        r_single  = snapshot(acc_single)
        # Identical push order → identical arithmetic → bitwise equality.
        @test isequal(r_chunked.dev, r_single.dev)
        @test r_chunked.neff == r_single.neff

        # And both match the batch kernel on the full record.
        rb = batch_fns[dev](PhaseData(x, tau0), ms; ci=false)
        for j in eachindex(ms)
            if isnan(rb.dev[j])
                @test isnan(r_chunked.dev[j])
            else
                @test isapprox(r_chunked.dev[j], rb.dev[j]; rtol=1e-10)
            end
        end
    end

    @testset "do-nothing accumulator: defined empty/NaN behavior" begin
        acc = StreamingStability(:adev, 1.0, [1, 2, 4])
        @test nsamples(acc) == 0
        r = snapshot(acc)
        @test r isa StabilityResult
        @test r.deviation_type == :adev
        @test r.tau == [1.0, 2.0, 4.0]
        @test all(isnan, r.dev)
        @test r.neff == [0, 0, 0]
        @test isempty(r.ci) && isempty(r.edf) && isempty(r.noise_type)

        # One sample in: still below every activation threshold.
        push!(acc, 1e-9)
        @test nsamples(acc) == 1
        @test all(isnan, snapshot(acc).dev)
        @test snapshot(acc).neff == [0, 0, 0]
    end

    @testset "memory bound: ring buffer sized to the kernel's deepest lag" begin
        for (dev, lag) in ((:adev, 2), (:mdev, 3), (:hdev, 3), (:mhdev, 4))
            acc = StreamingStability(dev, 1.0, [1, 5, 64])
            @test length(acc.buf) == lag * 64 + 1
        end
    end

    @testset "constructor validation" begin
        @test_throws ArgumentError StreamingStability(:totdev, 1.0, [1, 2])  # can't stream
        @test_throws ArgumentError StreamingStability(:theo1, 1.0, [2, 4])   # can't stream
        @test_throws ArgumentError StreamingStability(:adev, 0.0, [1, 2])    # tau0 ≤ 0
        @test_throws ArgumentError StreamingStability(:adev, -1.0, [1, 2])
        @test_throws ArgumentError StreamingStability(:adev, 1.0, Int[])     # empty m grid
        @test_throws ArgumentError StreamingStability(:adev, 1.0, [0, 2])    # m < 1
    end
end

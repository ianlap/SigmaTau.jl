# Dual-channel estimators: identities with the single-channel ones, recovery of
# a device below two noisier references, and the floor for channels that share
# nothing.

@testset "cross estimators" begin
    a = noise_gen(PhaseData, 4096, 1.0; h=Dict(0 => 1e-24, 2 => 1e-26), rng=Xoshiro(1))

    @testset "gaps are rejected" begin
        g = PhaseData([a.x[1:99]; NaN; a.x[101:end]], 1.0)
        @test_throws ArgumentError cross_adev(g, a)
        @test_throws ArgumentError cross_psd_x(a, g)
    end

    @testset "a × a is the single-channel estimate" begin
        m = [1, 2, 4, 16, 64]
        @test cross_adev(a, a, m).covar ≈ adev(a, m; ci=false).dev .^ 2 rtol=1e-12
        c = cross_psd_y(a, a)
        @test c.psd ≈ psd_y(a).psd rtol=1e-12
        @test c.floor ≈ c.psd ./ sqrt(SigmaTau._independent_averages(:hann, c.nperseg, c.nperseg - c.noverlap, c.nseg)) rtol=1e-12
        # Overlap adds segments, not independent averages: at 90 % the floor
        # stays near the one without overlap.
        k0 = cross_psd_y(a, a; nperseg=256, noverlap=0)
        k9 = cross_psd_y(a, a; nperseg=256, noverlap=230)
        @test k9.nseg > 8k0.nseg
        @test SigmaTau._independent_averages(:hann, 256, 26, k9.nseg) < 2.5k0.nseg
        @test SigmaTau._independent_averages(:hann, 256, 256, 16) == 16
        @test cross_psd_x(a, a).psd ≈ psd_x(a).psd rtol=1e-12
    end

    @testset "device recovered below two references" begin
        # Two channels: the same device (h0 = 1e-24) against independent
        # references ten times noisier (h0 = 1e-23 each).
        N = 1 << 16
        dut = noise_gen(FrequencyData, N, 1.0; h=Dict(0 => 1e-24), rng=Xoshiro(2)).y
        r1  = noise_gen(FrequencyData, N, 1.0; h=Dict(0 => 1e-23), rng=Xoshiro(3)).y
        r2  = noise_gen(FrequencyData, N, 1.0; h=Dict(0 => 1e-23), rng=Xoshiro(4)).y
        A, B = FrequencyData(dut .- r1, 1.0), FrequencyData(dut .- r2, 1.0)
        c = cross_psd_y(A, B; nperseg=512)
        band = 20:200
        @test isapprox(sum(c.psd[band]) / length(band), 1e-24; rtol=0.2)
        @test sum(c.psd_a[band]) / length(band) > 9e-24          # one channel sees the reference
        x = cross_adev(A, B, [1, 4, 16])
        @test all(abs.(x.covar .- adev(FrequencyData(dut, 1.0), [1, 4, 16]; ci=false).dev .^ 2) .< 3x.uncertainty)
        @test all(x.dev_a .> 3x.dev)
    end

    @testset "uncertainty for channels that share nothing" begin
        # The reported 1σ against the actual scatter over 60 independent pairs,
        # for white PM, white FM and random-walk FM.
        m = [1, 8, 64]
        for α in (2, 0, -2)
            gen(k) = noise_gen(PhaseData, 4096, 1.0; h=Dict(α => 1e-24), rng=Xoshiro(1000α + k))
            c = [cross_adev(gen(k), gen(500 + k), m) for k in 1:60]
            spread = [sqrt(sum(r.covar[j]^2 for r in c) / length(c)) for j in eachindex(m)]
            said = [sum(r.uncertainty[j] for r in c) / length(c) for j in eachindex(m)]
            @test all(0.7 .< spread ./ said .< 1.4)
        end
    end

    @testset "argument checks" begin
        @test_throws ArgumentError cross_adev(a, PhaseData(a.x, 2.0))
        @test_throws ArgumentError cross_adev(a, PhaseData(a.x[1:100], 1.0))
        @test_throws ArgumentError cross_psd_y(a, FrequencyData(a.x[2:end], 1.0))
        @test_throws ArgumentError cross_adev(a, a, [3000])
    end
end

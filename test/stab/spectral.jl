# test/stab/spectral.jl — Sy / Sx / L spectral-density coverage.
#
# Validates the Welch core's variance normalization, the S_x ↔ S_y conversion
# identity, ℒ(f) carrier scaling, and power-law slope recovery via noise_gen.

using Test
using Random
using Statistics
using FFTW                                  # AbstractFFTs backend
using SigmaTau
using SigmaTau: _welch_psd, _window, _phase_to_freq, _freq_to_phase

# Least-squares slope of log(psd) vs log(freq) over the given bin range.
function _loglog_slope(freq, psd, idx)
    lf = log.(freq[idx]); ls = log.(psd[idx])
    n  = length(lf)
    (n * sum(lf .* ls) - sum(lf) * sum(ls)) / (n * sum(lf .^ 2) - sum(lf)^2)
end

@testset "Spectral (Sy / Sx / L)" begin

    @testset "Welch window helper" begin
        @test _window(:rectangular, 8) == ones(8)
        @test _window(:hann, 4)[1] == 0.0           # periodic Hann starts at 0
        @test length(_window(:hamming, 16)) == 16
        @test_throws ArgumentError _window(:blackman, 8)
    end

    @testset "variance normalization (Parseval)" begin
        Random.seed!(1)
        z = randn(8192)
        for win in (:hann, :hamming, :rectangular), np in (256, 512)
            freq, psd = _welch_psd(z, 1.0; nperseg=np, noverlap=np ÷ 2, window=win)
            @test freq[1] == 0.0
            @test freq[end] ≈ 0.5                    # Nyquist at fs/2
            @test length(freq) == np ÷ 2 + 1
            df = freq[2] - freq[1]
            @test isapprox(sum(psd) * df, var(z); rtol=0.05)
        end
    end

    @testset "Welch argument validation" begin
        z = randn(64)
        @test_throws ArgumentError _welch_psd(z, 1.0; nperseg=1, noverlap=0, window=:hann)
        @test_throws ArgumentError _welch_psd(z, 1.0; nperseg=128, noverlap=0, window=:hann)
        @test_throws ArgumentError _welch_psd(z, 1.0; nperseg=32, noverlap=32, window=:hann)
        @test_throws ArgumentError _welch_psd(z, 0.0; nperseg=16, noverlap=8, window=:hann)
    end

    @testset "power-law slope recovery (S_y ∝ f^α)" begin
        Random.seed!(7)
        for α in (2, 0, -2)
            p = noise_gen(PhaseData, 16384, 1.0; h = Dict(α => 1.0))
            r = Sy(p; nperseg=2048)
            @test r.spectral_type === :Sy
            @test r.units === :per_Hz
            idx = filter(k -> r.freq[k] > 0 && r.psd[k] > 0, eachindex(r.freq))[2:end-3]
            @test isapprox(_loglog_slope(r.freq, r.psd, idx), α; atol=0.3)
        end
    end

    @testset "S_x ↔ S_y identity at low frequency" begin
        Random.seed!(11)
        p  = noise_gen(PhaseData, 16384, 1.0; sigma1 = Dict(0 => 1e-12))
        ry = Sy(p; nperseg=2048)
        rx = Sx(p; nperseg=2048)
        @test rx.units === :s2_per_Hz
        @test ry.freq == rx.freq
        lo   = 2:50
        pred = ry.psd[lo] ./ (2π .* ry.freq[lo]) .^ 2     # S_y/(2πf)²
        @test isapprox(median(rx.psd[lo] ./ pred), 1.0; rtol=0.05)
    end

    @testset "ℒ(f) carrier scaling" begin
        Random.seed!(13)
        p  = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(2 => 1e-11))
        fc = 1e7
        rx = Sx(p; nperseg=1024)
        rl = L(p; f_carrier=fc, nperseg=1024)
        @test rl.spectral_type === :L
        @test rl.units === :dBc_per_Hz
        @test length(rl.freq) == length(rx.freq) - 1        # DC dropped
        @test rl.freq == rx.freq[2:end]
        expected = 10 .* log10.(0.5 .* (2π * fc)^2 .* rx.psd[2:end])
        @test rl.psd ≈ expected
        @test_throws ArgumentError L(p; f_carrier=0.0)
    end

    @testset "PhaseData / FrequencyData entry points agree" begin
        Random.seed!(17)
        y  = randn(8192) .* 1e-10
        fd = FrequencyData(y, 1.0)
        pd = _freq_to_phase(fd)
        # Sy on the frequency record vs the integrated-then-differenced phase
        # record differ only by the dropped first sample; spectra should match.
        a = Sy(fd; nperseg=1024)
        b = Sy(pd; nperseg=1024)
        @test a.units === b.units === :per_Hz
        @test isapprox(median(a.psd[2:end] ./ b.psd[2:end]), 1.0; rtol=0.1)
        # Sx via either input is the phase spectrum.
        @test Sx(pd; nperseg=1024).units === :s2_per_Hz
        @test Sx(fd; nperseg=1024).spectral_type === :Sx
    end

    @testset "default Welch parameters" begin
        Random.seed!(19)
        p = noise_gen(PhaseData, 4096, 1.0; sigma1 = Dict(0 => 1e-12))
        r = Sy(p)
        @test r.nperseg == 256                  # min(N, 256)
        @test r.noverlap == 128                 # 50 % overlap
        @test r.window === :hann
    end
end

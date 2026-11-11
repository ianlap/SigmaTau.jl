# test/stab/spectral.jl — psd_y / psd_x / ssb_phase_noise coverage.
#
# Validates the Welch core's variance normalization, the S_x ↔ S_y conversion
# identity, ℒ(f) carrier scaling, and power-law slope recovery via noise_gen.

using Test
using Random
using Statistics
using SigmaTau
using SigmaTau: _welch_psd, _window, _freq_to_phase

# Least-squares slope of log(psd) vs log(freq) over the given bin range.
function _loglog_slope(freq, psd, idx)
    lf = log.(freq[idx]); ls = log.(psd[idx])
    n  = length(lf)
    (n * sum(lf .* ls) - sum(lf) * sum(ls)) / (n * sum(lf .^ 2) - sum(lf)^2)
end

@testset "Spectral (psd_y / psd_x / ssb_phase_noise)" begin

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
            r = psd_y(p; nperseg=2048)
            @test r.spectral_type === :Sy
            @test r.units === :per_Hz
            idx = filter(k -> r.freq[k] > 0 && r.psd[k] > 0, eachindex(r.freq))[2:end-3]
            @test isapprox(_loglog_slope(r.freq, r.psd, idx), α; atol=0.3)
        end
    end

    @testset "S_x is S_y over the first-difference response" begin
        p  = noise_gen(PhaseData, 16384, 1.0; h = Dict(-2 => 1.0), rng = Xoshiro(11))
        ry = psd_y(p; nperseg=2048)
        rx = psd_x(p; nperseg=2048)
        @test rx.units === :s2_per_Hz && ry.freq == rx.freq && isnan(rx.psd[1])
        @test rx.psd[2:end] ≈ ry.psd[2:end] ./ (4 .* sinpi.(ry.freq[2:end]) .^ 2)
        # Random-walk FM, S_x = f⁻⁴/(2π)²: the low bins read true (a Welch
        # estimate on the phase itself is 15 to 100 times high at bins 2 to 4).
        lo = 3:10
        @test isapprox(median(rx.psd[lo] .* (2π) .^ 2 .* rx.freq[lo] .^ 4), 1.0; rtol=0.25)
        # A frequency offset is the mean of y and leaves the spectrum alone.
        off = PhaseData(p.x .+ 50.0 .* (0:16383), 1.0)
        @test psd_x(off; nperseg=2048).psd[2:end] ≈ rx.psd[2:end] rtol=1e-6
        # White PM, S_x = 2τ₀σ²: the first bin stays true on a long segment
        # (with the plain segment mean it reads 1 + nperseg/(12π²) = 70 times high).
        w = PhaseData(randn(Xoshiro(3), 2^17), 1.0)
        @test 0.5 < psd_x(w; nperseg=8192).psd[2] / 2 < 2
        sy = psd_y(w; nperseg=8192)
        @test sy.psd[1] < 1e-20 * sy.psd[2]
    end

    @testset "ℒ(f) carrier scaling" begin
        Random.seed!(13)
        p  = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(2 => 1e-11))
        fc = 1e7
        rx = psd_x(p; nperseg=1024)
        rl = ssb_phase_noise(p; f_carrier=fc, nperseg=1024)
        @test rl.spectral_type === :L
        @test rl.units === :dBc_per_Hz
        @test length(rl.freq) == length(rx.freq) - 1        # DC dropped
        @test rl.freq == rx.freq[2:end]
        expected = 10 .* log10.(0.5 .* (2π * fc)^2 .* rx.psd[2:end])
        @test rl.psd ≈ expected
        @test_throws ArgumentError ssb_phase_noise(p; f_carrier=0.0)
    end

    @testset "PhaseData / FrequencyData entry points agree" begin
        Random.seed!(17)
        y  = randn(8192) .* 1e-10
        fd = FrequencyData(y, 1.0)
        pd = _freq_to_phase(fd)
        # psd_y on the frequency record vs the integrated-then-differenced phase
        # record differ only by the dropped first sample; spectra should match.
        a = psd_y(fd; nperseg=1024)
        b = psd_y(pd; nperseg=1024)
        @test a.units === b.units === :per_Hz
        @test isapprox(median(a.psd[2:end] ./ b.psd[2:end]), 1.0; rtol=0.1)
        # psd_x via either input is the phase spectrum.
        @test psd_x(pd; nperseg=1024).units === :s2_per_Hz
        @test psd_x(fd; nperseg=1024).spectral_type === :Sx
    end

    @testset "default Welch parameters" begin
        Random.seed!(19)
        p = noise_gen(PhaseData, 4096, 1.0; sigma1 = Dict(0 => 1e-12))
        r = psd_y(p)
        # psd_y differences phase to N − 1 samples: 4095 ÷ 4 = 1023 → 512.
        @test r.nperseg == 512                  # largest 2^k ≤ N ÷ 4
        @test r.noverlap == 256                 # 50 % overlap
        @test r.taper === :hann
        @test psd_x(p).nperseg == 1024          # from the 4096 phase samples
        @test psd_x(PhaseData(p.x[1:10], 1.0)).nperseg == 9
        @test SigmaTau._default_nperseg(4095) == 512
        @test SigmaTau._default_nperseg(40) == 16
        @test SigmaTau._default_nperseg(10) == 10
    end

    @testset "integrated_jitter (SP1065 §5.2.19)" begin
        # ℒ falling 10 dB a decade: S_φ = 2e-10·(f/1e3)^-1, so Δφ² = 2e-7·ln 100
        # and Δf² = 2e-7·(1e5² − 1e3²)/2, both exact on straight-line segments.
        f = [1e3, 1e4, 1e5]
        j = integrated_jitter(f, [-100.0, -110.0, -120.0]; f_carrier=10e6)
        @test j.phase ≈ √(2e-7 * log(100))
        @test j.time ≈ j.phase / (2π * 10e6)
        @test j.residual_fm ≈ √(2e-7 * (1e10 - 1e6) / 2)
        @test integrated_jitter(f, [-100.0, -110.0, -120.0]; f_carrier=10e6, band=(1e3, 1e4)).phase ≈ √(2e-7 * log(10))
        # White PM: the jitter over the whole band is the phase rms (Parseval),
        # and over a tenth of it √(1/10) of that, from any of the three spectra.
        Random.seed!(23)
        p = PhaseData(1e-9 .* randn(2^16), 1e-3)
        for s in (psd_x(p), psd_y(p), ssb_phase_noise(p; f_carrier=1e6))
            @test isapprox(integrated_jitter(s; f_carrier=1e6).time, 1e-9; rtol=0.02)
        end
        @test isapprox(integrated_jitter(psd_x(p); f_carrier=1e6, band=(400.0, 450.0)).time, √0.1 * 1e-9; rtol=0.05)
        @test_throws ArgumentError integrated_jitter(f, [-100.0, -110.0, -120.0]; f_carrier=10e6, band=(10.0, 1e4))
        @test_throws ArgumentError integrated_jitter(reverse(f), [-100.0, -110.0, -120.0]; f_carrier=10e6)
        @test integrated_jitter(f, [2e-10, 2e-11, 2e-12]; quantity=:Sphi, f_carrier=10e6) == j
    end

    @testset "spectrum_sigma against the closed forms" begin
        # A two-point table is one power law, S_y = f^α, and the band with its
        # two continuations covers every frequency: σ² = K·τ^(−α−1).
        f = [1e-3, 2e3]
        function K(dev, α; τ=1.0)
            r = spectrum_sigma(f, f .^ α, [τ]; quantity=:Sy, deviation=dev)
            return r.dev[1]^2 * (1 + r.below[1] + r.above[1]) * τ^(α + 1)
        end
        l2, l3 = log(2), log(3)
        # ADEV: SP1065 Eq. (67); IEEE Std 1139-1999 Table B.2
        @test K(:adev, 0)   ≈ 1 / 2                          rtol=1e-10
        @test K(:adev, -1)  ≈ 2l2                            rtol=1e-10
        @test K(:adev, -2)  ≈ 2π^2 / 3                       rtol=1e-10
        @test K(:adev, -2; τ=37.0) ≈ 2π^2 / 3                rtol=1e-10    # past the resolved lobes
        # MDEV: the large-n limits of NIST TN 1337 §A.6, Table 2
        @test K(:mdev, 2)   ≈ 3 / (8π^2)                     rtol=1e-10
        @test K(:mdev, 1)   ≈ (24l2 - 9l3) / (8π^2)          rtol=1e-10
        @test K(:mdev, 0)   ≈ 1 / 4                          rtol=1e-10
        @test K(:mdev, -1)  ≈ (27l3 - 32l2) / 8              rtol=1e-10
        @test K(:mdev, -2)  ≈ 11π^2 / 20                     rtol=1e-10
        # HDEV and MHDEV: the same integrals of their transfer functions
        @test K(:hdev, 0)   ≈ 1 / 2                          rtol=1e-10
        @test K(:hdev, -1)  ≈ (8l2 - 3l3) / 2                rtol=1e-10
        @test K(:hdev, -2)  ≈ π^2 / 3                        rtol=1e-10
        @test K(:hdev, -3)  ≈ π^2 * (27l3 - 32l2) / 6        rtol=1e-10
        @test K(:hdev, -4)  ≈ 11π^4 / 15                     rtol=1e-10
        @test K(:mhdev, 2)  ≈ 5 / (12π^2)                    rtol=1e-10
        @test K(:mhdev, 1)  ≈ (6l2 - 3l3) / π^2              rtol=1e-10
        @test K(:mhdev, 0)  ≈ 2 / 9                          rtol=1e-10
        @test K(:mhdev, -1) ≈ (27l3 - 40l2) / 3              rtol=1e-10
        @test K(:mhdev, -2) ≈ 2π^2 / 9                       rtol=1e-10
        @test K(:mhdev, -3) ≈ π^2 * (832l2 - 486l3) / 45     rtol=1e-10
        @test K(:mhdev, -4) ≈ 604π^4 / 945                   rtol=1e-10

        # PM noise under ADEV depends on the bandwidth f_h = 2 kHz, and its
        # continuation above the band does not converge. The 1.038 of the
        # standard's flicker PM term is 3γ − ln 2.
        wh = spectrum_sigma(f, f .^ 2, [1.0]; quantity=:Sy)
        fl = spectrum_sigma(f, f, [1.0]; quantity=:Sy)
        @test wh.dev[1]^2 ≈ 3 * 2e3 / (4π^2)                                        rtol=1e-10
        @test fl.dev[1]^2 ≈ (3MathConstants.γ - l2 + 3log(2π * 2e3)) / (4π^2)       rtol=1e-8
        @test wh.above == [Inf] && fl.above == [Inf]
        # Flicker-walk FM: the Allan variance does not converge, the Hadamard does.
        @test spectrum_sigma(f, f .^ -3, [1.0]; quantity=:Sy).below == [Inf]
        @test isfinite(spectrum_sigma(f, f .^ -3, [1.0]; quantity=:Sy, deviation=:hdev).below[1])
    end

    @testset "spectrum_sigma: quantities and band" begin
        # ℒ(f) of a 10 MHz carrier with white PM, white FM and flicker FM
        ν₀, τ = 10e6, [1e-3, 1e-2, 0.1]
        f = exp10.(0:0.25:4)
        Sy = powerlaw_psd(Dict(2 => 1e-28, 0 => 1e-24, -1 => 1e-22), f)
        Sx = Sy ./ (2π .* f) .^ 2
        Sφ = (2π * ν₀)^2 .* Sx
        r = spectrum_sigma(f, 10 .* log10.(Sφ ./ 2), τ; f_carrier=ν₀)
        @test r.dev ≈ spectrum_sigma(f, Sφ, τ; quantity=:Sphi, f_carrier=ν₀).dev
        @test r.dev ≈ spectrum_sigma(f, Sx, τ; quantity=:Sx).dev
        @test r.dev ≈ spectrum_sigma(f, Sy, τ; quantity=:Sy).dev
        # The table follows the model to 4 % between its points, and ADEV at
        # 10 ms is the model's with f_h = 10 kHz.
        @test r.dev[2]^2 ≈ 3e4 * 1e-28 / (4π^2 * 1e-4) + 1e-24 / 2e-2 + 2log(2) * 1e-22 rtol=0.02
        # Two bands that share an edge add in variance.
        v(band) = spectrum_sigma(f, Sy, τ; quantity=:Sy, deviation=:mdev, band=band).dev .^ 2
        @test v((1.0, 300.0)) .+ v((300.0, 1e4)) ≈ v(nothing)
        @test r.band == (1.0, 1e4)
        @test spectrum_sigma(f, Sy, τ; quantity=:Sy, band=(1, 300)).band == (1.0, 300.0)
        # The continuation below the band follows the least-squares line of
        # the lowest decade, here 1/f under the scatter, not the first two points.
        g = exp2.(0:10)
        clean = spectrum_sigma(g, 1 ./ g, [0.1]; quantity=:Sy)
        noisy = spectrum_sigma(g, [2; 0.5; 0.5; 2; ones(7)] ./ g, [0.1]; quantity=:Sy)
        @test noisy.below .* noisy.dev .^ 2 ≈ clean.below .* clean.dev .^ 2
        @test_throws ArgumentError spectrum_sigma(f, Sy, τ; quantity=:Sy, deviation=:totdev)
        @test_throws ArgumentError spectrum_sigma(f, Sy, τ; quantity=:Sy, band=(0.1, 10.0))
        @test_throws ArgumentError spectrum_sigma(f, Sy, [0.0]; quantity=:Sy)
        @test_throws ArgumentError spectrum_sigma(f, Sy, τ)              # ℒ(f) without a carrier
        @test_throws DimensionMismatch spectrum_sigma(f, Sy[2:end], τ; quantity=:Sy)
    end
end

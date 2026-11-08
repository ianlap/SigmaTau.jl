# noise_gen against theory: the Welch spectrum of the generated record (of
# its d-th difference for α ≤ −1, which makes it stationary) matches
# h_α (sin(πfτ₀)/(πτ₀))^α (2 sin(πfτ₀))^(2d) across the band, at every α.

@testset "noise_gen spectrum matches theory" begin
    τ₀, N, K, nperseg = 0.5, 4096, 40, 256
    band = 8:120                                  # bins well away from DC and Nyquist
    for α in 2:-1:-4
        d = α <= -1 ? cld(-α, 2) : 0
        hα = 1e-24
        rng = Xoshiro(7000 + α)
        S = zeros(nperseg ÷ 2 + 1)
        for _ in 1:K
            z = noise_gen(FrequencyData, N, τ₀; h=Dict(α => hα), rng=rng).y
            for _ in 1:d
                z = diff(z)
            end
            S .+= SigmaTau._welch_psd(z, 1 / τ₀; nperseg=nperseg, noverlap=nperseg ÷ 2,
                                      window=:hann)[2] ./ K
        end
        f = (0:nperseg ÷ 2) ./ (nperseg * τ₀)
        theory = @. hα * (sinpi(f * τ₀) / (π * τ₀))^α * (2 * sinpi(f * τ₀))^(2d)
        ratio = S[band] ./ theory[band]
        @test isapprox(sum(ratio) / length(ratio), 1; atol=0.03)
        @test all(r -> 0.8 < r < 1.25, ratio)
    end
end

@testset "noise_gen levels and keys" begin
    p = noise_gen(PhaseData, 4096, 2.0; sigma1=Dict(:RRFM => 3e-13), rng=Xoshiro(1))
    @test adev(p, [1]; ci=false).dev[1] ≈ 3e-13 rtol=1e-12
    a = noise_gen(FrequencyData, 512, 1.0; h=Dict(:FLFM => 1e-26), rng=Xoshiro(2))
    b = noise_gen(FrequencyData, 512, 1.0; h=Dict(-1 => 1e-26), rng=Xoshiro(2))
    @test a.y == b.y
    @test_throws ArgumentError noise_gen(PhaseData, 64, 1.0; h=Dict(-5 => 1.0))
    @test_throws ArgumentError noise_gen(PhaseData, 64, 1.0; h=Dict(:PINK => 1.0))
end

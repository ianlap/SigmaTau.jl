# test/stab/stated_noise.jl — the `noise` keyword: a stated noise type in place
# of the identified one.

using Test
using Random
using SigmaTau

@testset "Stated noise type" begin
    p = noise_gen(PhaseData, 2048, 1.0; h=Dict(0 => 1.0), rng=Xoshiro(11))
    ms = [1, 4, 16, 64]

    # The identified labels, stated, give the identified result.
    auto = htotdev(p, ms)
    stated = htotdev(p, ms; noise=auto.noise_type)
    @test (stated.dev, stated.edf, stated.ci) == (auto.dev, auto.edf, auto.ci)

    # One type holds at every τ and sets the EDF and the bias factor.
    r = totdev(p, ms; noise=:RWFM)
    @test r.noise_type == fill(:RWFM, 4)
    raw = totdev(p, ms; correct_bias=false).dev
    @test r.dev == raw ./ sqrt.(bias_correction(r.noise_type, :totvar, ms, 2048))
    @test r.edf == calculate_edf(:totdev, raw, r.noise_type, ms, r.tau, 2048, 2047.0)
    @test adev(p, ms; noise=:FLFM).edf != adev(p, ms).edf

    # A type redder than the kernel covers counts as its reddest.
    @test mdev(p, ms; noise=:RRFM).noise_type == fill(:RWFM, 4)
    @test mhdev(p, ms; noise=:RRFM).noise_type == fill(:RRFM, 4)

    # A model gives the term with the most variance at each τ: white PM and
    # white FM cross in the Allan variance at τ = 3 h₂ / (4π² τ₀ h₀) ≈ 7.6 s.
    model = PowerLawFit([2, 0], [1e-22, 1e-24], :adev, 1.0, Float64[], Float64[], Float64[])
    @test adev(p, ms; noise=model).noise_type == [:WHPM, :WHPM, :WHFM, :WHFM]

    # A record with gaps, a frequency record and the name-keyed calls take it too.
    x = copy(p.x); x[100:120] .= NaN
    @test adev(PhaseData(x, 1.0), ms; noise=:FLFM).noise_type == fill(:FLFM, 4)
    @test hdev(FrequencyData(diff(p.x), 1.0), ms; noise=:WHFM).noise_type == fill(:WHFM, 4)
    @test deviation(:tdev, p, ms; noise=:FLPM).noise_type == fill(:FLPM, 4)
    @test stability(p; devs=[:adev, :mtie], m=ms, noise=:FLFM)[:adev].edf == adev(p, ms; noise=:FLFM).edf

    @test_throws ArgumentError adev(p, ms; noise=:PINK)
    @test_throws ArgumentError adev(p, ms; noise=[:WHFM])
end

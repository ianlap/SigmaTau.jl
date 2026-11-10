# test/stab/noise_id.jl — lag-1 ACF / B1 noise identification and the
# B1(N, μ), R(n) closed forms.

using Test
using Random
using SigmaTau

@testset "Noise identification" begin
    # Mock data equivalent to NIST SP1065 sample
    x = collect(1.0:100.0) .^ 2  # simple quadratic
    m_values = [1, 2, 4]
    tau0 = 1.0

    @testset "flicker PM at long τ" begin
        # Decimated flicker-PM phase has too little lag-1 correlation left to
        # tell it from white PM; R(n) does, and leaves white PM alone.
        pm(alpha) = identify_noise(noise_gen(PhaseData, 4096, 1.0; h=Dict(alpha => 1.0),
                                             rng=Xoshiro(11)), [16, 64])
        @test pm(1) == [:FLPM, :FLPM]
        @test pm(2) == [:WHPM, :WHPM]
    end

    @testset "flicker FM at long τ" begin
        # On decimated phase flicker FM reads α ≈ −1.45 and random-walk FM
        # −2.4, so for m ≥ 2 the limit between them sits at −1.92. This
        # flicker FM record reads −1.73 at m = 64.
        fm(alpha) = identify_noise(noise_gen(PhaseData, 4096, 1.0; h=Dict(alpha => 1.0),
                                             rng=Xoshiro(11)), [16, 64])
        @test fm(-1) == [:FLFM, :FLFM]
        @test fm(-2) == [:RWFM, :RWFM]
        @test SigmaTau._lag1_label(-1.73, 1, -2) == -2    # not decimated: nearest integer
        # The readings follow from the generalized autocovariances the EDF uses.
        r1(alpha) = (s = t -> SigmaTau._compute_sw(Float64(t), alpha + 2);
                     (s(3) - 4s(2) + 7s(1)) / (2s(2) - 8s(1)))
        @test SigmaTau._LAG1_FLFM ≈ SigmaTau._lag1_reading(r1(-1), 2)
        @test SigmaTau._LAG1_RWFM ≈ SigmaTau._lag1_reading(r1(-2), 2) ≈ -2.4
    end

    @testset "flicker-walk and random-run FM" begin
        # A third difference (dmax = 3, the Hadamard family) tells the two
        # from random-walk FM; sampled, the three read α = −3.0, −3.93, −4.57.
        record(alpha) = noise_gen(PhaseData, 4096, 1.0; h=Dict(alpha => 1.0), rng=Xoshiro(11))
        h3(alpha) = identify_noise(record(alpha), [1, 16, 64]; dmax=3)
        @test h3(-2) == fill(:RWFM, 3)
        @test h3(-3) == fill(:FWFM, 3)
        @test h3(-4) == fill(:RRFM, 3)
        @test hdev(record(-3), [16]).noise_type == [:FWFM]
        r1(alpha) = (s = t -> SigmaTau._compute_sw(Float64(t), alpha + 2);
                     (s(4) - 6s(3) + 16s(2) - 26s(1)) / (2s(3) - 12s(2) + 30s(1)))
        @test SigmaTau._LAG1_RWFM3 ≈ SigmaTau._lag1_reading(r1(-2), 3) ≈ -3
        @test SigmaTau._LAG1_FWFM ≈ SigmaTau._lag1_reading(r1(-3), 3)
        @test SigmaTau._LAG1_RRFM ≈ SigmaTau._lag1_reading(r1(-4), 3)
        # HTOT has published bias factors down to α = −4 (Howe & Peppler 2001,
        # Table 3); MHTOT borrows random-walk FM's.
        @test bias_correction([:FWFM, :RRFM], :htot, [16.0, 16.0], 4095.0) ≈ [0.717, 0.679]
        @test allequal(bias_correction([:RWFM, :FWFM, :RRFM], :mhtot, fill(16.0, 3), 4095.0))
    end

    @testset "noise-ID scale invariance" begin
        # Regression for the eps*N threshold bug in `_lag1_acf` that
        # produced :unknown classifications on any phase record with std
        # below ~√eps ≈ 1.5e-8 (i.e. all real-world records, in seconds).
        # Same WPM+RWFM fixture used by the ADEV example: every
        # m must classify to a real power-law noise, never :unknown.
        using Random
        Random.seed!(20260509)
        N    = 4096
        tau0 = 1.0
        wpm  = randn(N) .* 1e-9                       # ~ns scale
        rwfm = cumsum(cumsum(randn(N) .* 1e-12))      # ~ps drift integrator
        x    = wpm .+ rwfm

        m_grid = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
        noises = identify_noise(x, m_grid; dmin=0, dmax=2)
        @test !any(==(:unknown), noises)
        # Sanity: every classification falls in the expected SP1065 alphabet.
        @test all(n -> n in (:WHPM, :FLPM, :WHFM, :FLFM, :RWFM), noises)

        # Same fixture rescaled by 1e6 (now in the µs band, std ~ 1e-3) —
        # classification must be invariant under positive linear rescaling.
        # This is the core scale-invariance property that the fixed
        # `_lag1_acf` guards.
        noises_scaled = identify_noise(x .* 1e6, m_grid; dmin=0, dmax=2)
        @test noises_scaled == noises
    end

    @testset "noise ID never returns :unknown" begin
        using Random
        known = (:WHPM, :FLPM, :WHFM, :FLFM, :RWFM)

        # Lag-1 estimates outside the labelled range clamp to the nearest end.
        # White PM scatters above α = 2.5 on some realizations, and
        # random-run FM (α = -4) lands below the Allan family's floor of -2;
        # both used to come back :unknown.
        rng = MersenneTwister(7)
        m_grid = [1, 2, 4, 8, 16, 32, 64]
        @test all(all(in(known), identify_noise(randn(rng, 2048), m_grid)) for _ in 1:50)
        rrfm = cumsum(cumsum(cumsum(randn(rng, 4096))))
        @test identify_noise(rrfm, m_grid) == fill(:RWFM, 7)

        # Averaging factors whose estimate fails take the successful one
        # nearest in log(m), including leading ones (previously :unknown),
        # with the shorter m winning a tie. No success at all → :WHFM.
        fill_syms = SigmaTau._fill_noise_symbols
        @test fill_syms([0, 0, -2, 0, 2], Bool[0, 0, 1, 0, 1], [1, 2, 4, 8, 16]) ==
              [:RWFM, :RWFM, :RWFM, :RWFM, :WHPM]
        @test fill_syms([1, 0, 0], Bool[1, 0, 0], [1, 2, 4]) == [:FLPM, :FLPM, :FLPM]
        @test fill_syms([0, 0], Bool[0, 0], [1, 2]) == [:WHFM, :WHFM]
        # Non-uniform grid: m = 10 is far closer to 11 than to 1.
        @test fill_syms([2, 0, -1], Bool[1, 0, 1], [1, 10, 11]) ==
              [:WHPM, :FLFM, :FLFM]

        # Degenerate records: constant phase and a record too short for any
        # estimate at the long averaging factors.
        @test identify_noise(fill(3.0, 256), [1, 2, 4]) == fill(:WHFM, 3)
        short = cumsum(randn(rng, 40))
        @test all(in(known), identify_noise(short, [1, 2, 4, 8, 13]))
    end

    @testset "noise-ID B1 length (b1_length)" begin
        # Long-τ points (N/m < NEFF_RELIABLE) go through the B1 fallback.
        # `b1_length=:averaged` evaluates B1 theory at the averaged-point
        # count and finds random-walk FM far more often than the
        # Stable32-compatible default `:record`.
        rng = MersenneTwister(20261108)
        ms = [40, 64, 128]
        hits = Dict(:record => 0, :averaged => 0)
        trials = 60
        for _ in 1:trials
            x = noise_gen(PhaseData, 1024, 1.0; h=Dict(-2 => 1.0), rng).x
            for mode in (:record, :averaged)
                hits[mode] += count(==(:RWFM), identify_noise(x, ms; b1_length=mode))
            end
        end
        total = trials * length(ms)
        @test hits[:averaged] / total > 0.4
        @test hits[:averaged] > 4 * hits[:record]
        # Default is :record, and the lag-1 ACF path (N/m ≥ 30) ignores it.
        x = noise_gen(PhaseData, 1024, 1.0; h=Dict(0 => 1.0), rng).x
        @test identify_noise(x, ms) == identify_noise(x, ms; b1_length=:record)
        @test identify_noise(x, [1, 4, 16]; b1_length=:averaged) == identify_noise(x, [1, 4, 16])
        @test identify_noise(PhaseData(x, 1.0), ms; b1_length=:averaged) ==
              identify_noise(x, ms; b1_length=:averaged)
        @test_throws ArgumentError identify_noise(x, ms; b1_length=:bogus)
    end

    @testset "B1(N, μ) closed-form parity vs allantools canonical" begin
        # Locks the closed-form B1 ratios used by the B1/R(n) noise-ID
        # fallback (`_noise_id_b1rn`) against the canonical allantools
        # `b1_theory(N, mu)` reference (Wallin 2018, citing Howe 2000).
        # Drift in any constant here changes which noise type the
        # B1 fallback picks at every m where lag-1 ACF can't run, so
        # the table is checked at machine precision.
        import SigmaTau as S

        # Canonical reference: allantools/ci.py b1_theory(N, mu).
        # mu ∈ {1, 0, -1, -2} are the values actually consumed by
        # _noise_id_b1rn (mu=2 is dead-code defensive, exercised below
        # for completeness).
        @test all([S._b1_theory(N, mu) for mu in (2, 1, 0, -1, -2)] ==
                  [N * (N + 1) / 6, N / 2, N * log(N) / (2 * (N - 1) * log(2)), 1.0,
                   (N^2 - 1) / (1.5 * N * (N - 1))]
                  for N in (50, 100, 1000, 10_000))

        # Generic Howe formula B1 = N(1 - N^μ) / (2(N-1)(1 - 2^μ))
        # — the closed forms above are special simplifications. Verify
        # the simplifications agree with the generic form at a sanity
        # μ value where both branches return finite values.
        generic(N) = N * (1 - Float64(N)^(-2)) / (2 * (N - 1) * (1 - 2.0^(-2)))
        @test all(isapprox(S._b1_theory(N, -2), generic(N); rtol=1e-14) for N in (100, 1000))
    end

    @testset "R(n)(af, b) closed-form parity vs allantools canonical" begin
        # Locks the R(n) = MVAR/AVAR closed-form used by the WPM/FLPM
        # disambiguation step in `_noise_id_b1rn`. b=0 is the
        # WHFM-shape (also covers WPM via slope ratio); b=-1 is the
        # FLFM-shape (also covers FLPM). All other b values fall
        # through to the af⁰ = 1 branch.
        import SigmaTau as S

        # b=0 branch: R(n) = 1/af for every af ≥ 1.
        afs = (1, 2, 4, 16, 64, 256, 1024)
        @test all(S._rn_theory(af, 0) == 1 / af for af in afs)

        # b=-1 branch: explicit Howe/Greenhall formula for the
        # FLFM/FLPM shape. avar = (1.038 + 3·log(2π·0.5·af)) / (4π²);
        # mvar = 3·log(256/27) / (8π²); R(n) = mvar/avar.
        avar(af) = (1.038 + 3 * log(2π * 0.5 * af)) / (4π^2)
        mvar = 3 * log(256 / 27) / (8π^2)
        @test all(isapprox(S._rn_theory(af, -1), mvar / avar(af); rtol=1e-15) for af in afs)

        # All other b values (including b=2, 1, -2) are not in the
        # closed-form table and fall through to the constant-1 branch.
        @test all(S._rn_theory(af, b) == 1.0 for af in (1, 4, 64), b in (1, 2, -2))
    end
end

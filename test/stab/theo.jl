# theo.jl — Thêo1 / ThêoBR / ThêoH tests (NIST SP1065 §5.2.15–5.2.16).
#
# The allantools reference values below were computed once with
# allantools 2024.06 (`theo1(x, rate=1.0, data_type="phase", taus=...)`)
# and frozen at %.17e. allantools implements the same eq. 30 statistic
# with the Howe & Peppler 0.75 normalization, but reports τ = m·τ0;
# SigmaTau reports the SP1065/Stable32 effective τ = 0.75·m·τ0, so the
# deviation values match at equal m while the τ columns differ by 0.75.

using Random
using SigmaTau: _theo1_core, _theobr_ratio, _theoh_theo_m, _theoh_is_adev,
                _theo1_edf, _neff_counts, _f64

@testset "Theo1 / ThêoBR / ThêoH" begin

    @testset "kernel sanity on the 5 power-law types" begin
        for alpha in (2, 1, 0, -1, -2)
            Random.seed!(100 + alpha)
            p = noise_gen(PhaseData, 2048, 1.0; sigma1 = Dict(alpha => 1e-12))
            ms = [10, 16, 32, 64, 128, 256]
            devs = _theo1_core(_f64(p.x), ms, 1.0)
            @test all(isfinite, devs)
            @test all(>(0.0), devs)
        end

        # Thêo1 has the same expected value as AVAR for white FM (SP1065
        # §5.2.15): raw Thêo1 at m = 1024 (τ_eff = 768) tracks overlapping
        # ADEV at m = 768 on the same record.
        Random.seed!(42)
        p = noise_gen(PhaseData, 4096, 1.0; sigma1 = Dict(0 => 1e-12))
        t1 = _theo1_core(_f64(p.x), [1024], 1.0)[1]
        ad = adev(p, [768]; ci=false).dev[1]
        @test isapprox(t1, ad; rtol=0.2)
    end

    @testset "even-m grid handling" begin
        ms = tau_values(Octave, 1000, :theo1)
        @test all(iseven, ms)
        @test ms[1] == 2                       # m = 1 rounds up to 2
        @test allunique(ms)
        @test issorted(ms)
        @test maximum(ms) <= 2 * ((1000 - 1) ÷ 2)
        @test SigmaTau._default_m_values(1000, :theo1) == ms

        # AllTaus: odd entries round to even and collapse.
        ms_all = tau_values(AllTaus, 21, :theo1)
        @test ms_all == [2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

        # :theoh grids are τ-grid factors — plain integers from 1 up to
        # ⌊0.75·(largest even ≤ N−1)⌋.
        mh = tau_values(Octave, 1000, :theoh)
        @test mh[1] == 1
        @test maximum(mh) <= floor(Int, 0.75 * (2 * ((1000 - 1) ÷ 2)))

        # Odd / out-of-range m yield NaN rows and zero neff from theo1.
        Random.seed!(1)
        p = PhaseData(randn(100), 1.0)
        r = theo1(p, [3, 4, 200]; ci=false, correct_bias=false)
        @test isnan(r.dev[1]) && r.neff[1] == 0     # odd
        @test isfinite(r.dev[2]) && r.neff[2] == 96 # N − m
        @test isnan(r.dev[3]) && r.neff[3] == 0     # m > N − 1
    end

    @testset "allantools cross-check (frozen, stable32gen.DAT)" begin
        dat_path = joinpath(@__DIR__, "..", "fixtures", "validation", "stable32gen.DAT")
        if !isfile(dat_path)
            @info "Stable32 phase fixture not present; skipping Theo1 allantools cross-check"
        else
            lines = readlines(dat_path)
            x = parse.(Float64, strip.(lines[11:end]))
            @test length(x) == 8192

            ms = [10, 16, 64, 256, 1024, 4096]
            # allantools 2024.06 theo1 devs at the same m (frozen 2026-06-10).
            ref = [2.09595443942596105e-01, 1.41170063482723346e-01,
                   4.24241065528329561e-02, 1.23509431294939673e-02,
                   3.88050646788194644e-03, 3.50353759490139739e-03]
            devs = _theo1_core(x, ms, 1.0)
            for (d, r) in zip(devs, ref)
                @test isapprox(d, r; rtol=1e-10)
            end

            # τ convention: SigmaTau reports the SP1065 effective τ = 0.75·m·τ0
            # (allantools reports m·τ0 for the same values).
            r = theo1(PhaseData(x, 1.0), ms; ci=false, correct_bias=false)
            @test r.tau ≈ 0.75 .* ms
            @test r.dev ≈ devs
        end
    end

    @testset "allantools cross-check (frozen, LCG synthetic)" begin
        # Deterministic LCG random-walk phase — exact integer recurrence,
        # bit-identical to the Python generator used to freeze the refs.
        N = 512
        s = 12345
        steps = Vector{Float64}(undef, N)
        for i in 1:N
            s = (1103515245 * s + 12345) % 2147483648
            steps[i] = s / 2147483648 - 0.5
        end
        x = cumsum(steps)

        ms = [2, 4, 8, 16, 64, 128, 256, 500]
        # allantools 2024.06 theo1 devs at the same m (frozen 2026-06-10).
        ref = [2.42475715196185104e-01, 1.68279473538979846e-01,
               1.13241924036000613e-01, 8.68106704528284201e-02,
               4.18605701646885126e-02, 2.73803955576344751e-02,
               1.91500560772033153e-02, 1.00433496067985177e-02]
        devs = _theo1_core(x, ms, 1.0)
        for (d, r) in zip(devs, ref)
            @test isapprox(d, r; rtol=1e-10)
        end
    end

    @testset "ThêoBR bias-correction behavior" begin
        # The eq. 33 correction is a single scalar variance ratio: the
        # corrected deviations are exactly raw·√ratio at every τ.
        Random.seed!(7)
        p = noise_gen(PhaseData, 2048, 1.0; sigma1 = Dict(0 => 1e-12))
        ms = [16, 32, 64, 128]
        raw = theo1(p, ms; ci=false, correct_bias=false)
        cor = theo1(p, ms; ci=false, correct_bias=true)
        ratio = _theobr_ratio(_f64(p.x), 1.0)
        @test cor.dev ≈ raw.dev .* sqrt(ratio)
        @test cor.tau == raw.tau
        # Thêo1 is unbiased for white FM, so the ratio sits near 1
        # (eq. 33's long-τ terms make it noisy; bounds are deliberately loose).
        @test 0.4 < ratio < 1.8

        # White PM: Table 2 bias a + b/m^c ≈ 0.09 + 0.74/m^0.4 — the ratio
        # average lands well below 1.
        Random.seed!(3)
        ppm = noise_gen(PhaseData, 2048, 1.0; sigma1 = Dict(2 => 1e-12))
        @test 0.05 < _theobr_ratio(_f64(ppm.x), 1.0) < 0.45

        # Records too short for the eq. 33 ladder (N < 19): NaN, and the
        # corrected wrapper propagates it.
        @test isnan(_theobr_ratio(randn(16), 1.0))
        rshort = theo1(PhaseData(randn(16), 1.0), [4]; ci=false)
        @test isnan(rshort.dev[1])
    end

    @testset "ThêoH assembly and crossover continuity" begin
        Random.seed!(42)
        N = 4096
        p = noise_gen(PhaseData, N, 1.0; sigma1 = Dict(0 => 1e-12))
        k = 0.2 * (N - 1)                      # eq. 34 crossover (≈ 819 τ0)

        # Branch bookkeeping on an octave grid: ADEV rows report τ = m·τ0,
        # Thêo1 rows report τ = 0.75·m_theo·τ0 with even m_theo ≈ 4m/3.
        mh = tau_values(Octave, N, :theoh)
        r = theoh(p, mh)
        @test r.deviation_type === :theoh
        @test all(isfinite, r.dev)
        for (j, m) in enumerate(mh)
            if _theoh_is_adev(m, N)
                @test r.tau[j] == m * 1.0
                @test r.neff[j] == N - 2m
            else
                mt = _theoh_theo_m(m)
                @test iseven(mt)
                @test r.tau[j] == 0.75 * mt
                @test r.neff[j] == N - mt
            end
        end

        # Continuity at the crossover: adjacent grid points either side of
        # k estimate the same AVAR on the same record, so the step between
        # them is small (seeded → deterministic; observed ≈ 2.5 %).
        rd = theoh(p, [floor(Int, k), floor(Int, k) + 1]; ci=false)
        @test _theoh_is_adev(floor(Int, k), N)
        @test !_theoh_is_adev(floor(Int, k) + 1, N)
        @test abs(rd.dev[2] - rd.dev[1]) / rd.dev[1] < 0.10

        # The ADEV segment is exactly overlapping ADEV.
        adev_ms = [m for m in mh if _theoh_is_adev(m, N)]
        ra = adev(p, adev_ms; ci=false)
        @test r.dev[1:length(adev_ms)] ≈ ra.dev
    end

    @testset "neff / EDF / CI population" begin
        Random.seed!(11)
        p = noise_gen(PhaseData, 1024, 1.0; sigma1 = Dict(0 => 1e-12))

        r1 = theo1(p, [16, 32, 64])
        @test r1.neff == _neff_counts(:theo1, 1024, [16, 32, 64])
        @test all(>(0), r1.neff)
        @test length(r1.ci) == 3 && length(r1.edf) == 3 && length(r1.noise_type) == 3
        @test all(isfinite, r1.edf)
        @test all(i -> r1.ci[i].lo <= r1.dev[i] <= r1.ci[i].hi, 1:3)

        rh = theoh(p, [1, 4, 16, 300, 600])
        @test rh.neff == _neff_counts(:theoh, 1024, [1, 4, 16, 300, 600])
        @test all(>(0), rh.neff)
        @test length(rh.ci) == 5 && length(rh.edf) == 5
        @test all(i -> rh.ci[i].lo <= rh.dev[i] <= rh.ci[i].hi, 1:5)

        # ci=false leaves noise/ci/edf empty, neff populated.
        r0 = theo1(p, [16]; ci=false)
        @test isempty(r0.ci) && isempty(r0.edf) && isempty(r0.noise_type)
        @test r0.neff == [1024 - 16]

        # Spot-check the SP1065 Table 3 W FM EDF formula by hand:
        # r = 0.75·64 = 48, N = 1024.
        rr = 48.0
        expected = ((4.1 * 1024 + 0.8) / rr - (3.1 * 1024 + 6.5) / 1024) *
                   (rr^1.5 / (rr^1.5 + 5.2))
        @test _theo1_edf(0, 64, 1024) ≈ expected
        # Theo1's selling point: more EDF than ADEV at the same effective τ.
        @test _theo1_edf(0, 64, 1024) > 0
    end

    @testset "dispatch chains and suite integration" begin
        Random.seed!(5)
        p = noise_gen(PhaseData, 512, 1.0; sigma1 = Dict(0 => 1e-12))
        fd = FrequencyData(diff(p.x) ./ 1.0, 1.0)
        pd_eq = PhaseData(cumsum(diff(p.x)), 1.0)

        # FrequencyData round-trips through _freq_to_phase.
        a = theo1(fd, [4, 8]; ci=false)
        b = theo1(pd_eq, [4, 8]; ci=false)
        @test a.dev ≈ b.dev
        ah = theoh(fd, [2, 4]; ci=false)
        bh = theoh(pd_eq, [2, 4]; ci=false)
        @test ah.dev ≈ bh.dev

        # TauMode and zero-arg conveniences.
        @test theo1(p, Octave; ci=false).deviation_type === :theo1
        @test theoh(p, HalfOctave; ci=false).deviation_type === :theoh
        @test theo1(p; ci=false).deviation_type === :theo1
        @test theoh(p; ci=false).deviation_type === :theoh

        # Suite integration with per-deviation grids.
        suite = stability(p; devs=(:adev, :theo1, :theoh), ci=false)
        @test keys(suite) == [:adev, :theo1, :theoh]
        @test all(iseven(round(Int, t / 0.75)) for t in suite[:theo1].tau)
    end
end

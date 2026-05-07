using Test
using SigmaTauBase
using SigmaTauStability
using SigmaTauStability: NEFF_RELIABLE

include("legacy_kernels.jl")
const LK = LegacyKernels

@testset "SigmaTauStability Math Core Parity" begin
    # Mock data equivalent to NIST SP1065 sample
    x = collect(1.0:100.0) .^ 2  # simple quadratic
    m_values = [1, 2, 4]
    tau0 = 1.0

    @testset "ADEV Core" begin
        devs = _adev_core(x, m_values, tau0)
        @test length(devs) == 3
        @test all(isfinite, devs)
        # 2nd diff of x^2 is constant 2. 
        # For m=1, d2 = x[i+2] - 2x[i+1] + x[i] = 2. sum_sq = sum(4) = 4 * 98. 
        # dev = sqrt(4 * 98 / (2 * 98 * 1^2 * 1^2)) = sqrt(2) ≈ 1.414
        @test isapprox(devs[1], sqrt(2.0))
    end

    @testset "MDEV Core" begin
        devs = _mdev_core(x, m_values, tau0)
        @test length(devs) == 3
        @test all(isfinite, devs)
    end
    
    @testset "TDEV Core" begin
        devs = _tdev_core(x, m_values, tau0)
        @test length(devs) == 3
        @test all(isfinite, devs)
    end
    
    @testset "HDEV & MHDEV Core" begin
        # 3rd difference of x^2 is 0. So dev should be very close to 0.
        devs_h = _hdev_core(x, m_values, tau0)
        devs_mh = _mhdev_core(x, m_values, tau0)
        @test length(devs_h) == 3
        @test all(isfinite, devs_h)
        @test all(d -> d < 1e-10, devs_h)
        @test length(devs_mh) == 3
        @test all(isfinite, devs_mh)
        @test all(d -> d < 1e-10, devs_mh)
    end

    @testset "Total Deviations Core" begin
        # Just ensure they run and produce finite values.
        devs_tot = _totdev_core(x, m_values, tau0)
        devs_mtot = _mtotdev_core(x, m_values, tau0)
        devs_htot = _htotdev_core(x, m_values, tau0)
        devs_mhtot = _mhtotdev_core(x, m_values, tau0)
        
        @test length(devs_tot) == 3
        @test all(isfinite, devs_tot)
        @test length(devs_mtot) == 3
        @test all(isfinite, devs_mtot)
        @test length(devs_htot) == 3
        @test all(isfinite, devs_htot)
        @test length(devs_mhtot) == 3
        @test all(isfinite, devs_mhtot)
    end
    
    @testset "Typed API Wrappers & Stats (Lag1 and B1-ratio)" begin
        using Random
        Random.seed!(42)
        # N=100. For m=1, N_eff = 100 >= NEFF_RELIABLE (30) -> Lag1 ACF used.
        # For m=4, N_eff = 25 < NEFF_RELIABLE             -> B1-Ratio used.
        pd = PhaseData(x .+ randn(100), tau0)
        
        res_adev = adev(pd, m_values)
        @test length(res_adev.dev) == 3
        @test res_adev.noise_type[1] != :unknown
        @test all(isfinite, res_adev.ci_lower)
        @test all(isfinite, res_adev.ci_upper)
        
        res_hdev = hdev(pd, m_values)
        @test length(res_hdev.dev) == 3
        @test all(isfinite, res_hdev.ci_upper)
        
        res_tot = totdev(pd, m_values)
        @test length(res_tot.dev) == 3
        @test all(isfinite, res_tot.ci_upper)
        
        res_ldev = ldev(pd, m_values)
        @test length(res_ldev.dev) == 3

        # tdev wraps mdev with the τ/√3 scaling — point estimates and CI bounds
        # must be consistent with that identity.
        res_tdev = tdev(pd, m_values)
        res_mdev_ref = mdev(pd, m_values)
        factor = res_tdev.tau ./ sqrt(3.0)
        @test length(res_tdev.dev) == 3
        @test res_tdev.deviation_type == :tdev
        @test res_tdev.dev ≈ res_mdev_ref.dev .* factor
        @test res_tdev.ci_lower ≈ res_mdev_ref.ci_lower .* factor
        @test res_tdev.ci_upper ≈ res_mdev_ref.ci_upper .* factor

        # calc_ci=false path returns empty CI and EDF vectors.
        res_tdev_nci = tdev(pd, m_values; calc_ci=false)
        @test isempty(res_tdev_nci.ci_lower)
        @test isempty(res_tdev_nci.ci_upper)
        @test isempty(res_tdev_nci.edf)

        # edf is populated when calc_ci=true.
        @test length(res_adev.edf) == length(m_values)
        @test all(isfinite, res_adev.edf)
        @test all(>=(0.0), res_adev.edf)
    end

    @testset "Legacy parity (extracted SP1065 kernels)" begin
        # Strict numerical parity against the legacy SigmaTau Julia reference
        # kernels, inlined verbatim under test/legacy_kernels.jl. Uses a seeded
        # power-law-flavoured fixture and a wide m-grid.
        using Random
        Random.seed!(20260507)

        N    = 4096
        tau0 = 1.0
        # Mixed power-law phase fixture: WPM + RWFM components.
        wpm  = randn(N) .* 1e-9
        rwfm = cumsum(cumsum(randn(N) .* 1e-12))
        x    = wpm .+ rwfm

        # cumsum prefix used by modified kernels.
        x_cs = pushfirst!(cumsum(x), 0.0)

        # Octave-spaced m-grid with extra coverage near the small-m regime.
        m_grid = [1, 2, 4, 8, 16, 32, 64, 128]

        rt = 1e-12   # rtol for direct kernel-vs-kernel parity
        at = 1e-25   # atol for near-zero comparisons

        # ADEV: 8 m-values × 1 kernel = 8 assertions.
        for m in m_grid
            new_dev = sqrt(LK.adev_var(x, m, tau0))
            @test SigmaTauStability._adev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # MDEV.
        for m in m_grid
            new_dev = sqrt(LK.mdev_var(x, m, tau0, x_cs))
            @test SigmaTauStability._mdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # HDEV.
        for m in m_grid
            new_dev = sqrt(LK.hdev_var(x, m, tau0))
            @test SigmaTauStability._hdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # MHDEV.
        for m in m_grid
            new_dev = sqrt(LK.mhdev_var(x, m, tau0, x_cs))
            @test SigmaTauStability._mhdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # TOTDEV. (Smaller grid — kernel is O(N) per m but allocates an extended
        # 3N-4 array each call.)
        for m in [1, 2, 4, 8, 16, 32]
            new_dev = sqrt(LK.totdev_var(x, m, tau0))
            @test SigmaTauStability._totdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # MTOTDEV.
        for m in [1, 2, 4, 8, 16]
            new_dev = sqrt(LK.mtotdev_var(x, m, tau0))
            @test SigmaTauStability._mtotdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # HTOTDEV.
        for m in [1, 2, 4, 8, 16]
            new_dev = sqrt(LK.htotdev_var(x, m, tau0))
            @test SigmaTauStability._htotdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # MHTOTDEV.
        for m in [1, 2, 4, 8]
            new_dev = sqrt(LK.mhtotdev_var(x, m, tau0))
            @test SigmaTauStability._mhtotdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end
    end

    @testset "FrequencyData ↔ PhaseData equivalence" begin
        # adev(FrequencyData(y, τ₀)) must equal adev(PhaseData(cumsum(y)·τ₀, τ₀)).
        # Spot-checked on adev (Allan family) and hdev (Hadamard family).
        using Random
        Random.seed!(7)

        tau0 = 1.0
        y = randn(200) .* 1e-9          # fractional frequency
        fd = FrequencyData(y, tau0)
        pd_equiv = PhaseData(cumsum(y) .* tau0, tau0)

        m_values_eq = [1, 2, 4, 8]

        @test adev(fd, m_values_eq; calc_ci=false).dev ≈
              adev(pd_equiv, m_values_eq; calc_ci=false).dev

        @test hdev(fd, m_values_eq; calc_ci=false).dev ≈
              hdev(pd_equiv, m_values_eq; calc_ci=false).dev
    end

    @testset "MTOTDEV across noise regimes" begin
        # Verify _mtotdev_core (and the :mtot bias correction) behaves on the
        # three power-law noise types whose synthesis is FFT-free: WPM, WHFM,
        # RWFM. (FLPM/FLFM need a 1/f filter and are exercised indirectly via
        # the mixed-noise fixture in the legacy parity testset above.)
        using Random
        Random.seed!(20260507)
        N    = 1024
        tau0 = 1.0
        ms   = [1, 2, 4, 8, 16]
        rt   = 1e-12

        noise_fixtures = (
            WPM  = randn(N) .* 1e-9,
            WHFM = cumsum(randn(N) .* 1e-9),
            RWFM = cumsum(cumsum(randn(N) .* 1e-12)),
        )

        for (label, x) in pairs(noise_fixtures)
            for m in ms
                ref = sqrt(LK.mtotdev_var(x, m, tau0))
                got = SigmaTauStability._mtotdev_core(x, [m], tau0)[1]
                @test got ≈ ref atol=1e-25 rtol=rt
            end

            # End-to-end pipeline (bias correction + CI bounds) should produce
            # finite, ordered outputs on all three noise types.
            res = mtotdev(PhaseData(x, tau0), ms)
            @test res.deviation_type == :mtotdev
            @test all(isfinite, res.dev)
            @test all(.>(0.0), res.dev)
            @test all(.<=(0.0), res.ci_lower .- res.dev)   # lower ≤ dev
            @test all(.>=(0.0), res.ci_upper .- res.dev)   # upper ≥ dev
        end
    end

    @testset "TOTDEV/HTOTDEV EDF fallback for WPM/FLPM" begin
        # _coeff_totvar returns (NaN, NaN) for α=2,1 (SP1065 Table 9 only
        # covers α∈{0,-1,-2}). The fallback dispatches to the ADEV-style
        # EDF formula so calculate_edf produces a finite EDF for every
        # noise type.
        m_values_eq = [1, 2, 4, 8, 16]
        taus_eq     = Float64.(m_values_eq)
        N_eq        = 4096
        T_eq        = (N_eq - 1) * 1.0
        # α=2 (WPM)
        noises_wpm = fill(:WHPM, length(m_values_eq))
        edfs = SigmaTauStability.calculate_edf(:totdev, ones(length(m_values_eq)),
                                               noises_wpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs)
        @test all(>(0.0), edfs)

        # α=1 (FLPM)
        noises_flpm = fill(:FLPM, length(m_values_eq))
        edfs = SigmaTauStability.calculate_edf(:totdev, ones(length(m_values_eq)),
                                               noises_flpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs)

        # htotdev gets the HDEV-style fallback at α=2,1 too.
        edfs_h = SigmaTauStability.calculate_edf(:htotdev, ones(length(m_values_eq)),
                                                 noises_wpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs_h)
    end

    @testset "NEFF_RELIABLE boundary" begin
        # NEFF_RELIABLE = 30 per legacy GEMINI.md §2 mandate. The boundary
        # determines whether identify_noise uses lag-1 ACF (N_eff ≥ threshold)
        # or the B1-ratio fallback. Construct two cases that straddle the
        # boundary and verify both produce a finite, classified noise type.
        @test NEFF_RELIABLE == 30

        using Random
        Random.seed!(2026)

        # m=1 with N=29 → N_eff=29 (one below) → B1-ratio path.
        pd_below = PhaseData(cumsum(randn(29)), 1.0)
        noises_below = identify_noise(pd_below.x, [1])
        @test noises_below[1] != :unknown

        # m=1 with N=31 → N_eff=31 (one above) → lag-1 ACF path.
        pd_above = PhaseData(cumsum(randn(31)), 1.0)
        noises_above = identify_noise(pd_above.x, [1])
        @test noises_above[1] != :unknown
    end
end

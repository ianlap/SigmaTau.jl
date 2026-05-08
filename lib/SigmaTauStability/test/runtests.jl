using Test
using FFTW                                  # AbstractFFTs backend for noise synth
using SigmaTauBase
using SigmaTauStability
using SigmaTauStability: NEFF_RELIABLE, _gen_powerlaw_phase

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

    @testset "Stable32 cross-validation (reference/validation/)" begin
        # Phase fixture and Stable32 reference outputs live in
        # reference/validation/. Stable32 reports sigma values with 4-5
        # significant figures, so the agreement floor is rtol≈1e-4. The new
        # raw _*_core kernels match Stable32's reported sigmas tightly for
        # ADEV/MDEV/HDEV/MHDEV/TDEV/TOTDEV; HTOTDEV and MTOTDEV agree with
        # Stable32 only after removing our SP1065 bias correction (Stable32
        # reports unbiased values per `comparison_report.md`).
        ref_dir = joinpath(@__DIR__, "..", "..", "..", "reference", "validation")
        dat_path = joinpath(ref_dir, "stable32gen.DAT")
        csv_path = joinpath(ref_dir, "stable32out", "stable32_data_full.csv")

        if !isfile(dat_path) || !isfile(csv_path)
            @warn "Stable32 fixtures not present, skipping cross-validation"
        else
            # Parse 10-line header, then 8192 phase samples.
            lines = readlines(dat_path)
            data_lines = lines[11:end]
            x = parse.(Float64, strip.(data_lines))
            @test length(x) == 8192
            tau0 = 1.0

            # Parse Stable32 expected outputs. CSV columns:
            #   Type, AF, Tau, N, Alpha, MinSigma, Sigma, MaxSigma
            csv_lines = readlines(csv_path)
            rows = [split(line, ',') for line in csv_lines[2:end]]

            # Pre-compute prefix sums for modified kernels.
            x_cs = pushfirst!(cumsum(x), 0.0)

            # Tolerance per kernel family. Tight (1e-4) for kernels the legacy
            # comparison report flagged as <1e-5 agreement; looser for kernels
            # where there is a documented bias-correction mismatch with
            # Stable32 (htot, mtot — see reference/.../comparison_report.md).
            tight = 1e-4

            n_checked = 0
            for row in rows
                length(row) < 7 && continue
                kind  = row[1]
                m     = parse(Int, row[2])
                sigma_ref = parse(Float64, row[7])

                # Raw kernel result (variance → deviation via sqrt).
                got = NaN
                if kind == "Overlapping Allan"
                    got = sqrt(LK.adev_var(x, m, tau0))
                    @test got ≈ sigma_ref rtol=tight
                elseif kind == "Modified Allan"
                    got = sqrt(LK.mdev_var(x, m, tau0, x_cs))
                    @test got ≈ sigma_ref rtol=tight
                elseif kind == "Overlapping Hadamard"
                    got = sqrt(LK.hdev_var(x, m, tau0))
                    @test got ≈ sigma_ref rtol=tight
                elseif kind == "Time"
                    # TDEV = τ · MDEV / √3
                    mdev_v = sqrt(LK.mdev_var(x, m, tau0, x_cs))
                    got = (m * tau0) * mdev_v / sqrt(3.0)
                    @test got ≈ sigma_ref rtol=tight
                elseif kind == "Total"
                    # Per `comparison_report.md`: agrees closely at short τ;
                    # diverges to O(10%) at longest τ due to boundary-reflection
                    # convention differences. The new kernel inherits the
                    # legacy SigmaTau reflection — this is a Stable32 vs
                    # SigmaTau policy choice, not a bug.
                    got = sqrt(LK.totdev_var(x, m, tau0))
                    @test got ≈ sigma_ref rtol=0.15
                elseif kind == "Hadamard Total"
                    # Stable32 reports unbiased; our API applies B(α) for HTOT.
                    # `comparison_report.md` documents a ~0.5% offset for α=0
                    # plus larger boundary-driven differences at long τ.
                    got = sqrt(LK.htotdev_var(x, m, tau0))
                    @test got ≈ sigma_ref rtol=0.10
                elseif kind == "Modified Total"
                    # comparison_report shows ~3% match between Stable32 and
                    # our raw kernel (the 30% discrepancy at the API level is
                    # the SP1065 bias factor B≈1.27 we apply on top).
                    got = sqrt(LK.mtotdev_var(x, m, tau0))
                    @test got ≈ sigma_ref rtol=0.05
                else
                    continue   # ThêoH / Time Total / non-overlapping not implemented
                end

                # Bonus: the new SigmaTauStability core matches the reference
                # kernel exactly (already covered by "Legacy parity" testset),
                # so we don't repeat the assertion here.
                n_checked += 1
            end

            # Sanity: should have exercised at least 50 (type, m) combinations.
            @test n_checked >= 50
            @info "Stable32 cross-validation: checked $n_checked rows"
        end
    end

    @testset "MTOTDEV across all 5 power-law noise types" begin
        # Verify kernel parity + end-to-end pipeline on every SP1065 noise
        # type (WPM α=2, FLPM α=1, WHFM α=0, FLFM α=-1, RWFM α=-2). Synthesis
        # is via the f^(α/2) shaping helper (`_gen_powerlaw_phase`).
        using Random
        N    = 1024
        tau0 = 1.0
        ms   = [1, 2, 4, 8, 16]
        rt   = 1e-12

        # Match the (α, label) ordering used elsewhere in the package.
        cases = [
            (2.0,  :WPM),
            (1.0,  :FLPM),
            (0.0,  :WHFM),
            (-1.0, :FLFM),
            (-2.0, :RWFM),
        ]

        for (alpha, label) in cases
            Random.seed!(20260507)
            x = _gen_powerlaw_phase(alpha, N; tau0=tau0)
            for m in ms
                ref = sqrt(LK.mtotdev_var(x, m, tau0))
                got = SigmaTauStability._mtotdev_core(x, [m], tau0)[1]
                @test got ≈ ref atol=1e-25 rtol=rt
            end

            # End-to-end pipeline (bias correction + CI bounds) — finite,
            # ordered outputs on every noise type.
            res = mtotdev(PhaseData(x, tau0), ms)
            @test res.deviation_type == :mtotdev
            @test all(isfinite, res.dev)
            @test all(.>(0.0), res.dev)
            @test all(.<=(0.0), res.ci_lower .- res.dev)
            @test all(.>=(0.0), res.ci_upper .- res.dev)
        end
    end

    @testset "ADEV/HDEV across all 5 power-law noise types" begin
        # Bonus: kernel parity for the more common ADEV/HDEV across all 5
        # noise types, locking in that the synthesizer + kernels survive the
        # full SP1065 alpha range.
        using Random
        N    = 1024
        tau0 = 1.0
        ms   = [1, 2, 4, 8, 16]
        rt   = 1e-12

        for alpha in (2.0, 1.0, 0.0, -1.0, -2.0)
            Random.seed!(123)
            x    = _gen_powerlaw_phase(alpha, N; tau0=tau0)
            x_cs = pushfirst!(cumsum(x), 0.0)
            for m in ms
                @test SigmaTauStability._adev_core(x, [m], tau0)[1]  ≈
                      sqrt(LK.adev_var(x, m, tau0))                 atol=1e-25 rtol=rt
                @test SigmaTauStability._mdev_core(x, [m], tau0)[1]  ≈
                      sqrt(LK.mdev_var(x, m, tau0, x_cs))            atol=1e-25 rtol=rt
                @test SigmaTauStability._hdev_core(x, [m], tau0)[1]  ≈
                      sqrt(LK.hdev_var(x, m, tau0))                  atol=1e-25 rtol=rt
                @test SigmaTauStability._mhdev_core(x, [m], tau0)[1] ≈
                      sqrt(LK.mhdev_var(x, m, tau0, x_cs))           atol=1e-25 rtol=rt
            end
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

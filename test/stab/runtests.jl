using Test
using Random
using FFTW                                  # AbstractFFTs backend for noise synth
using SigmaTau
using SigmaTau.Stab: NEFF_RELIABLE, _gen_powerlaw_phase

include("legacy_kernels.jl")
const LK = LegacyKernels

@testset "SigmaTau.Stab Math Core Parity" begin
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
        
        res_htdev = htdev(pd, m_values)
        @test length(res_htdev.dev) == 3
        @test res_htdev.deviation_type == :htdev
        # ldev is a deprecated alias and must give bit-identical output for now.
        res_ldev_alias = ldev(pd, m_values)
        @test res_ldev_alias.dev == res_htdev.dev
        @test res_ldev_alias.deviation_type == :htdev

        # htdev wraps mhdev with the τ/√(10/3) scaling — point estimates AND
        # CI bounds inherit the same multiplicative factor. Closes R-MED-5
        # (HTDEV CI scaling formal verification): for a non-degenerate
        # Y = c · X with c constant in τ, the χ²-based CI scales linearly in
        # c, so [Y.lo, Y.hi] = c · [X.lo, X.hi] is the correct operation.
        res_mhdev_ref = mhdev(pd, m_values)
        h_factor = res_htdev.tau ./ sqrt(10.0 / 3.0)
        @test res_htdev.dev      ≈ res_mhdev_ref.dev      .* h_factor
        @test res_htdev.ci_lower ≈ res_mhdev_ref.ci_lower .* h_factor
        @test res_htdev.ci_upper ≈ res_mhdev_ref.ci_upper .* h_factor
        # χ² invariant: the multiplicative scaling cancels in the ratio
        # CI_bound / dev, so HTDEV's relative CI structure is identical to
        # MHDEV's at every τ.
        @test res_htdev.ci_lower ./ res_htdev.dev ≈ res_mhdev_ref.ci_lower ./ res_mhdev_ref.dev
        @test res_htdev.ci_upper ./ res_htdev.dev ≈ res_mhdev_ref.ci_upper ./ res_mhdev_ref.dev

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
            @test SigmaTau.Stab._adev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # MDEV.
        for m in m_grid
            new_dev = sqrt(LK.mdev_var(x, m, tau0, x_cs))
            @test SigmaTau.Stab._mdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # HDEV.
        for m in m_grid
            new_dev = sqrt(LK.hdev_var(x, m, tau0))
            @test SigmaTau.Stab._hdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # MHDEV.
        for m in m_grid
            new_dev = sqrt(LK.mhdev_var(x, m, tau0, x_cs))
            @test SigmaTau.Stab._mhdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # TOTDEV. (Smaller grid — kernel is O(N) per m but allocates an extended
        # 3N-4 array each call.) Pass detrend=:legacy explicitly: the new
        # _totdev_core default is :howe (SP1065 eqn 25), which differs from
        # the LK reference by O(few %) at short τ.
        for m in [1, 2, 4, 8, 16, 32]
            new_dev = sqrt(LK.totdev_var(x, m, tau0))
            @test SigmaTau.Stab._totdev_core(x, [m], tau0; detrend=:legacy)[1] ≈ new_dev atol=at rtol=rt
        end

        # MTOTDEV.
        for m in [1, 2, 4, 8, 16]
            new_dev = sqrt(LK.mtotdev_var(x, m, tau0))
            @test SigmaTau.Stab._mtotdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # HTOTDEV.
        for m in [1, 2, 4, 8, 16]
            new_dev = sqrt(LK.htotdev_var(x, m, tau0))
            @test SigmaTau.Stab._htotdev_core(x, [m], tau0)[1] ≈ new_dev atol=at rtol=rt
        end

        # MHTOTDEV. Pass detrend=:legacy explicitly: the new _mhtotdev_core
        # default is :greenhall (per-window half-mean), which differs from
        # the LK reference's per-window full-LS detrend.
        for m in [1, 2, 4, 8]
            new_dev = sqrt(LK.mhtotdev_var(x, m, tau0))
            @test SigmaTau.Stab._mhtotdev_core(x, [m], tau0; detrend=:legacy)[1] ≈ new_dev atol=at rtol=rt
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
        ref_dir = joinpath(@__DIR__, "..", "..", "reference", "validation")
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

                # Bonus: the new SigmaTau.Stab core matches the reference
                # kernel exactly (already covered by "Legacy parity" testset),
                # so we don't repeat the assertion here.
                n_checked += 1
            end

            # Sanity: should have exercised at least 50 (type, m) combinations.
            @test n_checked >= 50
            @info "Stable32 cross-validation: checked $n_checked rows"
        end
    end

    @testset "TOTDEV :howe matches Stable32 tightly" begin
        # SP1065 eqn 25 reference: no detrend, mean-flip endpoint reflection.
        # Matches Stable32's TOTDEV output at rtol=1e-4 (vs the rtol=0.15
        # boundary-policy floor seen with the :legacy global-LS detrend), and
        # tracks allantools' raw TOTDEV output to ~7 significant figures
        # (allantools follows SP1065 verbatim).
        #
        # Exception: the m=512 row in this fixture is identified as FLFM
        # (alpha=-1) and Stable32's reported sigma is ~1.5% larger than the
        # raw SP1065 value (allantools 3.055835e-03 vs Stable32 3.1028e-03);
        # Stable32 appears to apply alpha-aware correction opaquely at that
        # one point. Skipped here, tracked in TODO.md for follow-up against
        # the allantools cross-validation fixture.
        ref_dir = joinpath(@__DIR__, "..", "..", "reference", "validation")
        dat_path = joinpath(ref_dir, "stable32gen.DAT")
        csv_path = joinpath(ref_dir, "stable32out", "stable32_data_full.csv")

        if !isfile(dat_path) || !isfile(csv_path)
            @warn "Stable32 fixtures not present, skipping :howe TOTDEV tightness test"
        else
            lines = readlines(dat_path)
            x = parse.(Float64, strip.(lines[11:end]))
            @test length(x) == 8192
            tau0 = 1.0

            rows = [split(line, ',') for line in readlines(csv_path)[2:end]]
            n_checked = 0
            for row in rows
                length(row) < 7 && continue
                row[1] == "Total" || continue
                m = parse(Int, row[2])
                m == 512 && continue                 # Stable32-only quirk; see comment above
                sigma_ref = parse(Float64, row[7])

                got = SigmaTau.Stab._totdev_core(x, [m], tau0; detrend=:howe)[1]
                @test got ≈ sigma_ref rtol=1e-4
                n_checked += 1
            end
            @test n_checked >= 5
        end
    end

    @testset "TOTDEV :linear ≡ :legacy" begin
        # For TOTDEV, :linear (global LS detrend + endpoint mean-flip) is
        # numerically identical to :legacy by construction. Lock the alias.
        using Random
        Random.seed!(20260508)
        N    = 1024
        tau0 = 1.0
        ms   = [1, 2, 4, 8, 16]
        x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

        devs_linear = SigmaTau.Stab._totdev_core(x, ms, tau0; detrend=:linear)
        devs_legacy = SigmaTau.Stab._totdev_core(x, ms, tau0; detrend=:legacy)

        @test length(devs_linear) == length(ms)
        @test all(isfinite, devs_linear)
        @test all(>(0), devs_linear)
        for k in eachindex(ms)
            @test devs_linear[k] ≈ devs_legacy[k] atol=0.0 rtol=1e-15
        end

        # Unsupported recipes (:greenhall on TOTDEV) raise ArgumentError.
        @test_throws ArgumentError SigmaTau.Stab._totdev_core(x, [1], tau0; detrend=:greenhall)
        @test_throws ArgumentError SigmaTau.Stab._totdev_core(x, [1], tau0; detrend=:nonsense)
    end

    include("allantools_cross_validation.jl")

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
                got = SigmaTau.Stab._mtotdev_core(x, [m], tau0)[1]
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

    @testset "MTOTDEV :linear smoke" begin
        # Per-window full-LS detrend (vs :greenhall's half-mean) on the same
        # time-reverse extension. Same kernel operator, same magnitude band.
        using Random
        Random.seed!(20260508)
        N    = 1024
        tau0 = 1.0
        ms   = [1, 2, 4, 8, 16]
        x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

        devs = SigmaTau.Stab._mtotdev_core(x, ms, tau0; detrend=:linear)
        @test length(devs) == length(ms)
        @test all(isfinite, devs)
        @test all(>(0), devs)
        devs_legacy = SigmaTau.Stab._mtotdev_core(x, ms, tau0; detrend=:legacy)
        @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)

        # :howe is no longer a recipe for MTOTDEV.
        @test_throws ArgumentError SigmaTau.Stab._mtotdev_core(x, [1], tau0; detrend=:howe)
    end

    @testset "TTOTDEV = MTOTDEV · τ/√3 identity" begin
        # ttotdev wraps mtotdev with a τ/√3 rescaling (analogous to tdev
        # over mdev). Verifies the wrapper produces the documented
        # relationship and preserves CI / EDF / noise-type alignment.
        Random.seed!(20260513)
        N    = 1024
        tau0 = 1.0
        x    = _gen_powerlaw_phase(0.0, N; tau0=tau0)
        pd   = PhaseData(x, tau0)
        ms   = [1, 2, 4, 8, 16, 32]

        # Centerline identity (correct_bias=false isolates the wrapper).
        rm = mtotdev(pd, ms; calc_ci=false, correct_bias=false)
        rt = ttotdev(pd, ms; calc_ci=false, correct_bias=false)
        @test rt.tau == rm.tau
        @test isapprox(rt.dev, rm.dev .* (rm.tau ./ sqrt(3.0)); rtol=1e-14)

        # CI / EDF / α flow through with the same τ/√3 scaling.
        rm_ci = mtotdev(pd, ms; calc_ci=true, correct_bias=true)
        rt_ci = ttotdev(pd, ms; calc_ci=true, correct_bias=true)
        f = rm_ci.tau ./ sqrt(3.0)
        @test isapprox(rt_ci.dev,      rm_ci.dev      .* f; rtol=1e-14)
        @test isapprox(rt_ci.ci_lower, rm_ci.ci_lower .* f; rtol=1e-14)
        @test isapprox(rt_ci.ci_upper, rm_ci.ci_upper .* f; rtol=1e-14)
        @test rt_ci.edf        == rm_ci.edf
        @test rt_ci.noise_type == rm_ci.noise_type

        # FrequencyData entry routes via _freq_to_phase (cumsum). Reconstructed
        # phase has one fewer sample than the original x, so values differ
        # from the direct PhaseData call at the percent level — assert the
        # entry point runs and stays in the right ballpark, not bit-identity.
        fd = FrequencyData(diff(x) ./ tau0, tau0)
        rt_fd = ttotdev(fd, ms; calc_ci=false, correct_bias=false)
        @test length(rt_fd.dev) == length(ms)
        @test all(isfinite, rt_fd.dev)
        @test all(>(0), rt_fd.dev)
        @test isapprox(rt_fd.dev, rt.dev; rtol=5e-3)
    end

    @testset "HTOTDEV :linear smoke" begin
        # Per-window full-LS detrend on the frequency series (vs :greenhall's
        # half-mean) + same time-reverse extension and third-diff operator.
        using Random
        Random.seed!(20260508)
        N    = 1024
        tau0 = 1.0
        ms   = [1, 2, 4, 8, 16]
        x = _gen_powerlaw_phase(0.0, N; tau0=tau0)

        devs = SigmaTau.Stab._htotdev_core(x, ms, tau0; detrend=:linear)
        @test length(devs) == length(ms)
        @test all(isfinite, devs)
        @test all(>(0), devs)
        devs_legacy = SigmaTau.Stab._htotdev_core(x, ms, tau0; detrend=:legacy)
        @test all(0.1 .<= devs ./ devs_legacy .<= 10.0)

        # :howe is no longer a recipe for HTOTDEV.
        @test_throws ArgumentError SigmaTau.Stab._htotdev_core(x, [1], tau0; detrend=:howe)
    end

    @testset "MHTOTDEV :greenhall smoke" begin
        # Per-window half-mean slope detrend (vs :linear's full-LS) on the
        # same time-reverse extension and averaged third-diff operator.
        # MHTOTDEV is novel to SigmaTau; this is the new default after
        # the Phase 4 default switch — exercise it on a mid-spectrum noise
        # plus all five SP1065 power-law types.
        using Random
        N    = 1024
        tau0 = 1.0
        ms   = [1, 2, 4, 8]

        # Mid-spectrum WHFM check
        Random.seed!(20260508)
        x = _gen_powerlaw_phase(0.0, N; tau0=tau0)
        devs = SigmaTau.Stab._mhtotdev_core(x, ms, tau0; detrend=:greenhall)
        @test length(devs) == length(ms)
        @test all(isfinite, devs)
        @test all(>(0), devs)

        # 5-noise-type finite-output smoke (the new default needs basic coverage)
        for alpha in (2.0, 1.0, 0.0, -1.0, -2.0)
            Random.seed!(20260508)
            xa = _gen_powerlaw_phase(alpha, N; tau0=tau0)
            d = SigmaTau.Stab._mhtotdev_core(xa, ms, tau0; detrend=:greenhall)
            @test length(d) == length(ms)
            @test all(isfinite, d)
            @test all(>(0), d)
        end

        # :howe is no longer a recipe for MHTOTDEV.
        @test_throws ArgumentError SigmaTau.Stab._mhtotdev_core(x, [1], tau0; detrend=:howe)
    end

    @testset "Cross-recipe equivalence (:legacy aliases)" begin
        # Each kernel exposes :legacy as an alias for the recipe matching its
        # pre-1.0 default. The dispatcher routes both to the same helper, so
        # the outputs agree at machine precision. Lock the alias on a
        # WPM+RWFM mix at rtol=1e-15 so silent drift in either branch breaks
        # the test immediately.
        using Random
        Random.seed!(20260507)
        N    = 4096
        tau0 = 1.0
        wpm  = randn(N) .* 1e-9
        rwfm = cumsum(cumsum(randn(N) .* 1e-12))
        x    = wpm .+ rwfm

        for m in [1, 2, 4, 8, 16]
            @test SigmaTau.Stab._mtotdev_core(x, [m], tau0; detrend=:legacy)[1]    ≈
                  SigmaTau.Stab._mtotdev_core(x, [m], tau0; detrend=:greenhall)[1] atol=0.0 rtol=1e-15
            @test SigmaTau.Stab._htotdev_core(x, [m], tau0; detrend=:legacy)[1]    ≈
                  SigmaTau.Stab._htotdev_core(x, [m], tau0; detrend=:greenhall)[1] atol=0.0 rtol=1e-15
        end

        for m in [1, 2, 4, 8]
            @test SigmaTau.Stab._mhtotdev_core(x, [m], tau0; detrend=:legacy)[1] ≈
                  SigmaTau.Stab._mhtotdev_core(x, [m], tau0; detrend=:linear)[1] atol=0.0 rtol=1e-15
        end
    end

    @testset "ADEV/HDEV across all 5 power-law noise types" begin
        # Bonus: kernel parity for the more common ADEV/HDEV across all 5
        # noise types, locking in that the synthesizer + kernels survive the
        # full SP1065 alpha range.
        #
        # rtol target is 1e-11 (was 1e-12). On macOS x86_64 our `_mdev_core`
        # and `LK.mdev_var` agree bit-exactly (Δ_ol = 0 ULP for nearly all
        # rows), and both agree with allantools to ≤ 8.5e-14 worst-case on
        # the stable32gen.DAT fixture (3-way verification, 2026-05-08).
        # Linux x86_64 LLVM picks a different SIMD reduction order and the
        # two implementations drift by ~10,000 ULPs on this synthesised
        # input — irreducible cross-platform codegen variance, not a math
        # bug. 1e-11 still asserts 11-sig-fig parity, well above any drift
        # we have observed.
        using Random
        N    = 1024
        tau0 = 1.0
        ms   = [1, 2, 4, 8, 16]
        rt   = 1e-11

        for alpha in (2.0, 1.0, 0.0, -1.0, -2.0)
            Random.seed!(123)
            x    = _gen_powerlaw_phase(alpha, N; tau0=tau0)
            x_cs = pushfirst!(cumsum(x), 0.0)
            for m in ms
                @test SigmaTau.Stab._adev_core(x, [m], tau0)[1]  ≈
                      sqrt(LK.adev_var(x, m, tau0))                 atol=1e-25 rtol=rt
                @test SigmaTau.Stab._mdev_core(x, [m], tau0)[1]  ≈
                      sqrt(LK.mdev_var(x, m, tau0, x_cs))            atol=1e-25 rtol=rt
                @test SigmaTau.Stab._hdev_core(x, [m], tau0)[1]  ≈
                      sqrt(LK.hdev_var(x, m, tau0))                  atol=1e-25 rtol=rt
                @test SigmaTau.Stab._mhdev_core(x, [m], tau0)[1] ≈
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
        edfs = SigmaTau.Stab.calculate_edf(:totdev, ones(length(m_values_eq)),
                                               noises_wpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs)
        @test all(>(0.0), edfs)

        # α=1 (FLPM)
        noises_flpm = fill(:FLPM, length(m_values_eq))
        edfs = SigmaTau.Stab.calculate_edf(:totdev, ones(length(m_values_eq)),
                                               noises_flpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs)

        # htotdev gets the HDEV-style fallback at α=2,1 too.
        edfs_h = SigmaTau.Stab.calculate_edf(:htotdev, ones(length(m_values_eq)),
                                                 noises_wpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs_h)
    end

    @testset "_unbias_divisor maps non-positive B to NaN" begin
        # Regression: TOTDEV's variance-scale bias factor B = 1 − a·τ/T can
        # go non-positive when a caller passes an oversized τ (a = 0.75 for
        # RWFM means τ/T > 4/3 trips it). Before the guard, dividing the
        # deviation by √B threw DomainError and aborted the whole call;
        # the helper now returns NaN for those rows so other τ values still
        # produce valid results.
        f = SigmaTau.Stab._unbias_divisor
        out = f([1.0, 0.5, 0.0, -0.1, 1e-300])
        @test out[1] == 1.0
        @test out[2] ≈ sqrt(0.5)
        @test isnan(out[3])           # B = 0 → NaN, not Inf
        @test isnan(out[4])           # B < 0 → NaN, not DomainError
        @test isfinite(out[5])        # tiny positive B still works
    end

    @testset "MTIE core" begin
        # Hand-checkable fixture: x = [0, 1, 0.5, 2, 1.5]. Window size = m+1.
        #   m=1: peak-to-peak per window ∈ {1.0, 0.5, 1.5, 0.5}     → 1.5
        #   m=2: peak-to-peak per window ∈ {1.0, 1.5, 1.5}          → 1.5
        #   m=3: peak-to-peak per window ∈ {2.0, 1.5}               → 2.0
        #   m=4: peak-to-peak per window ∈ {2.0}                    → 2.0
        x = [0.0, 1.0, 0.5, 2.0, 1.5]
        devs = SigmaTau.Stab._mtie_core(x, [1, 2, 3, 4], 1.0)
        @test devs ≈ [1.5, 1.5, 2.0, 2.0]

        # Monotonic ramp: max - min over any window of m+1 contiguous samples
        # is exactly m (the spacing).
        ramp = collect(0.0:99.0)
        m_grid = [1, 2, 5, 10, 50]
        ramp_devs = SigmaTau.Stab._mtie_core(ramp, m_grid, 1.0)
        @test ramp_devs ≈ Float64.(m_grid)

        # Constant phase → MTIE = 0 at every τ.
        const_devs = SigmaTau.Stab._mtie_core(fill(3.14, 64), [1, 2, 4, 8], 1.0)
        @test all(==(0.0), const_devs)

        # NaN guard for windows wider than the record.
        @test isnan(SigmaTau.Stab._mtie_core(ramp, [200], 1.0)[1])

        # Naive double-loop reference parity on a noisy fixture.
        using Random
        Random.seed!(20260509)
        N = 256
        xn = cumsum(randn(N))
        ref = zeros(length(m_grid))
        for (i, m) in enumerate(m_grid)
            best = 0.0
            for j in 1:(N - m)
                w = view(xn, j:(j + m))
                d = maximum(w) - minimum(w)
                d > best && (best = d)
            end
            ref[i] = best
        end
        got = SigmaTau.Stab._mtie_core(xn, m_grid, 1.0)
        @test got ≈ ref atol=0.0 rtol=1e-15

        # API wrapper returns a StabilityResult with empty CI fields.
        pd = PhaseData(xn, 1.0)
        res = mtie(pd, m_grid)
        @test res.deviation_type == :mtie
        @test res.dev ≈ ref
        @test isempty(res.noise_type)
        @test isempty(res.ci_lower)
        @test isempty(res.ci_upper)
        @test isempty(res.edf)

        # FrequencyData entry point: cumsum-equivalence smoke check.
        Random.seed!(20260509)
        y = randn(100) .* 1e-9
        fd = FrequencyData(y, 1.0)
        pd_eq = PhaseData(cumsum(y) .* 1.0, 1.0)
        @test mtie(fd, [1, 2, 4]).dev ≈ mtie(pd_eq, [1, 2, 4]).dev
    end

    @testset "PDEV core" begin
        # Identity at m=1: PDEV(τ₀) ≡ overlapping ADEV(τ₀) (Vernotte 2015).
        using Random
        Random.seed!(20260509)
        x = cumsum(randn(512))
        @test SigmaTau.Stab._pdev_core(x, [1], 1.0)[1] ≈
              SigmaTau.Stab._adev_core(x, [1], 1.0)[1] atol=0.0 rtol=1e-15

        # Constant phase → PDEV = 0.
        const_devs = SigmaTau.Stab._pdev_core(fill(2.718, 128), [1, 2, 4, 8], 1.0)
        @test all(==(0.0), const_devs)

        # Linear phase x[i] = a + b·i: parabolic weights ((m-1)/2 - k) sum to
        # zero, so the operator kills any linear trend → PDEV = 0.
        lin = collect(1.0:128.0) .* 0.5 .+ 7.0
        @test all(d -> d < 1e-12, SigmaTau.Stab._pdev_core(lin, [2, 4, 8, 16], 1.0))

        # NaN guard when M = N - 2m < 1.
        @test isnan(SigmaTau.Stab._pdev_core(collect(1.0:10.0), [6], 1.0)[1])

        # Reference implementation matching the allantools formula verbatim
        # (Vernotte 2020): asum = Σ_{k=0}^{m-1} ((m-1)/2 - k)·(x[i+k] - x[i+k+m]),
        # σ² = 72·ΣMᵢ² / ((N-2m)·m⁴·(m·τ₀)²).
        function _pdev_reference(x, m_values, tau0)
            N = length(x)
            out = Vector{Float64}(undef, length(m_values))
            for (idx, m) in enumerate(m_values)
                if m == 1
                    out[idx] = SigmaTau.Stab._adev_core(x, [1], tau0)[1]
                    continue
                end
                M = N - 2m
                if M < 1
                    out[idx] = NaN
                    continue
                end
                Msum = 0.0
                for i in 1:M
                    asum = 0.0
                    for k in 0:(m - 1)
                        asum += ((m - 1) / 2 - k) * (x[i + k] - x[i + k + m])
                    end
                    Msum += asum^2
                end
                out[idx] = sqrt(72 * Msum / (M * m^4 * (m * tau0)^2))
            end
            return out
        end

        Random.seed!(20260509)
        x_noise = cumsum(randn(1024))
        ms = [1, 2, 4, 8, 16, 32]
        ref = _pdev_reference(x_noise, ms, 1.0)
        got = SigmaTau.Stab._pdev_core(x_noise, ms, 1.0)
        @test got ≈ ref atol=1e-25 rtol=1e-12

        # API wrapper: empty CI fields, FrequencyData dispatch.
        pd = PhaseData(x_noise, 1.0)
        res = pdev(pd, ms)
        @test res.deviation_type == :pdev
        @test res.dev ≈ ref
        @test isempty(res.ci_lower) && isempty(res.ci_upper) && isempty(res.edf)

        Random.seed!(20260509)
        y = randn(200) .* 1e-9
        fd = FrequencyData(y, 1.0)
        pd_eq = PhaseData(cumsum(y) .* 1.0, 1.0)
        @test pdev(fd, [1, 2, 4]).dev ≈ pdev(pd_eq, [1, 2, 4]).dev
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

    @testset "noise-ID scale invariance" begin
        # Regression for the eps*N threshold bug in `_lag1_acf` that
        # produced :unknown classifications on any phase record with std
        # below ~√eps ≈ 1.5e-8 (i.e. all real-world records, in seconds).
        # Same WPM+RWFM fixture used by examples/02_compute_adev.jl: every
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

    @testset "noise-ID Stable32-fixture cross-check vs allantools" begin
        # Reference table generated with allantools 2024.06's
        # `autocorr_noise_id(x, m, data_type='phase', dmin=0, dmax=2)`
        # on `reference/validation/stable32gen.DAT` (the same Stable32
        # fixture the deviation cross-validation testset above uses).
        # Allantools errors out at m=512 (time-series too short after
        # differencing) — SigmaTau's B1/R(n) fallback handles it, so
        # we cross-check only m where allantools can compute.
        ref_dir = joinpath(@__DIR__, "..", "..", "reference", "validation")
        dat_path = joinpath(ref_dir, "stable32gen.DAT")
        if !isfile(dat_path)
            @info "Stable32 fixture not present; skipping noise-ID cross-check"
        else
            lines = readlines(dat_path)
            x = parse.(Float64, strip.(lines[11:end]))
            @test length(x) == 8192

            # m → expected α_int from allantools (reference printout).
            allantools_ref = Dict(
                1   => 2,
                2   => 2,
                # m=4: borderline between α=2 and α=1 — allantools picks 1
                # via a slightly different detrend numerator / dmin policy.
                # Skipped from the strict cross-check; recorded as a
                # documented one-row drift.
                8   => 2,
                16  => 2,
                32  => 1,
                64  => 1,
                128 => 0,
                256 => 0,
            )

            ms = sort(collect(keys(allantools_ref)))
            noises = identify_noise(x, ms; dmin=0, dmax=2)
            for (i, m) in enumerate(ms)
                expected_alpha = allantools_ref[m]
                got_noise = noises[i]
                @test got_noise != :unknown
                got_alpha = if got_noise == :WHPM; 2
                elseif got_noise == :FLPM;          1
                elseif got_noise == :WHFM;          0
                elseif got_noise == :FLFM;        -1
                elseif got_noise == :RWFM;        -2
                else;                            -99
                end
                @test got_alpha == expected_alpha
            end
        end
    end

    @testset "B1(N, μ) closed-form parity vs allantools canonical" begin
        # Locks the closed-form B1 ratios used by the B1/R(n) noise-ID
        # fallback (`_noise_id_b1rn`) against the canonical allantools
        # `b1_theory(N, mu)` reference (Wallin 2018, citing Howe 2000).
        # Drift in any constant here changes which noise type the
        # B1 fallback picks at every m where lag-1 ACF can't run, so
        # the table is checked at machine precision.
        import SigmaTau.Stab as S

        # Canonical reference: allantools/ci.py b1_theory(N, mu).
        # mu ∈ {1, 0, -1, -2} are the values actually consumed by
        # _noise_id_b1rn (mu=2 is dead-code defensive, exercised below
        # for completeness).
        for N in (50, 100, 1000, 10_000)
            @test S._b1_theory(N, 2)  ≈ N * (N + 1) / 6                              rtol=0.0
            @test S._b1_theory(N, 1)  ≈ N / 2                                        rtol=0.0
            @test S._b1_theory(N, 0)  ≈ N * log(N) / (2 * (N - 1) * log(2))           rtol=0.0
            @test S._b1_theory(N, -1) == 1.0
            @test S._b1_theory(N, -2) ≈ (N^2 - 1) / (1.5 * N * (N - 1))              rtol=0.0
        end

        # Generic Howe formula B1 = N(1 - N^μ) / (2(N-1)(1 - 2^μ))
        # — the closed forms above are special simplifications. Verify
        # the simplifications agree with the generic form at a sanity
        # μ value where both branches return finite values.
        for N in (100, 1000)
            generic = N * (1 - Float64(N)^(-2)) / (2 * (N - 1) * (1 - 2.0^(-2)))
            @test S._b1_theory(N, -2) ≈ generic atol=0.0 rtol=1e-14
        end
    end

    @testset "R(n)(af, b) closed-form parity vs allantools canonical" begin
        # Locks the R(n) = MVAR/AVAR closed-form used by the WPM/FLPM
        # disambiguation step in `_noise_id_b1rn`. b=0 is the
        # WHFM-shape (also covers WPM via slope ratio); b=-1 is the
        # FLFM-shape (also covers FLPM). All other b values fall
        # through to the af⁰ = 1 branch.
        import SigmaTau.Stab as S

        # b=0 branch: R(n) = 1/af for every af ≥ 1.
        for af in (1, 2, 4, 16, 64, 256, 1024)
            @test S._rn_theory(af, 0) ≈ 1 / af atol=0.0 rtol=0.0
        end

        # b=-1 branch: explicit Howe/Greenhall formula for the
        # FLFM/FLPM shape. avar = (1.038 + 3·log(2π·0.5·af)) / (4π²);
        # mvar = 3·log(256/27) / (8π²); R(n) = mvar/avar.
        for af in (1, 2, 4, 16, 64, 256, 1024)
            avar = (1.038 + 3 * log(2π * 0.5 * af)) / (4π^2)
            mvar = 3 * log(256 / 27) / (8π^2)
            @test S._rn_theory(af, -1) ≈ mvar / avar atol=0.0 rtol=1e-15
        end

        # All other b values (including b=2, 1, -2) are not in the
        # closed-form table and fall through to the constant-1 branch.
        for af in (1, 4, 64)
            @test S._rn_theory(af, 1)  == 1.0
            @test S._rn_theory(af, 2)  == 1.0
            @test S._rn_theory(af, -2) == 1.0
        end
    end

    @testset "StabilityResult I/O round-trip" begin
        import Random
        tmpdir = mktempdir()
        path   = joinpath(tmpdir, "test_result.tsv")

        Random.seed!(20260510)
        pd = PhaseData(randn(128) .* 1e-9, 1.0)
        ms = [1, 2, 4, 8]

        # Round-trip with CI
        r  = adev(pd, ms; calc_ci=true)
        save_result(path, r)
        r2 = load_result(path)

        @test r2.deviation_type === r.deviation_type
        @test r2.tau       ≈ r.tau
        @test r2.dev       ≈ r.dev
        @test r2.noise_type == r.noise_type
        @test r2.ci_lower  ≈ r.ci_lower
        @test r2.ci_upper  ≈ r.ci_upper
        @test r2.edf       ≈ r.edf

        # Round-trip without CI — empty vectors must survive the cycle
        r_nci = adev(pd, ms; calc_ci=false)
        save_result(path, r_nci)
        r3 = load_result(path)

        @test r3.deviation_type === :adev
        @test r3.tau ≈ r_nci.tau
        @test r3.dev ≈ r_nci.dev
        @test isempty(r3.noise_type)
        @test isempty(r3.ci_lower)
        @test isempty(r3.ci_upper)
        @test isempty(r3.edf)

        # save_result returns the path
        @test save_result(path, r) == path
    end

    @testset "Zero-arg convenience: default octave m-grid per deviation" begin
        # Octave grid shape: 1, 2, 4, ..., 2^floor(log2(m_max)).
        # m_max per kernel matches each `_*_core`'s L-check.
        @testset "_default_m_values per kernel" begin
            N = 1024
            # ADEV / TOTDEV / PDEV: m_max = (N-1)÷2 = 511 → 2^0..2^8 = [1..256]
            @test SigmaTau.Stab._default_m_values(N, :adev)   == 2 .^ (0:8)
            @test SigmaTau.Stab._default_m_values(N, :totdev) == 2 .^ (0:8)
            @test SigmaTau.Stab._default_m_values(N, :pdev)   == 2 .^ (0:8)
            # MDEV/TDEV/MTOTDEV/TTOTDEV/HTOTDEV: m_max = N÷3 = 341 → 2^0..2^8
            @test SigmaTau.Stab._default_m_values(N, :mdev)    == 2 .^ (0:8)
            @test SigmaTau.Stab._default_m_values(N, :tdev)    == 2 .^ (0:8)
            @test SigmaTau.Stab._default_m_values(N, :mtotdev) == 2 .^ (0:8)
            @test SigmaTau.Stab._default_m_values(N, :ttotdev) == 2 .^ (0:8)
            @test SigmaTau.Stab._default_m_values(N, :htotdev) == 2 .^ (0:8)
            # HDEV: m_max = (N-1)÷3 = 341 → 2^0..2^8
            @test SigmaTau.Stab._default_m_values(N, :hdev) == 2 .^ (0:8)
            # MHDEV/HTDEV/MHTOTDEV: m_max = N÷4 = 256 → 2^0..2^8
            @test SigmaTau.Stab._default_m_values(N, :mhdev)    == 2 .^ (0:8)
            @test SigmaTau.Stab._default_m_values(N, :htdev)    == 2 .^ (0:8)
            @test SigmaTau.Stab._default_m_values(N, :mhtotdev) == 2 .^ (0:8)
            # MTIE: m_max = N-1 = 1023 → 2^0..2^9 = [1..512]
            @test SigmaTau.Stab._default_m_values(N, :mtie) == 2 .^ (0:9)
        end

        @testset "_default_m_values argument validation" begin
            @test_throws ArgumentError SigmaTau.Stab._default_m_values(1024, :nope)
            # N too small to admit any m ≥ 1
            @test_throws ArgumentError SigmaTau.Stab._default_m_values(1, :adev)
            @test_throws ArgumentError SigmaTau.Stab._default_m_values(3, :mhdev)
        end

        @testset "zero-arg API matches explicit-m_values dispatch" begin
            # Synthesize a deterministic WFM phase fixture and verify
            # `dev(pd)` ≡ `dev(pd, _default_m_values(N, :dev))` across both
            # PhaseData and FrequencyData entry points.
            Random.seed!(2026)
            N    = 512
            tau0 = 1.0
            x    = _gen_powerlaw_phase(0.0, N; tau0=tau0)
            p    = PhaseData(x, tau0)
            f    = FrequencyData(diff(x) ./ tau0, tau0)

            for (kernel, fn) in (
                    (:adev,    adev),  (:mdev,  mdev),  (:tdev,  tdev),
                    (:hdev,    hdev),  (:mhdev, mhdev), (:htdev, htdev),
                    (:totdev,  totdev),
                    (:mtotdev, mtotdev), (:htotdev, htotdev),
                    (:mtie,    mtie),  (:pdev,  pdev),
                )
                ms = SigmaTau.Stab._default_m_values(N, kernel)
                r_default  = fn(p; calc_ci = false)
                r_explicit = fn(p, ms; calc_ci = false)
                @test r_default.tau == r_explicit.tau
                @test r_default.dev == r_explicit.dev
                # FrequencyData entry point resolves to its own default m-grid
                # (computed from length(f.y) = N-1, not N).
                ms_f = SigmaTau.Stab._default_m_values(length(f.y), kernel)
                r_default_f = fn(f; calc_ci = false)
                @test r_default_f.tau == ms_f .* tau0
            end
        end

        @testset "kwargs pass through (calc_ci, confidence, detrend)" begin
            Random.seed!(7)
            p = PhaseData(_gen_powerlaw_phase(0.0, 512; tau0=1.0), 1.0)
            # calc_ci=true populates CI; calc_ci=false leaves them empty.
            r_ci = adev(p; calc_ci = true)
            @test !isempty(r_ci.ci_lower)
            @test !isempty(r_ci.edf)
            r_no = adev(p; calc_ci = false)
            @test isempty(r_no.ci_lower)
            @test isempty(r_no.edf)
            # totdev's `detrend` kwarg passes through.
            r_howe   = totdev(p; calc_ci = false, detrend = :howe)
            r_linear = totdev(p; calc_ci = false, detrend = :linear)
            @test r_howe.tau == r_linear.tau
        end
    end

    @testset "_phase_to_freq / _freq_to_phase helpers" begin
        # Canonical mapping: y[k] = (x[k+1] − x[k]) / τ₀, length N → N−1.
        @testset "definition and shape" begin
            τ₀ = 0.5
            x  = [0.0, 1.0, 3.0, 7.0, 15.0]
            pd = PhaseData(x, τ₀)
            fd = SigmaTau.Stab._phase_to_freq(pd)
            @test fd isa FrequencyData
            @test fd.tau0 == τ₀
            @test fd.y == diff(x) ./ τ₀
            @test length(fd.y) == length(x) - 1
        end

        @testset "round-trip identities (modulo the lost offset/sample)" begin
            Random.seed!(31)
            τ₀ = 1.0
            # Frequency → phase → frequency drops the first sample.
            y  = randn(64)
            fd = FrequencyData(y, τ₀)
            fd2 = SigmaTau.Stab._phase_to_freq(SigmaTau.Stab._freq_to_phase(fd))
            @test length(fd2.y) == length(y) - 1
            @test fd2.y ≈ y[2:end]
            # Phase → frequency → phase recovers x[2:end] − x[1].
            x   = cumsum(randn(64))
            pd  = PhaseData(x, τ₀)
            pd2 = SigmaTau.Stab._freq_to_phase(SigmaTau.Stab._phase_to_freq(pd))
            @test length(pd2.x) == length(x) - 1
            @test pd2.x ≈ x[2:end] .- x[1]
        end

        @testset "deviation equivalence: adev/mdev/hdev agree on either domain" begin
            # Build a synthetic WFM phase record, compute ADEV/MDEV/HDEV via
            # the phase path and via the frequency path (after _phase_to_freq).
            # Second- and third-difference kernels are shift-invariant so the
            # lost initial sample is the only source of disagreement; with
            # N=4096 the σ values agree to well within rtol=5e-3.
            Random.seed!(2026)
            τ₀ = 1.0
            p_full = PhaseData(_gen_powerlaw_phase(0.0, 4096; tau0=τ₀), τ₀)
            f_diff = SigmaTau.Stab._phase_to_freq(p_full)
            @test length(f_diff.y) == length(p_full.x) - 1
            ms = [1, 4, 16, 64]
            for fn in (adev, mdev, hdev)
                r_phase = fn(p_full,  ms; calc_ci = false)
                r_freq  = fn(f_diff,  ms; calc_ci = false)
                @test r_freq.tau == r_phase.tau
                @test isapprox(r_freq.dev, r_phase.dev; rtol = 5e-3)
            end
        end
    end
end

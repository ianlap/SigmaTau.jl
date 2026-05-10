using Test
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
end

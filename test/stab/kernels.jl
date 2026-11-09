# test/stab/kernels.jl — the _*_core kernels against the reference kernels,
# window counts, MTIE / PDEV cores and the typed wrappers.

using Test
using Random
using SigmaTau

@testset "Deviation kernels" begin
    # Mock data equivalent to NIST SP1065 sample
    x = collect(1.0:100.0) .^ 2  # simple quadratic
    m_values = [1, 2, 4]
    tau0 = 1.0

    @testset "kernels on x = n²" begin
        # Second differences are constant 2, so ADEV(m=1) = √2; third
        # differences vanish, so the Hadamard family is zero.
        @test _adev_core(x, m_values, tau0)[1] ≈ sqrt(2.0)
        cores = (_adev_core, _mdev_core, _tdev_core, _hdev_core, _mhdev_core,
                 _totdev_core, _mtotdev_core, _htotdev_core, _mhtotdev_core)
        @test all(all(isfinite, c(x, m_values, tau0)) for c in cores)
        @test all(all(<(1e-10), c(x, m_values, tau0)) for c in (_hdev_core, _mhdev_core))
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
        @test all(isfinite, ci_lower(res_adev))
        @test all(isfinite, ci_upper(res_adev))
        # Tuple access reads naturally and agrees with the accessors.
        @test res_adev.ci[1].lo == ci_lower(res_adev)[1]
        @test res_adev.ci[1].hi == ci_upper(res_adev)[1]

        res_hdev = hdev(pd, m_values)
        @test length(res_hdev.dev) == 3
        @test all(isfinite, ci_upper(res_hdev))

        res_tot = totdev(pd, m_values)
        @test length(res_tot.dev) == 3
        @test all(isfinite, ci_upper(res_tot))
        
        res_htdev = htdev(pd, m_values)
        @test length(res_htdev.dev) == 3
        @test res_htdev.deviation_type == :htdev

        # htdev wraps mhdev with the τ/√(10/3) scaling — point estimates AND
        # CI bounds inherit the same multiplicative factor. For a non-degenerate
        # Y = c · X with c constant in τ, the χ²-based CI scales linearly in
        # c, so [Y.lo, Y.hi] = c · [X.lo, X.hi] is the correct operation.
        res_mhdev_ref = mhdev(pd, m_values)
        h_factor = res_htdev.tau ./ sqrt(10.0 / 3.0)
        @test res_htdev.dev      ≈ res_mhdev_ref.dev      .* h_factor
        @test ci_lower(res_htdev) ≈ ci_lower(res_mhdev_ref) .* h_factor
        @test ci_upper(res_htdev) ≈ ci_upper(res_mhdev_ref) .* h_factor
        # χ² invariant: the multiplicative scaling cancels in the ratio
        # CI_bound / dev, so HTDEV's relative CI structure is identical to
        # MHDEV's at every τ.
        @test ci_lower(res_htdev) ./ res_htdev.dev ≈ ci_lower(res_mhdev_ref) ./ res_mhdev_ref.dev
        @test ci_upper(res_htdev) ./ res_htdev.dev ≈ ci_upper(res_mhdev_ref) ./ res_mhdev_ref.dev
        # htdev inherits mhdev's window count.
        @test res_htdev.neff == res_mhdev_ref.neff

        # tdev wraps mdev with the τ/√3 scaling — point estimates and CI bounds
        # must be consistent with that identity.
        res_tdev = tdev(pd, m_values)
        res_mdev_ref = mdev(pd, m_values)
        factor = res_tdev.tau ./ sqrt(3.0)
        @test length(res_tdev.dev) == 3
        @test res_tdev.deviation_type == :tdev
        @test res_tdev.dev ≈ res_mdev_ref.dev .* factor
        @test ci_lower(res_tdev) ≈ ci_lower(res_mdev_ref) .* factor
        @test ci_upper(res_tdev) ≈ ci_upper(res_mdev_ref) .* factor
        @test res_tdev.neff == res_mdev_ref.neff

        # ci=false path returns empty CI and EDF vectors; neff stays populated.
        res_tdev_nci = tdev(pd, m_values; ci=false)
        @test isempty(res_tdev_nci.ci) && isempty(ci_lower(res_tdev_nci))
        @test isempty(res_tdev_nci.edf)
        @test res_tdev_nci.neff == res_tdev.neff

        # edf is populated when ci=true.
        @test length(res_adev.edf) == length(m_values)
        @test all(isfinite, res_adev.edf)
        @test all(>=(0.0), res_adev.edf)
    end

    @testset "neff window counts mirror the kernel loop bounds" begin
        using Random
        Random.seed!(20260610)
        N  = 256
        pd = PhaseData(cumsum(randn(N)) .* 1e-9, 1.0)
        ms = [1, 2, 4, 8, 16]

        # adev sums N−2m second differences (_adev_core: i in 1:N−2m).
        @test adev(pd, ms; ci=false).neff == N .- 2 .* ms
        # mhdev sums N−4m+1 fourth differences (_mhdev_core: i in 1:N−4m+1).
        @test mhdev(pd, ms; ci=false).neff == N .- 4 .* ms .+ 1
        # mhtotdev: same N−4m+1 subsequence count; the wrapper's global drift
        # removal does not change the record length, so remove_drift=true and
        # =false agree.
        @test mhtotdev(pd, ms; ci=false, remove_drift=true).neff == N .- 4 .* ms .+ 1
        @test mhtotdev(pd, ms; ci=false, remove_drift=false).neff == N .- 4 .* ms .+ 1
        # htotdev: m=1 falls back to the HDEV sum over N−3 terms; m≥2 runs
        # (N−1)−3m+1 subsequences of y = diff(x).
        r_htot = htotdev(pd, ms; ci=false)
        @test r_htot.neff[1] == N - 3
        @test r_htot.neff[2:end] == (N - 1) .- 3 .* ms[2:end] .+ 1
        # Other kernels, per their loop bounds.
        @test mdev(pd, ms; ci=false).neff    == N .- 3 .* ms .+ 1
        @test hdev(pd, ms; ci=false).neff    == N .- 3 .* ms
        @test totdev(pd, ms; ci=false).neff  == fill(N - 2, length(ms))
        @test mtotdev(pd, ms; ci=false).neff == N .- 3 .* ms .+ 1
        @test pdev(pd, ms; ci=false).neff    == N .- 2 .* ms

        # neff is populated identically with and without CI.
        @test adev(pd, ms; ci=true).neff == adev(pd, ms; ci=false).neff

        # NaN grid points (m too large for N) report neff = 0.
        r_edge = adev(pd, [N ÷ 2 - 1, N ÷ 2]; ci=false)
        @test isnan(r_edge.dev[2]) && r_edge.neff[2] == 0
        @test r_edge.neff[1] == 2
    end

    @testset "Reference-kernel parity" begin
        # Strict numerical parity against the textbook reference kernels in
        # test/stab/reference_kernels.jl. Uses a seeded
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

        # kernel => (reference variance at m, m grid). TOTDEV is the canonical
        # Howe form (SP1065 eq. 25, no detrend). MHTOTDEV (per-window
        # least-squares quadratic removal, against a Vandermonde-solve
        # reference) is novel to SigmaTau with no external reference, so this
        # machine-precision anchor is its primary correctness contract.
        refs = ((_adev_core,     m -> LK.adev_var(x, m, tau0),               m_grid),
                (_mdev_core,     m -> LK.mdev_var(x, m, tau0, x_cs),         m_grid),
                (_hdev_core,     m -> LK.hdev_var(x, m, tau0),               m_grid),
                (_mhdev_core,    m -> LK.mhdev_var(x, m, tau0, x_cs),        m_grid),
                (_totdev_core,   m -> LK.totdev_var_howe(x, m, tau0),        [1, 2, 4, 8, 16, 32]),
                (_mtotdev_core,  m -> LK.mtotdev_var(x, m, tau0),            [1, 2, 4, 8, 16]),
                (_htotdev_core,  m -> LK.htotdev_var(x, m, tau0),            [1, 2, 4, 8, 16]),
                (_mhtotdev_core, m -> LK.mhtotdev_var_greenhall(x, m, tau0), [1, 2, 4, 8]))
        @test isempty([(core, m) for (core, ref, ms) in refs for m in ms
                       if !isapprox(core(x, [m], tau0)[1], sqrt(ref(m)); atol=1e-25, rtol=1e-12)])
    end

    @testset "_unbias_divisor maps non-positive B to NaN" begin
        # Regression: TOTDEV's variance-scale bias factor B = 1 − a·τ/T can
        # go non-positive when a caller passes an oversized τ (a = 0.75 for
        # RWFM means τ/T > 4/3 trips it). Before the guard, dividing the
        # deviation by √B threw DomainError and aborted the whole call;
        # the helper now returns NaN for those rows so other τ values still
        # produce valid results.
        f = SigmaTau._unbias_divisor
        out = f([1.0, 0.5, 0.0, -0.1, 1e-300])
        @test out[1] == 1.0
        @test out[2] ≈ sqrt(0.5)
        @test isnan(out[3])           # B = 0 → NaN, not Inf
        @test isnan(out[4])           # B < 0 → NaN, not DomainError
        @test isfinite(out[5])        # tiny positive B still works
    end

    @testset "Non-Float64 PhaseData flows through the API" begin
        # Regression: PhaseData{Float32} (and any AbstractFloat eltype) used
        # to hit a MethodError on the Vector{Float64}-typed core kernels. The
        # API now promotes the sample vector at the boundary, so a Float32
        # record gives the same deviations as its Float64 twin.
        Random.seed!(20260529)
        x64 = cumsum(randn(512))
        x32 = Float32.(x64)
        ms  = [1, 2, 4, 8, 16]
        pd64 = PhaseData(x64, 1.0)
        pd32 = PhaseData(x32, 1.0)
        for f in (adev, mdev, tdev, hdev, mhdev, htdev,
                  totdev, mtotdev, ttotdev, htotdev, mhtotdev,
                  mtie, pdev)
            r32 = f(pd32, ms; ci=false)
            r64 = f(pd64, ms; ci=false)
            # Float32 input carries ~1e-7 relative precision; compare loosely.
            @test isapprox(r32.dev, r64.dev; rtol=1e-4, nans=true)
        end
        # CI path also runs without error on the promoted vector.
        @test adev(pd32, ms; ci=true).edf == adev(pd64, ms; ci=true).edf
    end

    @testset "confidence_intervals floors the lower bound at zero" begin
        # Regression: the Gaussian fallback (edf < 1) is a symmetric ±half
        # interval that could push the lower limit below zero — unphysical for
        # a deviation. It is now floored at 0.
        lo, hi = SigmaTau.confidence_intervals([1.0], [0.5], [:WHPM], 4, 0.99)
        @test lo[1] == 0.0            # would have been ≈ -0.275 before the floor
        @test hi[1] > 1.0
        @test all(lo .>= 0.0)
    end

    @testset "MTIE core" begin
        # Hand-checkable fixture: x = [0, 1, 0.5, 2, 1.5]. Window size = m+1.
        #   m=1: peak-to-peak per window ∈ {1.0, 0.5, 1.5, 0.5}     → 1.5
        #   m=2: peak-to-peak per window ∈ {1.0, 1.5, 1.5}          → 1.5
        #   m=3: peak-to-peak per window ∈ {2.0, 1.5}               → 2.0
        #   m=4: peak-to-peak per window ∈ {2.0}                    → 2.0
        x = [0.0, 1.0, 0.5, 2.0, 1.5]
        devs = SigmaTau._mtie_core(x, [1, 2, 3, 4], 1.0)
        @test devs ≈ [1.5, 1.5, 2.0, 2.0]

        # Monotonic ramp: max - min over any window of m+1 contiguous samples
        # is exactly m (the spacing).
        ramp = collect(0.0:99.0)
        m_grid = [1, 2, 5, 10, 50]
        ramp_devs = SigmaTau._mtie_core(ramp, m_grid, 1.0)
        @test ramp_devs ≈ Float64.(m_grid)

        # Constant phase → MTIE = 0 at every τ.
        const_devs = SigmaTau._mtie_core(fill(3.14, 64), [1, 2, 4, 8], 1.0)
        @test all(==(0.0), const_devs)

        # NaN guard for windows wider than the record.
        @test isnan(SigmaTau._mtie_core(ramp, [200], 1.0)[1])

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
        got = SigmaTau._mtie_core(xn, m_grid, 1.0)
        @test got ≈ ref atol=0.0 rtol=1e-15

        # API wrapper returns a StabilityResult with empty CI fields; neff
        # carries the N−m window count.
        pd = PhaseData(xn, 1.0)
        res = mtie(pd, m_grid)
        @test res.deviation_type == :mtie
        @test res.dev ≈ ref
        @test isempty(res.noise_type) && isempty(res.ci) && isempty(res.edf)
        @test res.neff == length(xn) .- m_grid

        # `confidence` is accepted for signature uniformity (no-op for MTIE):
        # passing it must not error and must not populate CI fields.
        res_conf = mtie(pd, m_grid; confidence=0.9)
        @test res_conf.dev ≈ ref
        @test isempty(res_conf.edf)
    end

    @testset "PDEV core" begin
        # Identity at m=1: PDEV(τ₀) ≡ overlapping ADEV(τ₀) (Vernotte 2015).
        using Random
        Random.seed!(20260509)
        x = cumsum(randn(512))
        @test SigmaTau._pdev_core(x, [1], 1.0)[1] ≈
              SigmaTau._adev_core(x, [1], 1.0)[1] atol=0.0 rtol=1e-15

        # Constant phase → PDEV = 0.
        const_devs = SigmaTau._pdev_core(fill(2.718, 128), [1, 2, 4, 8], 1.0)
        @test all(==(0.0), const_devs)

        # Linear phase x[i] = a + b·i: parabolic weights ((m-1)/2 - k) sum to
        # zero, so the operator kills any linear trend → PDEV = 0.
        lin = collect(1.0:128.0) .* 0.5 .+ 7.0
        @test all(d -> d < 1e-12, SigmaTau._pdev_core(lin, [2, 4, 8, 16], 1.0))

        # NaN guard when M = N - 2m < 2, plus invalid m.
        @test isnan(SigmaTau._pdev_core(collect(1.0:10.0), [6], 1.0)[1])
        @test isnan(SigmaTau._pdev_core(collect(1.0:10.0), [0], 1.0)[1])

        # Reference implementation matching the allantools formula verbatim
        # (Vernotte 2020): asum = Σ_{k=0}^{m-1} ((m-1)/2 - k)·(x[i+k] - x[i+k+m]),
        # σ² = 72·ΣMᵢ² / ((N-2m)·m⁴·(m·τ₀)²).
        function _pdev_reference(x, m_values, tau0)
            N = length(x)
            out = Vector{Float64}(undef, length(m_values))
            for (idx, m) in enumerate(m_values)
                if m == 1
                    out[idx] = SigmaTau._adev_core(x, [1], tau0)[1]
                    continue
                end
                M = N - 2m
                if M < 2
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
        got = SigmaTau._pdev_core(x_noise, ms, 1.0)
        @test got ≈ ref atol=1e-25 rtol=1e-12

        # Rolling-sum regression: non-octave factors and near-limit large m still
        # match the direct weighted-sum formula while avoiding O(Nm) runtime.
        ms_wide = [3, 7, 31, 127, 255, 511]
        ref_wide = _pdev_reference(x_noise, ms_wide, 1.0)
        got_wide = SigmaTau._pdev_core(x_noise, ms_wide, 1.0)
        @test got_wide ≈ ref_wide atol=1e-25 rtol=1e-12

        Random.seed!(20260607)
        x_refresh = cumsum(randn(10_000))
        ms_refresh = [2, 5]
        @test SigmaTau._pdev_core(x_refresh, ms_refresh, 1.0) ≈
              _pdev_reference(x_refresh, ms_refresh, 1.0) atol=1e-25 rtol=1e-12

        # API wrapper: raw devs unchanged; CI populated by default (Vernotte
        # PVAR EDF model), empty when ci=false. FrequencyData dispatch below.
        pd = PhaseData(x_noise, 1.0)
        res = pdev(pd, ms)
        @test res.deviation_type == :pdev
        @test res.dev ≈ ref
        @test length(res.edf) == length(ms) && all(isfinite, res.edf)
        @test !isempty(res.noise_type)
        @test all(ci_lower(res) .<= res.dev .+ 1e-12) && all(res.dev .<= ci_upper(res) .+ 1e-12)

        res_noci = pdev(pd, ms; ci=false)
        @test res_noci.dev ≈ ref
        @test isempty(res_noci.ci) && isempty(res_noci.edf)

        Random.seed!(20260509)
        y = randn(200) .* 1e-9
        fd = FrequencyData(y, 1.0)
        pd_eq = PhaseData([0.0; cumsum(y)] .* 1.0, 1.0)
        @test pdev(fd, [1, 2, 4]).dev ≈ pdev(pd_eq, [1, 2, 4]).dev
    end

    @testset "NEFF_RELIABLE boundary" begin
        # The threshold determines whether identify_noise uses lag-1 ACF
        # (N_eff ≥ threshold) or the B1-ratio fallback. Construct two cases
        # that straddle the boundary and verify both produce a finite,
        # classified noise type.
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

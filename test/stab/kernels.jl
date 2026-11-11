# test/stab/kernels.jl — the _*_core kernels against the reference kernels,
# window counts, MTIE / PDEV cores and the typed wrappers.

using Test
using Random
using SigmaTau

@testset "Deviation kernels" begin
    x = collect(1.0:100.0) .^ 2          # quadratic phase
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

    @testset "typed wrappers and the scaled families" begin
        Random.seed!(42)
        # m=1 has N_eff = 100 ≥ NEFF_RELIABLE (lag-1 ACF); m=4 has 25 (B1 ratio).
        pd = PhaseData(x .+ randn(100), tau0)
        r = adev(pd, m_values)
        @test length(r.dev) == length(r.edf) == 3
        @test r.noise_type[1] != :unknown
        @test all(isfinite, ci_lower(r)) && all(isfinite, ci_upper(r))
        @test (r.ci[1].lo, r.ci[1].hi) == (ci_lower(r)[1], ci_upper(r)[1])
        @test all(isfinite, ci_upper(hdev(pd, m_values))) &&
              all(isfinite, ci_upper(totdev(pd, m_values)))

        # htdev = mhdev · τ/√(10/3) and tdev = mdev · τ/√3: the point
        # estimates, both χ² bounds and the window counts scale together.
        for (f, g, c) in ((htdev, mhdev, sqrt(10 / 3)), (tdev, mdev, sqrt(3.0)))
            rf = f(pd, m_values); rg = g(pd, m_values); k = rf.tau ./ c
            @test rf.deviation_type === Symbol(f)
            @test rf.dev ≈ rg.dev .* k
            @test ci_lower(rf) ≈ ci_lower(rg) .* k && ci_upper(rf) ≈ ci_upper(rg) .* k
            @test rf.neff == rg.neff
        end

        # ci=false leaves CI and EDF empty; neff stays populated.
        rn = tdev(pd, m_values; ci=false)
        @test isempty(rn.ci) && isempty(ci_lower(rn)) && isempty(rn.edf)
        @test rn.neff == tdev(pd, m_values).neff
    end

    @testset "neff mirrors the kernel loop bounds" begin
        Random.seed!(20260610)
        N  = 256
        pd = PhaseData(cumsum(randn(N)) .* 1e-9, 1.0)
        ms = [1, 2, 4, 8, 16]
        count_of(f; kw...) = f(pd, ms; ci=false, kw...).neff
        @test count_of(adev)  == N .- 2 .* ms
        @test count_of(pdev)  == N .- 2 .* ms
        @test count_of(mdev)  == N .- 3 .* ms .+ 1
        @test count_of(hdev)  == N .- 3 .* ms
        @test count_of(mhdev) == N .- 4 .* ms .+ 1
        @test count_of(totdev)  == fill(N - 2, length(ms))
        @test count_of(mtotdev) == N .- 3 .* ms .+ 1
        # mhtotdev's global drift removal keeps the record length.
        @test count_of(mhtotdev; remove_drift=true) == count_of(mhtotdev; remove_drift=false) ==
              N .- 4 .* ms .+ 1
        # htotdev: m=1 is the HDEV sum over N−3 terms; m≥2 runs on y = diff(x).
        @test count_of(htotdev) == [N - 3; (N - 1) .- 3 .* ms[2:end] .+ 1]
        @test adev(pd, ms; ci=true).neff == count_of(adev)
        # m too large for N: NaN and neff = 0.
        r_edge = adev(pd, [N ÷ 2 - 1, N ÷ 2]; ci=false)
        @test isnan(r_edge.dev[2]) && r_edge.neff == [2, 0]
        # totdev: the reflected extension reaches m = N − 1 and no further.
        r_edge = totdev(pd, [N - 1, N, 2N]; ci=false, correct_bias=false)
        @test r_edge.dev[1] ≈ sqrt(LK.totdev_var_howe(pd.x, N - 1, pd.tau0)) rtol = 1e-12
        @test all(isnan, r_edge.dev[2:3]) && r_edge.neff == [N - 2, 0, 0]
    end

    @testset "reference-kernel parity" begin
        # Machine-precision parity with test/stab/reference_kernels.jl on a
        # WPM + RWFM fixture. MHTOTDEV has no external reference, so this
        # anchor is its primary correctness contract.
        Random.seed!(20260507)
        N = 4096
        x = randn(N) .* 1e-9 .+ cumsum(cumsum(randn(N) .* 1e-12))
        x_cs = pushfirst!(cumsum(x), 0.0)
        oct(k) = [2^i for i in 0:k]
        refs = ((_adev_core,     m -> LK.adev_var(x, m, tau0),               oct(7)),
                (_mdev_core,     m -> LK.mdev_var(x, m, tau0, x_cs),         oct(7)),
                (_hdev_core,     m -> LK.hdev_var(x, m, tau0),               oct(7)),
                (_mhdev_core,    m -> LK.mhdev_var(x, m, tau0, x_cs),        oct(7)),
                (_totdev_core,   m -> LK.totdev_var_howe(x, m, tau0),        oct(5)),
                (_mtotdev_core,  m -> LK.mtotdev_var(x, m, tau0),            oct(4)),
                (_htotdev_core,  m -> LK.htotdev_var(x, m, tau0),            oct(4)),
                (_mhtotdev_core, m -> LK.mhtotdev_var_greenhall(x, m, tau0), oct(3)))
        @test isempty([(core, m) for (core, ref, ms) in refs for m in ms
                       if !isapprox(core(x, [m], tau0)[1], sqrt(ref(m)); atol=1e-25, rtol=1e-12)])
    end

    @testset "_unbias_divisor maps non-positive B to NaN" begin
        # TOTDEV's bias factor B = 1 − a·τ/T goes non-positive for an
        # oversized τ; those rows become NaN instead of a DomainError.
        out = SigmaTau._unbias_divisor([1.0, 0.5, 0.0, -0.1, 1e-300])
        @test out[1] == 1.0 && out[2] ≈ sqrt(0.5) && isfinite(out[5])
        @test isnan(out[3]) && isnan(out[4])
    end

    @testset "non-Float64 PhaseData flows through the API" begin
        Random.seed!(20260529)
        x64 = cumsum(randn(512))
        ms  = [1, 2, 4, 8, 16]
        pd64 = PhaseData(x64, 1.0); pd32 = PhaseData(Float32.(x64), 1.0)
        for f in (adev, mdev, tdev, hdev, mhdev, htdev, totdev, mtotdev, ttotdev,
                  htotdev, mhtotdev, mtie, pdev)
            @test isapprox(f(pd32, ms; ci=false).dev, f(pd64, ms; ci=false).dev; rtol=1e-4, nans=true)
        end
        @test adev(pd32, ms; ci=true).edf == adev(pd64, ms; ci=true).edf
    end

    @testset "confidence_intervals: χ² for any ν > 0, none otherwise" begin
        # ν = 1 at 68.3 %: the χ² quantiles are the squared normal quantiles at
        # 0.92075 and 0.57925, 1.9885 and 0.039990.
        lo, hi = confidence_intervals([1.0, 1.0, 1.0, 1.0], [1.0, 0.9, 0.0, -0.3], 0.683)
        @test lo[1] ≈ 1 / sqrt(1.9885) rtol=1e-4
        @test hi[1] ≈ 1 / sqrt(0.039990) rtol=1e-4
        @test lo[2] < lo[1] && hi[2] > hi[1]          # fewer degrees of freedom, wider
        @test all(isnan, lo[3:4]) && all(isnan, hi[3:4])
        for c in (-0.1, 0.0, 1.0, 1.5, NaN)    # a confidence is a probability, between 0 and 1
            @test_throws "between 0 and 1" confidence_intervals([1.0], [10.0], c)
            @test_throws "between 0 and 1" adev(PhaseData(randn(64), 1.0), [1]; confidence=c)
        end
    end

    @testset "MTIE core" begin
        # Window size m+1: x = [0, 1, 0.5, 2, 1.5] gives 1.5, 1.5, 2, 2.
        @test SigmaTau._mtie_core([0.0, 1.0, 0.5, 2.0, 1.5], [1, 2, 3, 4], 1.0) ≈ [1.5, 1.5, 2.0, 2.0]
        ramp = collect(0.0:99.0)
        m_grid = [1, 2, 5, 10, 50]
        @test SigmaTau._mtie_core(ramp, m_grid, 1.0) ≈ Float64.(m_grid)
        @test all(==(0.0), SigmaTau._mtie_core(fill(3.14, 64), [1, 2, 4, 8], 1.0))
        @test isnan(SigmaTau._mtie_core(ramp, [200], 1.0)[1])

        # Sliding-window parity with the direct max − min over every window.
        Random.seed!(20260509)
        N = 256
        xn = cumsum(randn(N))
        ref = [maximum(maximum(view(xn, j:j+m)) - minimum(view(xn, j:j+m)) for j in 1:N-m) for m in m_grid]
        @test SigmaTau._mtie_core(xn, m_grid, 1.0) ≈ ref atol=0.0 rtol=1e-15

        # Wrapper: no CI fields, neff = N − m; `confidence` is accepted and ignored.
        res = mtie(PhaseData(xn, 1.0), m_grid; confidence=0.9)
        @test res.deviation_type == :mtie && res.dev ≈ ref
        @test isempty(res.noise_type) && isempty(res.ci) && isempty(res.edf)
        @test res.neff == N .- m_grid
    end

    @testset "PDEV core" begin
        # Direct weighted-sum formula (Vernotte 2015 / allantools):
        # Mᵢ = Σₖ ((m−1)/2 − k)(x[i+k] − x[i+k+m]), σ² = 72 ΣMᵢ² / ((N−2m) m⁴ τ²).
        function pdev_reference(x, ms, tau0)
            N = length(x)
            map(ms) do m
                m == 1 && return SigmaTau._adev_core(x, [1], tau0)[1]
                M = N - 2m
                M < 2 && return NaN
                s = sum(abs2(sum(((m - 1) / 2 - k) * (x[i+k] - x[i+k+m]) for k in 0:m-1)) for i in 1:M)
                sqrt(72 * s / (M * m^4 * (m * tau0)^2))
            end
        end

        Random.seed!(20260509)
        x = cumsum(randn(1024))
        @test SigmaTau._pdev_core(x, [1], 1.0)[1] ≈ SigmaTau._adev_core(x, [1], 1.0)[1] rtol=1e-15
        @test all(==(0.0), SigmaTau._pdev_core(fill(2.718, 128), [1, 2, 4, 8], 1.0))
        # The parabolic weights sum to zero, so a linear phase gives 0.
        @test all(<(1e-12), SigmaTau._pdev_core(collect(1.0:128.0) .* 0.5 .+ 7.0, [2, 4, 8, 16], 1.0))
        @test isnan(SigmaTau._pdev_core(collect(1.0:10.0), [6], 1.0)[1])
        @test isnan(SigmaTau._pdev_core(collect(1.0:10.0), [0], 1.0)[1])
        # Octave grid plus non-octave and near-limit m exercise the rolling sums.
        ms = [1, 2, 3, 4, 7, 8, 16, 31, 32, 127, 255, 511]
        ref = pdev_reference(x, ms, 1.0)
        @test SigmaTau._pdev_core(x, ms, 1.0) ≈ ref atol=1e-25 rtol=1e-12

        pd = PhaseData(x, 1.0)
        res = pdev(pd, ms)
        @test res.deviation_type == :pdev && res.dev ≈ ref
        @test length(res.edf) == length(ms) && all(isfinite, res.edf) && !isempty(res.noise_type)
        @test all(ci_lower(res) .<= res.dev .+ 1e-12) && all(res.dev .<= ci_upper(res) .+ 1e-12)
        @test isempty(pdev(pd, ms; ci=false).ci)
        y = randn(200) .* 1e-9
        @test pdev(FrequencyData(y, 1.0), [1, 2, 4]).dev ≈ pdev(PhaseData([0.0; cumsum(y)], 1.0), [1, 2, 4]).dev
    end

    @testset "NEFF_RELIABLE boundary" begin
        # N_eff ≥ 30 uses the lag-1 ACF, below it the B1 ratio; both classify.
        @test NEFF_RELIABLE == 30
        Random.seed!(2026)
        @test identify_noise(cumsum(randn(29)), [1])[1] != :unknown
        @test identify_noise(cumsum(randn(31)), [1])[1] != :unknown
    end
end

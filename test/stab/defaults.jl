# test/stab/defaults.jl — zero-argument default grids and the
# phase/frequency conversion helpers.

using Test
using Random
using SigmaTau

@testset "Default grids and conversions" begin
    # Mock data equivalent to NIST SP1065 sample
    x = collect(1.0:100.0) .^ 2  # simple quadratic
    m_values = [1, 2, 4]
    tau0 = 1.0

    @testset "Zero-arg convenience: default octave m-grid per deviation" begin
        # Octave grid shape: 1, 2, 4, ..., 2^floor(log2(m_max)).
        # m_max per kernel matches each `_*_core`'s L-check.
        @testset "_default_m_values per kernel" begin
            # Octave up to each kernel's m_max at N = 1024. Ordinary and modified
            # estimators need ≥ 2 windows at the largest τ; total estimators keep
            # the single-subsequence point (see `_kernel_m_max`). MHDEV/HTDEV
            # (m_max = 255) stop at 128, MTIE (1023) reaches 512.
            top = Dict(:mhdev => 7, :htdev => 7, :mtie => 9)
            ks = (:adev, :pdev, :totdev, :mdev, :tdev, :mtotdev, :ttotdev, :hdev,
                  :htotdev, :mhdev, :htdev, :mhtotdev, :mtie)
            @test isempty([k for k in ks
                           if SigmaTau._default_m_values(1024, k) != 2 .^ (0:get(top, k, 8))])
        end

        @testset "_default_m_values argument validation" begin
            @test_throws ArgumentError SigmaTau._default_m_values(1024, :nope)
            # N too small to admit any m ≥ 1
            @test_throws ArgumentError SigmaTau._default_m_values(1, :adev)
            @test_throws ArgumentError SigmaTau._default_m_values(3, :mhdev)
        end

        @testset "HTOTDEV at N=1536 stays clear of the diff(x) m-max" begin
            # Regression: HTOTDEV operates on
            # y = diff(x), so its valid m-range is `m ≤ (N−1) ÷ 3`. At
            # N=1536, naïvely using `N ÷ 3 = 512` would push past
            # `(1535) ÷ 3 = 511`, putting an invalid `m = 512` into the
            # octave grid and producing a trailing NaN.
            ms = SigmaTau._default_m_values(1536, :htotdev)
            @test maximum(ms) <= (1536 - 1) ÷ 3
            @test 512 ∉ ms
            Random.seed!(2026)
            p = PhaseData(noise_gen(PhaseData, 1536, 1.0; h=Dict(0 => 1.0)).x, 1.0)
            r = htotdev(p; ci = false)
            @test all(isfinite, r.dev)
        end

        @testset "zero-arg API matches explicit-m_values dispatch" begin
            # Synthesize a deterministic WFM phase fixture and verify
            # `dev(pd)` ≡ `dev(pd, _default_m_values(N, :dev))` across both
            # PhaseData and FrequencyData entry points.
            Random.seed!(2026)
            N    = 512
            tau0 = 1.0
            x    = noise_gen(PhaseData, N, tau0; h=Dict(0 => 1.0)).x
            p    = PhaseData(x, tau0)
            f    = FrequencyData(diff(x) ./ tau0, tau0)

            fns = (adev, mdev, tdev, hdev, mhdev, htdev, totdev, mtotdev, ttotdev,
                   htotdev, mhtotdev, mtie, pdev)
            # FrequencyData resolves its grid from the length(f.y) + 1 phase points.
            @test isempty([fn for fn in fns
                           if fn(p; ci=false).dev != fn(p, SigmaTau._default_m_values(N, Symbol(fn)); ci=false).dev ||
                              fn(f; ci=false).tau != SigmaTau._default_m_values(N, Symbol(fn)) .* tau0])
        end

        @testset "kwargs pass through (ci, confidence)" begin
            Random.seed!(7)
            p = PhaseData(noise_gen(PhaseData, 512, 1.0; h=Dict(0 => 1.0)).x, 1.0)
            # ci=true populates CI; ci=false leaves them empty (neff stays).
            r_ci = adev(p; ci = true)
            r_no = adev(p; ci = false)
            @test !isempty(r_ci.ci) && !isempty(r_ci.edf) && isempty(r_no.ci) && isempty(r_no.edf)
            @test r_no.neff == r_ci.neff && !isempty(r_no.neff)
            # totdev's `correct_bias` kwarg passes through.
            r_raw  = totdev(p; ci = false, correct_bias = false)
            r_corr = totdev(p; ci = false, correct_bias = true)
            @test r_raw.tau == r_corr.tau
        end
    end

    @testset "phase ↔ frequency conversion" begin
        # Canonical mapping: y[k] = (x[k+1] − x[k]) / τ₀, length N → N−1.
        @testset "definition and shape" begin
            τ₀ = 0.5
            x  = [0.0, 1.0, 3.0, 7.0, 15.0]
            pd = PhaseData(x, τ₀; source="clk", start=86400.0)
            fd = FrequencyData(pd)
            @test fd.y == diff(x) ./ τ₀
            @test (fd.tau0, fd.source, fd.start) == (τ₀, "clk", 86400.0)
            @test PhaseData(fd).start == 86400.0
            # A missing phase sample loses the frequency on both sides; a
            # missing frequency sample leaves the later phase unknown.
            fg = FrequencyData(PhaseData([0.0, 1.0, NaN, 3.0, 4.0]))
            @test isequal(fg.y, [1.0, NaN, NaN, 1.0])
            @test_throws ArgumentError PhaseData(fg)
        end

        @testset "round-trip identities (modulo the phase offset)" begin
            Random.seed!(31)
            τ₀ = 1.0
            # Frequency → phase → frequency is the identity: N samples give
            # N + 1 phase points starting at 0 (SP1065, Stable32, allantools).
            y  = randn(64)
            fd = FrequencyData(y, τ₀)
            pd1 = PhaseData(fd)
            @test length(pd1.x) == length(y) + 1
            @test pd1.x[1] == 0.0
            fd2 = FrequencyData(pd1)
            @test fd2.y ≈ y
            # Phase → frequency → phase recovers x − x[1].
            x   = cumsum(randn(64))
            pd  = PhaseData(x, τ₀)
            pd2 = PhaseData(FrequencyData(pd))
            @test length(pd2.x) == length(x)
            @test pd2.x ≈ x .- x[1]
        end

        @testset "deviation equivalence: adev/mdev/hdev agree on either domain" begin
            # Build a synthetic WFM phase record, compute ADEV/MDEV/HDEV via
            # the phase path and via the frequency path (after FrequencyData(p)).
            # Second- and third-difference kernels are shift-invariant and the
            # conversion is lossless up to the offset, so both paths agree to
            # rounding.
            Random.seed!(2026)
            τ₀ = 1.0
            p_full = PhaseData(noise_gen(PhaseData, 4096, τ₀; h=Dict(0 => 1.0)).x, τ₀)
            f_diff = FrequencyData(p_full)
            @test length(f_diff.y) == length(p_full.x) - 1
            ms = [1, 4, 16, 64]
            for fn in (adev, mdev, hdev)
                r_phase = fn(p_full,  ms; ci = false)
                r_freq  = fn(f_diff,  ms; ci = false)
                @test r_freq.tau == r_phase.tau
                @test isapprox(r_freq.dev, r_phase.dev; rtol = 1e-10)
            end
        end
    end
end

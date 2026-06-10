# tierms.jl — RMS Time Interval Error tests (NIST SP1065 §5.2.18 eq. 37).
#
# The allantools reference values below were computed once with
# allantools 2024.06 (`tierms(x, rate=1.0, data_type="phase", taus=...)`)
# and frozen at %.17e. allantools implements the same eq. 37 statistic
# (its per-window max−min over the pair (x[i], x[i+m]) is |x[i+m]−x[i]|)
# and reports τ = m·τ0 with n = N − m spans, exactly as SigmaTau does.

using Test
using Random
using SigmaTau
using SigmaTau: _tierms_core, _f64

@testset "TIE rms" begin

    @testset "kernel hand-checkable fixture" begin
        # x = [0, 1, 0.5, 2, 1.5]:
        #   m=1: diffs (1, −0.5, 1.5, −0.5)  → mean sq 3.75/4   → √0.9375
        #   m=2: diffs (0.5, 1, 1)           → mean sq 2.25/3   → √0.75
        #   m=3: diffs (2, 0.5)              → mean sq 4.25/2   → √2.125
        #   m=4: diffs (1.5,)                → mean sq 2.25     → 1.5
        x = [0.0, 1.0, 0.5, 2.0, 1.5]
        devs = _tierms_core(x, [1, 2, 3, 4], 1.0)
        @test devs ≈ [sqrt(0.9375), sqrt(0.75), sqrt(2.125), 1.5]

        # Constant phase → TIE rms = 0 at every τ.
        @test all(==(0.0), _tierms_core(fill(2.5, 64), [1, 2, 4, 8], 1.0))

        # Linear ramp with unit slope: every span of m·τ₀ changes by exactly
        # m, so the rms is m itself.
        ramp = collect(0.0:99.0)
        m_grid = [1, 2, 5, 10, 50]
        @test _tierms_core(ramp, m_grid, 1.0) ≈ Float64.(m_grid)

        # NaN guards: m ≥ N and m < 1.
        @test isnan(_tierms_core(ramp, [100], 1.0)[1])
        @test isnan(_tierms_core(ramp, [200], 1.0)[1])
        @test isnan(_tierms_core(ramp, [0], 1.0)[1])
    end

    @testset "naive reference parity" begin
        # Direct transcription of SP1065 eq. 37 on a noisy fixture.
        Random.seed!(20260610)
        N = 256
        xn = cumsum(randn(N))
        m_grid = [1, 2, 5, 10, 50, 200]
        ref = [sqrt(sum(abs2, xn[(1 + m):N] .- xn[1:(N - m)]) / (N - m))
               for m in m_grid]
        @test _tierms_core(xn, m_grid, 1.0) ≈ ref atol=0.0 rtol=1e-15
    end

    @testset "allantools cross-check (frozen, stable32gen.DAT)" begin
        dat_path = joinpath(@__DIR__, "..", "fixtures", "validation", "stable32gen.DAT")
        if !isfile(dat_path)
            @info "Stable32 phase fixture not present; skipping TIE rms allantools cross-check"
        else
            lines = readlines(dat_path)
            x = parse.(Float64, strip.(lines[11:end]))
            @test length(x) == 8192

            ms = [1, 2, 4, 8, 16, 64, 256, 1024, 4096]
            # allantools 2024.06 tierms devs at the same m (frozen 2026-06-10).
            ref = [8.24298917534318165e-01, 8.24007385339782128e-01,
                   8.24797519785257616e-01, 8.47165434097589709e-01,
                   8.29207543359109289e-01, 9.32617883330743247e-01,
                   1.81598469805780205e+00, 6.03698500637962354e+00,
                   1.67273413520060075e+01]
            devs = _tierms_core(x, ms, 1.0)
            for (d, r) in zip(devs, ref)
                @test isapprox(d, r; rtol=1e-12)
            end

            # Wrapper: τ = m·τ0 (same convention as allantools) and the
            # n = N − m span count.
            r = tierms(PhaseData(x, 1.0), ms)
            @test r.tau ≈ Float64.(ms)
            @test r.dev ≈ devs
            @test r.neff == length(x) .- ms
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

        ms = [1, 2, 4, 8, 16, 64, 128, 256, 500]
        # allantools 2024.06 tierms devs at the same m (frozen 2026-06-10).
        ref = [2.91376788072597259e-01, 4.04261931677273978e-01,
               5.69965604523575253e-01, 8.43948710052766171e-01,
               1.19297966750043383e+00, 1.87504525789628285e+00,
               2.54740653677371887e+00, 2.89748300023606342e+00,
               3.95001034791551353e+00]
        devs = _tierms_core(x, ms, 1.0)
        for (d, r) in zip(devs, ref)
            @test isapprox(d, r; rtol=1e-12)
        end
    end

    @testset "API wrapper conventions" begin
        Random.seed!(20260610)
        N = 256
        xn = cumsum(randn(N))
        pd = PhaseData(xn, 2.0)
        ms = [1, 2, 4, 8]

        res = tierms(pd, ms)
        @test res.deviation_type == :tierms
        @test res.tau ≈ ms .* 2.0
        # No published EDF model: CI fields stay empty even with ci=true.
        @test isempty(res.noise_type)
        @test isempty(res.ci)
        @test isempty(res.edf)
        @test res.neff == N .- ms

        # `confidence` is accepted for signature uniformity (no-op):
        res_conf = tierms(pd, ms; confidence=0.9)
        @test res_conf.dev ≈ res.dev
        @test isempty(res_conf.edf)

        # NaN grid points report neff = 0.
        r_edge = tierms(pd, [N - 1, N])
        @test isfinite(r_edge.dev[1]) && r_edge.neff[1] == 1
        @test isnan(r_edge.dev[2]) && r_edge.neff[2] == 0

        # FrequencyData entry point: cumsum-equivalence.
        Random.seed!(20260610)
        y = randn(100) .* 1e-9
        fd = FrequencyData(y, 1.0)
        pd_eq = PhaseData(cumsum(y) .* 1.0, 1.0)
        @test tierms(fd, [1, 2, 4]).dev ≈ tierms(pd_eq, [1, 2, 4]).dev

        # TauMode and zero-arg dispatch: octave grid up to m_max = N − 1.
        @test SigmaTau._default_m_values(1024, :tierms) == 2 .^ (0:9)
        r_oct = tierms(PhaseData(xn, 1.0), Octave)
        r_default = tierms(PhaseData(xn, 1.0))
        @test r_oct.tau == r_default.tau
        @test r_oct.dev ≈ r_default.dev
        @test r_default.tau[end] == 128.0    # 2^7 ≤ 255 < 2^8

        # stability suite integration.
        suite = stability(PhaseData(xn, 1.0); devs=(:adev, :tierms), ci=false)
        @test haskey(suite, :tierms)
        @test suite[:tierms].dev ≈ tierms(PhaseData(xn, 1.0), Octave).dev
    end
end

# test/stab/noise_types.jl — total, Allan and Hadamard families across the
# power-law noise types, identities and drift handling.

using Test
using Random
using SigmaTau

@testset "Estimators across noise types" begin
    # Mock data equivalent to NIST SP1065 sample
    x = collect(1.0:100.0) .^ 2  # simple quadratic
    m_values = [1, 2, 4]
    tau0 = 1.0

    @testset "MTOTDEV across all 5 power-law noise types" begin
        # Verify kernel parity + end-to-end pipeline on every SP1065 noise
        # type (WPM α=2, FLPM α=1, WHFM α=0, FLFM α=-1, RWFM α=-2). Synthesis
        # is via the Kasdin–Walter helper (`_gen_powerlaw_phase`).
        ms = [1, 2, 4, 8, 16]
        bad = []
        for alpha in (2.0, 1.0, 0.0, -1.0, -2.0)
            Random.seed!(20260507)
            x = _gen_powerlaw_phase(alpha, 1024; tau0=1.0)
            for m in ms
                isapprox(SigmaTau._mtotdev_core(x, [m], 1.0)[1], sqrt(LK.mtotdev_var(x, m, 1.0));
                         atol=1e-25, rtol=1e-12) || push!(bad, (alpha, m, :kernel))
            end
            # Bias correction and CI bounds: finite, positive and ordered.
            r = mtotdev(PhaseData(x, 1.0), ms)
            all(i -> 0 < ci_lower(r)[i] <= r.dev[i] <= ci_upper(r)[i] < Inf, eachindex(ms)) ||
                push!(bad, (alpha, :ci))
        end
        @test isempty(bad)
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
        rm = mtotdev(pd, ms; ci=false, correct_bias=false)
        rt = ttotdev(pd, ms; ci=false, correct_bias=false)
        @test rt.tau == rm.tau
        @test isapprox(rt.dev, rm.dev .* (rm.tau ./ sqrt(3.0)); rtol=1e-14)

        # CI / EDF / α flow through with the same τ/√3 scaling.
        rm_ci = mtotdev(pd, ms; ci=true, correct_bias=true)
        rt_ci = ttotdev(pd, ms; ci=true, correct_bias=true)
        f = rm_ci.tau ./ sqrt(3.0)
        @test isapprox(rt_ci.dev,          rm_ci.dev          .* f; rtol=1e-14)
        @test isapprox(ci_lower(rt_ci), ci_lower(rm_ci) .* f; rtol=1e-14)
        @test isapprox(ci_upper(rt_ci), ci_upper(rm_ci) .* f; rtol=1e-14)
        @test rt_ci.edf        == rm_ci.edf
        @test rt_ci.noise_type == rm_ci.noise_type
        @test rt_ci.neff       == rm_ci.neff

        # FrequencyData entry routes via _freq_to_phase (cumsum). Reconstructed
        # phase has one fewer sample than the original x, so values differ
        # from the direct PhaseData call at the percent level — assert the
        # entry point runs and stays in the right ballpark, not bit-identity.
        fd = FrequencyData(diff(x) ./ tau0, tau0)
        rt_fd = ttotdev(fd, ms; ci=false, correct_bias=false)
        @test isapprox(rt_fd.dev, rt.dev; rtol=5e-3)
    end

    @testset "MHTOTDEV finite and positive on every noise type" begin
        # MHTOTDEV has no external reference; this is its smoke test.
        @test all(all(d -> isfinite(d) && d > 0,
                      _mhtotdev_core(_gen_powerlaw_phase(alpha, 1024; tau0=1.0), [1, 2, 4, 8], 1.0))
                  for alpha in (2.0, 1.0, 0.0, -1.0, -2.0))
    end

    @testset "MHTOTDEV deterministic-trend invariance (wrapper)" begin
        # The mhtotdev wrapper removes the record's least-squares quadratic
        # (frequency drift) before the kernel by default. The fit is a linear
        # projection, so adding a deterministic offset/frequency-offset/drift
        # polynomial to a record cannot change the result — exactly, up to
        # roundoff. The raw kernel itself is NOT drift-immune (its per-window
        # detrend removes only the local frequency offset), which is the
        # reason the wrapper-level removal exists; lock in both facts.
        using Random
        N    = 512
        tau0 = 1.0
        ms   = [1, 2, 4, 8, 16, 32]

        Random.seed!(20260509)
        x = _gen_powerlaw_phase(0.0, N; tau0=tau0)
        t = collect(0:N-1) .* tau0

        # The synthetic noise is unit-scale (y has unit variance), so the
        # trend magnitudes are chosen relative to that — large enough that
        # the drift dominates the statistic at long τ, which is what makes
        # the remove_drift=false assertion below bite.
        offset = 0.1                      # frequency offset → linear phase
        drift  = 0.1                      # frequency drift  → quadratic phase
        x_trend = x .+ offset .* t .+ 0.5 .* drift .* t .^ 2

        base       = mhtotdev(PhaseData(x, tau0),       ms; ci=false, correct_bias=false)
        with_trend = mhtotdev(PhaseData(x_trend, tau0), ms; ci=false, correct_bias=false)
        @test all(isfinite, with_trend.dev)
        @test with_trend.dev ≈ base.dev rtol = 1e-8

        # remove_drift=false on the drifting record must differ at long τ —
        # the kernel alone re-admits drift through the window boundaries.
        raw = mhtotdev(PhaseData(x_trend, tau0), ms; ci=false, correct_bias=false,
                       remove_drift=false)
        @test raw.dev[end] > 2 * base.dev[end]
    end

    @testset "ADEV/HDEV across all 5 power-law noise types" begin
        # Kernel accuracy across the full SP1065 α range: each kernel against
        # the same formula evaluated in BigFloat, so the tolerance measures the
        # Float64 kernel's own rounding and not a Float64 reference's. ADEV and
        # HDEV are direct differences and sit near machine precision; MDEV and
        # MHDEV go through prefix sums, whose conditioning on a random-walk
        # phase (α = −2, excursions ~10³ against unit second differences)
        # bounds them at ~1e-11 over many seeds, hence 1e-10.
        bad = []
        for alpha in (2.0, 1.0, 0.0, -1.0, -2.0)
            Random.seed!(123)
            x  = _gen_powerlaw_phase(alpha, 1024; tau0=1.0)
            xb = BigFloat.(x)
            cs = pushfirst!(cumsum(xb), big(0.0))
            pairs = ((_adev_core, m -> LK.adev_var(xb, m, big(1.0))),
                     (_mdev_core, m -> LK.mdev_var(xb, m, big(1.0), cs)),
                     (_hdev_core, m -> LK.hdev_var(xb, m, big(1.0))),
                     (_mhdev_core, m -> LK.mhdev_var(xb, m, big(1.0), cs)))
            for (core, ref) in pairs, m in (1, 2, 4, 8, 16)
                isapprox(core(x, [m], 1.0)[1], Float64(sqrt(ref(m))); atol=1e-25, rtol=1e-10) ||
                    push!(bad, (alpha, core, m))
            end
        end
        @test isempty(bad)
    end

    @testset "MDEV/MHDEV precision on a clock with a frequency offset" begin
        # A phase ramp (uncorrected frequency offset) contributes nothing to the
        # modified estimators but makes naive prefix sums grow like N²·offset;
        # the kernels remove the least-squares line first. Compare against the
        # same formula evaluated in BigFloat; without the detrend the error is
        # ~3e-7 here, with it the ~1e-11 prefix-sum conditioning floor remains.
        rng  = Xoshiro(20261026)
        N    = 20000
        tau0 = 1.0
        p    = noise_gen(PhaseData, N, tau0; sigma1 = Dict(0 => 1e-12), rng = rng)
        x    = p.x .+ 1e-9 .* (0:N-1)
        xb   = BigFloat.(x)
        xb_cs = pushfirst!(cumsum(xb), big(0.0))
        ms = [1, 16, 256]
        exact(v) = [Float64(sqrt(v(xb, m, big(tau0), xb_cs))) for m in ms]
        @test all(isapprox.(_mdev_core(x, ms, tau0), exact(LK.mdev_var); rtol=1e-10))
        @test all(isapprox.(_mhdev_core(x, ms, tau0), exact(LK.mhdev_var); rtol=1e-10))
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
        edfs = SigmaTau.calculate_edf(:totdev, ones(length(m_values_eq)),
                                               noises_wpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs)
        @test all(>(0.0), edfs)

        # α=1 (FLPM)
        noises_flpm = fill(:FLPM, length(m_values_eq))
        edfs = SigmaTau.calculate_edf(:totdev, ones(length(m_values_eq)),
                                               noises_flpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs)

        # htotdev gets the HDEV-style fallback at α=2,1 too.
        edfs_h = SigmaTau.calculate_edf(:htotdev, ones(length(m_values_eq)),
                                                 noises_wpm, m_values_eq, taus_eq, N_eq, T_eq)
        @test all(isfinite, edfs_h)
    end
end

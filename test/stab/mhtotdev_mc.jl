# test/stab/mhtotdev_mc.jl — fast regression tripwire for the MHTOTDEV bias and
# EDF coefficients measured by tools/mc_mhtotdev.jl.
#
# This is NOT a precision check — the authoritative coefficients come from the
# full workstation sweep (N up to 32768, R=3000). Here a tiny fresh Monte Carlo
# (N=513, R=100) re-measures the bias B and EDF at a couple of (α, m) cells and
# asserts they fall within a loose multiplicative tolerance of the shipped
# `bias_correction(:mhtot)` / `_coeff_mhtot` values, so a silent drift in the
# kernels or the synthesis is caught without a multi-hour run. With R=100 the
# EDF sampling SE is ≈ √(2/100) ≈ 14 %, hence the deliberately loose 40 %
# tolerance; the bias is far tighter thanks to the paired V/W draw.

using Test
using Random
using Statistics
using FFTW                                  # AbstractFFTs backend for synthesis
using SigmaTau
using SigmaTau: _gen_powerlaw_y, _coeff_mhtot, bias_correction

@testset "MHTOTDEV bias/EDF regression (fast MC tripwire)" begin
    N    = 513
    tau0 = 1.0
    T    = (N - 1) * tau0
    R    = 100
    ms   = [16, 32]                          # both ≥ 16 (the τ/τ0 ≥ 16 fit window)
    seed = 20260531
    noise_sym = Dict(0 => :WHFM, -2 => :RWFM)
    # Process-stable seed (plain integer arithmetic, unlike Base `hash`), so the
    # tripwire is deterministic across Julia versions. α+3 maps {-2..2}→{1..5}.
    sd(alpha, r) = seed + 1_000_003 * r + 7919 * (alpha + 3) + N

    for alpha in (0, -2)
        V = Matrix{Float64}(undef, length(ms), R)   # raw MHTOTVAR
        W = Matrix{Float64}(undef, length(ms), R)    # MHVAR (mhdev²)
        for r in 1:R
            rng = Xoshiro(sd(alpha, r))
            pd  = PhaseData(cumsum(_gen_powerlaw_y(float(alpha), N; rng=rng)), tau0)
            V[:, r] = mhtotdev(pd, ms; ci=false, correct_bias=false).dev .^ 2
            W[:, r] = mhdev(pd,    ms; ci=false).dev .^ 2
        end

        b, c = _coeff_mhtot(alpha)

        for k in eachindex(ms)
            τ      = ms[k] * tau0
            meanV  = mean(@view V[k, :]); varV = var(@view V[k, :])
            edf    = 2 * meanV^2 / varV
            B      = meanV / mean(@view W[k, :])
            edf_pred = b * (T / τ) - c
            # B is τ/T-linear, so evaluate the prediction at this cell's τ.
            B_pred   = bias_correction([noise_sym[alpha]], :mhtot, [τ], T)[1]

            @test isfinite(edf) && edf > 0
            @test isapprox(edf, edf_pred; rtol = 0.40)   # ~2.8σ tripwire at R=100
            @test isapprox(B,   B_pred;   rtol = 0.15)
        end
    end
end

@testset "MTOT bias direction tripwire (B < 1 for red noise)" begin
    # bias_correction(:mtot) shipped inverted values (B up to 1.31) until
    # 2026-08; the correct direction is B < 1 — the half-average detrend
    # strips genuine low-frequency power, so MTOT understates MVAR (Riley's
    # Stable32 CI document, confirmed by MC). This tripwire re-measures one
    # red-noise cell so a future regression of the sign cannot ship silently.
    N = 513; tau0 = 1.0; T = (N - 1) * tau0
    R = 100; m = 16
    sd(r) = 77_000_003 + 1_000_003 * r
    V = zeros(R); W = zeros(R)
    for r in 1:R
        rng = Xoshiro(sd(r))
        pd  = PhaseData(cumsum(_gen_powerlaw_y(-2.0, N; rng=rng)), tau0)
        V[r] = mtotdev(pd, [m]; ci=false, correct_bias=false).dev[1]^2
        W[r] = mdev(pd, [m]; ci=false).dev[1]^2
    end
    B_meas = mean(V) / mean(W)
    B_ship = bias_correction([:RWFM], :mtot, [m * tau0], T)[1]
    @test B_meas < 0.9                       # direction: biased low
    @test B_ship < 1.0                       # shipped table on the right side
    @test isapprox(B_meas, B_ship; rtol = 0.15)
end

@testset "MHTOTDEV EDF extrapolation floor (CI inversion guard)" begin
    # Historical context: the pre-quadratic kernel's EDF model was calibrated
    # only to T/τ ≥ 8 and dipped below 1 at the default grid's last point,
    # silently tripping the Gaussian fallback and collapsing the CI ~13×.
    # The current calibration measures the edge cells, so the model stays
    # above 1 everywhere the default grid reaches — but the ν ≥ 1 floor in
    # calculate_edf stays as a dormant guard, and this testset pins both
    # facts so neither can regress silently.
    N = 1024; tau0 = 1.0; T = (N - 1) * tau0
    ms     = [64, 128, 256]              # 256 = default-grid maximum (T/τ ≈ 4)
    taus   = ms .* tau0
    devs   = ones(length(ms))
    noises = fill(:WHFM, length(ms))
    edfs = SigmaTau.calculate_edf(:mhtotdev, devs, noises, ms, taus, N, T)
    @test all(edfs .>= 1.0)
    @test edfs[end] > 1.5                # edge zone measured, not floored
    lo, hi = SigmaTau.confidence_intervals(devs, edfs, noises, N, 0.683)
    # CI width must widen monotonically toward the reach limit — no inversion.
    @test hi[end] - lo[end] > hi[end-1] - lo[end-1]

    # The floor itself, exercised beyond the estimator's reach (T/τ = 2 gives
    # a raw model value of ≈ 0.19): χ² stays in use, never the Gaussian branch.
    edf_far = SigmaTau.calculate_edf(:mhtotdev, [1.0], [:WHFM], [512], [512.0], N, T)
    @test edf_far[1] == 1.0
end

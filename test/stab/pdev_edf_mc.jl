# test/stab/pdev_edf_mc.jl — fast regression tripwire for the PVAR (PDEV) EDF
# model wired in `src/edf.jl::_pvar_edf`, after Vernotte–Chen–Rubiola
# 2020 (arXiv:2005.13631).
#
# This is NOT a precision check — the model is published, not measured here. A
# tiny fresh Monte Carlo (N=1025, R=200) re-measures the EDF as ν = 2·E[V]²/Var[V]
# (the paper's Eq. 19) at a few (α, m) cells inside the m ≤ N/4 validity window
# and asserts it falls within a loose multiplicative tolerance of `_pvar_edf`, so
# a silent drift in the kernel or the synthesis is caught. With R=200 the EDF
# sampling SE is ≈ √(2/200) ≈ 10 %, hence the deliberately loose 25 % tolerance.

using Test
using Random
using Statistics
using FFTW                                  # AbstractFFTs backend for synthesis
using SigmaTau
using SigmaTau: _gen_powerlaw_y, _pvar_edf

@testset "PDEV/PVAR EDF regression (fast MC tripwire)" begin
    N    = 1025
    tau0 = 1.0
    R    = 200
    ms   = [4, 16, 64]                       # all in the m ≤ N/4 (=256) window
    seed = 20260531
    # Process-stable seed (plain integer arithmetic), deterministic across Julia
    # versions. α+3 maps {-2..2} → {1..5}.
    sd(alpha, r) = seed + 1_000_003 * r + 7919 * (alpha + 3) + N

    for alpha in (2, 0, -2)
        V = Matrix{Float64}(undef, length(ms), R)   # PVAR estimates (dev²)
        for r in 1:R
            rng = Xoshiro(sd(alpha, r))
            pd  = PhaseData(cumsum(_gen_powerlaw_y(float(alpha), N; rng=rng)), tau0)
            V[:, r] = pdev(pd, ms; ci=false).dev .^ 2
        end

        for k in eachindex(ms)
            m       = ms[k]
            meanV   = mean(@view V[k, :]); varV = var(@view V[k, :])
            edf     = 2 * meanV^2 / varV
            edf_pred = _pvar_edf(alpha, m, N)

            @test isfinite(edf) && edf > 0
            @test isfinite(edf_pred) && edf_pred > 0
            @test isapprox(edf, edf_pred; rtol = 0.25)
        end
    end

    # Cheap structural checks on the public API (one representative record).
    rng = Xoshiro(sd(0, 1))
    pd  = PhaseData(cumsum(_gen_powerlaw_y(0.0, N; rng=rng)), tau0)
    grid = [1, 2, 4, 8, 16, 32, 64, 128]

    r = pdev(pd, grid)                        # ci=true (default)
    @test length(r.edf) == length(grid)
    @test all(isfinite, r.edf) && all(>(0), r.edf)
    @test all(isfinite, ci_lower(r)) && all(isfinite, ci_upper(r))
    @test all(ci_lower(r) .<= r.dev .+ 1e-12) && all(r.dev .<= ci_upper(r) .+ 1e-12)
    @test !isempty(r.noise_type)
    # EDF decreases with τ across the formula window (m ≤ N/4).
    @test issorted(r.edf[1:6]; rev = true)

    rf = pdev(pd, grid; ci = false)      # empty-CI contract preserved
    @test isempty(rf.edf) && isempty(rf.noise_type)
    @test isempty(rf.ci)
    @test rf.neff == N .- 2 .* grid      # pdev: N − 2m windows, even with ci=false

    # PVAR(τ₀) ≡ AVAR(τ₀): kernel identity at m = 1, and EDF uses ADEV's value.
    @test pdev(pd, [1]; ci = false).dev[1] == adev(pd, [1]; ci = false).dev[1]
    @test _pvar_edf(0, 1, N) == SigmaTau._calc_edf_core(0, 2, 1, 1, 1, N)
end

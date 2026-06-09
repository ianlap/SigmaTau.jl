# tools/mc_mhtotdev.jl — Monte Carlo study of MHTOTDEV bias B(α) and EDF.
#
# MHTOTDEV is novel to SigmaTau; no published reference gives its bias factor or
# its equivalent degrees of freedom. This harness measures both from synthesized
# known-α power-law noise and fits the coefficient tables that ship in
# `src/edf.jl` (`_coeff_mhtot`) and `bias_correction(:mhtot, …)`.
#
# Method (see docs/src/theory/mhtotdev_bias_edf.md):
#   • For each (α, N, m) cell, over R realizations of pure-α noise drawn from
#     independent RNG streams, compute V̂ = MHTOTVAR and Ŵ = MHVAR (= mhdev²) on
#     the SAME realization at matched m.
#   • EDF:  edf = 2·E[V̂]² / Var[V̂]      (χ²_edf scaling of a variance estimator)
#   • Bias: B(α) = E[V̂] / E[Ŵ]           (ratio-of-means; MHVAR is the reference)
#   • Cross-check: the measured MHVAR slope vs τ must match the analytic
#     power-law μ(α) before any ratio is trusted.
#   • Fit edf = b·(T/τ) − c per α by weighted least squares over τ ≤ T/10.
#
# Usage:
#   julia --project=. --threads=auto tools/mc_mhtotdev.jl            # full sweep
#   julia --project=. --threads=auto tools/mc_mhtotdev.jl quick      # laptop check
#
# Reproducible: every realization uses Xoshiro(hash((seed, α, N, r))), so results
# are independent of thread count. Emits a CSV (per-cell) and a JSON (metadata +
# fitted coefficients) under tools/artifacts/. Generated artifacts are not
# tracked in the repository.

using SigmaTau
using SigmaTau: _gen_powerlaw_y, _default_m_values
using Statistics
using Random
using Printf
using Dates

# ── Configuration ───────────────────────────────────────────────────────────
const QUICK = "quick" in ARGS

const MASTER_SEED = 20260531
const ALPHAS      = (2, 1, 0, -1, -2)               # WHPM … RWFM
# Record length N as the number of samples (powers of two). T = (N-1)·τ0.
const N_GRID      = QUICK ? (1024, 2048) :
                           (1024, 2048, 4096, 8192, 16384, 32768)
const R           = QUICK ? 200  : 3000             # realizations per cell
const R_ANCHOR    = QUICK ? 400  : 3000             # top-N tier (uniform with R)
const N_BOOT      = QUICK ? 500  : 2000             # bootstrap resamples
const MIN_NSUBS   = 16                              # drop tiny-window cells
# Total-estimator EDF model validity: Howe et al. 2000 (TotHvar, eqn 7) state
# the fit "should be used only if … τ/τ0 ≥ 16" — below that the data-extension
# gives no advantage over the plain estimator. Fit only m ≥ FIT_M_MIN.
const FIT_M_MIN   = 16

# Analytic MHVAR variance slope σ²_MH(τ) ∝ τ^μ for each α (modified family).
_mu(alpha::Int) = alpha <= 0 ? -alpha - 1 : (alpha == 1 ? -2 : -3)

# Process-stable seed derivation (splitmix64 finalizer). Base `hash` is NOT
# guaranteed stable across Julia versions/processes, so it cannot anchor a
# reproducible artifact; this mixes the integer components deterministically.
# `reinterpret` handles negative α via two's-complement bits.
function _seed(parts::Integer...)
    h = 0x243f6a8885a308d3
    for p in parts
        h ⊻= reinterpret(UInt64, Int64(p))
        h = (h ⊻ (h >> 30)) * 0xbf58476d1ce4e5b9
        h = (h ⊻ (h >> 27)) * 0x94d049bb133111eb
        h ⊻= (h >> 31)
    end
    return h
end

# ── Per-cell Monte Carlo ──────────────────────────────────────────────────────

"""
    cell_stats(alpha, N, ms, R) → (meanV, varV, meanW, edf, B, edf_se, B_se)

Run R independent realizations of pure-α noise at length N; return per-m
statistics of V̂=MHTOTVAR and Ŵ=MHVAR. Realizations parallelize over independent
RNG streams; the kernels' own threading composes underneath.
"""
function cell_stats(alpha::Int, N::Int, ms::Vector{Int}, nreal::Int)
    nc = length(ms)
    V = Matrix{Float64}(undef, nc, nreal)   # MHTOTVAR samples
    W = Matrix{Float64}(undef, nc, nreal)   # MHVAR samples
    Threads.@threads :dynamic for r in 1:nreal
        rng = Xoshiro(_seed(MASTER_SEED, alpha, N, r))
        y  = _gen_powerlaw_y(float(alpha), N; rng=rng)
        pd = PhaseData(cumsum(y), 1.0)
        # correct_bias=false: measure the RAW MHTOTVAR — the bias B we are
        # estimating must not be pre-applied (it would make this circular now
        # that mhtotdev applies the measured B by default).
        v  = mhtotdev(pd, ms; ci=false, correct_bias=false).dev
        w  = mhdev(pd,    ms; ci=false).dev
        @inbounds for k in 1:nc
            V[k, r] = v[k]^2
            W[k, r] = w[k]^2
        end
    end

    meanV = vec(mean(V; dims=2));  varV = vec(var(V; dims=2))   # Bessel
    meanW = vec(mean(W; dims=2))
    edf = 2 .* meanV .^ 2 ./ varV
    B   = meanV ./ meanW

    # Paired bootstrap SE (resample realization indices, preserving V–W pairing).
    edf_se = zeros(nc);  B_se = zeros(nc)
    boot_edf = Vector{Float64}(undef, N_BOOT)
    boot_B   = Vector{Float64}(undef, N_BOOT)
    brng = Xoshiro(_seed(MASTER_SEED, 0x_b001, alpha, N))   # 0xb001 = "boot" tag
    idx = Vector{Int}(undef, nreal)
    for k in 1:nc
        Vk = @view V[k, :];  Wk = @view W[k, :]
        for bidx in 1:N_BOOT
            rand!(brng, idx, 1:nreal)
            mV = 0.0; mV2 = 0.0; mW = 0.0
            @inbounds for j in idx
                vj = Vk[j]; mV += vj; mV2 += vj*vj; mW += Wk[j]
            end
            mV /= nreal; mV2 /= nreal; mW /= nreal
            vV = (mV2 - mV^2) * nreal / (nreal - 1)
            boot_edf[bidx] = 2 * mV^2 / vV
            boot_B[bidx]   = mV / mW
        end
        edf_se[k] = std(boot_edf)
        B_se[k]   = std(boot_B)
    end

    return meanV, varV, meanW, edf, B, edf_se, B_se
end

# Weighted 2-parameter least squares for edf = b·x − c (x = T/τ). Weights w.
# Returns (b, c, R²).
function wls_edf(x::Vector{Float64}, y::Vector{Float64}, w::Vector{Float64})
    Sw  = sum(w);  Sx = sum(w .* x);  Sy = sum(w .* y)
    Sxx = sum(w .* x .^ 2);  Sxy = sum(w .* x .* y)
    Δ   = Sw * Sxx - Sx^2
    b   = (Sw * Sxy - Sx * Sy) / Δ
    a   = (Sy - b * Sx) / Sw          # intercept = −c
    c   = -a
    ŷ   = b .* x .+ a
    ss_res = sum(w .* (y .- ŷ) .^ 2)
    ss_tot = sum(w .* (y .- sum(w .* y) / Sw) .^ 2)
    R2  = 1 - ss_res / ss_tot
    return b, c, R2
end

# TotHvar EDF form (Howe et al. 2000 eqn 7): edf = (T/τ) / (b0 + b1·u), u = τ/T.
# Linearize as 1/edf = b0·u + b1·u² (no intercept) and fit by WLS in edf-space:
# a point with edf y and SE σ_y maps to z = 1/y with σ_z ≈ σ_y/y², so its weight
# is y⁴/σ_y². Returns (b0, b1, R²) with R² measured back in edf-space.
function wls_htot(u::Vector{Float64}, y::Vector{Float64}, se::Vector{Float64})
    z  = 1.0 ./ y
    w  = (y .^ 4) ./ (se .^ 2)
    Su2 = sum(w .* u .^ 2);  Su3 = sum(w .* u .^ 3);  Su4 = sum(w .* u .^ 4)
    Suz = sum(w .* u .* z);  Su2z = sum(w .* u .^ 2 .* z)
    Δ   = Su2 * Su4 - Su3^2
    b0  = (Suz * Su4 - Su2z * Su3) / Δ
    b1  = (Su2 * Su2z - Su3 * Suz) / Δ
    ŷ   = 1.0 ./ (b0 .* u .+ b1 .* u .^ 2)   # back to edf
    wy  = 1.0 ./ se .^ 2
    ss_res = sum(wy .* (y .- ŷ) .^ 2)
    ss_tot = sum(wy .* (y .- sum(wy .* y) / sum(wy)) .^ 2)
    return b0, b1, 1 - ss_res / ss_tot
end

# Log-log slope of y vs x (unweighted), interior points already selected.
function loglog_slope(x::Vector{Float64}, y::Vector{Float64})
    lx = log.(x); ly = log.(y); n = length(lx)
    (n * sum(lx .* ly) - sum(lx) * sum(ly)) / (n * sum(lx .^ 2) - sum(lx)^2)
end

# ── Minimal JSON writer (no external dep under --project=.) ───────────────────
_jval(x::Real)    = (isfinite(x) ? @sprintf("%.8g", x) : "null")
_jval(x::Integer) = string(x)
_jval(x::AbstractString) = "\"$x\""
_jval(x::Symbol)  = "\"$x\""
_jval(x::AbstractVector) = "[" * join(_jval.(x), ", ") * "]"

# ── Driver ────────────────────────────────────────────────────────────────────

function main()
    git_sha = try
        strip(read(`git rev-parse HEAD`, String))
    catch
        "unknown"
    end
    stamp = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    outdir = joinpath(@__DIR__, "artifacts")
    mkpath(outdir)
    tag = QUICK ? "quick" : "full"
    csv_path  = joinpath(outdir, "mhtotdev_mc_$(tag).csv")
    json_path = joinpath(outdir, "mhtotdev_mc_$(tag).json")

    total_reals = length(ALPHAS) * ((length(N_GRID) - 1) * R + R_ANCHOR)
    println("MHTOTDEV Monte Carlo  [$(tag)]  threads=$(Threads.nthreads())  ",
            "R=$R  R_anchor=$R_ANCHOR  N=$(N_GRID)  seed=$MASTER_SEED")
    println("Total realizations: $total_reals (each = 1 MHTOTDEV + 1 MHDEV over the m-grid)")
    flush(stdout)

    # Accumulate per-cell rows and per-α pooled fit inputs.
    rows = NamedTuple[]
    fit_x  = Dict(a => Float64[] for a in ALPHAS)  # T/τ (fit window, m ≥ FIT_M_MIN)
    fit_y  = Dict(a => Float64[] for a in ALPHAS)  # edf
    fit_w  = Dict(a => Float64[] for a in ALPHAS)  # 1/SE²
    fit_se = Dict(a => Float64[] for a in ALPHAS)  # SE(edf)
    B_by_alpha   = Dict(a => Tuple{Float64,Float64,Float64}[] for a in ALPHAS)  # (τ/T, B, B_se)
    slope_ok = true

    for alpha in ALPHAS
        for N in N_GRID
            T = float(N - 1)
            ms = _default_m_values(N, :mhtotdev)
            ms = filter(m -> (N - 4m + 1) >= MIN_NSUBS, ms)
            isempty(ms) && continue
            nreal = (N == N_GRID[end]) ? R_ANCHOR : R   # anchor at the top N
            @printf("  → α=%+d  N=%5d  R=%d  cells=%d  …\n", alpha, N, nreal, length(ms))
            flush(stdout)
            local stats
            elapsed = @elapsed stats = cell_stats(alpha, N, ms, nreal)
            meanV, varV, meanW, edf, B, edf_se, B_se = stats

            # μ(α) cross-check on the MHVAR estimate (interior m only).
            if length(ms) >= 5
                lo, hi = 2, length(ms) - 1
                s = loglog_slope(Float64.(ms[lo:hi]), meanW[lo:hi])
                if abs(s - _mu(alpha)) > 0.25
                    @warn "MHVAR slope off" alpha N measured=round(s, digits=3) expected=_mu(alpha)
                    slope_ok = false
                end
            end

            for k in eachindex(ms)
                m = ms[k]; τ = float(m); tot = T / τ; nsubs = N - 4m + 1
                push!(rows, (; alpha, N, m, T_over_tau=tot, nsubs,
                             R=nreal, meanV=meanV[k], varV=varV[k],
                             edf=edf[k], edf_se=edf_se[k],
                             meanW=meanW[k], B=B[k], B_se=B_se[k]))
                # Both the bias and the EDF fit use the τ/τ0 ≥ 16 validity window
                # (Howe et al. 2000). The small-m cells — especially m=1, which is
                # near-degenerate and carries the tightest SE (largest weight) —
                # otherwise skew the bias toward a non-physical value.
                if m >= FIT_M_MIN && isfinite(edf[k]) && edf_se[k] > 0
                    push!(B_by_alpha[alpha], (τ / T, B[k], B_se[k]))
                    push!(fit_x[alpha],  tot)
                    push!(fit_y[alpha],  edf[k])
                    push!(fit_w[alpha],  1 / edf_se[k]^2)
                    push!(fit_se[alpha], edf_se[k])
                end
            end
            @printf("  ✓ α=%+d  N=%5d  cells=%2d  R=%d  (%.1f s)\n",
                    alpha, N, length(ms), nreal, elapsed)
            flush(stdout)
        end
    end

    # ── Write per-cell CSV ────────────────────────────────────────────────────
    open(csv_path, "w") do io
        println(io, "alpha,N,m,T_over_tau,nsubs,R,meanV,varV,edf,edf_se,meanW,B,B_se")
        for r in rows
            @printf(io, "%d,%d,%d,%.6e,%d,%d,%.10e,%.10e,%.6f,%.6f,%.10e,%.6f,%.6f\n",
                    r.alpha, r.N, r.m, r.T_over_tau, r.nsubs, r.R,
                    r.meanV, r.varV, r.edf, r.edf_se, r.meanW, r.B, r.B_se)
        end
    end

    # ── Per-α fits ─────────────────────────────────────────────────────────────
    # Two EDF parameterizations from the literature, fit over the τ/τ0 ≥ 16
    # validity window and selected by R² in edf-space:
    #   form A (Mod-Totvar, Vernotte–Howe eqn 2):  edf = b·(T/τ) − c
    #   form B (TotHvar, Howe et al. 2000 eqn 7):   edf = (T/τ)/(b0 + b1·τ/T)
    edf_fit = Dict{Int,NamedTuple}()
    B_fit   = Dict{Int,NamedTuple}()
    for alpha in ALPHAS
        if length(fit_x[alpha]) >= 3
            x = fit_x[alpha]; y = fit_y[alpha]; w = fit_w[alpha]; se = fit_se[alpha]
            b, c, R2a = wls_edf(x, y, w)
            b0, b1, R2b = wls_htot(1.0 ./ x, y, se)   # u = τ/T = 1/(T/τ)
            best = R2b > R2a ? "tothvar" : "modtot"
            edf_fit[alpha] = (; b, c, R2a, b0, b1, R2b, best)
        else
            edf_fit[alpha] = (; b=NaN, c=NaN, R2a=NaN, b0=NaN, b1=NaN, R2b=NaN, best="none")
        end

        # Bias: weighted-mean (α-only) and linear-in-τ/T fit; compare RSS.
        pts = B_by_alpha[alpha]
        tt  = Float64[p[1] for p in pts]
        bb  = Float64[p[2] for p in pts]
        ww  = Float64[1 / max(p[3], eps())^2 for p in pts]
        Bmean = sum(ww .* bb) / sum(ww)
        rss_const = sum(ww .* (bb .- Bmean) .^ 2)
        # B = b0 + b1·(τ/T). wls_edf models y = b·x + a and returns c = −a, so
        # the intercept b0 is −c_returned.
        b1, negint, _ = wls_edf(tt, bb, ww)
        b0 = -negint
        rss_lin = sum(ww .* (bb .- (b1 .* tt .+ b0)) .^ 2)
        B_fit[alpha] = (; Bmean, nbias = Bmean - 1, b0, b1, rss_const, rss_lin,
                        improve = 1 - rss_lin / rss_const)
    end

    # ── Write JSON (metadata + fits) ──────────────────────────────────────────
    open(json_path, "w") do io
        println(io, "{")
        println(io, "  \"meta\": {")
        println(io, "    \"git_sha\": ", _jval(git_sha), ",")
        println(io, "    \"timestamp\": ", _jval(stamp), ",")
        println(io, "    \"julia_version\": ", _jval(string(VERSION)), ",")
        println(io, "    \"nthreads\": ", Threads.nthreads(), ",")
        println(io, "    \"master_seed\": ", MASTER_SEED, ",")
        println(io, "    \"seed_scheme\": ", _jval("Xoshiro(splitmix64(seed, alpha, N, r))"), ",")
        println(io, "    \"N_grid\": ", _jval(collect(N_GRID)), ",")
        println(io, "    \"R\": ", R, ", \"R_anchor\": ", R_ANCHOR, ",")
        println(io, "    \"n_boot\": ", N_BOOT, ", \"min_nsubs\": ", MIN_NSUBS, ",")
        println(io, "    \"fit_m_min\": ", FIT_M_MIN, ",")
        println(io, "    \"mhvar_slope_check_passed\": ", slope_ok)
        println(io, "  },")
        println(io, "  \"edf_fit\": {")
        for (i, a) in enumerate(ALPHAS)
            f = edf_fit[a]
            comma = i < length(ALPHAS) ? "," : ""
            println(io, "    \"$a\": {",
                    "\"modtot\": {\"b\": ", _jval(f.b), ", \"c\": ", _jval(f.c), ", \"R2\": ", _jval(f.R2a), "}, ",
                    "\"tothvar\": {\"b0\": ", _jval(f.b0), ", \"b1\": ", _jval(f.b1), ", \"R2\": ", _jval(f.R2b), "}, ",
                    "\"best\": ", _jval(f.best), "}", comma)
        end
        println(io, "  },")
        println(io, "  \"bias_fit\": {")
        for (i, a) in enumerate(ALPHAS)
            f = B_fit[a]
            comma = i < length(ALPHAS) ? "," : ""
            println(io, "    \"$a\": {\"B_mean\": ", _jval(f.Bmean),
                    ", \"nbias\": ", _jval(f.nbias),
                    ", \"b0\": ", _jval(f.b0), ", \"b1\": ", _jval(f.b1),
                    ", \"linear_improve\": ", _jval(f.improve), "}", comma)
        end
        println(io, "  }")
        println(io, "}")
    end

    # ── Console summary (the numbers that go into edf.jl) ─────────────────────
    println("\n── Fitted MHTOTDEV coefficients ──  (fit window m ≥ $FIT_M_MIN, i.e. τ/τ0 ≥ 16)")
    println("EDF form A (Mod-Totvar): edf = b·(T/τ) − c   |   form B (TotHvar): edf = (T/τ)/(b0 + b1·τ/T)")
    for a in ALPHAS
        e = edf_fit[a]; f = B_fit[a]
        @printf("  α=%+d:  A b=%7.4f c=%8.4f (R²=%.4f)   B b0=%.4f b1=%.4f (R²=%.4f)  → best=%s\n",
                a, e.b, e.c, e.R2a, e.b0, e.b1, e.R2b, e.best)
        @printf("          bias B=%.4f  (nbias=%+.4f)  [lin τ/T improve %.1f%%]\n",
                f.Bmean, f.nbias, 100 * f.improve)
    end
    println("\nWrote:\n  $csv_path\n  $json_path")
    # Artifacts are written first (diagnostics available), but fail loudly so a
    # failed-slope run cannot be silently consumed by a coefficient update.
    slope_ok || error("MHVAR μ(α) slope check FAILED for ≥1 cell — the synthesis " *
                      "is miscalibrated and the bias/EDF fits are NOT trustworthy. " *
                      "Artifacts were written for inspection but must not be used.")
    return nothing
end

main()

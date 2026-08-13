# tools/nnls_noise_id.jl — weighted-NNLS power-law noise identification study.
#
# Question: given ONE composite record, can a nonnegative least-squares fit of the
# measured variance profile σ̂²(τ) against a per-noise-type basis recover the four
# component amplitudes — and does the estimator choice (MDEV / MHDEV / HDEV /
# MHTOTDEV) decide whether the α = −4 (RRFM) component is identifiable at all?
#
# Method:
#   • Unit generators for WPM, WFM, RWFM, RRFM emit fractional frequency y with
#     NO per-realization rescale, so the ensemble mean of D² is a clean power law
#     (a per-realization rescale would peg every draw to the same τ₀ level and
#     destroy the sampling statistics this study is about).
#   • Basis: K̄ᵢ^D(τ) = ensemble mean of D²(τ) over R_basis unit realizations and
#     Varᵢ^D(τ) = ensemble variance, one basis per estimator D. The basis is
#     EMPIRICAL, not analytic: MHTOTDEV's reflected/extended windows make the
#     response deviate from τ^μ near τ → T, and MDEV's α = −4 column does not
#     converge at all (that is the exhibit — MDEV alone). MHTOTDEV's per-window
#     QUADRATIC LS detrend is matched to the third difference's annihilation
#     order, so its kernel converges over the full Hadamard range α > −5 and its
#     α = −4 column is a genuine τ³ response, like MHDEV's and HDEV's.
#   • Composite: y = Σᵢ cᵢ yᵢ with cᵢ = σᵢ(τ₀)/√(K̄ᵢ^mhdev(τ₀)). Ground-truth
#     amplitudes are aᵢ = cᵢ²: every estimator is a quadratic form in the data and
#     the components are independent, so E[D²_composite] = Σᵢ cᵢ² K̄ᵢ^D exactly
#     (drift removal is a linear projection, so it does not break additivity).
#   • Fit: min ‖W^½(Ka − v̂)‖² subject to a ≥ 0, solved EXACTLY by enumerating all
#     2ⁿ − 1 nonempty supports (plus a = 0) and keeping the feasible minimum-
#     residual solution; cross-checked against Lawson–Hanson on every solve.
#     (NonNegLeastSquares.jl would be the drop-in escape hatch if n ever grew
#     past exhaustive range; at n = 4 enumeration is exact and instant.)
#
# Usage:
#   julia --project=. --threads=auto tools/nnls_noise_id.jl          # full study
#   julia --project=. --threads=auto tools/nnls_noise_id.jl quick    # laptop check
#
# Reproducible: every realization draws from Xoshiro(splitmix64(seed, tag, i, r)),
# so results are independent of thread count. Emits basis / trials / summary CSVs
# and a JSON (metadata + gate outcomes + per-estimator metrics) under
# tools/artifacts/.

using SigmaTau
using SigmaTau: _gen_powerlaw_y, calculate_edf, _nnls_exact
using Statistics
using Random
using Distributions
using LinearAlgebra
using Printf
using Dates

# ── Configuration ───────────────────────────────────────────────────────────
const QUICK = "quick" in ARGS

const MASTER_SEED = 20260805
const TAU0        = 1.0
# Full-mode record length: T ≈ 1.3e5 s puts the last usable octave at
# τ = 32768 s and gives RRFM ≥ 1.3 decades of dominance past the 1600 s
# crossover. R is halved relative to the N=65536 design to hold the sweep
# near the same wall-clock; basis/trial SEs grow only by √2.
const N_SAMPLES   = QUICK ? 8192 : 131072
const R_BASIS     = QUICK ? 200 : 1000      # unit realizations per noise type
const R_TRIAL     = QUICK ? 100 : 500       # composite realizations
const R_PURE      = QUICK ? 100 : 500       # pure-noise recovery gate
const N_CHI2      = QUICK ? 500 : 2000      # parametric χ²(ν) bootstrap draws
const MIN_NSUBS   = 16                      # keep only m with N − 4m + 1 ≥ 16

const TYPES      = (:WPM, :WFM, :RWFM, :RRFM)
const TYPE_ALPHA = Dict(:WPM => 2, :WFM => 0, :RWFM => -2, :RRFM => -4)
const ESTIMATORS = (:mdev, :mhdev, :hdev, :mhtotdev)

# Component MHDEV(τ₀) levels. The three crossovers are where consecutive
# components have equal MHDEV; the σ table is derived from them, so the composite
# spans all four regimes inside the τ grid.
#   WPM/WFM  at τ_a:  σ_WPM  = σ_WFM · τ_a       (slopes −3/2 and −1/2)
#   WFM/RWFM at τ_b:  σ_RWFM = σ_WFM / τ_b       (slopes −1/2 and +1/2)
#   RWFM/RRFM at τ_c: σ_RRFM = σ_RWFM / τ_c      (slopes +1/2 and +3/2)
const CROSSOVERS = QUICK ? (8.0, 64.0, 512.0) : (10.0, 125.0, 1600.0)
const SIGMA1 = QUICK ?
    Dict(:WPM => 8e-12, :WFM => 1e-12, :RWFM => 1.5625e-14, :RRFM => 3.0517578125e-17) :
    Dict(:WPM => 1e-11, :WFM => 1e-12, :RWFM => 8e-15,      :RRFM => 5e-18)

# MHDEV log-log slope of a pure α noise: MHVAR ∝ τ^(−α−1) ⇒ MHDEV ∝ τ^(−α−1)/2.
_mhdev_slope(alpha::Int) = -(alpha + 1) / 2
const SLOPE_TOL = 0.15

# A τ point counts toward a noise type's identification score only where the
# ANALYTIC design makes that type dominant by a clear margin. Octaves sitting on
# a crossover are excluded: there the two neighbouring components are equal by
# construction, so "which one is dominant" is a coin flip and scoring it would
# manufacture failures. 1.2 in variance ≈ 1.10 in deviation.
const DOMINANCE_MARGIN = 1.2

# Disjoint seed tags keep the basis, trial, pure-noise and bootstrap ensembles
# statistically independent (the basis must not be reused as its own test set).
const TAG_BASIS = 0x0ba5
const TAG_TRIAL = 0x7241
const TAG_PURE  = 0x9021
const TAG_CHI2  = 0xc812
const TAG_FFT   = 0x4ff7
const TAG_LEN   = 0x1e07

# The shipped EDF tables stop at α = −2 and `_alpha_from_noise` maps any
# unrecognized symbol to α = 0, so the α = −4 rows borrow the RWFM model. This is
# a documented approximation of the EDF-weighted fit, not a claim about RRFM EDF.
const EDF_SYMBOL = Dict(:WPM => :WHPM, :WFM => :WHFM, :RWFM => :RWFM, :RRFM => :RWFM)

# Gate thresholds for pure-noise recovery. Every third-difference estimator is
# gated on every type INCLUDING RRFM: MHTOTDEV's quadratic-detrend kernel
# converges at α = −4, so it is held to the same standard as MHDEV and HDEV.
# MDEV is the only estimator left out — its α = −4 column does not converge, and
# that failure is the paper's exhibit, not a gate.
const GATE_ESTIMATORS = (:mhdev, :hdev, :mhtotdev)
const GATE_DOMINANCE  = 0.90     # correct type's share of fitted variance
const GATE_CROSSTALK  = 0.10     # any wrong type's share of fitted variance
const GATE_SOLVER_TOL = 1e-10    # relative objective agreement, exact vs LH

# Process-stable seed derivation (splitmix64 finalizer), identical in spirit to
# tools/mc_mhtotdev.jl: base `hash` is not stable across Julia versions, so it
# cannot anchor a reproducible artifact.
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

# ── Unit noise generators (fractional frequency, unit scale) ─────────────────

"""
    unit_y(kind, N, rng) → Vector{Float64}

One realization of the named noise type as a fractional-frequency series at
whatever absolute scale the generator naturally produces. WPM and WFM come from
the package's internal f^(α/2) spectral shaper; RWFM and RRFM are the raw
single/double cumulative sums of white noise.

The random walks are NOT drawn from the FFT shaper: the shaper zeroes the DC bin
of a periodic record, which caps the very lowest Fourier component and flattens
the α = −4 response at long τ (quantified in the artifacts' `fft_rrfm` block).
The cumsum construction has the correct 1/(2 sin(πfτ₀))^{2,4} response all the
way down to the record's fundamental.
"""
function unit_y(kind::Symbol, N::Int, rng)
    y = if kind === :WPM
        _gen_powerlaw_y(2.0, N; rng = rng)
    elseif kind === :WFM
        _gen_powerlaw_y(0.0, N; rng = rng)
    elseif kind === :RWFM
        cumsum(randn(rng, N))
    elseif kind === :RRFM
        cumsum(cumsum(randn(rng, N)))
    else
        throw(ArgumentError("unit_y: unknown noise type $kind"))
    end
    return y .- mean(y)
end

# Every estimator is called raw and uniformly: no CI, no bias correction. The
# MHTOTDEV default `remove_drift=true` is kept, because that is what users
# compute and what the shipped EDF/bias tables were calibrated against.
function dev2(est::Symbol, pd::PhaseData, ms::Vector{Int})
    r = if est === :mdev
        mdev(pd, ms; ci = false)
    elseif est === :mhdev
        mhdev(pd, ms; ci = false)
    elseif est === :hdev
        hdev(pd, ms; ci = false)
    elseif est === :mhtotdev
        mhtotdev(pd, ms; ci = false, correct_bias = false)
    else
        throw(ArgumentError("dev2: unknown estimator $est"))
    end
    return r.dev .^ 2
end

phase_of(y::Vector{Float64}) = PhaseData(cumsum(y) .* TAU0, TAU0)

# Log-log slope of y vs x (unweighted).
function loglog_slope(x::Vector{Float64}, y::Vector{Float64})
    lx = log.(x); ly = log.(y); n = length(lx)
    (n * sum(lx .* ly) - sum(lx) * sum(ly)) / (n * sum(lx .^ 2) - sum(lx)^2)
end

iqr(v::AbstractVector{Float64}) = quantile(v, 0.75) - quantile(v, 0.25)

# Variance response exponent σ²_D(τ) ∝ τ^μ. The modified kernels separate the
# two PM noises (μ = −α−1 throughout); the unmodified ones collapse WPM and FPM
# onto μ = −2. Only α = 2 differs across the four types used here, but it
# differs by a full power, which moves the WPM/WFM crossover by ~8×.
_response_mu(est::Symbol, alpha::Int) =
    (est === :hdev && alpha >= 1) ? -2.0 : -Float64(alpha) - 1.0

"""
    truth_regions(est, taus) → (dominant, regions)

Ground-truth dominance for estimator `est`, purely from the analytic crossover
design: component `i` sits at `σᵢ(τ₀)²·(τ/τ₀)^μ(αᵢ)` with the TEXTBOOK response
exponent for that estimator's family.

Nothing measured enters. That is deliberate: anchoring on the estimator's own
measured τ₀ levels would let a divergent column define its own truth — MDEV's
α = −4 response at τ₀ is already inflated by the record's low-frequency content,
so an empirical anchor makes RRFM "truth-dominant" from τ ≈ 32 s and quietly
erases the RWFM region the estimator is failing to identify.

Family correction: the modified kernels separate the PM noises (μ = −α−1
throughout) while the unmodified ones put WPM at μ = −2, which moves the
WPM/WFM crossover by ~8× — big enough to change which octaves belong to which
type, so it must be applied. The cross-family O(1) difference in the τ₀ *level*
(HVAR/MHVAR at τ₀) is not corrected; the σ table is MHDEV-calibrated throughout
and the exponent is the order-of-magnitude effect.

`regions[t]` lists the τ indices where `t` leads the runner-up by at least
`DOMINANCE_MARGIN` in variance — the points at which a fit is scored for
identifying `t`. Octaves sitting on a crossover are excluded: there the two
neighbours are equal by construction, so scoring them manufactures failures. A
type with no such octave has its misID reported as null rather than as a pass.
"""
function truth_regions(est::Symbol, taus::Vector{Float64})
    mu = [_response_mu(est, TYPE_ALPHA[t]) for t in TYPES]
    dominant = Vector{Symbol}(undef, length(taus))
    regions  = Dict(t => Int[] for t in TYPES)
    for (k, tau) in enumerate(taus)
        levels = [SIGMA1[TYPES[i]] ^ 2 * (tau / TAU0) ^ mu[i] for i in eachindex(TYPES)]
        i = argmax(levels)
        dominant[k] = TYPES[i]
        runner = maximum(levels[j] for j in eachindex(levels) if j != i)
        levels[i] >= DOMINANCE_MARGIN * runner && push!(regions[TYPES[i]], k)
    end
    return dominant, regions
end

# ── NNLS solvers ─────────────────────────────────────────────────────────────

# The exact solver is the SHIPPED `SigmaTau._nnls_exact` (exhaustive active-set
# enumeration over all 2ⁿ − 1 nonempty supports, seeded with a = 0). Using the
# package's own routine here means every solve in this study — tens of thousands
# of them — is also a regression test of the code users call, and the
# Lawson–Hanson implementation below is a genuinely independent cross-check.
# (NonNegLeastSquares.jl is the escape hatch if the basis ever outgrows
# exhaustive range; at n = 4 enumeration is exact and instant.)

"""
    nnls_lawson_hanson(A, b) → a

Classic Lawson–Hanson (1974, Algorithm NNLS) active-set solver, kept as an
independent cross-check on `SigmaTau._nnls_exact`. Every solve in this study is run
through both and the objectives must agree to `GATE_SOLVER_TOL`.
"""
function nnls_lawson_hanson(A::AbstractMatrix{Float64}, b::Vector{Float64};
                            tol::Float64 = 1e-12, maxiter::Int = 500)
    n = size(A, 2)
    x = zeros(n)
    passive = falses(n)
    grad = A' * b
    gscale = max(maximum(abs, grad), eps())
    outer = 0
    while any(!, passive) && outer < maxiter
        j = 0; best = -Inf
        for i in 1:n
            if !passive[i] && grad[i] > best
                best = grad[i]; j = i
            end
        end
        (j == 0 || best <= tol * gscale) && break
        passive[j] = true
        outer += 1
        inner = 0
        while inner < maxiter
            inner += 1
            idx = findall(passive)
            z = zeros(n)
            zs = try
                A[:, idx] \ b
            catch
                passive[j] = false
                break
            end
            z[idx] = zs
            if minimum(z[idx]) > 0
                x = z
                break
            end
            α = Inf
            for i in idx
                z[i] <= 0 && (α = min(α, x[i] / (x[i] - z[i])))
            end
            isfinite(α) || (α = 0.0)
            x = x .+ α .* (z .- x)
            xscale = max(maximum(abs, x), eps())
            for i in idx
                x[i] <= tol * xscale && (passive[i] = false; x[i] = 0.0)
            end
            any(passive) || break
        end
        grad = A' * (b .- A * x)
    end
    return x
end

"""
    solve_nnls(K, v, w) → NamedTuple

Weighted NNLS: minimize ‖W^½(Ka − v)‖² with a ≥ 0. Columns of W^½K are
equilibrated to unit 2-norm before the solve (the four power-law columns differ
by ~20 orders of magnitude otherwise) and the amplitudes are unscaled after.
κ₂ is reported for the equilibrated matrix — that is the conditioning the solver
actually sees, and the estimator-to-estimator comparison the study reports.
"""
function solve_nnls(K::Matrix{Float64}, v::Vector{Float64}, w::Vector{Float64})
    sw = sqrt.(w)
    A0 = sw .* K
    b  = sw .* v
    d  = [norm(@view A0[:, j]) for j in axes(A0, 2)]
    all(dj -> isfinite(dj) && dj > 0, d) ||
        error("solve_nnls: degenerate basis column (norms = $d)")
    A = A0 ./ d'
    a_ex, rss = _nnls_exact(A, b)
    a_lh      = nnls_lawson_hanson(A, b)
    rss_lh    = sum(abs2, A * a_lh .- b)
    obj_diff  = abs(rss_lh - rss) / max(rss, eps())
    coef_diff = norm(a_lh .- a_ex) / max(norm(a_ex), eps())
    return (a = a_ex ./ d, rss = rss, kappa2 = cond(A, 2),
            obj_diff = obj_diff, coef_diff = coef_diff,
            active = [a_ex[j] > 0 for j in eachindex(a_ex)])
end

# ── Ensembles ────────────────────────────────────────────────────────────────

"""
    basis_ensemble(ms, R) → (Kbar, Kvar, Kh1, Kh2)

Per (estimator, noise type): ensemble mean and variance of D²(τ) over R unit
realizations, plus the two half-ensemble means used for the split-half column
stability diagnostic.
"""
function basis_ensemble(ms::Vector{Int}, R::Int)
    nk = length(ms)
    raw = Dict((e, t) => Matrix{Float64}(undef, nk, R)
               for e in ESTIMATORS, t in TYPES)
    Threads.@threads :dynamic for r in 1:R
        for (ti, t) in enumerate(TYPES)
            rng = Xoshiro(_seed(MASTER_SEED, TAG_BASIS, ti, r))
            pd  = phase_of(unit_y(t, N_SAMPLES, rng))
            for e in ESTIMATORS
                raw[(e, t)][:, r] = dev2(e, pd, ms)
            end
        end
    end
    half = R ÷ 2
    Kbar = Dict(k => vec(mean(v; dims = 2)) for (k, v) in raw)
    Kvar = Dict(k => vec(var(v; dims = 2))  for (k, v) in raw)
    Kh1  = Dict(k => vec(mean(@view(v[:, 1:half]);     dims = 2)) for (k, v) in raw)
    Kh2  = Dict(k => vec(mean(@view(v[:, half+1:end]); dims = 2)) for (k, v) in raw)
    return Kbar, Kvar, Kh1, Kh2
end

"""
    trial_ensemble(ms, R, c) → Dict(estimator => nk × R matrix of D²)

R composite realizations, each y = Σᵢ cᵢ yᵢ with independent per-component
streams drawn from the trial seed tag (disjoint from the basis tag).
"""
function trial_ensemble(ms::Vector{Int}, R::Int, c::Dict{Symbol,Float64})
    nk = length(ms)
    V = Dict(e => Matrix{Float64}(undef, nk, R) for e in ESTIMATORS)
    Threads.@threads :dynamic for r in 1:R
        y = zeros(Float64, N_SAMPLES)
        for (ti, t) in enumerate(TYPES)
            rng = Xoshiro(_seed(MASTER_SEED, TAG_TRIAL, ti, r))
            y .+= c[t] .* unit_y(t, N_SAMPLES, rng)
        end
        pd = phase_of(y)
        for e in ESTIMATORS
            V[e][:, r] = dev2(e, pd, ms)
        end
    end
    return V
end

"""
    pure_ensemble(ms, R, c) → Dict((estimator, type) => nk × R matrix of D²)

R realizations of each single component at its composite amplitude, drawn from
the pure-noise seed tag. Feeds the pure-noise recovery gate.
"""
function pure_ensemble(ms::Vector{Int}, R::Int, c::Dict{Symbol,Float64})
    nk = length(ms)
    P = Dict((e, t) => Matrix{Float64}(undef, nk, R)
             for e in ESTIMATORS, t in TYPES)
    Threads.@threads :dynamic for r in 1:R
        for (ti, t) in enumerate(TYPES)
            rng = Xoshiro(_seed(MASTER_SEED, TAG_PURE, ti, r))
            pd  = phase_of(c[t] .* unit_y(t, N_SAMPLES, rng))
            for e in ESTIMATORS
                P[(e, t)][:, r] = dev2(e, pd, ms)
            end
        end
    end
    return P
end

# ── Fit drivers ──────────────────────────────────────────────────────────────

# Solve every column of V against basis Kmat with the supplied weight rule.
# `wrule(v) -> w` so the empirical (fixed) and EDF-model (v-dependent) schemes
# share one code path.
function solve_ensemble(Kmat::Matrix{Float64}, V::Matrix{Float64}, wrule)
    nk, R = size(V)
    n = size(Kmat, 2)
    A   = Matrix{Float64}(undef, n, R)
    rss = Vector{Float64}(undef, R)
    kap = Vector{Float64}(undef, R)
    od  = Vector{Float64}(undef, R)
    cd  = Vector{Float64}(undef, R)
    act = Vector{String}(undef, R)
    for r in 1:R
        v = V[:, r]
        s = solve_nnls(Kmat, v, wrule(v))
        A[:, r]  = s.a
        rss[r]   = s.rss
        kap[r]   = s.kappa2
        od[r]    = s.obj_diff
        cd[r]    = s.coef_diff
        act[r]   = join((string(TYPES[j]) for j in 1:n if s.active[j]), "+")
    end
    return (; A, rss, kappa2 = median(kap), obj_diff = maximum(od),
            coef_diff = maximum(cd), active = act)
end

"""
    misid_rates(Kmat, A, regions) → (misid, misid_pointwise, dominant_fraction)

Per-noise-type identification scoring over the trial ensemble.

- `misid[t]` — fraction of trials in which type `t` fails: its coefficient is
  zeroed, or it is not the dominant fitted component at *some* τ where the
  analytic design says it dominates. All-or-nothing per trial.
- `misid_pointwise[t]` — fraction of (trial, τ-in-region) pairs where `t` is not
  the dominant fitted component. Finer-grained, and the more informative number
  when a type owns only one or two octaves of the grid.
- `dominant_fraction[t]` — mean fraction of the WHOLE τ grid over which `t` is
  the dominant fitted component. This is the masking diagnostic: an estimator
  whose α = −4 column is inflated can make RRFM dominant everywhere, so RWFM is
  never dominant and its (unconstrained) amplitude error stops being independent
  evidence about that estimator.
"""
function misid_rates(Kmat::Matrix{Float64}, A::Matrix{Float64},
                     regions::Dict{Symbol,Vector{Int}})
    n, R = size(A)
    nk = size(Kmat, 1)
    misid   = Dict(t => 0.0 for t in TYPES)
    pointw  = Dict(t => 0.0 for t in TYPES)
    domfrac = Dict(t => 0.0 for t in TYPES)
    dom = Vector{Int}(undef, nk)
    for r in 1:R
        a = @view A[:, r]
        for k in 1:nk
            dom[k] = argmax(vec(Kmat[k, :]) .* a)
        end
        for k in 1:nk
            domfrac[TYPES[dom[k]]] += 1 / (R * nk)
        end
        for (i, t) in enumerate(TYPES)
            ks = regions[t]
            isempty(ks) && continue
            wrong = count(k -> dom[k] != i, ks)
            pointw[t] += wrong / (R * length(ks))
            (a[i] == 0 || wrong > 0) && (misid[t] += 1 / R)
        end
    end
    for t in TYPES
        isempty(regions[t]) && (misid[t] = NaN; pointw[t] = NaN)
    end
    return misid, pointw, domfrac
end

# Fraction of the fitted variance carried by each type at each τ.
function fitted_fractions(Kmat::Matrix{Float64}, a::Vector{Float64})
    contrib = Kmat .* a'                       # nk × n
    tot = vec(sum(contrib; dims = 2))
    frac = similar(contrib)
    for k in axes(contrib, 1)
        frac[k, :] = tot[k] > 0 ? contrib[k, :] ./ tot[k] : zeros(size(contrib, 2))
    end
    return frac
end

# ── Minimal JSON writer (no external dep under --project=.) ──────────────────
_jval(x::Real)    = (isfinite(x) ? @sprintf("%.8g", x) : "null")
_jval(x::Integer) = string(x)
_jval(x::Bool)    = x ? "true" : "false"
_jval(x::AbstractString) = "\"$x\""
_jval(x::Symbol)  = "\"$x\""
_jval(x::AbstractVector) = "[" * join(_jval.(x), ", ") * "]"

# {"k": v, …} from any iterable of (key, value) pairs.
_jobj(pairs) = "{" * join(("\"$k\": " * _jval(v) for (k, v) in pairs), ", ") * "}"

# ── Driver ───────────────────────────────────────────────────────────────────

function main()
    t_start = time()
    git_sha = try
        strip(read(`git rev-parse HEAD`, String))
    catch
        "unknown"
    end
    stamp  = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")
    outdir = joinpath(@__DIR__, "artifacts")
    mkpath(outdir)
    tag = QUICK ? "quick" : "full"
    basis_csv   = joinpath(outdir, "nnls_noise_id_basis_$(tag).csv")
    trials_csv  = joinpath(outdir, "nnls_noise_id_trials_$(tag).csv")
    summary_csv = joinpath(outdir, "nnls_noise_id_summary_$(tag).csv")
    json_path   = joinpath(outdir, "nnls_noise_id_$(tag).json")

    ms   = [2^k for k in 0:ceil(Int, log2(N_SAMPLES)) if N_SAMPLES - 4 * 2^k + 1 >= MIN_NSUBS]
    taus = Float64.(ms) .* TAU0
    nk   = length(ms)
    T    = (N_SAMPLES - 1) * TAU0

    println("NNLS noise identification  [$(tag)]  threads=$(Threads.nthreads())  ",
            "seed=$MASTER_SEED")
    println("N=$N_SAMPLES  τ₀=$TAU0  m-grid=$(ms)  R_basis=$R_BASIS  ",
            "R_trial=$R_TRIAL  R_pure=$R_PURE  n_chi2=$N_CHI2")
    println("crossovers (s): $(collect(CROSSOVERS))   σ MHDEV(τ₀): ",
            join(("$t=$(SIGMA1[t])" for t in TYPES), "  "))
    flush(stdout)

    # Cost projection from one record-evaluation (all four estimators on one
    # record); the whole study is 4·R_basis + R_trial + 4·R_pure of them.
    let rng = Xoshiro(_seed(MASTER_SEED, TAG_BASIS, 0, 0))
        pd = phase_of(unit_y(:WFM, N_SAMPLES, rng))
        for e in ESTIMATORS; dev2(e, pd, ms); end          # warm up the JIT
        t1 = @elapsed for e in ESTIMATORS; dev2(e, pd, ms); end
        nrec = 4 * R_BASIS + R_TRIAL + 4 * R_PURE
        @printf("First record-evaluation: %.3f s → %d evaluations ≈ %.1f min wall at %d threads\n",
                t1, nrec, t1 * nrec / 60 / Threads.nthreads(), Threads.nthreads())
        flush(stdout)
    end

    # ── 1. Basis ensemble ────────────────────────────────────────────────────
    print("  → basis ensemble … "); flush(stdout)
    t_basis = @elapsed ((Kbar, Kvar, Kh1, Kh2) = basis_ensemble(ms, R_BASIS))
    @printf("%.1f s\n", t_basis); flush(stdout)

    # ── 2. Hard slope self-check on the MHDEV basis columns ──────────────────
    # MHDEV is the reference estimator for calibration, so its columns must be
    # the textbook power laws before anything downstream is trusted.
    lo, hi = 2, nk - 1
    slopes = Dict{Symbol,Float64}()
    slope_ok = true
    for t in TYPES
        s = loglog_slope(taus[lo:hi], sqrt.(Kbar[(:mhdev, t)][lo:hi]))
        slopes[t] = s
        expected = _mhdev_slope(TYPE_ALPHA[t])
        if abs(s - expected) > SLOPE_TOL
            @warn "MHDEV slope off" noise=t measured=round(s, digits=4) expected
            slope_ok = false
        end
    end
    println("  slope check (MHDEV, interior τ): ",
            join((@sprintf("%s %.3f/%.1f", t, slopes[t], _mhdev_slope(TYPE_ALPHA[t]))
                  for t in TYPES), "   "),
            slope_ok ? "  ✓" : "  ✗")
    flush(stdout)

    # ── 3. FFT-shaper RRFM control (documents why cumsum is used) ────────────
    # `_gen_powerlaw_y(-4, N)` zeroes the DC bin of a periodic record, capping the
    # lowest Fourier component. Measure the resulting long-τ deficit against the
    # τ^(+3/2) MHDEV law extrapolated from the interior.
    fft_slope, fft_flat = let R = min(64, R_BASIS)
        acc = zeros(nk)
        for r in 1:R
            rng = Xoshiro(_seed(MASTER_SEED, TAG_FFT, r))
            y   = _gen_powerlaw_y(-4.0, N_SAMPLES; rng = rng)
            acc .+= dev2(:mhdev, phase_of(y .- mean(y)), ms)
        end
        acc ./= R
        s   = loglog_slope(taus[lo:hi], sqrt.(acc[lo:hi]))
        # extrapolate the interior fit to the last τ and compare
        ref = sqrt(acc[hi]) * (taus[end] / taus[hi])^_mhdev_slope(-4)
        (s, sqrt(acc[end]) / ref)
    end
    @printf("  FFT RRFM control: MHDEV slope %.3f (cumsum %.3f, ideal %.1f), long-τ %.3f× extrapolation\n",
            fft_slope, slopes[:RRFM], _mhdev_slope(-4), fft_flat)
    flush(stdout)

    # ── 3b. RRFM column diagnostics ──────────────────────────────────────────
    # (a) Record-length dependence. MVAR does not converge for α ≤ −3, so MDEV's
    #     α = −4 column is set by the record's lowest Fourier component: it grows
    #     with T instead of following τ³, and the same basis is therefore invalid
    #     at any other N. Repeat the column at N/4 and take the level ratio —
    #     ≈ 1 for a convergent estimator, ≈ 4 for a T-proportional divergence.
    # (b) MHTOTDEV's α = −4 column, with the wrapper's default global drift
    #     removal switched OFF. The quadratic-detrend kernel is drift-immune on
    #     its own, so the wrapper projection is no longer load-bearing for the
    #     column's shape: drift-off must give the same τ³ variance slope as
    #     drift-on. Kept as a standing check that the two agree.
    Rd   = max(R_BASIS ÷ 2, 20)
    Nq   = N_SAMPLES ÷ 4
    ms_q = filter(m -> Nq - 4m + 1 >= MIN_NSUBS, ms)
    kq   = length(ms_q)
    accq  = Dict(e => zeros(kq) for e in ESTIMATORS)
    accnd = zeros(nk)
    diaglock = ReentrantLock()
    t_diag = @elapsed Threads.@threads :dynamic for r in 1:Rd
        pdq  = phase_of(unit_y(:RRFM, Nq, Xoshiro(_seed(MASTER_SEED, TAG_LEN, r))))
        vals = Dict(e => dev2(e, pdq, ms_q) for e in ESTIMATORS)
        # Same seed as the basis ensemble → paired drift-on/drift-off comparison.
        pdf  = phase_of(unit_y(:RRFM, N_SAMPLES,
                               Xoshiro(_seed(MASTER_SEED, TAG_BASIS,
                                             findfirst(==(:RRFM), TYPES), r))))
        nd   = mhtotdev(pdf, ms; ci = false, correct_bias = false,
                        remove_drift = false).dev .^ 2
        lock(diaglock) do
            for e in ESTIMATORS
                accq[e] .+= vals[e]
            end
            accnd .+= nd
        end
    end
    len_ratio = Dict(e => median(Kbar[(e, :RRFM)][1:kq] ./ (accq[e] ./ Rd))
                     for e in ESTIMATORS)
    nodrift_slope = loglog_slope(taus[lo:hi], (accnd ./ Rd)[lo:hi])
    @printf("  RRFM column vs record length (K̄(N)/K̄(N/4), 1 = convergent): %s   [%.1f s]\n",
            join((@sprintf("%s %.2f", e, len_ratio[e]) for e in ESTIMATORS), "  "), t_diag)
    @printf("  MHTOTDEV RRFM column slope with remove_drift=false: %.3f (ideal +3)\n",
            nodrift_slope)
    flush(stdout)

    # ── 4. Calibration and ground truth ──────────────────────────────────────
    c      = Dict(t => SIGMA1[t] / sqrt(Kbar[(:mhdev, t)][1]) for t in TYPES)
    a_true = Dict(t => c[t]^2 for t in TYPES)
    a_vec  = [a_true[t] for t in TYPES]
    println("  calibration: ", join((@sprintf("c_%s=%.4g", t, c[t]) for t in TYPES), "  "))
    flush(stdout)

    Kmat = Dict(e => hcat((Kbar[(e, t)] for t in TYPES)...) for e in ESTIMATORS)

    # ── 5. Trials ────────────────────────────────────────────────────────────
    print("  → trial ensemble … "); flush(stdout)
    t_trial = @elapsed (V = trial_ensemble(ms, R_TRIAL, c))
    @printf("%.1f s\n", t_trial); flush(stdout)

    # ── 6. Pure-noise ensemble (recovery gate) ───────────────────────────────
    print("  → pure-noise ensemble … "); flush(stdout)
    t_pure = @elapsed (P = pure_ensemble(ms, R_PURE, c))
    @printf("%.1f s\n", t_pure); flush(stdout)

    # ── 7. Per-estimator fits and metrics ────────────────────────────────────
    println("  truth regions (τ where a type leads by ≥$(DOMINANCE_MARGIN)× in variance), per estimator:")
    results   = Dict{Symbol,Any}()
    trial_rows = NamedTuple[]
    solver_obj = 0.0
    solver_coef = 0.0
    edf_clamped = 0

    for e in ESTIMATORS
        K = Kmat[e]

        # Ground-truth dominance for this estimator, from the analytic design
        # with family-correct exponents. Also feeds the EDF model's noise types.
        dom_true, regions = truth_regions(e, taus)
        @printf("    %-9s %s\n", e,
                join((@sprintf("%s %s", t, Int.(taus[regions[t]])) for t in TYPES), "  "))
        any(isempty, (regions[t] for t in TYPES)) &&
            @warn "a noise type has no scoring τ for $e — its misID rate is reported as null"
        noises = [EDF_SYMBOL[t] for t in dom_true]
        nu_raw   = calculate_edf(e, ones(nk), noises, ms, taus, N_SAMPLES, T)
        # The b·(T/τ) − c total-family form can return a non-physical EDF at the
        # longest τ; floor it at one degree of freedom and count the repairs.
        nu = [(isfinite(v) && v >= 1.0) ? v : 1.0 for v in nu_raw]
        edf_clamped += count(v -> !(isfinite(v) && v >= 1.0), nu_raw)

        # Empirical weights: inverse of the measured composite variance of σ̂².
        w_emp = 1.0 ./ vec(var(V[e]; dims = 2))

        fits = Dict{Symbol,Any}()
        for (wname, wrule) in ((:empirical, _ -> w_emp),
                               (:edf,       v -> nu ./ (2 .* v .^ 2)))
            f = solve_ensemble(K, V[e], wrule)
            solver_obj  = max(solver_obj, f.obj_diff)
            solver_coef = max(solver_coef, f.coef_diff)

            relerr = Dict(t => abs.(f.A[i, :] .- a_true[t]) ./ a_true[t]
                          for (i, t) in enumerate(TYPES))
            spread = Dict(t => begin
                              med = median(f.A[i, :])
                              med > 0 ? iqr(f.A[i, :]) / med : Inf
                          end for (i, t) in enumerate(TYPES))
            zero_rate = Dict(t => count(iszero, f.A[i, :]) / R_TRIAL
                             for (i, t) in enumerate(TYPES))

            mis, misp, domfrac = misid_rates(K, f.A, regions)

            fits[wname] = (; relerr, spread, zero_rate, kappa2 = f.kappa2,
                           rss = median(f.rss), misid = mis, misid_pointwise = misp,
                           dominant_fraction = domfrac,
                           obj_diff = f.obj_diff, coef_diff = f.coef_diff)

            for r in 1:R_TRIAL
                push!(trial_rows, (; estimator = e, weights = wname, trial = r,
                                   a = f.A[:, r], rss = f.rss[r], active = f.active[r]))
            end
        end

        # ── Parametric χ²(ν) bootstrap at fixed truth ────────────────────────
        # v̂_k = μ_k · χ²_{ν_k}/ν_k with μ = K·a_true and the EDF model's ν. The
        # basis is exact by construction, so whatever error survives is the EDF /
        # χ²-sampling contribution alone.
        mu = K * a_vec
        w_oracle = nu ./ (2 .* mu .^ 2)
        Vb = Matrix{Float64}(undef, nk, N_CHI2)
        brng = Xoshiro(_seed(MASTER_SEED, TAG_CHI2, findfirst(==(e), ESTIMATORS)))
        chisq = [Chisq(nu[k]) for k in 1:nk]
        for r in 1:N_CHI2, k in 1:nk
            Vb[k, r] = mu[k] * rand(brng, chisq[k]) / nu[k]
        end
        fb = solve_ensemble(K, Vb, _ -> w_oracle)
        solver_obj  = max(solver_obj, fb.obj_diff)
        solver_coef = max(solver_coef, fb.coef_diff)
        chi2_relerr = Dict(t => abs.(fb.A[i, :] .- a_true[t]) ./ a_true[t]
                           for (i, t) in enumerate(TYPES))
        chi2_mis, chi2_misp, _ = misid_rates(K, fb.A, regions)

        # ── Split-half basis-column stability (MDEV's α = −4 divergence exhibit)
        split = Dict(t => begin
                         h1 = Kh1[(e, t)]; h2 = Kh2[(e, t)]
                         d = abs.(h1 .- h2) ./ ((h1 .+ h2) ./ 2)
                         (max = maximum(d), last = d[end])
                     end for t in TYPES)
        col_slope = Dict(t => loglog_slope(taus[lo:hi], Kbar[(e, t)][lo:hi]) for t in TYPES)

        # ── Pure-noise recovery ─────────────────────────────────────────────
        pure = Dict{Symbol,Any}()
        for (i, t) in enumerate(TYPES)
            wp = 1.0 ./ (c[t]^4 .* Kvar[(e, t)])
            domin = Vector{Float64}(undef, R_PURE)
            cross = Vector{Float64}(undef, R_PURE)
            recov = Vector{Float64}(undef, R_PURE)
            for r in 1:R_PURE
                s = solve_nnls(K, P[(e, t)][:, r], wp)
                solver_obj  = max(solver_obj, s.obj_diff)
                solver_coef = max(solver_coef, s.coef_diff)
                fr = fitted_fractions(K, s.a)
                # Interior τ only: the grid endpoints carry the largest estimator
                # variance and the least basis leverage, so a per-τ dominance floor
                # there measures window edge effects, not identification.
                domin[r] = minimum(@view fr[lo:hi, i])
                cross[r] = maximum(maximum(@view fr[lo:hi, j]) for j in 1:length(TYPES) if j != i)
                recov[r] = abs(s.a[i] - a_true[t]) / a_true[t]
            end
            pure[t] = (dominance = median(domin), crosstalk = median(cross),
                       relerr = median(recov))
        end

        results[e] = (; fits, chi2_relerr, chi2_mis, chi2_misp, split, col_slope,
                      pure, nu, dom_true, regions)
        @printf("  ✓ %-9s κ₂(emp)=%.3g  RRFM misID emp=%.3f edf=%.3f χ²=%.3f  relerr(RRFM)=%.3f\n",
                e, fits[:empirical].kappa2, fits[:empirical].misid[:RRFM],
                fits[:edf].misid[:RRFM], chi2_mis[:RRFM],
                median(fits[:empirical].relerr[:RRFM]))
        flush(stdout)
    end

    # ── 8. Gate evaluation ───────────────────────────────────────────────────
    solver_ok = solver_obj <= GATE_SOLVER_TOL
    gate_rows = NamedTuple[]
    recovery_ok = true
    for e in GATE_ESTIMATORS, t in TYPES
        p = results[e].pure[t]
        ok = p.dominance >= GATE_DOMINANCE && p.crosstalk < GATE_CROSSTALK
        ok || (recovery_ok = false)
        push!(gate_rows, (; estimator = e, noise = t, dominance = p.dominance,
                          crosstalk = p.crosstalk, ok))
    end

    # ── 9. Artifacts ─────────────────────────────────────────────────────────
    open(basis_csv, "w") do io
        println(io, "estimator,noise,alpha,m,tau,Kbar,Kvar,Kbar_h1,Kbar_h2,R,a_true,sigma1")
        for e in ESTIMATORS, t in TYPES, k in 1:nk
            @printf(io, "%s,%s,%d,%d,%.6e,%.10e,%.10e,%.10e,%.10e,%d,%.10e,%.6e\n",
                    e, t, TYPE_ALPHA[t], ms[k], taus[k],
                    Kbar[(e, t)][k], Kvar[(e, t)][k], Kh1[(e, t)][k], Kh2[(e, t)][k],
                    R_BASIS, a_true[t], SIGMA1[t])
        end
    end

    open(trials_csv, "w") do io
        println(io, "estimator,weights,trial,", join(("a_$t" for t in TYPES), ","),
                ",rss,active")
        for r in trial_rows
            @printf(io, "%s,%s,%d,%s,%.10e,%s\n", r.estimator, r.weights, r.trial,
                    join((@sprintf("%.10e", v) for v in r.a), ","), r.rss, r.active)
        end
    end

    open(summary_csv, "w") do io
        println(io, "estimator,weights,noise,metric,value")
        for e in ESTIMATORS
            R = results[e]
            for wname in (:empirical, :edf)
                f = R.fits[wname]
                @printf(io, "%s,%s,-,kappa2,%.10e\n", e, wname, f.kappa2)
                @printf(io, "%s,%s,-,rss_median,%.10e\n", e, wname, f.rss)
                for t in TYPES
                    @printf(io, "%s,%s,%s,relerr_median,%.10e\n", e, wname, t, median(f.relerr[t]))
                    @printf(io, "%s,%s,%s,relerr_iqr,%.10e\n",    e, wname, t, iqr(f.relerr[t]))
                    @printf(io, "%s,%s,%s,spread_iqr_over_median,%.10e\n", e, wname, t, f.spread[t])
                    @printf(io, "%s,%s,%s,zero_rate,%.10e\n",      e, wname, t, f.zero_rate[t])
                    @printf(io, "%s,%s,%s,misid,%.10e\n",          e, wname, t, f.misid[t])
                    @printf(io, "%s,%s,%s,misid_pointwise,%.10e\n", e, wname, t, f.misid_pointwise[t])
                    @printf(io, "%s,%s,%s,dominant_fraction,%.10e\n", e, wname, t, f.dominant_fraction[t])
                end
            end
            @printf(io, "%s,basis,RRFM,length_ratio_N_over_Nq,%.10e\n", e, len_ratio[e])
            for t in TYPES
                @printf(io, "%s,chi2,%s,relerr_median,%.10e\n", e, t, median(R.chi2_relerr[t]))
                @printf(io, "%s,chi2,%s,misid,%.10e\n", e, t, R.chi2_mis[t])
                @printf(io, "%s,chi2,%s,misid_pointwise,%.10e\n", e, t, R.chi2_misp[t])
                @printf(io, "%s,basis,%s,split_half_max,%.10e\n", e, t, R.split[t].max)
                @printf(io, "%s,basis,%s,split_half_last,%.10e\n", e, t, R.split[t].last)
                @printf(io, "%s,basis,%s,column_slope,%.10e\n", e, t, R.col_slope[t])
                @printf(io, "%s,pure,%s,dominance,%.10e\n", e, t, R.pure[t].dominance)
                @printf(io, "%s,pure,%s,crosstalk,%.10e\n", e, t, R.pure[t].crosstalk)
                @printf(io, "%s,pure,%s,relerr_median,%.10e\n", e, t, R.pure[t].relerr)
            end
        end
    end

    runtime = time() - t_start
    open(json_path, "w") do io
        println(io, "{")
        println(io, "  \"meta\": {")
        println(io, "    \"git_sha\": ", _jval(git_sha), ",")
        println(io, "    \"timestamp\": ", _jval(stamp), ",")
        println(io, "    \"julia_version\": ", _jval(string(VERSION)), ",")
        println(io, "    \"nthreads\": ", Threads.nthreads(), ",")
        println(io, "    \"quick\": ", _jval(QUICK), ",")
        println(io, "    \"master_seed\": ", MASTER_SEED, ",")
        println(io, "    \"seed_scheme\": ", _jval("Xoshiro(splitmix64(seed, tag, type, r))"), ",")
        println(io, "    \"N\": ", N_SAMPLES, ", \"tau0\": ", _jval(TAU0), ",")
        println(io, "    \"m_grid\": ", _jval(ms), ",")
        println(io, "    \"R_basis\": ", R_BASIS, ", \"R_trial\": ", R_TRIAL,
                ", \"R_pure\": ", R_PURE, ", \"n_chi2\": ", N_CHI2, ",")
        println(io, "    \"min_nsubs\": ", MIN_NSUBS, ",")
        println(io, "    \"crossovers\": ", _jval(collect(CROSSOVERS)), ",")
        println(io, "    \"dominance_margin\": ", _jval(DOMINANCE_MARGIN), ",")
        println(io, "    \"sigma1\": ", _jobj((t, SIGMA1[t]) for t in TYPES), ",")
        println(io, "    \"scale_c\": ", _jobj((t, c[t]) for t in TYPES), ",")
        println(io, "    \"a_true\": ", _jobj((t, a_true[t]) for t in TYPES), ",")
        println(io, "    \"edf_symbol_note\": ",
                _jval("alpha=-4 borrows the RWFM EDF model; shipped tables stop at alpha=-2"), ",")
        println(io, "    \"edf_clamped_points\": ", edf_clamped, ",")
        println(io, "    \"runtime_s\": ", _jobj([("basis", t_basis), ("trials", t_trial),
                                                 ("pure", t_pure), ("total", runtime)]))
        println(io, "  },")
        println(io, "  \"slope_check\": {")
        println(io, "    \"passed\": ", _jval(slope_ok), ", \"tol\": ", _jval(SLOPE_TOL), ",")
        println(io, "    \"mhdev_slope\": ", _jobj((t, slopes[t]) for t in TYPES), ",")
        println(io, "    \"expected\": ",
                _jobj((t, _mhdev_slope(TYPE_ALPHA[t])) for t in TYPES), ",")
        println(io, "    \"fft_rrfm\": ",
                _jobj([("mhdev_slope_interior", fft_slope),
                       ("longtau_ratio_vs_power_law", fft_flat),
                       ("note", "FFT shaper zeroes DC on a periodic record; the cumsum-of-cumsum generator is used instead")]))
        println(io, "  },")
        println(io, "  \"rrfm_diagnostics\": ",
                _jobj([("N_reference", N_SAMPLES), ("N_short", Nq), ("R", Rd),
                       ("mhtotdev_column_slope_no_drift_removal", nodrift_slope),
                       ("note", "K(N)/K(N/4) ~ 1 for a convergent column; MVAR diverges for alpha <= -3 so MDEV's scales with T")]), ",")
    println(io, "  \"solver_crosscheck\": ",
                _jobj([("max_rel_objective_diff", solver_obj),
                       ("max_rel_coefficient_diff", solver_coef),
                       ("tol", GATE_SOLVER_TOL), ("passed", solver_ok)]), ",")
        println(io, "  \"estimators\": {")
        for (ei, e) in enumerate(ESTIMATORS)
            R = results[e]
            comma = ei < length(ESTIMATORS) ? "," : ""
            println(io, "    \"$e\": {")
            for wname in (:empirical, :edf)
                f = R.fits[wname]
                println(io, "      \"$wname\": {")
                println(io, "        \"kappa2\": ", _jval(f.kappa2),
                        ", \"rss_median\": ", _jval(f.rss), ",")
                println(io, "        \"misid\": ",
                        _jobj((t, f.misid[t]) for t in TYPES), ",")
                println(io, "        \"misid_pointwise\": ",
                        _jobj((t, f.misid_pointwise[t]) for t in TYPES), ",")
                println(io, "        \"dominant_fraction\": ",
                        _jobj((t, f.dominant_fraction[t]) for t in TYPES), ",")
                println(io, "        \"relerr_median\": ",
                        _jobj((t, median(f.relerr[t])) for t in TYPES), ",")
                println(io, "        \"relerr_iqr\": ",
                        _jobj((t, iqr(f.relerr[t])) for t in TYPES), ",")
                println(io, "        \"spread_iqr_over_median\": ",
                        _jobj((t, f.spread[t]) for t in TYPES), ",")
                println(io, "        \"zero_rate\": ",
                        _jobj((t, f.zero_rate[t]) for t in TYPES))
                println(io, "      },")
            end
            println(io, "      \"chi2_bootstrap\": {")
            println(io, "        \"misid\": ", _jobj((t, R.chi2_mis[t]) for t in TYPES), ",")
            println(io, "        \"misid_pointwise\": ",
                    _jobj((t, R.chi2_misp[t]) for t in TYPES), ",")
            println(io, "        \"relerr_median\": ",
                    _jobj((t, median(R.chi2_relerr[t])) for t in TYPES), ",")
            println(io, "        \"edf\": ", _jval(R.nu))
            println(io, "      },")
            println(io, "      \"truth_regions_tau\": ",
                    _jobj((t, taus[R.regions[t]]) for t in TYPES), ",")
            println(io, "      \"basis\": {")
            println(io, "        \"split_half_max\": ",
                    _jobj((t, R.split[t].max) for t in TYPES), ",")
            println(io, "        \"split_half_last_tau\": ",
                    _jobj((t, R.split[t].last) for t in TYPES), ",")
            println(io, "        \"column_slope\": ",
                    _jobj((t, R.col_slope[t]) for t in TYPES), ",")
            println(io, "        \"rrfm_length_ratio_N_over_Nquarter\": ",
                    _jval(len_ratio[e]))
            println(io, "      },")
            println(io, "      \"pure_noise\": {")
            println(io, "        \"dominance\": ",
                    _jobj((t, R.pure[t].dominance) for t in TYPES), ",")
            println(io, "        \"crosstalk\": ",
                    _jobj((t, R.pure[t].crosstalk) for t in TYPES), ",")
            println(io, "        \"relerr_median\": ",
                    _jobj((t, R.pure[t].relerr) for t in TYPES))
            println(io, "      }")
            println(io, "    }", comma)
        end
        println(io, "  },")
        println(io, "  \"gates\": {")
        println(io, "    \"slope_check\": ", _jval(slope_ok), ",")
        println(io, "    \"solver_agreement\": ", _jval(solver_ok), ",")
        println(io, "    \"pure_noise_recovery\": ", _jval(recovery_ok), ",")
        println(io, "    \"dominance_floor\": ", _jval(GATE_DOMINANCE),
                ", \"crosstalk_ceiling\": ", _jval(GATE_CROSSTALK), ",")
        println(io, "    \"rows\": [")
        for (i, g) in enumerate(gate_rows)
            println(io, "      ", _jobj([("estimator", g.estimator), ("noise", g.noise),
                                         ("dominance", g.dominance),
                                         ("crosstalk", g.crosstalk), ("ok", g.ok)]),
                    i < length(gate_rows) ? "," : "")
        end
        println(io, "    ]")
        println(io, "  }")
        println(io, "}")
    end

    # ── 10. Console summary ──────────────────────────────────────────────────
    println("\n── Median relative error of â (empirical weights) ──")
    println("estimator   ", join((rpad(string(t), 12) for t in TYPES)), " κ₂")
    for e in ESTIMATORS
        f = results[e].fits[:empirical]
        @printf("%-11s %s %.3g\n", e,
                join((rpad(@sprintf("%.3g", median(f.relerr[t])), 12) for t in TYPES)),
                f.kappa2)
    end
    println("\n── Per-type misidentification (empirical weights): per-trial / pointwise ──")
    println("scored only where the analytic design makes the type dominant by ≥$(DOMINANCE_MARGIN)× in variance")
    println("estimator   ", join((rpad(string(t), 14) for t in TYPES)))
    for e in ESTIMATORS
        f = results[e].fits[:empirical]
        @printf("%-11s %s\n", e,
                join((rpad(@sprintf("%.2f/%.2f", f.misid[t], f.misid_pointwise[t]), 14)
                      for t in TYPES)))
    end
    println("\n── Fraction of the τ grid each type dominates in the fit (masking check) ──")
    println("estimator   ", join((rpad(string(t), 12) for t in TYPES)))
    for e in ESTIMATORS
        f = results[e].fits[:empirical]
        @printf("%-11s %s\n", e,
                join((rpad(@sprintf("%.2f", f.dominant_fraction[t]), 12) for t in TYPES)))
    end
    println("  (a type that is never dominant has its relative error masked — the fit is")
    println("   free to absorb its amplitude into a neighbouring column with no visible cost)")
    println("\n── RRFM basis-column stability (split-half |Δ|/mean) and fitted slope ──")
    for e in ESTIMATORS
        R = results[e]
        @printf("%-11s split-half max %.3g   column slope %.3f (ideal %+.1f)   K(N)/K(N/4) %.2f\n",
                e, R.split[:RRFM].max, R.col_slope[:RRFM],
                Float64(-TYPE_ALPHA[:RRFM] - 1), len_ratio[e])
    end
    println("\n── Pure-noise recovery gate (median over trials) ──")
    for g in gate_rows
        @printf("  %-9s %-5s dominance %.3f  crosstalk %.3f  %s\n",
                g.estimator, g.noise, g.dominance, g.crosstalk, g.ok ? "ok" : "FAIL")
    end
    @printf("\nSolver cross-check: max rel objective diff %.3g, max rel coefficient diff %.3g\n",
            solver_obj, solver_coef)
    println("\nWrote:\n  $basis_csv\n  $trials_csv\n  $summary_csv\n  $json_path")
    @printf("Total runtime: %.1f s\n", runtime)

    # Artifacts land first so a failed run is still inspectable, then fail loudly.
    failures = String[]
    slope_ok    || push!(failures, "pure-noise MHDEV slope check")
    solver_ok   || push!(failures, "exact-vs-Lawson-Hanson agreement " *
                                   "(max rel objective diff $(solver_obj) > $GATE_SOLVER_TOL)")
    recovery_ok || push!(failures, "pure-noise recovery for " *
                                   join(string.(unique(g.estimator for g in gate_rows if !g.ok)), ", "))
    isempty(failures) || error("GATE FAILURE: " * join(failures, "; ") *
                               ". Artifacts were written for inspection but must not be used.")
    return nothing
end

main()

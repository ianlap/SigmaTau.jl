# Ensembles of three or more clocks measured against one another: the
# N-cornered hat and the Groslambert covariance, each separating the pairwise
# measurements into one deviation per clock, with confidence intervals.

# Measured pair (i, j) from a possibly part-filled matrix; the transposed
# entry stands in for it (s²(x_i − x_j) = s²(x_j − x_i)).
function _nch_pair(pairwise::AbstractMatrix{StabilityResult}, i::Int, j::Int)
    isassigned(pairwise, i, j) && return pairwise[i, j]
    isassigned(pairwise, j, i) && return pairwise[j, i]
    throw(ArgumentError("ncornered_hat: no pairwise result for clocks ($i, $j); fill pairwise[$i, $j] (or its transpose)"))
end

"""
    ncornered_hat(pairwise::AbstractMatrix{StabilityResult}) → Vector{StabilityResult}
    ncornered_hat(r_ab, r_bc, r_ca) → Vector{StabilityResult}

N-cornered-hat noise separation (the classic Gray–Allan three-cornered
hat at `N = 3`; see NIST SP1065 §10.7 [riley-2008-sp1065](@cite)).
Given the pairwise deviations between every pair of `N ≥ 3` clocks,
estimates each clock's *individual* deviation under the assumption that
the clock noises are statistically independent, so the measured pairwise
variances satisfy `s²ᵢⱼ(τ) = σ²ᵢ(τ) + σ²ⱼ(τ)`.

With `Rᵢ = Σ_{j≠i} s²ᵢⱼ` (clock `i`'s row sum) and `S = Σ_{j<k} s²ⱼₖ`
(the sum over all pairs), the separation is

```math
\\hat\\sigma^2_i(\\tau) = \\frac{1}{N-2}\\left( R_i - \\frac{S}{N-1} \\right)
```

At `N = 3` this is the classic three-cornered hat,
`σ̂²_A = ½(s²_AB + s²_AC − s²_BC)`.

`pairwise[i, j]` (for `i < j`; the transposed entry is also accepted, so
fill either triangle, or both) is the deviation computed on the
difference record between clocks `i` and `j`. All entries must share the
same `deviation_type` and the same τ grid; otherwise an `ArgumentError`
is thrown. The 3-argument convenience form takes the three pairwise
results in the cyclic order A−B, B−C, C−A and returns `[A, B, C]`
(the orientation of each difference record does not matter).

Returned conventions (one `StabilityResult` per clock, in clock order):

- `deviation_type` is kept from the inputs.
- Non-positive variance estimates (common when one clock dominates or
  the independence assumption fails) are reported as `NaN` in `dev`
  rather than clamped to zero.
- When every input carries `edf` (computed with `ci=true`), `edf` and `ci`
  are filled at the inputs' confidence level. Each pairwise variance is
  χ² with its own `edf`; two pairs that share clock `c` covary by
  `2σ_c⁴/√(ν₁ν₂)`. Propagating these through the separation gives the
  variance of each `σ̂²ᵢ`, reported as the equivalent `edf = 2σ̂⁴ᵢ / Var(σ̂²ᵢ)`
  (the covariance model of [vernotte-2019-tch-gcov-uncertainty](@cite)).
  The intervals widen sharply for a clock much quieter than the others.
- `noise_type` is empty; `neff` is the elementwise minimum over the inputs.

See also [`groslambert`](@ref), which needs the difference records rather
than their deviations, never goes negative on average, and tolerates noise
common to two measurement channels.

# Examples

```jldoctest
julia> using SigmaTau

julia> x1, x2, x3 = randn(256), randn(256), randn(256);

julia> r(a, b) = adev(PhaseData(a .- b, 1.0), [1, 2, 4]);

julia> clocks = ncornered_hat(r(x1, x2), r(x2, x3), r(x3, x1));

julia> length(clocks), length(clocks[1].ci)
(3, 3)
```
"""
function ncornered_hat(pairwise::AbstractMatrix{StabilityResult})
    n_clocks, n_cols = size(pairwise)
    n_clocks == n_cols || throw(ArgumentError(
        "ncornered_hat: pairwise matrix must be square, got $(n_clocks)×$(n_cols)"))
    n_clocks >= 3 || throw(ArgumentError(
        "ncornered_hat: need at least 3 clocks for a cornered hat, got $n_clocks"))

    ref = _nch_pair(pairwise, 1, 2)
    n_tau = length(ref.tau)
    pairs = [(i, j) for i in 1:(n_clocks - 1) for j in (i + 1):n_clocks]
    s2 = zeros(length(pairs), n_tau)            # pairwise variances
    ν = fill(NaN, length(pairs), n_tau)         # their edf, if every input has one
    neff_min = fill(typemax(Int), n_tau)
    for (p, (i, j)) in enumerate(pairs)
        r = _nch_pair(pairwise, i, j)
        r.deviation_type === ref.deviation_type || throw(ArgumentError(
            "ncornered_hat: pairwise ($i, $j) is :$(r.deviation_type) but (1, 2) is :$(ref.deviation_type); all pairs must share one deviation type"))
        r.tau == ref.tau || throw(ArgumentError(
            "ncornered_hat: pairwise ($i, $j) was computed on a different τ grid; all pairs must share one grid"))
        s2[p, :] .= r.dev .^ 2
        length(r.edf) == n_tau && (ν[p, :] .= r.edf)
        neff_min .= min.(neff_min, r.neff)
    end

    # σ̂²ᵢ = Σ_pairs wᵢ(p) s²_p with wᵢ(p) = (1[i ∈ p] − 1/(N−1)) / (N−2).
    w(i, (j, k)) = ((i == j || i == k) - 1 / (n_clocks - 1)) / (n_clocks - 2)
    var = [sum(w(i, pairs[p]) * s2[p, k] for p in eachindex(pairs)) for i in 1:n_clocks, k in 1:n_tau]
    with_ci = !any(isnan, ν)
    # Pairs covary only through a shared clock: list, for each pair, the
    # other pairs that share one with it and that clock.
    nbrs = [[(q, only(intersect(pairs[p], pairs[q]))) for q in eachindex(pairs)
             if q != p && !isempty(intersect(pairs[p], pairs[q]))] for p in eachindex(pairs)]
    confidence = isnan(ref.confidence) ? DEFAULT_CONFIDENCE : ref.confidence
    results = Vector{StabilityResult}(undef, n_clocks)
    for i in 1:n_clocks
        dev = [v > 0 ? sqrt(v) : NaN for v in var[i, :]]
        edf = Float64[]
        if with_ci
            edf = [_nch_edf(i, k, pairs, nbrs, s2, ν, var, w) for k in 1:n_tau]
        end
        results[i] = StabilityResult(ref.deviation_type, copy(ref.tau), dev, Symbol[],
                                     with_ci ? _zip_ci(confidence_intervals(dev, edf, confidence)...) : _CIBound[], edf,
                                     copy(neff_min), ref.tau0, ref.n, with_ci ? confidence : NaN)
    end
    return results
end

# Equivalent edf of clock i's separated variance at τ index k.
function _nch_edf(i, k, pairs, nbrs, s2, ν, var, w)
    v = var[i, k]
    v > 0 || return NaN
    total = 0.0
    for p in eachindex(pairs)
        wp = w(i, pairs[p])
        total += wp^2 * 2 * s2[p, k]^2 / ν[p, k]
        for (q, c) in nbrs[p]
            σ2 = max(var[c, k], 0.0)
            total += wp * w(i, pairs[q]) * 2 * σ2^2 / sqrt(ν[p, k] * ν[q, k])
        end
    end
    return 2 * v^2 / total
end

# Three-cornered hat from A−B, B−C, C−A; C−A fills pairwise[1, 3] since
# orientation is irrelevant.
function ncornered_hat(r_ab::StabilityResult, r_bc::StabilityResult, r_ca::StabilityResult)
    pairwise = Matrix{StabilityResult}(undef, 3, 3)
    pairwise[1, 2] = r_ab
    pairwise[2, 3] = r_bc
    pairwise[1, 3] = r_ca
    return ncornered_hat(pairwise)
end

"""
    groslambert(pairwise::AbstractMatrix{<:AbstractTimingData}, m_values=Octave;
                confidence=DEFAULT_CONFIDENCE, interval=:chi2) → Vector{StabilityResult}
    groslambert(x_ab, x_bc, x_ca, m_values=Octave; confidence, interval) → Vector{StabilityResult}

Groslambert covariance (Vernotte et al., *IEEE Trans. UFFC* 63, 2016): each
clock's Allan deviation from the cross Allan variance of two difference
records that share it. For clock `i` and two others `j`, `k`, the second
differences of `x_i − x_j` and `x_i − x_k` have only clock `i` in common, so

    σ²_i(τ) = ⟨Δ²(x_i − x_j) · Δ²(x_i − x_k)⟩ / (2τ²),

averaged over every pair `(j, k)` of the other clocks. Unlike the
[`ncornered_hat`](@ref) it does not subtract variances, so its estimates
scatter about the true value rather than turning negative for a clock much
quieter than the rest, and noise shared by two measurement channels enters
only the clock both channels see.

`pairwise[i, j]` is the difference record `x_i − x_j` (`PhaseData` or
`FrequencyData`, same τ₀ and length throughout); the transposed entry,
negated, stands in where one is missing. The 3-argument form takes the
records A−B, B−C, C−A and returns `[A, B, C]`.

Each result is an `:adev` `StabilityResult` with an EDF from the covariance
model of [vernotte-2019-tch-gcov-uncertainty](@cite): each record's Allan
variance carries its EDF for the identified noise type
[greenhall-2003-edf-stability](@cite) (as in [`adev`](@ref)), the products of
Gaussian second differences give the variance of every covariance estimate in
closed form, and `edf = 2σ̂⁴/Var(σ̂²)` sets a χ² interval `ci` (none below
1 EDF, where the χ² model degenerates). A non-positive covariance gives a
`NaN` estimate.

The χ² interval holds for a clock that dominates its two records. For a
clock quieter than the others it is too narrow once the clock is no longer
resolved, with its lower bound above the truth: in a Monte Carlo with one
clock at a third of the others' level, its 95 % interval covered that
clock's true value 94 % of the time at 10 τ₀ and 60 % at 100 τ₀.

With three clocks, `interval = :klts` gives the KLTS interval of
[lantz-2019-klts](@cite) instead: the quantiles of the posterior of the three
clock variances given the sample covariance of the three records, under a
log-uniform prior from 10⁻⁵ times the smallest to 10³ times the largest
estimate, with the measurement noise read from the closure of the records
and the geometric mean of their EDF. It is given also where the estimate is
negative, and a lower bound of 0 means the data do not exclude zero (the
bound would otherwise follow the prior's lower limit). In the same Monte
Carlo it covered 92 % at 10 τ₀, 86 % at 20 τ₀, where the clock is barely
resolved and the upper bound follows the prior, and 99–100 % from 50 τ₀ on.
Each τ costs about 0.1 s.
"""
function groslambert(pairwise::AbstractMatrix{<:AbstractTimingData}, m_values::Vector{Int};
                     confidence::Real=DEFAULT_CONFIDENCE, interval::Symbol=:chi2)
    n_clocks = size(pairwise, 1)
    size(pairwise, 2) == n_clocks ||
        throw(ArgumentError("groslambert: pairwise matrix must be square"))
    n_clocks >= 3 || throw(ArgumentError("groslambert: need at least 3 clocks, got $n_clocks"))
    interval in (:chi2, :klts) ||
        throw(ArgumentError("groslambert: interval must be :chi2 or :klts, got :$interval"))
    interval === :klts && n_clocks != 3 &&
        throw(ArgumentError("groslambert: interval = :klts is defined for 3 clocks, got $n_clocks"))
    x = Matrix{Vector{Float64}}(undef, n_clocks, n_clocks)
    ref = nothing
    for i in 1:n_clocks, j in 1:n_clocks
        i == j && continue
        d = isassigned(pairwise, i, j) ? pairwise[i, j] :
            isassigned(pairwise, j, i) ? pairwise[j, i] :
            throw(ArgumentError("groslambert: no record for clocks ($i, $j)"))
        ref === nothing && (ref = d)
        _cross_pair(ref, d)
        xi = _f64(_as_phase(d).x)
        x[i, j] = isassigned(pairwise, i, j) ? xi : -xi
    end
    τ₀ = Float64(ref.tau0); N = length(x[1, 2])
    all(m -> 1 <= m && 2m < N, m_values) ||
        throw(ArgumentError("groslambert: every m must satisfy 1 ≤ m < N/2 (N = $N)"))
    K = length(m_values)
    # Each record's Allan variance and its EDF for the identified noise type.
    s2 = Matrix{Vector{Float64}}(undef, n_clocks, n_clocks)
    ν = Matrix{Vector{Float64}}(undef, n_clocks, n_clocks)
    for i in 1:n_clocks, j in (i + 1):n_clocks
        r = adev(PhaseData(x[i, j], τ₀), m_values)
        s2[i, j] = s2[j, i] = r.dev .^ 2
        ν[i, j] = ν[j, i] = r.edf
    end
    c = Float64(confidence)
    klts = interval === :klts ? _gl_klts(x, s2, ν, m_values, τ₀, c) : nothing
    results = Vector{StabilityResult}(undef, n_clocks)
    for i in 1:n_clocks
        others = [j for j in 1:n_clocks if j != i]
        couples = [(j, k) for j in others for k in others if j < k]
        covar = zeros(K)
        for (q, m) in enumerate(m_values)
            p = zeros(N - 2m)
            for (j, k) in couples
                p .+= _second_diff_products(x[i, j], x[i, k], m)[1]
            end
            covar[q] = sum(p) / (length(p) * length(couples) * 2 * (m * τ₀)^2)
        end
        edf = [_gl_edf(i, q, couples, covar[q], s2, ν) for q in 1:K]
        dev = [v > 0 ? sqrt(v) : NaN for v in covar]
        ci = klts === nothing ? _gl_chi2_ci(dev, edf, c) : klts[i]
        results[i] = StabilityResult(:adev, m_values .* τ₀, dev, Symbol[], ci, edf,
                                     N .- 2 .* m_values, τ₀, N, c)
    end
    return results
end

# Equivalent EDF of clock i's Groslambert covariance at τ index q, the mean
# over `couples` of Ĝ_jk = ⟨z_ij z_ik⟩. For Gaussian records
# Cov(z_ij z_ik, z_il z_in) = E[z_ij z_il] E[z_ik z_in] + E[z_ij z_in] E[z_ik z_il],
# where E[z_ij z_il] is record (i, j)'s Allan variance for l = j and clock i's
# own for l ≠ j; each term is divided by the geometric mean of the EDF of the
# records involved (Vernotte & Lantz 2019, appendix B).
function _gl_edf(i, q, couples, v, s2, ν)
    v > 0 || return NaN
    e(j, l) = j == l ? s2[i, j][q] : v
    total = 0.0
    for (j, k) in couples, (l, n) in couples
        c = e(j, l) * e(k, n) + e(j, n) * e(k, l)
        total += c / sqrt(sqrt(ν[i, j][q] * ν[i, k][q] * ν[i, l][q] * ν[i, n][q]))
    end
    return 2 * v^2 * length(couples)^2 / total
end

# χ² interval on the equivalent EDF. Below 1 EDF the χ² model degenerates (its
# lower bound can pass the estimate), so no interval is given there.
function _gl_chi2_ci(dev, edf, confidence)
    lo, hi = confidence_intervals(dev, edf, confidence)
    for q in eachindex(edf)
        edf[q] >= 1 || (lo[q] = hi[q] = NaN)
    end
    return _zip_ci(lo, hi)
end

# KLTS intervals [lantz-2019-klts](@cite) for clocks 1, 2, 3 at each τ, from the
# sample covariance of the records 1−2, 2−3, 3−1, the measurement noise read
# from their closure (σ²_N = σ²_closure / 3, eq 9) and the geometric mean of
# the three records' EDF. A lower bound that depends on the prior is 0.
function _gl_klts(x, s2, ν, m_values, τ₀, confidence)
    rec = (x[1, 2], x[2, 3], x[3, 1])
    cls = rec[1] .+ rec[2] .+ rec[3]
    probs = ((1 - confidence) / 2, (1 + confidence) / 2)
    ci = [Vector{_CIBound}(undef, length(m_values)) for _ in 1:3]
    for (q, m) in enumerate(m_values)
        w = 2 * (m * τ₀)^2 * (length(cls) - 2m)
        S = [j == k ? s2[_gl_pair(j)...][q] : sum(_second_diff_products(rec[j], rec[k], m)[1]) / w
             for j in 1:3, k in 1:3]
        vN = sum(_second_diff_products(cls, cls, m)[1]) / w / 3
        n = cbrt(ν[1, 2][q] * ν[2, 3][q] * ν[1, 3][q])
        if !(isfinite(n) && all(isfinite, S))
            for i in 1:3; ci[i][q] = (lo = NaN, hi = NaN); end
            continue
        end
        vN <= 1e-9 * sum(S[k, k] for k in 1:3) / 3 && (S = S[1:2, 1:2])
        qs, zero_lower = _klts(S, vN, n, probs)
        for i in 1:3
            ci[i][q] = (lo = zero_lower[i] ? 0.0 : sqrt(qs[i, 1]), hi = sqrt(qs[i, 2]))
        end
    end
    return ci
end

# Clock indices of record j among (1−2, 2−3, 3−1).
_gl_pair(j) = ((1, 2), (2, 3), (3, 1))[j]

groslambert(pairwise::AbstractMatrix{<:AbstractTimingData}, taus::TauMode; kw...) =
    groslambert(pairwise, tau_values(taus, _gl_n(pairwise), :adev); kw...)
groslambert(pairwise::AbstractMatrix{<:AbstractTimingData}; kw...) =
    groslambert(pairwise, _default_m_values(_gl_n(pairwise), :adev); kw...)

# Phase samples per record: those of the first off-diagonal record.
function _gl_n(pairwise)
    for I in CartesianIndices(pairwise)
        I[1] != I[2] && isassigned(pairwise, I) && return length(_as_phase(pairwise[I]).x)
    end
    throw(ArgumentError("groslambert: the pairwise matrix is empty"))
end

# Three clocks from A−B, B−C, C−A; C−A fills pairwise[3, 1].
function groslambert(x_ab::AbstractTimingData, x_bc::AbstractTimingData,
                     x_ca::AbstractTimingData, args...; kw...)
    pairwise = Matrix{AbstractTimingData}(undef, 3, 3)
    pairwise[1, 2] = x_ab
    pairwise[2, 3] = x_bc
    pairwise[3, 1] = x_ca
    return groslambert(pairwise, args...; kw...)
end

# Three or more clocks measured against one another: the N-cornered hat and
# the Groslambert covariance separate the pairwise measurements into one
# deviation per clock, with a confidence interval.

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
    ncornered_hat(cc::ClockComparisons, m_values=Octave; deviation=:adev, reference=true,
                  confidence=DEFAULT_CONFIDENCE) → Vector{StabilityResult}

Deviation of each of ``N ≥ 3`` clocks from the deviations measured between
every pair of them: the three-cornered hat of
[Gray and Allan (1974)](@cite gray-1974-three-cornered-hat), and its
extension to ``N`` clocks.

For independent clocks the variance measured between clocks ``i`` and
``j`` is ``s_{ij}^2(τ) = σ_i^2(τ) + σ_j^2(τ)``, and the variance of clock
``i`` follows from the ``N(N-1)/2`` measured ones as

```math
\\hat σ_i^2(τ) = \\frac{1}{N-2} \\left( \\sum_{j \\ne i} s_{ij}^2
    - \\frac{1}{N-1} \\sum_{j<k} s_{jk}^2 \\right)
```

([SP1065 §10.14, Eqs. (71) and (72)](@cite riley-2008-sp1065)). For three
clocks this is
``\\hat σ_a^2 = \\tfrac{1}{2}(s_{ab}^2 + s_{ac}^2 - s_{bc}^2)``.

The estimate can be negative for a clock much more stable than the others,
for clocks that are correlated, and where few terms are averaged, most
often at long ``τ``; `dev` is `NaN` there. Noise of the instrument that
measured a pair counts as noise of its two clocks.
[`groslambert`](@ref) keeps it out, given the difference records.

# Arguments

- `pairwise`: an ``N × N`` matrix whose entry `[i, j]` is the result of an
  estimator on the difference record of clocks ``i`` and ``j``. One entry
  per pair is enough, in either triangle. All entries must come from one
  estimator on one ``τ`` grid. The estimator must be one whose variance
  adds for independent clocks, such as [`adev`](@ref), [`mdev`](@ref) or
  [`hdev`](@ref).
- `r_ab`, `r_bc`, `r_ca`: for three clocks, the results for A−B, B−C and
  C−A. The sign of each difference does not matter. The clocks are
  returned in the order A, B, C.
- `cc`: clocks read against one reference at the same instants. The
  estimator is run on the difference record of every pair, in which the
  reference cancels, and the reference is one more clock, whose record
  against clock ``i`` is the reading of clock ``i``. Two clocks and
  their reference make a three-cornered hat. The clocks are returned in
  the order of `cc`, with the reference last. Where `cc` has a column of
  zeros, that clock is the reference and no clock is added.
- `m_values`: with `cc`, averaging factors or a [`TauMode`](@ref) grid.

# Keywords

These go with `cc`.

- `deviation=:adev`: the estimator, one of `:adev`, `:mdev`, `:hdev`,
  `:mhdev`, `:pdev`, `:tdev` and `:htdev`.
- `reference=true`: whether the reference counts as a clock. With `false`
  the hat is that of the clocks' differences alone, three or more, and a
  noisy reference stays out of it.
- `confidence=0.683`: confidence level of the intervals.

# Returns

One [`StabilityResult`](@ref) per clock, of the estimator of the inputs.

- `edf` and `ci` are filled when every input has an EDF, at the confidence
  level of the inputs. The variance of ``\\hat σ_i^2`` follows from that of
  each pairwise variance, ``2s_{ij}^4/ν_{ij}``, and from the covariance
  ``2σ_c^4/\\sqrt{ν_1 ν_2}`` of two pairs that share clock ``c``
  ([Vernotte and Lantz (2019), Appendix B](@cite vernotte-2019-tch-gcov-uncertainty),
  with the EDF ``ν`` of each record in place of their number of independent
  measurements). `edf` is the equivalent
  ``2\\hat σ_i^4/\\mathrm{Var}(\\hat σ_i^2)``, and `ci` its ``χ^2`` interval.
  For three clocks with a common EDF this is the loss of degrees of freedom
  of [SP1065, Eq. (70)](@cite riley-2008-sp1065). The interval is too
  narrow for a clock that is no longer resolved under the others (see
  [`groslambert`](@ref)), and none is given below 1 EDF: the limits are
  `NaN` there.
- `noise_type` is empty, and `neff` is the smallest of the inputs' at each
  ``τ``.

# Examples

```jldoctest
julia> using SigmaTau

julia> x1, x2, x3 = randn(256), randn(256), randn(256);

julia> r(a, b) = adev(PhaseData(a .- b, 1.0), [1, 2, 4]);

julia> clocks = ncornered_hat(r(x1, x2), r(x2, x3), r(x3, x1));

julia> length(clocks), length(clocks[1].ci)
(3, 3)
```

Two clocks read against a third, the reference, which comes last:

```jldoctest
julia> using SigmaTau

julia> cc = ClockComparisons(randn(256, 2), 1.0);

julia> length(ncornered_hat(cc, [1, 2, 4]; deviation = :mdev))
3
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
                                     with_ci ? _sep_chi2_ci(dev, edf, confidence) : _CIBound[], edf,
                                     copy(neff_min), ref.tau0, ref.n, with_ci ? confidence : NaN)
    end
    return results
end

# Equivalent EDF of clock i's variance at τ index k, 2σ̂⁴/Var(σ̂²). Var(σ̂²)
# sums w(p) w(q) Cov(s²_p, s²_q) over the pairs, with Var(s²_p) = 2s⁴_p/ν_p
# and, for two pairs that share clock c, Cov = 2σ_c⁴/√(ν_p ν_q) (Vernotte and
# Lantz (2019), Appendix B, with the EDF in place of the number of
# independent measurements).
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

# Estimators a hat can separate: quadratic in the record with no factor that
# depends on its noise type, so the variance of a difference of independent
# clocks is the sum of theirs.
const _HAT_DEVIATIONS = (:adev, :mdev, :hdev, :mhdev, :pdev, :tdev, :htdev)

# The difference records clock i − clock j of `cc`, i < j, in the upper
# triangle. With `reference`, the reference of the readings is one more clock,
# the last, whose record against clock i is the reading of clock i; a clock
# with a column of zeros is the reference itself.
function _hat_records(cc::ClockComparisons, reference::Bool, name::String)
    k = size(cc.x, 2)
    n = k + (reference && !any(iszero, eachcol(cc.x)))
    n >= 3 || throw(ArgumentError(
        "$name: need at least 3 clocks, got $k without the reference"))
    records = Matrix{PhaseData{Float64}}(undef, n, n)
    for i in 1:n, j in (i + 1):n
        records[i, j] = j > k ? cc[i] : cc[i, j]
    end
    return records
end

function ncornered_hat(cc::ClockComparisons, m_values::Vector{Int}; deviation::Symbol=:adev,
                       reference::Bool=true, confidence::Real=DEFAULT_CONFIDENCE)
    deviation in _HAT_DEVIATIONS || throw(ArgumentError(
        "ncornered_hat: deviation must be one of $(join(_HAT_DEVIATIONS, ", ")), got :$deviation"))
    records = _hat_records(cc, reference, "ncornered_hat")
    pairwise = Matrix{StabilityResult}(undef, size(records))
    c = Float64(confidence)
    for i in axes(records, 1), j in (i + 1):size(records, 2)
        pairwise[i, j] = SigmaTau.deviation(deviation, records[i, j], m_values; confidence=c)
    end
    return ncornered_hat(pairwise)
end

ncornered_hat(cc::ClockComparisons, taus::TauMode=Octave; deviation::Symbol=:adev, kw...) =
    ncornered_hat(cc, tau_values(taus, size(cc.x, 1), deviation); deviation, kw...)

"""
    groslambert(pairwise::AbstractMatrix{<:AbstractTimingData}, m_values=Octave;
                confidence=DEFAULT_CONFIDENCE, interval=:chi2) → Vector{StabilityResult}
    groslambert(x_ab, x_bc, x_ca, m_values=Octave; kwargs...) → Vector{StabilityResult}
    groslambert(cc::ClockComparisons, m_values=Octave; reference=true,
                kwargs...) → Vector{StabilityResult}

Allan deviation of each of ``N ≥ 3`` clocks from the difference records
between them, by the Groslambert covariance: the cross Allan variance of
two records that share the clock
([Vernotte et al. (2016), Eq. (7)](@cite vernotte-2016-gcov), after the
two-sample cross-variance of
[Groslambert et al. (1981), Eqs. (12) to (16)](@cite groslambert-1981-crosscorrelations)).

The second differences of ``x_i - x_j`` and ``x_i - x_k`` have only clock
``i`` in common, so for independent clocks

```math
\\hat σ_i^2(τ) = \\frac{1}{2τ^2} \\left\\langle Δ^2(x_i - x_j)\\; Δ^2(x_i - x_k) \\right\\rangle ,
\\qquad Δ^2 x_n = x_{n+2m} - 2x_{n+m} + x_n ,
```

averaged over the record and, with more than three clocks, over every pair
``(j, k)`` of the others.

On records that close, ``(x_i - x_j) + (x_j - x_k) + (x_k - x_i) = 0``, the
estimate equals that of [`ncornered_hat`](@ref) on the Allan variances of
the records. The two differ when each record carries noise of its own,
such as that of the counter that measured it: that noise biases the hat and
averages out of the covariance
([Groslambert et al. (1981), Eqs. (9) and (17) to (19)](@cite groslambert-1981-crosscorrelations);
[Vernotte et al. (2016), Eqs. (13) and (15)](@cite vernotte-2016-gcov)).
Like the hat, the estimate can be negative for a clock much more stable
than the others; `dev` is `NaN` there.

A phase record may have missing samples (`NaN`). A term that reads one is
left out, of the products as of the Allan variance of each record
([SP1065 §10.4](@cite riley-2008-sp1065)), and the average runs over the
terms that remain, of all pairs ``(j, k)`` together. The handbook gives the
rule for the Allan variance of one record; its use on the products is the
package's own. With three clocks a sample missing in one record costs every
clock its terms there. In simulation (three clocks, 20 000 samples; a clock
unread for a tenth or for half of the record, a clock that joins late, a
tenth of all readings missing at random), the mean of the estimates of a
resolved clock stayed within 2 % of that without gaps up to ``80τ_0``.

# Arguments

- `pairwise`: an ``N × N`` matrix whose entry `[i, j]` is the record
  ``x_i - x_j``. One entry per pair is enough, in either triangle. The
  records are all `PhaseData` or all `FrequencyData`, with one sample
  interval and length. Missing samples are taken in phase records, with
  `interval=:chi2`.
- `x_ab`, `x_bc`, `x_ca`: for three clocks, the records A−B, B−C and C−A.
  The clocks are returned in the order A, B, C.
- `cc`: clocks read against one reference, as for [`ncornered_hat`](@ref):
  the records are the differences of every pair, and the reference is one
  more clock, returned last. These records close, so the estimate is that
  of the hat on [`adev`](@ref).
- `m_values`: averaging factors, or a [`TauMode`](@ref) grid.

# Keywords

- `reference=true`: with `cc`, whether the reference counts as a clock.

- `confidence=0.683`: confidence level of the interval.
- `interval=:chi2`: the interval reported, `:chi2` or, for three clocks,
  `:klts`.

# Returns

One [`StabilityResult`](@ref) of type `:adev` per clock, without noise
types. `neff` is the number of terms averaged for each pair ``(j, k)``.
`edf` is the equivalent
``2\\hat σ_i^4/\\mathrm{Var}(\\hat σ_i^2)``. The variance of the estimate is
that of a mean of products of Gaussian second differences
([Vernotte and Lantz (2019), Appendix B](@cite vernotte-2019-tch-gcov-uncertainty)).
Each record's Allan variance enters it, and the geometric mean of the
records' EDF for their identified noise types stands for the number of
independent measurements, a choice of this package. Where missing samples
cost terms, the EDF of each record is scaled to the terms that remain; in
the simulation above the 95 % interval then covered the true value of a
resolved clock in 90 to 100 % of the records, and in 93 to 96 % without
gaps.

# Confidence intervals

With `interval=:chi2` the interval is the ``χ^2`` interval on `edf`; none
is given below 1 EDF. It holds for a clock that dominates its records. For
a clock more stable than the others it becomes too narrow once the clock
is no longer resolved, and its lower bound then lies above the true value.

With `interval=:klts`, for three clocks, the interval is that of
[Lantz et al. (2019)](@cite lantz-2019-klts): quantiles of the posterior
distribution of the three clock variances, given the sample covariance of
the three records. The prior is uniform in the logarithm of each variance,
from ``10^{-5}`` times the smallest estimate to ``10^3`` times the largest.
The noise of the measuring instruments is estimated from the closure of
the records (their Eq. (9)), and the number of degrees of freedom is the
geometric mean of the records' EDF. A lower bound that would depend on the
lower limit of the prior is reported as 0: the record does not exclude
zero. The interval is given also where the estimate is negative. It is
computed on a grid in place of the random draws of the paper, in about
0.1 s for each ``τ``, and agrees within 2 % with the KLTS bounds of
[SigmaTheta](@cite vernotte-2025-sigmatheta) on the same inputs.

In simulation (4096 samples; a white FM clock beside two clocks two and
three times noisier, with instrument noise in each record; 400 records),
the 95 % ``χ^2`` interval covered the true value of the most stable clock
in 94 % of the records at ``10τ_0``. At ``100τ_0`` four records in ten
gave a negative estimate, and the interval covered the true value in 60 %
of the others. The KLTS interval covered it in 92 % of all records at
``10τ_0``, in 86 % at ``20τ_0``, where the clock is barely resolved and
the prior draws the upper bound down, and in 99 % or more from ``50τ_0``.

# Examples

```julia
ab, bc, ca = read_phase("ab.txt"), read_phase("bc.txt"), read_phase("ca.txt")
a, b, c = groslambert(ab, bc, ca; interval=:klts, confidence=0.95)
```
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
        _cross_pair(ref, d; gaps=d isa PhaseData && interval === :chi2)
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
    own = Matrix{Vector{Int}}(undef, n_clocks, n_clocks)        # terms of each record
    for i in 1:n_clocks, j in (i + 1):n_clocks
        r = adev(PhaseData(x[i, j], τ₀), m_values)
        s2[i, j] = s2[j, i] = r.dev .^ 2
        ν[i, j] = ν[j, i] = r.edf
        own[i, j] = own[j, i] = r.neff
    end
    c = Float64(confidence)
    klts = interval === :klts ? _gl_klts(x, s2, ν, m_values, τ₀, c) : nothing
    results = Vector{StabilityResult}(undef, n_clocks)
    for i in 1:n_clocks
        others = [j for j in 1:n_clocks if j != i]
        couples = [(j, k) for j in others for k in others if j < k]
        covar, terms = zeros(K), zeros(Int, K)
        for (q, m) in enumerate(m_values)
            p = zeros(N - 2m)
            n = 0
            for (j, k) in couples, (t, v) in enumerate(_second_diff_products(x[i, j], x[i, k], m)[1])
                isnan(v) || (p[t] += v; n += 1)     # NaN: a sample of the term is missing
            end
            covar[q] = sum(p) / (n * 2 * (m * τ₀)^2)
            terms[q] = n ÷ length(couples)
        end
        # Each record's EDF, for the terms the products kept.
        νi = [j == i ? Float64[] : ν[i, j] .* (terms ./ own[i, j]) for j in 1:n_clocks]
        edf = [_gl_edf(i, q, couples, covar[q], s2, νi) for q in 1:K]
        dev = [v > 0 ? sqrt(v) : NaN for v in covar]
        ci = klts === nothing ? _sep_chi2_ci(dev, edf, c) : klts[i]
        results[i] = StabilityResult(:adev, m_values .* τ₀, dev, Symbol[], ci, edf,
                                     terms, τ₀, N, c)
    end
    return results
end

# Equivalent EDF of clock i's Groslambert covariance at τ index q, the mean
# over `couples` of Ĝ_jk = ⟨z_ij z_ik⟩. For Gaussian records
# Cov(z_ij z_ik, z_il z_in) = E[z_ij z_il] E[z_ik z_in] + E[z_ij z_in] E[z_ik z_il],
# where E[z_ij z_il] is record (i, j)'s Allan variance for l = j and clock i's
# own for l ≠ j (Vernotte and Lantz (2019), Appendix B). Each term is divided
# by the geometric mean of the EDF of the records involved.
function _gl_edf(i, q, couples, v, s2, ν)
    v > 0 || return NaN
    e(j, l) = j == l ? s2[i, j][q] : v
    total = 0.0
    for (j, k) in couples, (l, n) in couples
        c = e(j, l) * e(k, n) + e(j, n) * e(k, l)
        total += c / sqrt(sqrt(ν[j][q] * ν[k][q] * ν[l][q] * ν[n][q]))
    end
    return 2 * v^2 * length(couples)^2 / total
end

# χ² interval of a separated variance on its equivalent EDF. Below 1 EDF the χ²
# model degenerates (its lower bound can pass the estimate), so no interval is
# given there.
function _sep_chi2_ci(dev, edf, confidence)
    lo, hi = confidence_intervals(dev, edf, confidence)
    for q in eachindex(edf)
        edf[q] >= 1 || (lo[q] = hi[q] = NaN)
    end
    return _zip_ci(lo, hi)
end

# KLTS intervals (Lantz et al. (2019)) for clocks 1, 2, 3 at each τ, from the
# sample covariance of the records 1−2, 2−3, 3−1, the measurement noise read
# from their closure (σ²_N = σ²_closure / 3, Eq. (9)) and the geometric mean
# of the three records' EDF. A lower bound that depends on the prior is 0.
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
        Sk = vN <= 1e-9 * sum(S[k, k] for k in 1:3) / 3 ? S[1:2, 1:2] : S
        qs, zero_lower = _klts(Sk, vN, n, probs)
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

groslambert(cc::ClockComparisons, args...; reference::Bool=true, kw...) =
    groslambert(_hat_records(cc, reference, "groslambert"), args...; kw...)

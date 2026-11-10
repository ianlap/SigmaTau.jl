# KLTS confidence intervals for the Groslambert covariance of three clocks
# [lantz-2019-klts](@cite): the Bayesian posterior of the three clock variances
# given the sample covariance of the three difference records.

# Prior range of each clock variance, relative to the smallest and largest
# |estimate|: log-uniform ("total ignorance") over [1e-5 min, 1e3 max]. The
# paper takes 1e-5 for the lower limit; the upper one is SigmaTheta's.
const _KLTS_PRIOR = (1e-5, 1e3)

# Log-likelihood of clock variances (a, b, c) for records (ab, bc, ca) with
# sample covariance S, measurement noise vN per record and ν degrees of
# freedom: −(ν/2)(log det Σ + tr(Σ⁻¹ S)), the logarithm of eqs (13)–(14).
# Σ has Var(ab) = a + b + vN, … and Cov(ab, bc) = −b, Cov(bc, ca) = −c,
# Cov(ab, ca) = −a. Without measurement noise Σ is singular (the records sum
# to zero) and the first two records carry everything (the paper's 2 × 2 case).
function _klts_loglik(a, b, c, vN, S, ν)
    if size(S, 1) == 2
        s11, s22, s12 = a + b + vN, b + c + vN, -b
        d = s11 * s22 - s12^2
        d > 0 || return -Inf
        return -ν / 2 * (log(d) + (s22 * S[1, 1] + s11 * S[2, 2] - 2s12 * S[1, 2]) / d)
    end
    s11, s22, s33, s12, s23, s13 = a + b + vN, b + c + vN, c + a + vN, -b, -c, -a
    a11 = s22 * s33 - s23^2; a22 = s11 * s33 - s13^2; a33 = s11 * s22 - s12^2
    a12 = s13 * s23 - s12 * s33; a13 = s12 * s23 - s13 * s22; a23 = s12 * s13 - s11 * s23
    d = s11 * a11 + s12 * a12 + s13 * a13
    d > 0 || return -Inf
    tr = a11 * S[1, 1] + a22 * S[2, 2] + a33 * S[3, 3] + 2 * (a12 * S[1, 2] + a13 * S[1, 3] + a23 * S[2, 3])
    return -ν / 2 * (log(d) + tr / d)
end

"""
    _klts(S, vN, ν, probs) → (q, zero_lower)

Quantiles `probs` of the posterior of clock variances (a, b, c) (rows of the
3 × length(probs) matrix `q`), and for each clock whether its lower bound is
to be read as 0. `S` is the 3 × 3 sample covariance of the records ab, bc, ca
(or the 2 × 2 one of ab, bc when there is no measurement noise).

The paper draws 10⁷ random triplets; here the posterior is integrated on a
log grid instead: a coarse grid over the whole prior box, shrunk onto the
posterior mass five times, then a fine 128³ grid, with trapezoidal CDFs.
The lower bound reads 0 when the lowest decade of the prior holds more than
the 0.135 % (−3σ) tail, because the bound then depends on the prior's lower
limit: the paper's rule, which replaces such a bound by 0.
"""
function _klts(S::AbstractMatrix{Float64}, vN::Float64, ν::Float64, probs)
    est = size(S, 1) == 3 ? abs.((S[1, 3], S[1, 2], S[2, 3])) :
          abs.((S[1, 1] + S[1, 2], S[1, 2], S[2, 2] + S[1, 2]))
    pos = filter(>(0), collect(est))
    lo = log(_KLTS_PRIOR[1] * minimum(pos)); hi = log(_KLTS_PRIOR[2] * maximum(pos))
    box = fill((lo, hi), 3)
    u = [range(b[1], b[2]; length=40) for b in box]; m = Vector{Vector{Float64}}(undef, 3)
    for pass in 1:6
        n = pass == 6 ? 128 : 40
        u .= [range(b[1], b[2]; length=n) for b in box]
        e = [exp.(ui) for ui in u]
        L = [_klts_loglik(e[1][i], e[2][j], e[3][k], vN, S, ν) for i in 1:n, j in 1:n, k in 1:n]
        L .= exp.(L .- maximum(L))
        m .= [vec(sum(L; dims=(2, 3))), vec(sum(L; dims=(1, 3))), vec(sum(L; dims=(1, 2)))]
        pass == 6 && break
        for d in 1:3
            c = _trapz_cdf(m[d]); h = step(u[d])
            i1 = something(findfirst(>(1e-10), c), 1); i2 = something(findfirst(>=(1 - 1e-10), c), n)
            box[d] = (max(lo, u[d][i1] - 2h), min(hi, u[d][i2] + 2h))
        end
    end
    q = Matrix{Float64}(undef, 3, length(probs)); zero_lower = falses(3)
    for d in 1:3
        c = _trapz_cdf(m[d])
        q[d, :] .= [exp(_cdf_inverse(c, u[d], p)) for p in probs]
        zero_lower[d] = _cdf_at(c, u[d], lo + log(10)) > 0.00135
    end
    return q, zero_lower
end

# Normalized CDF at the grid points of a density sampled on a uniform grid (trapezoids).
function _trapz_cdf(m)
    c = [0.0; cumsum((m[1:end-1] .+ m[2:end]) ./ 2)]
    return c ./ c[end]
end

# CDF `c` over grid `u` at the point `x`, linear between points.
function _cdf_at(c, u, x)
    x <= u[1] && return 0.0
    x >= u[end] && return 1.0
    i = searchsortedlast(u, x)
    return c[i] + (x - u[i]) / (u[i+1] - u[i]) * (c[i+1] - c[i])
end

# Point where the CDF `c` over grid `u` crosses `p`, linear between points.
function _cdf_inverse(c, u, p)
    i = findfirst(>=(p), c)
    i === nothing && return u[end]
    i == 1 && return u[1]
    return u[i-1] + (p - c[i-1]) / (c[i] - c[i-1]) * (u[i] - u[i-1])
end

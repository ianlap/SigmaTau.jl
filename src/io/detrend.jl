# Mean, endpoint, linear and quadratic detrending of PhaseData / FrequencyData,
# and the least-squares fit behind the polynomial modes, shared with
# `drift_fit`, the noise identifier and the CLI's drift removal. Dispatching
# on the record types lets the bare name `detrend` coexist with
# `detrend(::Vector)` from other packages.

const _DETREND_MODES = (:none, :mean, :endpoint, :linear, :quadratic)
# Precomputed so error paths do not print a tuple at runtime (keeps juliac --trim happy).
const _DETREND_MODES_STR = string(_DETREND_MODES)

# ── Least squares ────────────────────────────────────────────────────────

# Column norms of X over the rows where `ok`, to scale the columns to unit
# length so a t² column in seconds does not wreck the conditioning.
_col_scale(X::Matrix{Float64}, ok::BitVector) =
    [max(sqrt(sum(abs2(X[i, j]) for i in eachindex(ok) if ok[i])), floatmin()) for j in 1:size(X, 2)]

# Inverse of a symmetric 2×2 or 3×3 matrix, by cofactors.
function _sym_inv(A::Matrix{Float64})
    if size(A, 1) == 2
        a, b, d = A[1, 1], A[1, 2], A[2, 2]
        return [d -b; -b a] ./ (a * d - b^2)
    end
    a, b, c, d, e, f = A[1, 1], A[1, 2], A[1, 3], A[2, 2], A[2, 3], A[3, 3]
    C11 = d * f - e^2; C12 = c * e - b * f; C13 = b * e - c * d
    C22 = a * f - c^2; C23 = b * c - a * e; C33 = a * d - b^2
    return [C11 C12 C13; C12 C22 C23; C13 C23 C33] ./ (a * C11 + b * C12 + c * C13)
end

# Inverse of the scaled normal matrix XᵀX (two or three columns) over the
# rows where `ok`.
function _normal_inv(X::Matrix{Float64}, ok::BitVector, sc::Vector{Float64})
    p = size(X, 2)
    return _sym_inv([sum(X[i, j] * X[i, k] for i in eachindex(ok) if ok[i]) / (sc[j] * sc[k])
                     for j in 1:p, k in 1:p])
end

# Least squares on the columns of X over the rows where `ok`; returns the
# coefficients and the residual sum of squares. A second pass on the
# residual recovers the digits the normal equations lose.
function _lsq(X::Matrix{Float64}, y::Vector{Float64}, ok::BitVector)
    p = size(X, 2)
    sc = _col_scale(X, ok)
    Ainv = _normal_inv(X, ok, sc)
    beta = zeros(p)
    r = [ok[i] ? y[i] : 0.0 for i in eachindex(y)]
    for _ in 1:2
        b = [sum(X[i, j] * r[i] for i in eachindex(ok) if ok[i]) / sc[j] for j in 1:p]
        step = [sum(Ainv[j, k] * b[k] for k in 1:p) for j in 1:p]
        beta .+= step
        for i in eachindex(r)
            ok[i] && (r[i] -= sum(X[i, j] * step[j] / sc[j] for j in 1:p))
        end
    end
    return beta ./ sc, sum(abs2(r[i]) for i in eachindex(r) if ok[i])
end

# Least-squares polynomial of `degree` 1 or 2 in the centred sample index
# t = (i − 1) − (N − 1)/2 over the rows where `ok`, as [c0, c1, c2] (c2 = 0
# for a line). Centring makes c1 the slope at the record's middle.
function _index_polyfit(x::Vector{Float64}, ok::BitVector, degree::Int)
    N = length(x)
    t = [(i - 1) - (N - 1) / 2 for i in 1:N]
    X = degree == 1 ? hcat(ones(N), t) : hcat(ones(N), t, t .^ 2)
    c, _ = _lsq(X, x, ok)
    return degree == 1 ? [c[1], c[2], 0.0] : c
end

# ── Detrending ───────────────────────────────────────────────────────────

# In-place detrend of a real vector over its finite samples; NaN samples stay
# NaN and do not enter the fit. `:none`, `:mean`, `:endpoint` (the line
# through the first and last finite samples), `:linear` and `:quadratic`
# (least squares in the sample index; on phase data the quadratic term is a
# frequency drift, the removal SP1065 prescribes before the total estimators).
function _detrend_core!(x::Vector{Float64}, mode::Symbol)
    mode === :none && return x
    mode in _DETREND_MODES ||
        throw(ArgumentError("detrend: unknown method $mode (expected one of $_DETREND_MODES_STR)"))
    length(x) <= 1 && return x
    ok = isfinite.(x)
    nok = count(ok)
    if mode === :mean
        nok == 0 && return x
        x .-= sum(x[i] for i in eachindex(x) if ok[i]) / nok
    elseif mode === :endpoint
        nok < 2 && return x
        i1, i2 = findfirst(ok), findlast(ok)
        x1, slope = x[i1], (x[i2] - x[i1]) / (i2 - i1)
        @inbounds for i in eachindex(x)
            x[i] -= x1 + slope * (i - i1)
        end
    else
        degree = mode === :linear ? 1 : 2
        nok > degree || return x
        c = _index_polyfit(x, ok, degree)
        mu = (length(x) - 1) / 2
        @inbounds for i in eachindex(x)
            t = (i - 1) - mu
            x[i] -= c[1] + (c[2] + c[3] * t) * t
        end
    end
    return x
end

# Out-of-place `:quadratic` detrend: the drift removal before `mhtotdev` and
# inside the noise identifier.
_remove_quadratic(x::Vector{Float64}) = _detrend_core!(copy(x), :quadratic)

"""
    detrend(pd::PhaseData; method::Symbol=:linear) → PhaseData
    detrend(fd::FrequencyData; method::Symbol=:linear) → FrequencyData

Return a new `PhaseData` / `FrequencyData` with a detrended sample vector.
The original is left untouched.

`method` is one of `:linear` (least-squares straight line, default),
`:quadratic` (least-squares parabola; on phase data this removes a frequency
drift, the preprocessing SP1065 prescribes before the total estimators),
`:endpoint` (line through first and last samples), `:mean` (zero-mean), or
`:none` (returns a copy unchanged). On phase data `:linear` removes a
frequency offset, not a drift. Gaps (`NaN`) stay where they are and do not
enter the fit.

# Examples

```julia
pd = PhaseData(cumsum(randn(1024)) .+ 1e-3 .* (0:1023), 1.0)
pd_clean = detrend(pd)                    # :linear default
pd_zero  = detrend(pd; method=:mean)
```
"""
detrend(pd::PhaseData; method::Symbol=:linear) =
    PhaseData(_detrend_core!(Vector{Float64}(pd.x), method), pd.tau0; source=pd.source)

detrend(fd::FrequencyData; method::Symbol=:linear) =
    FrequencyData(_detrend_core!(Vector{Float64}(fd.y), method), fd.tau0; source=fd.source)

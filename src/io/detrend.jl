# io/detrend.jl — Linear / endpoint / mean detrending for PhaseData and FrequencyData.
#
# Multiple dispatch on our timing-data types lets us export the bare name
# `detrend` without colliding with `detrend(::Vector)` from other packages —
# they're different methods on different types.

const _DETREND_MODES = (:none, :mean, :endpoint, :linear, :quadratic)

"""
    _detrend_core!(x::Vector{Float64}, mode::Symbol) → x

In-place detrend of a real vector. `mode` is one of `:none`, `:mean`,
`:endpoint`, `:linear`, `:quadratic`. Internal helper for
`detrend(::PhaseData)` and `detrend(::FrequencyData)`.

- `:none`      — no-op.
- `:mean`      — subtract the sample mean.
- `:endpoint`  — subtract the line connecting the first and last samples.
- `:linear`    — subtract the least-squares-fit straight line in index `n=0..N-1`.
- `:quadratic` — subtract the least-squares-fit quadratic in index `n=0..N-1`.
  On phase data the quadratic term is frequency drift, so this is the drift
  removal SP1065 prescribes before the total estimators.
"""
function _detrend_core!(x::Vector{Float64}, mode::Symbol)
    mode === :none && return x
    N = length(x)
    N <= 1 && return x

    if mode === :mean
        m = sum(x) / N
        @inbounds @simd for i in eachindex(x)
            x[i] -= m
        end
        return x
    elseif mode === :endpoint
        x1    = x[1]
        slope = (x[end] - x1) / (N - 1)
        @inbounds @simd for i in eachindex(x)
            x[i] -= x1 + slope * (i - 1)
        end
        return x
    elseif mode === :linear
        # closed-form OLS for y = a + b·n on n = 0..N-1
        sx  = (N - 1) * N / 2                 # Σ n
        sxx = (N - 1) * N * (2N - 1) / 6      # Σ n²
        sy  = 0.0
        sxy = 0.0
        @inbounds @simd for i in 1:N
            xi   = x[i]
            sy  += xi
            sxy += (i - 1) * xi
        end
        denom = N * sxx - sx * sx
        b     = (N * sxy - sx * sy) / denom
        a     = (sy - b * sx) / N
        @inbounds @simd for i in eachindex(x)
            x[i] -= a + b * (i - 1)
        end
        return x
    elseif mode === :quadratic
        # Closed-form OLS for y = c0 + c1·t + c2·t² on the centered index
        # t = (i−1) − (N−1)/2. Centering zeroes the odd power sums
        # (Σt = Σt³ = 0), so c1 decouples and [c0, c2] solve a 2×2 system.
        # The fit is a linear projection, so a deterministic quadratic added
        # to the data is removed exactly.
        mu = (N - 1) / 2.0
        S2 = 0.0; S4 = 0.0
        Sy = 0.0; Sty = 0.0; St2y = 0.0
        @inbounds for i in 1:N
            t  = (i - 1) - mu
            t2 = t * t
            xi = x[i]
            S2   += t2
            S4   += t2 * t2
            Sy   += xi
            Sty  += t * xi
            St2y += t2 * xi
        end
        c1  = Sty / S2
        det = N * S4 - S2 * S2
        c2  = (N * St2y - S2 * Sy) / det
        c0  = (Sy - c2 * S2) / N
        @inbounds @simd for i in eachindex(x)
            t = (i - 1) - mu
            x[i] -= c0 + c1 * t + c2 * t * t
        end
        return x
    else
        throw(ArgumentError("detrend: unknown method $mode (expected one of $_DETREND_MODES)"))
    end
end

"""
    _remove_quadratic(x::Vector{Float64}) → Vector{Float64}

Out-of-place least-squares quadratic removal (the `:quadratic` mode of
[`_detrend_core!`](@ref) on a copy). Used by `mhtotdev`'s default
`remove_drift=true` path.
"""
_remove_quadratic(x::Vector{Float64}) = _detrend_core!(copy(x), :quadratic)

"""
    detrend(pd::PhaseData; method::Symbol=:linear) → PhaseData
    detrend(fd::FrequencyData; method::Symbol=:linear) → FrequencyData

Return a new `PhaseData` / `FrequencyData` with a detrended sample vector.
The original is left untouched.

`method` is one of `:linear` (least-squares straight line, default),
`:endpoint` (line through first and last samples), `:mean` (zero-mean), or
`:none` (returns a copy unchanged).

# Examples

```julia
pd = PhaseData(cumsum(randn(1024)) .+ 1e-3 .* (0:1023), 1.0)
pd_clean = detrend(pd)                    # :linear default
pd_zero  = detrend(pd; method=:mean)
```
"""
function detrend(pd::PhaseData; method::Symbol=:linear)
    method in _DETREND_MODES ||
        throw(ArgumentError("detrend: unknown method $method (expected one of $_DETREND_MODES)"))
    x = Vector{Float64}(pd.x)
    _detrend_core!(x, method)
    return PhaseData(x, pd.tau0)
end

function detrend(fd::FrequencyData; method::Symbol=:linear)
    method in _DETREND_MODES ||
        throw(ArgumentError("detrend: unknown method $method (expected one of $_DETREND_MODES)"))
    y = Vector{Float64}(fd.y)
    _detrend_core!(y, method)
    return FrequencyData(y, fd.tau0)
end

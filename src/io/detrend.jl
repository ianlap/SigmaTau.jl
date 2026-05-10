# io/detrend.jl — Linear / endpoint / mean detrending for PhaseData and FrequencyData.
#
# Multiple dispatch on our timing-data types lets us export the bare name
# `detrend` without colliding with `detrend(::Vector)` from other packages —
# they're different methods on different types.

const _DETREND_MODES = (:none, :mean, :endpoint, :linear)

"""
    _detrend_core!(x::Vector{Float64}, mode::Symbol) → x

In-place detrend of a real vector. `mode` is one of `:none`, `:mean`,
`:endpoint`, `:linear`. Internal helper for `detrend(::PhaseData)` and
`detrend(::FrequencyData)`.

- `:none`     — no-op.
- `:mean`     — subtract the sample mean.
- `:endpoint` — subtract the line connecting the first and last samples.
- `:linear`   — subtract the least-squares-fit straight line in index `n=0..N-1`.
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
    else
        throw(ArgumentError("detrend: unknown method $mode (expected one of $_DETREND_MODES)"))
    end
end

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

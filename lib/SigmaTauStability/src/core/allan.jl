# core/allan.jl — Core Stability Kernels

"""
    _adev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Overlapping Allan Deviation (ADEV) for a set of averaging factors `m`.
This is a highly optimized, SIMD-vectorized loop with no allocations in the inner loop.
"""
function _adev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))
    
    for (k, m) in enumerate(m_values)
        L = N - 2m
        if L <= 0
            devs[k] = NaN
            continue
        end
        
        sum_sq = 0.0
        @inbounds @simd for i in 1:L
            d2 = x[i+2m] - 2x[i+m] + x[i]
            sum_sq += d2^2
        end
        
        devs[k] = sqrt(sum_sq / (2.0 * L * Float64(m)^2 * tau0^2))
    end
    
    return devs
end

"""
    _mdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Modified Allan Deviation (MDEV) using prefix sums for O(N) performance per `m`.
"""
function _mdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # Precompute cumulative sum
    x_cs = cumsum(pushfirst!(copy(x), 0.0))

    for (k, m) in enumerate(m_values)
        Ne = N - 3m + 1
        if Ne <= 0
            devs[k] = NaN
            continue
        end

        # Pairwise sum (Julia's default for `sum`) is bit-stable across
        # platforms; the previous `@simd` accumulator allowed CPU-dependent
        # reordering that drifted ~1 ULP from the legacy reference on Linux.
        sum_sq = @inbounds sum(1:Ne) do i
            d = x_cs[i+3m] - 3x_cs[i+2m] + 3x_cs[i+m] - x_cs[i]
            d * d
        end

        devs[k] = sqrt(sum_sq / (2.0 * Ne * Float64(m)^4 * tau0^2))
    end

    return devs
end

"""
    _tdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Time Deviation (TDEV).
"""
function _tdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    mdevs = _mdev_core(x, m_values, tau0)
    taus = m_values .* tau0
    return taus .* mdevs ./ sqrt(3.0)
end

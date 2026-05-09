# core/hadamard.jl — Core Hadamard Stability Kernels

"""
    _hdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Overlapping Hadamard Deviation (HDEV) for a set of averaging factors `m`.
This uses a highly optimized, SIMD-vectorized loop with no allocations in the inner loop.
"""
function _hdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))
    
    for (k, m) in enumerate(m_values)
        L = N - 3m
        if L <= 0
            devs[k] = NaN
            continue
        end
        
        sum_sq = 0.0
        @inbounds @simd for i in 1:L
            d3 = x[i+3m] - 3.0 * x[i+2m] + 3.0 * x[i+m] - x[i]
            sum_sq += d3^2
        end
        
        devs[k] = sqrt(sum_sq / (6.0 * L * Float64(m)^2 * tau0^2))
    end
    
    return devs
end

"""
    _mhdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Modified Hadamard Deviation (MHDEV) using prefix sums for O(N) performance per `m`.
"""
function _mhdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # Precompute cumulative sum
    x_cs = cumsum(pushfirst!(copy(x), 0.0))

    for (k, m) in enumerate(m_values)
        Ne = N - 4m + 1
        if Ne <= 0
            devs[k] = NaN
            continue
        end

        sum_sq = 0.0
        @inbounds @simd for i in 1:Ne
            d = x_cs[i+4m] - 4.0 * x_cs[i+3m] + 6.0 * x_cs[i+2m] - 4.0 * x_cs[i+m] + x_cs[i]
            sum_sq += d^2
        end

        devs[k] = sqrt(sum_sq / (6.0 * Ne * Float64(m)^4 * tau0^2))
    end

    return devs
end

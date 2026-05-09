# core/pdev.jl — Parabolic deviation (Vernotte 2015 / 2020)

"""
    _pdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Parabolic deviation σ_PDEV(τ) per Vernotte–Lenczner–Bourgeois–Rubiola
(IEEE T-UFFC 63(4), 2016) and Vernotte 2020. Built from a least-squares
parabolic fit to the phase record over each window:

```math
\\sigma^2_{\\text{PDEV}}(m\\tau_0) = \\frac{72}{(N-2m) \\, m^4 \\, (m\\tau_0)^2}
   \\sum_{i=1}^{N-2m} \\left[
       \\sum_{k=0}^{m-1} \\left(\\frac{m-1}{2} - k\\right)
       \\bigl(x_{i+k} - x_{i+k+m}\\bigr)
   \\right]^2
```

For `m = 1` the parabolic weight collapses to zero, so we fall back to
overlapping ADEV (the canonical PDEV(τ₀) ≡ ADEV(τ₀) identity from
Vernotte 2015).
"""
function _pdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    for (k, m) in enumerate(m_values)
        if m == 1
            devs[k] = _adev_core(x, [m], tau0)[1]
            continue
        end

        M = N - 2m
        if M < 1
            devs[k] = NaN
            continue
        end

        half = (m - 1) / 2.0
        Msum = 0.0
        @inbounds for i in 1:M
            asum = 0.0
            @simd for kk in 0:(m - 1)
                w = half - kk
                asum += w * (x[i + kk] - x[i + kk + m])
            end
            Msum += asum * asum
        end

        var = 72.0 * Msum / (M * Float64(m)^6 * tau0^2)
        devs[k] = sqrt(var)
    end

    return devs
end

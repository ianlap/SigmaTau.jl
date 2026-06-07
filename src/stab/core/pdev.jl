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
Vernotte 2015). For `m > 1` the weighted parabolic sum is evaluated by a
rolling two-sum recurrence, reducing each averaging factor from O(Nm) to O(N).
"""
function _pdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    for (k, m) in enumerate(m_values)
        if m < 1
            devs[k] = NaN
            continue
        end

        if m == 1
            devs[k] = _adev_core(x, [m], tau0)[1]
            continue
        end

        M = N - 2m
        if M < 2          # need ≥2 windows
            devs[k] = NaN
            continue
        end

        half = (m - 1) / 2.0
        A = 0.0
        B = 0.0
        @inbounds for kk in 0:(m - 1)
            y = x[1 + kk] - x[1 + kk + m]
            A += y
            B += kk * y
        end

        # Periodically rebuild the rolling sums exactly. The interval scales with
        # m so the refresh work remains amortized O(N) while bounding accumulated
        # roundoff on very long records.
        refresh_every = max(4096, m)
        next_refresh = refresh_every

        Msum = 0.0
        @inbounds for i in 1:M
            asum = half * A - B
            Msum += asum * asum

            if i < M
                if i == next_refresh
                    A = 0.0
                    B = 0.0
                    i_next = i + 1
                    @simd for kk in 0:(m - 1)
                        y = x[i_next + kk] - x[i_next + kk + m]
                        A += y
                        B += kk * y
                    end
                    next_refresh += refresh_every
                else
                    y_old = x[i] - x[i + m]
                    y_new = x[i + m] - x[i + 2m]
                    old_A = A
                    A += y_new - y_old
                    B += (m - 1) * y_new - old_A + y_old
                end
            end
        end

        var = 72.0 * Msum / (M * Float64(m)^6 * tau0^2)
        devs[k] = sqrt(var)
    end

    return devs
end

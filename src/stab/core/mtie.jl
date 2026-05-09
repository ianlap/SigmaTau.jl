# core/mtie.jl — Maximum Time Interval Error kernel

"""
    _mtie_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Maximum Time Interval Error per ITU-T G.810. For each averaging factor
`m`, slides a window of `m+1` phase samples (spanning `τ = m·τ₀`) across
the record and returns the largest peak-to-peak phase excursion observed
in any window.

Units are seconds (a σ_x quantity, like TDEV) — MTIE is itself a phase
measure, not a relative-frequency measure, so no τ rescaling is applied.

The implementation is the straightforward O(N·m) sweep — sufficient for
typical clock records. A monotonic-deque or sparse-table optimisation
to O(N log m) is tracked in `TODO.md`.
"""
function _mtie_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    for (k, m) in enumerate(m_values)
        L = N - m
        if L <= 0
            devs[k] = NaN
            continue
        end

        max_excursion = 0.0
        @inbounds for i in 1:L
            lo = x[i]
            hi = x[i]
            @simd for j in (i + 1):(i + m)
                xj = x[j]
                lo = ifelse(xj < lo, xj, lo)
                hi = ifelse(xj > hi, xj, hi)
            end
            d = hi - lo
            if d > max_excursion
                max_excursion = d
            end
        end

        devs[k] = max_excursion
    end

    return devs
end

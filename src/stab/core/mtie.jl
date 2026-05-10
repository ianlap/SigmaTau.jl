# core/mtie.jl — Maximum Time Interval Error kernel

"""
    _mtie_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Maximum Time Interval Error per ITU-T G.810. For each averaging factor
`m`, slides a window of `m+1` phase samples (spanning `τ = m·τ₀`) across
the record and returns the largest peak-to-peak phase excursion observed
in any window.

Units are seconds (a σ_x quantity, like TDEV) — MTIE is itself a phase
measure, not a relative-frequency measure, so no τ rescaling is applied.

# Implementation

Monotonic-deque sliding window: each m value runs in O(N) total work via
two parallel index deques (one tracking the running window maximum, one
the running minimum). Each phase sample enters and leaves each deque at
most once, so overall complexity is O(N · |m_values|).

The deque is materialised as a single pre-allocated `Vector{Int}` of
length N with explicit head/tail cursors per m, avoiding any
DataStructures.jl dependency or per-step allocation.
"""
function _mtie_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # Pre-allocate deque storage once; the deque holds at most `win`
    # indices but cursors can advance up to N before wrapping considerations.
    # Sizing to N avoids any modular indexing.
    max_dq = Vector{Int}(undef, N)
    min_dq = Vector{Int}(undef, N)

    for (k, m) in enumerate(m_values)
        L = N - m
        if L <= 0
            devs[k] = NaN
            continue
        end

        win = m + 1
        max_h, max_t = 1, 0   # empty when max_t < max_h
        min_h, min_t = 1, 0
        max_excursion = 0.0

        @inbounds for j in 1:N
            # Drop indices that have left the trailing edge of the window.
            while max_t >= max_h && max_dq[max_h] <= j - win
                max_h += 1
            end
            while min_t >= min_h && min_dq[min_h] <= j - win
                min_h += 1
            end

            xj = x[j]

            # Maintain monotonic-decreasing deque (running maximum).
            while max_t >= max_h && x[max_dq[max_t]] <= xj
                max_t -= 1
            end
            max_t += 1
            max_dq[max_t] = j

            # Maintain monotonic-increasing deque (running minimum).
            while min_t >= min_h && x[min_dq[min_t]] >= xj
                min_t -= 1
            end
            min_t += 1
            min_dq[min_t] = j

            # Once the window is fully populated, record the excursion.
            if j >= win
                d = x[max_dq[max_h]] - x[min_dq[min_h]]
                if d > max_excursion
                    max_excursion = d
                end
            end
        end

        devs[k] = max_excursion
    end

    return devs
end

# core/total.jl — Core Total Stability Kernels

"""
    _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:howe) → Vector{Float64}

Computes the Total Deviation (TOTDEV) for a set of averaging factors `m`.

TOTDEV uses the Howe 1995 / NIST SP1065 eqn 25 endpoint mean-flip extension.
`detrend` selects the detrending applied before the extension:

- `:howe` — no detrend (canonical SP1065 eqn 25, matches allantools).
  Default.
- `:linear` — global least-squares detrend over the whole vector; alias
  for `:legacy` on this kernel.
- `:legacy` — pre-1.0 SigmaTau behavior; identical to `:linear` here.
"""
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:howe)
    detrend === :howe && return _totdev_howe(x, m_values, tau0)
    if detrend === :legacy || detrend === :linear
        return _totdev_legacy(x, m_values, tau0)
    end
    throw(ArgumentError("unknown detrend recipe: $detrend; valid for TOTDEV: :howe, :linear, :legacy"))
end

function _totdev_legacy(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # Linear detrend of the whole vector (analytic LS sums)
    N_float = Float64(N)
    sum_i = (N_float * (N_float + 1.0)) / 2.0
    sum_i2 = (N_float * (N_float + 1.0) * (2.0*N_float + 1.0)) / 6.0
    delta = N_float * sum_i2 - sum_i^2

    sum_x = sum(x)
    sum_ix = 0.0
    @inbounds @simd for i in 1:N
        sum_ix += i * x[i]
    end

    a = (sum_x * sum_i2 - sum_ix * sum_i) / delta
    b = (N_float * sum_ix - sum_x * sum_i) / delta

    xd = Vector{Float64}(undef, N)
    @inbounds @simd for i in 1:N
        xd[i] = x[i] - (a + b * i)
    end

    # Mean-flip endpoint reflection: x_star of length 3N-4
    x_star = Vector{Float64}(undef, 3N - 4)
    @inbounds for i in 1:N-2
        x_star[i] = 2.0*xd[1] - xd[i+1]
    end
    @inbounds for i in 1:N
        x_star[N-2+i] = xd[i]
    end
    @inbounds for i in 1:N-2
        x_star[2N-2+i] = 2.0*xd[N] - xd[N-i]
    end

    off = N - 2
    Threads.@threads :dynamic for k in eachindex(m_values)
        m = m_values[k]
        D = 0.0
        count = 0
        @inbounds @simd for i in 1:N
            lo = off + i
            hi = off + i + 2m
            if hi <= length(x_star)
                d2 = x_star[hi] - 2.0*x_star[off + i + m] + x_star[lo]
                D += d2^2
                count += 1
            end
        end

        if count == 0
            devs[k] = NaN
        else
            devs[k] = sqrt(D / (2.0 * (N - 2) * Float64(m)^2 * tau0^2))
        end
    end

    return devs
end

function _totdev_howe(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # NIST SP1065 eqn 25 / Greenhall-Howe-Percival 1998 eq (3):
    #   Totvar = (1 / (2(m·τ₀)²(N−2))) · Σ_{n=2}^{N-1} (x*_{n-m} − 2x*_n + x*_{n+m})²
    # Mean-flip endpoint reflection (SP1065 eq 2):
    #   x*_k =  x[k]                   for 1 ≤ k ≤ N
    #   x*_k =  2·x[1] − x[2 − k]      for k ≤ 0       (left reflect)
    #   x*_k =  2·x[N] − x[2N − k]     for k ≥ N + 1   (right reflect)
    #
    # Stream the extension on the fly rather than materialising a (3N−4)-Float64
    # buffer. The inner n-loop splits into three ranges so the central majority
    # (which only reads x[n±m] without reflection) stays SIMD-vectorisable;
    # the two short tail loops resolve the reflection inline. Saves a 24·N-byte
    # allocation per call (~72 MiB at N = 3 × 10⁶).
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    if N <= 2
        fill!(devs, NaN)
        return devs
    end

    x1 = x[1]
    xN = x[N]

    Threads.@threads :dynamic for k in eachindex(m_values)
        m = m_values[k]
        D = 0.0

        # Left tail: lo_k = n − m < 1. hi_k may or may not reflect.
        n_left_end = min(N - 1, m)
        @inbounds for n in 2:n_left_end
            lo_k = n - m
            hi_k = n + m
            lo_val = lo_k >= 1 ? x[lo_k] : 2.0*x1 - x[2 - lo_k]
            hi_val = hi_k <= N ? x[hi_k] : 2.0*xN - x[2N - hi_k]
            d2 = hi_val - 2.0*x[n] + lo_val
            D += d2^2
        end

        # Central: lo_k ≥ 1 AND hi_k ≤ N. Plain three-point stencil, vectorises.
        n_central_start = max(2, m + 1)
        n_central_end   = min(N - 1, N - m)
        if n_central_start <= n_central_end
            @inbounds @simd for n in n_central_start:n_central_end
                d2 = x[n+m] - 2.0*x[n] + x[n-m]
                D += d2^2
            end
        end

        # Right tail: hi_k = n + m > N. lo_k is central for m ≤ (N−1)/2; the
        # branch below covers the rare m > (N−1)/2 case where lo_k still reflects.
        # `max(n_left_end + 1, …)` prevents double-counting n values already
        # visited by the left loop when the two regions would otherwise overlap.
        n_right_start = max(n_left_end + 1, N - m + 1)
        @inbounds for n in n_right_start:N-1
            lo_k = n - m
            hi_k = n + m
            lo_val = lo_k >= 1 ? x[lo_k] : 2.0*x1 - x[2 - lo_k]
            hi_val = hi_k <= N ? x[hi_k] : 2.0*xN - x[2N - hi_k]
            d2 = hi_val - 2.0*x[n] + lo_val
            D += d2^2
        end

        devs[k] = sqrt(D / (2.0 * (N - 2) * Float64(m)^2 * tau0^2))
    end

    return devs
end

"""
    _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall) → Vector{Float64}

Computes the Modified Total Deviation (MTOTDEV).

MTOTDEV uses Greenhall 2003's per-window time-reverse extension. `detrend`
selects the per-window detrending applied before the extension:

- `:greenhall` — half-mean slope removal (Greenhall 2003 canonical).
- `:linear` — full least-squares (slope + intercept) per window; tighter
  detrend at the cost of slightly higher per-window variance.
- `:legacy` — pre-1.0 SigmaTau behavior; alias for `:greenhall` here.
"""
function _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _mtotdev_greenhall(x, m_values, tau0)
    end
    detrend === :linear && return _mtotdev_linear(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid for MTOTDEV: :greenhall, :linear, :legacy"))
end

function _mtotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # Sequential m-loop, parallel over subsequences. The subsequence list
    # (length nsubs ≈ N − 3m) is long and uniform; chunking it across
    # threads gives near-perfect speedup per m. Per-m parallelism was
    # tried and lost — uneven m-work + nested @threads overhead outweighed
    # the gain. Reduction order changes vs sequential: a few-ULP drift is
    # documented in CHANGELOG.
    #
    # ext-buffer pool: one buffer per chunk slot, sized once for the
    # largest m. Reused across all m's; cheap m's only touch the first
    # 3·seg_len entries. Eliminates per-m, per-chunk allocation churn
    # (was ~600 MB of churn on the 3M-sample file, max_m=2^19).
    nthreads = Threads.nthreads()
    max_m = maximum(m_values; init=0)
    if max_m < 1 || N - 3*max_m + 1 < 1
        # All m's invalid or no subsequences possible: fall through and
        # let the per-k NaN branch handle each.
    end
    ext_pool = [Vector{Float64}(undef, 3 * 3 * max(max_m, 1)) for _ in 1:nthreads]

    for k in eachindex(m_values)
        m = m_values[k]
        nsubs = N - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        nchunks = max(1, min(nthreads, nsubs))
        chunk_size = cld(nsubs, nchunks)
        chunk_sums = zeros(nchunks)

        half_n = seg_len / 2.0
        hi_half = floor(Int, half_n)
        lo_half_len = seg_len - hi_half
        inv_hi = 1.0 / hi_half
        inv_lo = 1.0 / lo_half_len
        inv_half_n_tau0 = 1.0 / (half_n * tau0)

        Threads.@threads :dynamic for c in 1:nchunks
            n_lo = (c - 1) * chunk_size + 1
            n_hi = min(c * chunk_size, nsubs)
            ext = ext_pool[c]
            local_sum = 0.0

            # Seed the running half-window sums for the chunk's first window.
            # For m ≥ 2 we then update s1/s2 in O(1) per n inside the loop
            # rather than re-summing hi/lo_half_len samples each window.
            # Seed the running half-sums for the chunk's first window (m ≥ 2).
            s1 = 0.0
            s2 = 0.0
            if m >= 2
                @inbounds @simd for i in 1:hi_half
                    s1 += x[n_lo - 1 + i]
                end
                @inbounds @simd for i in hi_half+1:seg_len
                    s2 += x[n_lo - 1 + i]
                end
            end

            for n in n_lo:n_hi
                if m == 1
                    slope = (x[n+2] - x[n]) / (2.0 * tau0)
                else
                    slope = (s2 * inv_lo - s1 * inv_hi) * inv_half_n_tau0
                end

                @inbounds for j in 1:seg_len
                    val = x[n-1+j] - slope * tau0 * (j - 1)
                    rev_val = x[n-1 + seg_len - j + 1] - slope * tau0 * (seg_len - j)

                    ext[j] = rev_val
                    ext[seg_len + j] = val
                    ext[2seg_len + j] = rev_val
                end

                # Sliding-window inner reduction. a1, a2, a3 are running
                # sums of m consecutive ext values, equivalent to
                # cs[j+km+1] - cs[j+(k-1)m+1] from the cumulative-sum
                # variant. Eliminating cs saves the build pass (9m memory
                # ops per subseq) and one per-chunk allocation. The carry
                # chain on a1/a2/a3 prevents @simd vectorisation of the
                # slide loop, but the savings on cs traffic outweigh the
                # lost SIMD.
                a1 = 0.0
                a2 = 0.0
                a3 = 0.0
                @inbounds @simd for i in 1:m
                    a1 += ext[i]
                    a2 += ext[i+m]
                    a3 += ext[i+2m]
                end

                d2 = (a3 - 2.0*a2 + a1) / m
                block_sum = d2^2

                @inbounds for j in 1:(6m - 1)
                    a1 += ext[j+m]  - ext[j]
                    a2 += ext[j+2m] - ext[j+m]
                    a3 += ext[j+3m] - ext[j+2m]
                    d2 = (a3 - 2.0*a2 + a1) / m
                    block_sum += d2^2
                end
                local_sum += block_sum / (6.0 * m)

                # O(1) half-sum update for the next window: drop the head of
                # each half and pick up the new tail. s1 covers x[n..n+hi_half-1],
                # s2 covers x[n+hi_half..n+seg_len-1]; shifting by 1 swaps x[n]
                # out of s1, x[n+hi_half] from s2 into s1, and x[n+seg_len] into s2.
                if m >= 2 && n < n_hi
                    @inbounds begin
                        s1 += x[n + hi_half] - x[n]
                        s2 += x[n + seg_len] - x[n + hi_half]
                    end
                end
            end
            chunk_sums[c] = local_sum
        end

        outer_sum = sum(chunk_sums)
        devs[k] = sqrt(outer_sum / (2.0 * Float64(m)^2 * tau0^2 * nsubs))
    end

    return devs
end

function _mtotdev_linear(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window full LS detrend (analytic slope+intercept) + per-window
    # time-reverse extension + modified 2nd-difference operator. Same
    # extension/operator shape as `_mtotdev_greenhall`; only the slope
    # estimate differs (full LS instead of half-mean).
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # See `_mtotdev_greenhall` for the sequential-outer / parallel-inner
    # rationale, sliding-window inner reduction, and ext-buffer pool.
    nthreads = Threads.nthreads()
    max_m = maximum(m_values; init=0)
    ext_pool = [Vector{Float64}(undef, 3 * 3 * max(max_m, 1)) for _ in 1:nthreads]

    for k in eachindex(m_values)
        m = m_values[k]
        nsubs = N - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        L_float = Float64(seg_len)
        sum_i = (L_float * (L_float + 1.0)) / 2.0
        sum_i2 = (L_float * (L_float + 1.0) * (2.0*L_float + 1.0)) / 6.0
        delta = L_float * sum_i2 - sum_i^2

        nchunks = max(1, min(nthreads, nsubs))
        chunk_size = cld(nsubs, nchunks)
        chunk_sums = zeros(nchunks)

        Threads.@threads :dynamic for c in 1:nchunks
            n_lo = (c - 1) * chunk_size + 1
            n_hi = min(c * chunk_size, nsubs)
            ext = ext_pool[c]
            local_sum = 0.0
            for n in n_lo:n_hi
                sum_x = 0.0
                sum_ix = 0.0
                @inbounds @simd for j in 1:seg_len
                    v = x[n-1+j]
                    sum_x += v
                    sum_ix += j * v
                end

                a = (sum_x * sum_i2 - sum_ix * sum_i) / delta
                b = (L_float * sum_ix - sum_x * sum_i) / delta

                @inbounds for j in 1:seg_len
                    val = x[n-1+j] - (a + b * j)
                    rev_val = x[n-1 + seg_len - j + 1] - (a + b * (seg_len - j + 1))

                    ext[j] = rev_val
                    ext[seg_len + j] = val
                    ext[2seg_len + j] = rev_val
                end

                a1 = 0.0
                a2 = 0.0
                a3 = 0.0
                @inbounds @simd for i in 1:m
                    a1 += ext[i]
                    a2 += ext[i+m]
                    a3 += ext[i+2m]
                end

                d2 = (a3 - 2.0*a2 + a1) / m
                block_sum = d2^2

                @inbounds for j in 1:(6m - 1)
                    a1 += ext[j+m]  - ext[j]
                    a2 += ext[j+2m] - ext[j+m]
                    a3 += ext[j+3m] - ext[j+2m]
                    d2 = (a3 - 2.0*a2 + a1) / m
                    block_sum += d2^2
                end
                local_sum += block_sum / (6.0 * m)
            end
            chunk_sums[c] = local_sum
        end

        outer_sum = sum(chunk_sums)
        devs[k] = sqrt(outer_sum / (2.0 * Float64(m)^2 * tau0^2 * nsubs))
    end

    return devs
end

"""
    _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall) → Vector{Float64}

Computes the Hadamard Total Deviation (HTOTDEV) on the frequency series
`y = diff(x) / tau0`.

HTOTDEV uses Greenhall 2003's per-window time-reverse extension. `detrend`
selects the per-window detrending applied to `y` before the extension:

- `:greenhall` — half-mean slope removal on `y` (Greenhall 2003 canonical).
- `:linear` — full least-squares (slope + intercept) per window on `y`.
- `:legacy` — pre-1.0 SigmaTau behavior; alias for `:greenhall` here.
"""
function _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _htotdev_greenhall(x, m_values, tau0)
    end
    detrend === :linear && return _htotdev_linear(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid for HTOTDEV: :greenhall, :linear, :legacy"))
end

function _htotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    y = Vector{Float64}(undef, N-1)
    @inbounds @simd for i in 1:N-1
        y[i] = (x[i+1] - x[i]) / tau0
    end
    Ny = length(y)

    # See `_mtotdev_greenhall` for the sequential-outer / parallel-inner
    # rationale, sliding-window inner reduction, and ext-buffer pool.
    nthreads = Threads.nthreads()
    max_m = maximum(m_values; init=0)
    ext_pool = [Vector{Float64}(undef, 3 * 3 * max(max_m, 1)) for _ in 1:nthreads]

    for k in eachindex(m_values)
        m = m_values[k]
        if m == 1
            L = N - 3
            if L <= 0
                devs[k] = NaN
                continue
            end
            sum_sq = 0.0
            @inbounds @simd for i in 1:L
                d3 = x[i+3] - 3.0*x[i+2] + 3.0*x[i+1] - x[i]
                sum_sq += d3^2
            end
            devs[k] = sqrt(sum_sq / (6.0 * L * tau0^2))
            continue
        end

        n_iter = Ny - 3m + 1
        if n_iter < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        nchunks = max(1, min(nthreads, n_iter))
        chunk_size = cld(n_iter, nchunks)
        chunk_sums = zeros(nchunks)

        Threads.@threads :dynamic for c in 1:nchunks
            i_lo = (c - 1) * chunk_size
            i_hi = min(c * chunk_size, n_iter) - 1
            ext = ext_pool[c]
            local_sum = 0.0
            for i in i_lo:i_hi
                hi = floor(Int, seg_len / 2)
                lo_start = ceil(Int, seg_len / 2) + 1

                s1 = 0.0
                @inbounds @simd for j in 1:hi
                    s1 += y[i+j]
                end
                m1 = s1 / hi

                s2 = 0.0
                @inbounds @simd for j in lo_start:seg_len
                    s2 += y[i+j]
                end
                m2 = s2 / (seg_len - lo_start + 1)

                slope = isodd(seg_len) ? (m2 - m1) / (0.5*(seg_len - 1) + 1.0) : (m2 - m1) / (0.5*seg_len)
                mid = floor(seg_len / 2)

                @inbounds for j in 1:seg_len
                    val = y[i+j] - slope * (j - 1 - mid)
                    rev_val = y[i + seg_len - j + 1] - slope * (seg_len - j - mid)

                    ext[j] = rev_val
                    ext[seg_len + j] = val
                    ext[2seg_len + j] = rev_val
                end

                a1 = 0.0
                a2 = 0.0
                a3 = 0.0
                @inbounds @simd for j in 1:m
                    a1 += ext[j]
                    a2 += ext[j+m]
                    a3 += ext[j+2m]
                end

                d3 = (a3 - 2.0*a2 + a1) / m
                sq = d3^2

                @inbounds for j in 1:(6m - 1)
                    a1 += ext[j+m]  - ext[j]
                    a2 += ext[j+2m] - ext[j+m]
                    a3 += ext[j+3m] - ext[j+2m]
                    d3 = (a3 - 2.0*a2 + a1) / m
                    sq += d3^2
                end
                local_sum += sq / (6.0 * m)
            end
            chunk_sums[c] = local_sum
        end

        dev_sum = sum(chunk_sums)
        devs[k] = sqrt(dev_sum / (6.0 * n_iter))
    end

    return devs
end

function _htotdev_linear(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window full LS detrend on the frequency series y = diff(x)/tau0
    # + per-window time-reverse extension + third-difference operator.
    # Same extension/operator shape as `_htotdev_greenhall`; only the slope
    # estimate differs (full LS slope + intercept instead of half-mean slope).
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    y = Vector{Float64}(undef, N-1)
    @inbounds @simd for i in 1:N-1
        y[i] = (x[i+1] - x[i]) / tau0
    end
    Ny = length(y)

    # See `_mtotdev_greenhall` for the sequential-outer / parallel-inner
    # rationale, sliding-window inner reduction, and ext-buffer pool.
    nthreads = Threads.nthreads()
    max_m = maximum(m_values; init=0)
    ext_pool = [Vector{Float64}(undef, 3 * 3 * max(max_m, 1)) for _ in 1:nthreads]

    for k in eachindex(m_values)
        m = m_values[k]
        if m == 1
            L = N - 3
            if L <= 0
                devs[k] = NaN
                continue
            end
            sum_sq = 0.0
            @inbounds @simd for i in 1:L
                d3 = x[i+3] - 3.0*x[i+2] + 3.0*x[i+1] - x[i]
                sum_sq += d3^2
            end
            devs[k] = sqrt(sum_sq / (6.0 * L * tau0^2))
            continue
        end

        n_iter = Ny - 3m + 1
        if n_iter < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        L_float = Float64(seg_len)
        sum_i = (L_float * (L_float + 1.0)) / 2.0
        sum_i2 = (L_float * (L_float + 1.0) * (2.0*L_float + 1.0)) / 6.0
        delta = L_float * sum_i2 - sum_i^2

        nchunks = max(1, min(nthreads, n_iter))
        chunk_size = cld(n_iter, nchunks)
        chunk_sums = zeros(nchunks)

        Threads.@threads :dynamic for c in 1:nchunks
            i_lo = (c - 1) * chunk_size
            i_hi = min(c * chunk_size, n_iter) - 1
            ext = ext_pool[c]
            local_sum = 0.0
            for i in i_lo:i_hi
                sum_y = 0.0
                sum_iy = 0.0
                @inbounds @simd for j in 1:seg_len
                    v = y[i+j]
                    sum_y += v
                    sum_iy += j * v
                end

                a = (sum_y * sum_i2 - sum_iy * sum_i) / delta
                b = (L_float * sum_iy - sum_y * sum_i) / delta

                @inbounds for j in 1:seg_len
                    val = y[i+j] - (a + b * j)
                    rev_val = y[i + seg_len - j + 1] - (a + b * (seg_len - j + 1))

                    ext[j] = rev_val
                    ext[seg_len + j] = val
                    ext[2seg_len + j] = rev_val
                end

                a1 = 0.0
                a2 = 0.0
                a3 = 0.0
                @inbounds @simd for j in 1:m
                    a1 += ext[j]
                    a2 += ext[j+m]
                    a3 += ext[j+2m]
                end

                d3 = (a3 - 2.0*a2 + a1) / m
                sq = d3^2

                @inbounds for j in 1:(6m - 1)
                    a1 += ext[j+m]  - ext[j]
                    a2 += ext[j+2m] - ext[j+m]
                    a3 += ext[j+3m] - ext[j+2m]
                    d3 = (a3 - 2.0*a2 + a1) / m
                    sq += d3^2
                end
                local_sum += sq / (6.0 * m)
            end
            chunk_sums[c] = local_sum
        end

        dev_sum = sum(chunk_sums)
        devs[k] = sqrt(dev_sum / (6.0 * n_iter))
    end

    return devs
end

"""
    _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall) → Vector{Float64}

Computes the Modified Hadamard Total Deviation (MHTOTDEV).

MHTOTDEV is novel to SigmaTau and uses a per-window time-reverse extension
on phase. `detrend` selects the per-window detrending applied before the
extension:

- `:greenhall` — half-mean slope removal (matches the MTOTDEV/HTOTDEV
  Greenhall convention). Default.
- `:linear` — full least-squares (slope + intercept) per window. Alias
  for `:legacy` on this kernel.
- `:legacy` — pre-1.0 SigmaTau behavior; identical to `:linear` here.
"""
function _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    detrend === :greenhall && return _mhtotdev_greenhall(x, m_values, tau0)
    if detrend === :linear || detrend === :legacy
        return _mhtotdev_linear(x, m_values, tau0)
    end
    throw(ArgumentError("unknown detrend recipe: $detrend; valid for MHTOTDEV: :greenhall, :linear, :legacy"))
end

function _mhtotdev_linear(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # See `_mtotdev_greenhall` for the sequential-outer / parallel-inner
    # rationale, sliding-window inner reduction, and buffer pool.
    nthreads = Threads.nthreads()
    max_m = maximum(m_values; init=0)
    max_Lp = 3 * max(max_m, 1) + 1
    max_L3 = 3 * max_Lp - 3 * max(max_m, 1)
    ext_pool = [Vector{Float64}(undef, 3 * max_Lp) for _ in 1:nthreads]
    d3_pool  = [Vector{Float64}(undef, max_L3)     for _ in 1:nthreads]

    for k in eachindex(m_values)
        m = m_values[k]
        if m < 1
            devs[k] = NaN
            continue
        end
        nsubs = N - 4m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        Lp = 3m + 1
        L3 = 3Lp - 3m
        Lp_float = Float64(Lp)
        sum_i = (Lp_float * (Lp_float + 1.0)) / 2.0
        sum_i2 = (Lp_float * (Lp_float + 1.0) * (2.0*Lp_float + 1.0)) / 6.0
        delta = Lp_float * sum_i2 - sum_i^2

        nchunks = max(1, min(nthreads, nsubs))
        chunk_size = cld(nsubs, nchunks)
        chunk_sums = zeros(nchunks)

        Threads.@threads :dynamic for c in 1:nchunks
            n_lo = (c - 1) * chunk_size + 1
            n_hi = min(c * chunk_size, nsubs)
            ext = ext_pool[c]
            d3_vec = d3_pool[c]
            local_sum = 0.0
            for n in n_lo:n_hi
                sum_x = 0.0
                sum_ix = 0.0
                @inbounds @simd for j in 1:Lp
                    val = x[n-1+j]
                    sum_x += val
                    sum_ix += j * val
                end

                a = (sum_x * sum_i2 - sum_ix * sum_i) / delta
                b = (Lp_float * sum_ix - sum_x * sum_i) / delta

                @inbounds for j in 1:Lp
                    val = x[n-1+j] - (a + b * j)
                    rev_val = x[n-1 + Lp - j + 1] - (a + b * (Lp - j + 1))

                    ext[j] = rev_val
                    ext[Lp + j] = val
                    ext[2Lp + j] = rev_val
                end

                @inbounds for j in 1:L3
                    d3_vec[j] = ext[j] - 3.0*ext[j+m] + 3.0*ext[j+2m] - ext[j+3m]
                end

                n_avg = L3 + 1 - m
                if n_avg > 0
                    # Sliding-window sum over m-windows of d3_vec, replacing
                    # the cumulative-sum buffer S used previously. A_j is
                    # the running m-window sum starting at index j.
                    A = 0.0
                    @inbounds @simd for j in 1:m
                        A += d3_vec[j]
                    end
                    block_var = A^2

                    @inbounds for j in 1:(n_avg - 1)
                        A += d3_vec[j+m] - d3_vec[j]
                        block_var += A^2
                    end
                    block_var /= (n_avg * 6.0 * Float64(m)^2)
                else
                    block_var = 0.0
                end

                local_sum += block_var
            end
            chunk_sums[c] = local_sum
        end

        total_sum = sum(chunk_sums)
        devs[k] = sqrt(total_sum / (nsubs * Float64(m)^2 * tau0^2))
    end

    return devs
end

function _mhtotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window half-mean slope removal + per-window time-reverse extension
    # + averaged third-difference operator. Same window/operator structure as
    # `_mhtotdev_linear`; only the slope estimate differs (half-mean instead
    # of full LS).
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    # See `_mtotdev_greenhall` for the sequential-outer / parallel-inner
    # rationale, sliding-window inner reduction, and buffer pool.
    nthreads = Threads.nthreads()
    max_m = maximum(m_values; init=0)
    max_Lp = 3 * max(max_m, 1) + 1
    max_L3 = 3 * max_Lp - 3 * max(max_m, 1)
    ext_pool = [Vector{Float64}(undef, 3 * max_Lp) for _ in 1:nthreads]
    d3_pool  = [Vector{Float64}(undef, max_L3)     for _ in 1:nthreads]

    for k in eachindex(m_values)
        m = m_values[k]
        if m < 1
            devs[k] = NaN
            continue
        end
        nsubs = N - 4m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        Lp = 3m + 1
        L3 = 3Lp - 3m

        nchunks = max(1, min(nthreads, nsubs))
        chunk_size = cld(nsubs, nchunks)
        chunk_sums = zeros(nchunks)

        Threads.@threads :dynamic for c in 1:nchunks
            n_lo = (c - 1) * chunk_size + 1
            n_hi = min(c * chunk_size, nsubs)
            ext = ext_pool[c]
            d3_vec = d3_pool[c]
            local_sum = 0.0
            for n in n_lo:n_hi
                half = floor(Int, Lp / 2)
                s1 = 0.0
                @inbounds @simd for j in 1:half
                    s1 += x[n-1+j]
                end
                s1 /= half

                s2 = 0.0
                @inbounds @simd for j in (half+1):Lp
                    s2 += x[n-1+j]
                end
                s2 /= (Lp - half)

                slope = (s2 - s1) / ((Lp / 2.0) * tau0)

                @inbounds for j in 1:Lp
                    val = x[n-1+j] - slope * tau0 * (j - 1)
                    rev_val = x[n-1 + Lp - j + 1] - slope * tau0 * (Lp - j)

                    ext[j] = rev_val
                    ext[Lp + j] = val
                    ext[2Lp + j] = rev_val
                end

                @inbounds for j in 1:L3
                    d3_vec[j] = ext[j] - 3.0*ext[j+m] + 3.0*ext[j+2m] - ext[j+3m]
                end

                n_avg = L3 + 1 - m
                if n_avg > 0
                    A = 0.0
                    @inbounds @simd for j in 1:m
                        A += d3_vec[j]
                    end
                    block_var = A^2

                    @inbounds for j in 1:(n_avg - 1)
                        A += d3_vec[j+m] - d3_vec[j]
                        block_var += A^2
                    end
                    block_var /= (n_avg * 6.0 * Float64(m)^2)
                else
                    block_var = 0.0
                end

                local_sum += block_var
            end
            chunk_sums[c] = local_sum
        end

        total_sum = sum(chunk_sums)
        devs[k] = sqrt(total_sum / (nsubs * Float64(m)^2 * tau0^2))
    end

    return devs
end

# kernels.jl — raw deviation kernels (Vector{Float64} in, arrays out).
#
# Merged from stab/core/{allan, hadamard, total, mtie, pdev}.jl. These are the
# internal _*_core functions; the public PhaseData/FrequencyData API lives in
# deviations.jl.


# ──────────────────────────────────────────────────────────────────────
# ── stab/core/allan.jl ──────────────────────────────────────────────────

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
        if L < 2          # need ≥2 analysis windows; one window is a single difference
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

    # Prefix sum x_cs[1] = 0, x_cs[i+1] = Σⱼ₌₁ⁱ x[j]. One length-(N+1)
    # allocation; the carry chain blocks SIMD, so no @simd here.
    x_cs = Vector{Float64}(undef, N + 1)
    x_cs[1] = 0.0
    acc = 0.0
    @inbounds for i in 1:N
        acc += x[i]
        x_cs[i+1] = acc
    end

    for (k, m) in enumerate(m_values)
        Ne = N - 3m + 1
        if Ne < 2          # need ≥2 estimates
            devs[k] = NaN
            continue
        end

        sum_sq = 0.0
        @inbounds @simd for i in 1:Ne
            d = x_cs[i+3m] - 3x_cs[i+2m] + 3x_cs[i+m] - x_cs[i]
            sum_sq += d^2
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


# ──────────────────────────────────────────────────────────────────────
# ── stab/core/hadamard.jl ───────────────────────────────────────────────

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
        if L < 2          # need ≥2 analysis windows
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

    # Prefix sum x_cs[1] = 0, x_cs[i+1] = Σⱼ₌₁ⁱ x[j]. One length-(N+1)
    # allocation; the carry chain blocks SIMD, so no @simd here.
    x_cs = Vector{Float64}(undef, N + 1)
    x_cs[1] = 0.0
    acc = 0.0
    @inbounds for i in 1:N
        acc += x[i]
        x_cs[i+1] = acc
    end

    for (k, m) in enumerate(m_values)
        Ne = N - 4m + 1
        if Ne < 2          # need ≥2 estimates; MHDEV is not a total estimator (raw 4th-difference)
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


# ──────────────────────────────────────────────────────────────────────
# ── stab/core/total.jl ──────────────────────────────────────────────────

# core/total.jl — Core Total Stability Kernels

"""
    _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Total Deviation (TOTDEV) for a set of averaging factors `m`,
using the Howe 1995 / NIST SP1065 eqn 25 endpoint mean-flip extension with no
per-window detrend — the canonical TOTDEV form (matches allantools).
"""
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    return _totdev_howe(x, m_values, tau0)
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
    _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Modified Total Deviation (MTOTDEV), using Greenhall 2003's
per-window time-reverse extension with half-mean slope removal — the canonical
MTOTDEV form.
"""
function _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    return _mtotdev_greenhall(x, m_values, tau0)
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

"""
    _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Hadamard Total Deviation (HTOTDEV) on the frequency series
`y = diff(x) / tau0`, using Greenhall 2003's per-window time-reverse extension
with half-mean slope removal on `y` — the canonical HTOTDEV form.
"""
function _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    return _htotdev_greenhall(x, m_values, tau0)
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

"""
    _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Modified Hadamard Total Deviation (MHTOTDEV).

MHTOTDEV is novel to SigmaTau; no canonical/published form exists. SigmaTau
adopts the Greenhall methodology — per-window detrend + a 3× time-reverse
extension on phase — by consistency with MTOTDEV and HTOTDEV. The detrend is
degree-matched to the third-difference kernel: the window's frequency drift
(TOTHVAR's half-average estimate on the frequency increments) and frequency
offset (MTOTDEV's half-mean slope estimate on phase) are both removed, so the
extension cannot re-admit either polynomial the kernel annihilates. With only
the degree-1 removal, residual quadratic phase makes the reflected copies
non-smooth at the boundaries and a drifting record reads several times high
at long τ.
"""
function _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    return _mhtotdev_greenhall(x, m_values, tau0)
end

function _mhtotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per-window drift + half-mean slope removal + per-window time-reverse
    # extension + averaged third-difference operator.
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

        # Per-window drift removal (degree-2 phase detrend). The third
        # difference annihilates quadratic phase inside each copy of the
        # window, but the time-reverse extension re-admits it at the copy
        # boundaries unless it is removed first — the degree-2 analog of the
        # degree-1 half-mean slope removal MTOTDEV needs for its
        # second-difference kernel. The drift estimate is TOTHVAR's
        # half-average of the window's frequency increments (exact for a
        # deterministic drift), computed O(1) from telescoping phase sums.
        Ly     = 3m                          # frequency increments per window
        hy     = fld(Ly, 2)
        lo0    = cld(Ly, 2) + 1              # skip the middle sample when Ly is odd
        ddenom = isodd(Ly) ? 0.5 * (Ly - 1) + 1.0 : 0.5 * Ly
        dmid   = fld(Ly, 2)

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
                # τ0·(half-averages of the window's frequency increments);
                # the increment sums telescope to phase differences.
                ty1 = (x[n-1+hy+1] - x[n]) / hy
                ty2 = (x[n-1+Lp] - x[n-1+lo0]) / (Ly - lo0 + 1)
                sd  = (ty2 - ty1) / ddenom   # τ0 · per-sample drift in y
                # Quadratic phase implied by that drift, subtracted below:
                # q(j) = sd · (½(j−1)(j−2) − dmid·(j−1)).

                half = floor(Int, Lp / 2)
                s1 = 0.0
                @inbounds @simd for j in 1:half
                    s1 += x[n-1+j] - sd * (0.5 * (j - 1) * (j - 2) - dmid * (j - 1))
                end
                s1 /= half

                s2 = 0.0
                @inbounds @simd for j in (half+1):Lp
                    s2 += x[n-1+j] - sd * (0.5 * (j - 1) * (j - 2) - dmid * (j - 1))
                end
                s2 /= (Lp - half)

                slope = (s2 - s1) / ((Lp / 2.0) * tau0)

                @inbounds for j in 1:Lp
                    jr = Lp - j + 1
                    val     = x[n-1+j]  - sd * (0.5 * (j - 1) * (j - 2) - dmid * (j - 1)) -
                              slope * tau0 * (j - 1)
                    rev_val = x[n-1+jr] - sd * (0.5 * (jr - 1) * (jr - 2) - dmid * (jr - 1)) -
                              slope * tau0 * (jr - 1)

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


# ──────────────────────────────────────────────────────────────────────
# ── stab/core/mtie.jl ───────────────────────────────────────────────────

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


# ──────────────────────────────────────────────────────────────────────
# ── stab/core/pdev.jl ───────────────────────────────────────────────────

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

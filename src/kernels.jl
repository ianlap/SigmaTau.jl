# Deviation kernels: Vector{Float64} in, deviations out. The public API in
# deviations.jl wraps these.

"""
    _adev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; overlapping=true) → Vector{Float64}

Allan deviation at each averaging factor `m`: overlapping, or with
`overlapping=false` the original non-overlapped estimator, whose second
differences start every `m` samples (SP1065 §5.2.2). Allocation-free SIMD
inner loop.
"""
function _adev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; overlapping::Bool=true)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    for (k, m) in enumerate(m_values)
        idx = 1:(overlapping ? 1 : m):(N - 2m)
        n = length(idx)
        if n < 2          # need ≥2 analysis windows; one window is a single difference
            devs[k] = NaN
            continue
        end

        sum_sq = 0.0
        @inbounds @simd for i in idx
            d2 = x[i+2m] - 2x[i+m] + x[i]
            sum_sq += d2^2
        end

        devs[k] = sqrt(sum_sq / (2.0 * n * Float64(m)^2 * tau0^2))
    end

    return devs
end

# Prefix sums x_cs[1] = 0, x_cs[i+1] = Σⱼ₌₁ⁱ x[j] of the phase with its
# least-squares line removed. The modified estimators are exact second
# (MDEV) and third (MHDEV) differences of x, so a linear phase ramp, i.e. a
# frequency offset, contributes nothing to them; removing it before the
# prefix sums costs nothing in exactness and keeps the sums from growing like
# N²·offset. Without this, a record with an uncorrected 1e-9 frequency offset
# loses about seven significant digits to cancellation in the prefix sums.
# One length-(N+1) allocation; the carry chain blocks SIMD, so no @simd.
function _prefix_sums_detrended(x::Vector{Float64})
    N = length(x)
    # Least-squares line through (i, x[i]), i = 1..N, in one pass.
    sx = 0.0; sxi = 0.0
    @inbounds @simd for i in 1:N
        sx  += x[i]
        sxi += x[i] * i
    end
    mean_i = (N + 1) / 2
    var_i  = (N^2 - 1) / 12 * N          # Σ (i - mean_i)²
    slope  = N > 1 ? (sxi - mean_i * sx) / var_i : 0.0
    icpt   = sx / N - slope * mean_i
    x_cs = Vector{Float64}(undef, N + 1)
    x_cs[1] = 0.0
    acc = 0.0
    @inbounds for i in 1:N
        acc += x[i] - (icpt + slope * i)
        x_cs[i+1] = acc
    end
    return x_cs
end

"""
    _mdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Modified Allan deviation at each `m`, as a third difference of prefix sums
(O(N) per `m`).
"""
function _mdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    x_cs = _prefix_sums_detrended(x)

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

Time deviation (τ/√3) · MDEV.
"""
function _tdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    mdevs = _mdev_core(x, m_values, tau0)
    taus = m_values .* tau0
    return taus .* mdevs ./ sqrt(3.0)
end


"""
    _hdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; overlapping=true) → Vector{Float64}

Hadamard deviation at each averaging factor `m`: overlapping, or with
`overlapping=false` non-overlapped, third differences every `m` samples
(SP1065 §5.2.8). Allocation-free SIMD inner loop.
"""
function _hdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; overlapping::Bool=true)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    for (k, m) in enumerate(m_values)
        idx = 1:(overlapping ? 1 : m):(N - 3m)
        n = length(idx)
        if n < 2          # need ≥2 analysis windows
            devs[k] = NaN
            continue
        end

        sum_sq = 0.0
        @inbounds @simd for i in idx
            d3 = x[i+3m] - 3.0 * x[i+2m] + 3.0 * x[i+m] - x[i]
            sum_sq += d3^2
        end

        devs[k] = sqrt(sum_sq / (6.0 * n * Float64(m)^2 * tau0^2))
    end

    return devs
end

"""
    _mhdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Modified Hadamard deviation at each `m`, as a fourth difference of prefix
sums (O(N) per `m`).
"""
function _mhdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    x_cs = _prefix_sums_detrended(x)

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


"""
    _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Total deviation at each averaging factor `m`, with the Howe 1995 / NIST SP1065
eq. 25 endpoint mean-flip extension and no per-window detrend (the canonical
form; matches allantools).
"""
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
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

# ── Shared machinery for the Greenhall total deviations ──────────────────
#
# MTOTDEV, HTOTDEV and MHTOTDEV all take, for each window of L samples, the
# window minus a fitted ramp (`val`), extend it to [reverse(val); val;
# reverse(val)], and sum the squares of an m-spaced difference of m-sample
# box sums over that extension. A box sum is a difference of the cumulative
# sum E of the extension, and E is piecewise a cumulative sum V of `val`
# read forward or backward:
#
#   E(q) = V(L) − V(L−q)        0  ≤ q ≤ L
#   E(q) = V(L) + V(q−L)        L  ≤ q ≤ 2L
#   E(q) = 3V(L) − V(3L−q)      2L ≤ q ≤ 3L
#
# with V(k) = W(k) − c·Q(k), W the window's own cumulative sum and Q the
# ramp's. W comes from a compensated (double-double) prefix sum of the whole
# record, so each window costs one pass to build V and its reverse and one
# vectorised pass over the box differences, where the per-sample running box
# sums this replaces carried a serial dependency through every step. The
# extension's mirror symmetry halves the second pass.

# `v` minus its least-squares line. MTOTDEV, HTOTDEV and MHTOTDEV are
# unchanged by adding a line to the phase (the per-window frequency removal
# takes out its slope and the differences take out the rest), so the
# kernels work on the residual, which keeps the prefix sums below small.
function _minus_line(v::Vector{Float64})
    n = length(v)
    n < 2 && return copy(v)
    c = (n + 1) / 2.0
    sy = 0.0; sty = 0.0
    @inbounds for i in 1:n
        sy += v[i]; sty += (i - c) * v[i]
    end
    b = sty / (n * (Float64(n)^2 - 1) / 12.0)
    a = sy / n
    return Float64[v[i] - (a + b * (i - c)) for i in 1:n]
end

# Prefix sums S[k+1] = Σ_{i≤k} v[i] as hi + lo (TwoSum compensation), so a
# window sum S[b] − S[a] keeps full precision on long records.
function _prefix_dd(v::Vector{Float64})
    n = length(v)
    hi = Vector{Float64}(undef, n + 1)
    lo = Vector{Float64}(undef, n + 1)
    hi[1] = 0.0; lo[1] = 0.0
    s = 0.0; c = 0.0
    @inbounds for i in 1:n
        a = v[i]
        t = s + a
        bp = t - s
        c += (s - (t - bp)) + (a - bp)
        s = t
        hi[i+1] = s; lo[i+1] = c
    end
    return hi, lo
end

# V(k), k = 0..L, into Vb[1..L+1] for the window starting after sample `o`:
# W(k) = S[o+1+k] − S[o+1], minus c·(k(k−1)/2 − mid·k). The window mean is
# taken out as well (every statistic here differences it away), which keeps
# E small next to the differences taken of it when the record has an offset.
@inline function _fill_window_cumsum!(Vb::Vector{Float64}, hi::Vector{Float64}, lo::Vector{Float64},
                                      o::Int, L::Int, c::Float64, mid::Float64)
    h0 = hi[o+1]; l0 = lo[o+1]
    Lf = Float64(L)
    mean = (((hi[o+1+L] - h0) + (lo[o+1+L] - l0)) - c * (0.5 * Lf * (Lf - 1.0) - mid * Lf)) / Lf
    @inbounds @simd for k in 0:L
        kf = Float64(k)
        Vb[k+1] = ((hi[o+1+k] - h0) + (lo[o+1+k] - l0)) - (c * (0.5 * kf * (kf - 1.0) - mid * kf) + mean * kf)
    end
    return Vb
end

# The extension's cumulative sum E is never stored. With Vr the reverse of
# V (Vr[t+1] = V(L−t)), each piece of E is a forward read:
#
#   piece 1 (q ≤ L):       E(q) = V(L) − Vr[q+1]
#   piece 2 (L ≤ q ≤ 2L):  E(q) = V(L) + V[q−L+1]
#   piece 3 (2L ≤ q):      E(q) = 3V(L) − Vr[q−2L+1]
#
# so a box difference Σ_p w_p E(j−1+p·m) over a run of starts j where each
# term stays in one piece is a unit-stride loop over two arrays.
@inline _extension_piece(q::Int, L::Int) = q <= L ? 1 : (q <= 2L ? 2 : 3)

# (constant, coefficient, use Vr?, index offset from j) for term p.
@inline function _extension_term(piece::Int, VL::Float64, L::Int, p::Int, m::Int, w::Float64)
    piece == 1 && return (w * VL, -w, true, p * m)
    piece == 2 && return (w * VL, w, false, p * m - L)
    return (3.0 * w * VL, -w, true, p * m - 2L)
end

# Σ_{j=ja}^{jb} (Σ_p w_p E(j−1+p·m))² for a run of starts in which no term
# changes piece (each term's piece is read at jb).
function _sumsq_box_run(w::NTuple{K,Float64}, Vb::Vector{Float64}, Vr::Vector{Float64},
                        L::Int, m::Int, ja::Int, jb::Int) where {K}
    ja > jb && return 0.0
    VL = Vb[L+1]
    t = ntuple(p -> _extension_term(_extension_piece(jb - 1 + (p - 1) * m, L), VL, L, p - 1, m, w[p]), Val(K))
    cst = 0.0
    for p in 1:K
        cst += t[p][1]
    end
    a = ntuple(p -> t[p][3] ? Vr : Vb, Val(K))
    c = ntuple(p -> t[p][2], Val(K))
    o = ntuple(p -> t[p][4], Val(K))
    s = 0.0
    @inbounds @simd for j in ja:jb
        d = cst
        for p in 1:K
            d += c[p] * a[p][j+o[p]]
        end
        s += d * d
    end
    return s
end

# The same over any run of starts, split where a term crosses a seam.
function _sumsq_boxes(w::NTuple{K,Float64}, Vb::Vector{Float64}, Vr::Vector{Float64},
                      L::Int, m::Int, ja::Int, jb::Int) where {K}
    s = 0.0
    lo = ja
    while lo <= jb
        hi = jb
        for p in 0:K-1
            for seam in (L, 2L)
                jc = seam + 2 - p * m        # first start with q = seam + 1
                (jc > lo && jc <= hi) && (hi = jc - 1)
            end
        end
        s += _sumsq_box_run(w, Vb, Vr, L, m, lo, hi)
        lo = hi + 1
    end
    return s
end

# The extension is mirror-symmetric about both seams, so a box difference
# starting at j and one at its mirror image are equal up to sign. Over a run
# of starts that the mirror maps onto itself, sum the first half twice and
# the middle start once.
function _sumsq_boxes_mirrored(w::NTuple{K,Float64}, Vb::Vector{Float64}, Vr::Vector{Float64},
                               L::Int, m::Int, ja::Int, jb::Int) where {K}
    n = jb - ja + 1
    n <= 0 && return 0.0
    h = n ÷ 2
    s = 2.0 * _sumsq_boxes(w, Vb, Vr, L, m, ja, ja + h - 1)
    isodd(n) && (s += _sumsq_boxes(w, Vb, Vr, L, m, ja + h, ja + h))
    return s
end

@inline function _reverse_into!(Vr::Vector{Float64}, Vb::Vector{Float64}, L::Int)
    @inbounds @simd for t in 0:L
        Vr[t+1] = Vb[L-t+1]
    end
    return Vr
end

const _W_D2 = (-1.0, 3.0, -3.0, 1.0)          # E(j−1), E(j−1+m), …, E(j−1+3m)
const _W_D3 = (1.0, -4.0, 6.0, -4.0, 1.0)     # E(j−1), …, E(j−1+4m)

# Run `window_sum(n_lo, n_hi, Vb, Vr)` over `nwin` windows split into
# per-thread chunks, each with its own V buffers for windows of Lmax samples.
function _sum_over_windows(window_sum::F, nwin::Int, Lmax::Int) where {F}
    nchunks = max(1, min(Threads.nthreads(), nwin))
    chunk_size = cld(nwin, nchunks)
    chunk_sums = zeros(nchunks)
    Threads.@threads :dynamic for c in 1:nchunks
        n_lo = (c - 1) * chunk_size + 1
        n_hi = min(c * chunk_size, nwin)
        Vb = Vector{Float64}(undef, Lmax + 1)
        Vr = Vector{Float64}(undef, Lmax + 1)
        chunk_sums[c] = n_lo <= n_hi ? window_sum(n_lo, n_hi, Vb, Vr) : 0.0
    end
    return sum(chunk_sums)
end

"""
    _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Modified total deviation, using Greenhall 2003's
per-window time-reverse extension with half-mean slope removal, the canonical
MTOTDEV form.
"""
function _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Greenhall 2003: per window of 3m phase samples, remove the half-mean
    # frequency, extend by reflection, and average the squared second
    # difference of m-sample means over the 6m positions of the extension.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))
    hi, lo = _prefix_dd(_minus_line(x))
    for k in eachindex(m_values)
        m = m_values[k]
        nsubs = N - 3m + 1
        if m < 1 || nsubs < 1
            devs[k] = NaN
            continue
        end
        L = 3m
        # Half-average frequency (Stable32, allantools): the means of the
        # first and last ⌊3m/2⌋ samples, a distance 3m − ⌊3m/2⌋ apart; for
        # odd m the middle sample is in neither half.
        kh = L ÷ 2
        inv_span = 1.0 / (Float64(kh) * (L - kh))
        outer = _sum_over_windows(nsubs, L) do n_lo, n_hi, Vb, Vr
            acc = 0.0
            for n in n_lo:n_hi
                o = n - 1
                c = @inbounds begin
                    s1 = (hi[o+1+kh] - hi[o+1]) + (lo[o+1+kh] - lo[o+1])
                    s2 = (hi[o+1+L] - hi[o+1+L-kh]) + (lo[o+1+L] - lo[o+1+L-kh])
                    (s2 - s1) * inv_span
                end
                _fill_window_cumsum!(Vb, hi, lo, o, L, c, 0.0)
                _reverse_into!(Vr, Vb, L)
                # Starts 1..3m+1 mirror about the first seam, 3m+2..6m about the second.
                sq = _sumsq_boxes_mirrored(_W_D2, Vb, Vr, L, m, 1, 3m + 1) +
                     _sumsq_boxes_mirrored(_W_D2, Vb, Vr, L, m, 3m + 2, 6m)
                acc += sq / (Float64(m)^2 * 6.0 * m)
            end
            acc
        end
        devs[k] = sqrt(outer / (2.0 * Float64(m)^2 * tau0^2 * nsubs))
    end
    return devs
end

"""
    _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Hadamard total deviation on the frequency series
`y = diff(x) / tau0`, using Greenhall 2003's per-window time-reverse extension
with half-mean slope removal on `y`, the canonical HTOTDEV form.
"""
function _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))
    Ny = N - 1
    # Cumulative sums of y = diff(x)/τ0 are phase differences, so the
    # window sums need no prefix array: Σ y[i+1..i+k] = (x[i+1+k] − x[i+1])/τ0.
    # The compensated prefix machinery is shared, so pass x/τ0 as `hi`
    # and zeros as `lo`; x minus its line, as y's mean drops out too.
    xs = _minus_line(x) ./ tau0
    zlo = zeros(N)
    for k in eachindex(m_values)
        m = m_values[k]
        if m == 1
            L1 = N - 3
            if L1 <= 0
                devs[k] = NaN
                continue
            end
            sum_sq = 0.0
            @inbounds @simd for i in 1:L1
                d3 = x[i+3] - 3.0*x[i+2] + 3.0*x[i+1] - x[i]
                sum_sq += d3^2
            end
            devs[k] = sqrt(sum_sq / (6.0 * L1 * tau0^2))
            continue
        end
        n_iter = Ny - 3m + 1
        if m < 1 || n_iter < 1
            devs[k] = NaN
            continue
        end
        L = 3m
        hi_half = L ÷ 2
        lo_start = cld(L, 2) + 1
        nlo = L - lo_start + 1
        denom = isodd(L) ? 0.5 * (L - 1) + 1.0 : 0.5 * L
        mid = Float64(L ÷ 2)
        total = _sum_over_windows(n_iter, L) do n_lo, n_hi, Vb, Vr
            acc = 0.0
            for n in n_lo:n_hi
                i = n - 1                       # window y[i+1 .. i+L]
                @inbounds begin
                    m1 = (xs[i+1+hi_half] - xs[i+1]) / hi_half
                    m2 = (xs[i+1+L] - xs[i+lo_start]) / nlo
                end
                slope = (m2 - m1) / denom
                _fill_window_cumsum!(Vb, xs, zlo, i, L, slope, mid)
                _reverse_into!(Vr, Vb, L)
                # Starts 1..3m+1 mirror about the first seam, 3m+2..6m about the second.
                sq = _sumsq_boxes_mirrored(_W_D2, Vb, Vr, L, m, 1, 3m + 1) +
                     _sumsq_boxes_mirrored(_W_D2, Vb, Vr, L, m, 3m + 2, 6m)
                acc += sq / (Float64(m)^2 * 6.0 * m)
            end
            acc
        end
        devs[k] = sqrt(total / (6.0 * n_iter))
    end
    return devs
end

"""
    _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Computes the Modified Hadamard Total Deviation (MHTOTDEV).

MHTOTDEV is novel to SigmaTau; no canonical/published form exists. SigmaTau
adopts the Greenhall methodology — per-window detrend + a 3× time-reverse
extension on phase — by consistency with MTOTDEV and HTOTDEV, with one
deliberate departure: the per-window detrend is a QUADRATIC least-squares
fit, not the half-mean slope removal the second-difference totals use.

The detrend degree must match the kernel's annihilation degree (the third
difference annihilates phase polynomials of degree ≤ 2, as MTOT's second
difference annihilates degree ≤ 1 and its linear detrend matches). With
only a slope removal, the stochastic quadratic carried by noise at
α ≤ −3 survives into the time-reverse extension, the seam kink feeds the
boundary-spanning third differences, and the estimator inherits MVAR's
divergence (measured: log-log slope 2 instead of 3 at α = −4, level ∝ T).
The matched quadratic detrend restores convergence over the full Hadamard
range α > −5 (slope 3.00, level N-stable) while leaving the PM statistics
healthy (B(WHPM) = 0.90, flat in m and N) — verified by two independent
implementations. It also removes any deterministic frequency offset and
drift exactly within each window, so this kernel IS drift-immune on its
own; the wrapper's record-level removal (`remove_drift=true` default)
remains for noise-identification hygiene, not for the statistic.

Point (odd) reflection is not used: with it B(WHPM) grows ∝ m (≈ 3.8) and
the χ² EDF collapses; the even (mirror) reflection shows neither.
"""
function _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    # Per window of 4m phase samples: least-squares quadratic removal,
    # reflection extension, and the squared m-sample box sum of the third
    # difference at each of the 8m + 1 positions of the extension.
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))
    r = _minus_line(x)          # a line is absorbed by every window's quadratic
    hi, lo = _prefix_dd(r)
    for k in eachindex(m_values)
        m = m_values[k]
        nsubs = N - 4m + 1
        if m < 1 || nsubs < 1
            devs[k] = NaN
            continue
        end
        L = 4m
        n_avg = 8m + 1
        p1, p2, P1, P2, s11, s22 = _quad_basis(L)
        total = _sum_over_windows(nsubs, L) do n_lo, n_hi, Vb, Vr
            acc = 0.0
            for n in n_lo:n_hi
                o = n - 1
                _fill_window_cumsum_quad!(Vb, hi, lo, r, o, L, p1, p2, P1, P2, s11, s22)
                _reverse_into!(Vr, Vb, L)
                # Starts 1..4m+1 and 4m+2..8m mirror onto themselves; 8m+1 is the
                # mirror of 4m+1 across the second seam but is taken as is.
                sq = _sumsq_boxes_mirrored(_W_D3, Vb, Vr, L, m, 1, 4m + 1) +
                     _sumsq_boxes_mirrored(_W_D3, Vb, Vr, L, m, 4m + 2, 8m) +
                     _sumsq_boxes(_W_D3, Vb, Vr, L, m, 8m + 1, 8m + 1)
                acc += sq / (n_avg * 6.0 * Float64(m)^2)
            end
            acc
        end
        devs[k] = sqrt(total / (nsubs * Float64(m)^2 * tau0^2))
    end
    return devs
end

# Orthogonal quadratic basis on the window grid t = 0..L−1: p1 = t − t̄ and
# p2 = p1² − ⟨p1²⟩, with their cumulative sums P1, P2 (P[k+1] = Σ_{j≤k} p[j])
# and norms. Orthogonality turns each window's least-squares fit into three
# independent dot products.
function _quad_basis(L::Int)
    tbar = (L - 1) / 2.0
    p1 = Float64[(j - 1) - tbar for j in 1:L]
    s11 = sum(abs2, p1)
    c2 = s11 / L
    p2 = Float64[p1[j]^2 - c2 for j in 1:L]
    s22 = sum(abs2, p2)
    P1 = zeros(L + 1); P2 = zeros(L + 1)
    for j in 1:L
        P1[j+1] = P1[j] + p1[j]
        P2[j+1] = P2[j] + p2[j]
    end
    return p1, p2, P1, P2, s11, s22
end

# V(k), k = 0..L, for the window starting after sample `o`: the cumulative
# sum of the window minus its least-squares quadratic a0 + a1·p1 + a2·p2.
@inline function _fill_window_cumsum_quad!(Vb::Vector{Float64}, hi::Vector{Float64}, lo::Vector{Float64},
                                           r::Vector{Float64}, o::Int, L::Int,
                                           p1::Vector{Float64}, p2::Vector{Float64},
                                           P1::Vector{Float64}, P2::Vector{Float64},
                                           s11::Float64, s22::Float64)
    s1 = 0.0; s2 = 0.0
    @inbounds @simd for j in 1:L
        w = r[o+j]
        s1 += w * p1[j]
        s2 += w * p2[j]
    end
    h0 = hi[o+1]; l0 = lo[o+1]
    a0 = ((hi[o+1+L] - h0) + (lo[o+1+L] - l0)) / L
    a1 = s1 / s11
    a2 = s2 / s22
    @inbounds @simd for k in 0:L
        Vb[k+1] = ((hi[o+1+k] - h0) + (lo[o+1+k] - l0)) - (a0 * Float64(k) + a1 * P1[k+1] + a2 * P2[k+1])
    end
    return Vb
end


# ── Maximum time interval error ──────────────────────────────────────────

"""
    _mtie_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Maximum Time Interval Error per ITU-T G.810. For each averaging factor
`m`, slides a window of `m+1` phase samples (spanning `τ = m·τ₀`) across
the record and returns the largest peak-to-peak phase excursion observed
in any window.

Units are seconds (a σ_x quantity, like TDEV); MTIE is itself a phase
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


"""
    _tierms_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

RMS Time Interval Error per NIST SP1065 §5.2.18 eq. 37. For each
averaging factor `m`, the root-mean-square of the phase change over all
`N − m` spans of `τ = m·τ₀`:

```math
\\mathrm{TIE\\,rms}(m\\tau_0) = \\sqrt{\\frac{1}{N-m}
    \\sum_{i=1}^{N-m} \\bigl(x_{i+m} - x_i\\bigr)^2}
```

Units are seconds (a σ_x quantity, like MTIE/TDEV); TIE rms is a phase
measure, so no τ rescaling is applied. Matches allantools' `tierms`
(its per-window max−min of the two samples `x[i]`, `x[i+m]` is exactly
`|x[i+m] − x[i]|`). Returns NaN where `m ≥ N`.
"""
function _tierms_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    for (k, m) in enumerate(m_values)
        L = N - m
        if L <= 0 || m < 1
            devs[k] = NaN
            continue
        end
        s = 0.0
        @inbounds for i in 1:L
            d = x[i + m] - x[i]
            s += d * d
        end
        devs[k] = sqrt(s / L)
    end

    return devs
end


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

For `m = 1` the parabolic weight collapses to zero, so the kernel falls back
to overlapping ADEV (the canonical PDEV(τ₀) ≡ ADEV(τ₀) identity from
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
            L1 = N - 2
            if L1 < 2
                devs[k] = NaN
                continue
            end
            sum_sq1 = 0.0
            @inbounds @simd for i in 1:L1
                d2 = x[i+2] - 2x[i+1] + x[i]
                sum_sq1 += d2^2
            end
            devs[k] = sqrt(sum_sq1 / (2.0 * L1 * tau0^2))
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


"""
    _theo1_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64) → Vector{Float64}

Raw Thêo1 deviation per NIST SP1065 (Riley 2008) §5.2.15
eq. 30, with the 0.75 normalization of the original Howe & Peppler 2003
definition (used by Stable32 and allantools; SP1065's printed eq. 30
omits the 0.75, but its own Table 2 bias value a = 1.00 for white FM
only holds with the factor included):

```math
\\mathrm{Theo1}(m, \\tau_0, N) = \\frac{1}{0.75\\,(N-m)\\,(m\\tau_0)^2}
   \\sum_{i=1}^{N-m} \\sum_{\\delta=0}^{m/2-1} \\frac{1}{m/2-\\delta}
   \\bigl[(x_i - x_{i-\\delta+m/2}) + (x_{i+m} - x_{i+\\delta+m/2})\\bigr]^2
```

Defined for even `m` with `2 ≤ m ≤ N − 1` (SP1065 states the statistic
for `m ≥ 10`; smaller even factors are computed but carry little Theo1
benefit over ADEV). Odd or out-of-range `m` yield NaN. The effective
averaging time of a Thêo1 value is `τ = 0.75·m·τ0` (SP1065 §5.2.15).

The double sum is O((N−m)·m/2) per averaging factor, that is O(N·m), unlike
the O(N) Allan kernels; the largest factors dominate, so a full even-m
sweep on a long record costs O(N³) overall. The loop is δ-major with a
SIMD inner i-loop and no allocations; entries of `m_values` are
processed in parallel across threads.
"""
function _theo1_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    devs = Vector{Float64}(undef, length(m_values))
    # Thêo1 is unchanged by a phase offset. Measuring from the first sample
    # (exact for nearby values) keeps the pair sums small next to their
    # differences on records with a large offset.
    x1 = isempty(x) ? 0.0 : x[1]
    xr = x .- x1
    Threads.@threads :dynamic for k in eachindex(m_values)
        devs[k] = _theo1_dev_at(xr, m_values[k], tau0)
    end
    return devs
end

# Single-m Thêo1 deviation; NaN for odd m or m outside [2, N−1].
#
# Each term is (x_i + x_{i+m}) − (x_{c−δ} + x_{c+δ}) with c = i + m/2, so the
# pair sum A_i = x_i + x_{i+m} is formed once per m, and two δ are taken per
# pass over i, which reads 2.5 samples per term instead of 4.
function _theo1_dev_at(x::Vector{Float64}, m::Int, tau0::Float64)
    N = length(x)
    (m < 2 || isodd(m) || m > N - 1) && return NaN

    half = m ÷ 2
    L = N - m
    A = Vector{Float64}(undef, L)
    @inbounds @simd for i in 1:L
        A[i] = x[i] + x[i + m]
    end
    total = 0.0
    delta = 0
    @inbounds while delta + 1 <= half - 1
        s0 = 0.0
        s1 = 0.0
        @simd for i in 1:L
            c = i + half
            d0 = A[i] - (x[c - delta] + x[c + delta])
            d1 = A[i] - (x[c - delta - 1] + x[c + delta + 1])
            s0 += d0 * d0
            s1 += d1 * d1
        end
        total += s0 / (half - delta) + s1 / (half - delta - 1)
        delta += 2
    end
    @inbounds if delta <= half - 1
        s0 = 0.0
        @simd for i in 1:L
            c = i + half
            d0 = A[i] - (x[c - delta] + x[c + delta])
            s0 += d0 * d0
        end
        total += s0 / (half - delta)
    end

    return sqrt(total / (0.75 * L * (m * tau0)^2))
end

# Overlapping Allan *variance* at a single averaging factor, allowing a
# single analysis window (L ≥ 1). Used only for the ThêoBR ratio terms,
# where SP1065 eq. 33 reaches AVAR factors as large as m ≈ N/2, one
# window is acceptable inside an average over many ratio terms, while
# the public `_adev_core` keeps its stricter L ≥ 2 guard.
function _avar_var_at(x::Vector{Float64}, m::Int, tau0::Float64)
    N = length(x)
    L = N - 2m
    (m < 1 || L < 1) && return NaN
    s = 0.0
    @inbounds @simd for i in 1:L
        d2 = x[i + 2m] - 2.0 * x[i + m] + x[i]
        s += d2 * d2
    end
    return s / (2.0 * L * (m * tau0)^2)
end

# Thêo1 factors the fast ThêoBR ratio computes exactly: the first
# _THEOBR_DENSE of the ladder and _THEOBR_SAMPLES spread over the rest.
const _THEOBR_DENSE = 16
const _THEOBR_SAMPLES = 32

"""
    _theobr_ratio(x::Vector{Float64}, tau0::Float64; exact=false, sp1065=false) → Float64

ThêoBR automatic bias-removal factor [howe-2006-theoh](@cite), eq. 6: the
average over `i = 0…n`, `n = ⌊N/30 − 3⌋`, of the variance ratios
`AVAR(m = 9+3i) / Theo1(m = 12+4i)`. The two factor ladders sample the same
effective τ (`(9+3i)·τ0 = 0.75·(12+4i)·τ0`), so each term is an AVAR/Theo1
comparison at matched averaging time, out to a tenth of the record where
AVAR is still well determined. Multiply a Thêo1 *variance* by this factor (a
deviation by its square root) to get the bias-removed ThêoBR estimate.

[mcgee-2007-theoh](@cite) and [taylor-2010-fast-theobr](@cite) give the same
`n`. SP1065 §5.2.16 eq. 33 prints `⌊N/6 − 3⌋`, which runs AVAR out to half
the record, where it has one degree of freedom: on simulated white, flicker
and random-walk FM that ladder scatters the ratio by 25 to 30 % and biases
the corrected variance by −9 to +11 %, against 8 to 11 % and −5 to +3 % here
(N = 1024). `sp1065=true` takes that ladder, which is the one Stable32 runs.
The paper's formula needs `N ≥ 90` for `n ≥ 0`; shorter records, outside its
stated range, use the first rung alone.

Unless `exact=true` or the ladder is short, Thêo1 is computed only at the
first 16 factors and at 32 evenly spaced along the rest, and interpolated
between them linearly in log m–log σ² (exact for a power law); AVAR is
computed at every factor. Against the exact ratio on simulated power-law
noises (α = +2…−2, N = 8192 and 16 384) this stays within 0.5 % (0.25 % on
the deviation), far inside ThêoBR's own confidence interval. The
block-averaging form of [taylor-2010-fast-theobr](@cite) was tried and
rejected here because it moved the ratio by 10–40 % for the PM noises.

Ratio terms whose Thêo1 is not finite/positive are skipped and the average
is taken over the computable terms. Returns NaN when no term is computable
(`N < 19`).
"""
function _theobr_ratio(x::Vector{Float64}, tau0::Float64; exact::Bool=false, sp1065::Bool=false)
    N = length(x)
    N < 19 && return NaN
    n = max(N ÷ (sp1065 ? 6 : 30) - 3, 0)

    # Factors i at which Thêo1 is computed: all of them, or the first 16
    # plus evenly spaced samples (always including the last).
    idx = if exact || n + 1 <= _THEOBR_DENSE + _THEOBR_SAMPLES
        collect(0:n)
    else
        sort!(unique!(vcat(collect(0:_THEOBR_DENSE-1),
                           round.(Int, range(0, n; length=_THEOBR_SAMPLES)))))
    end
    tvars = _theo1_core(x, [12 + 4i for i in idx], tau0) .^ 2
    ok = [isfinite(t) && t > 0 for t in tvars]

    avars = Vector{Float64}(undef, n + 1)
    Threads.@threads :dynamic for i in 0:n
        avars[i+1] = _avar_var_at(x, 9 + 3i, tau0)
    end

    total = 0.0
    cnt = 0
    j = 1                                  # idx[j] ≤ i < idx[j+1]
    for i in 0:n
        while j < length(idx) && idx[j+1] <= i
            j += 1
        end
        tvar = if idx[j] == i
            ok[j] ? tvars[j] : NaN
        elseif j < length(idx) && ok[j] && ok[j+1]
            l0 = log(12.0 + 4idx[j]); l1 = log(12.0 + 4idx[j+1])
            exp(log(tvars[j]) + (log(tvars[j+1]) - log(tvars[j])) * (log(12.0 + 4i) - l0) / (l1 - l0))
        else
            NaN
        end
        avar = avars[i+1]
        if isfinite(avar) && isfinite(tvar) && tvar > 0
            total += avar / tvar
            cnt += 1
        end
    end

    return cnt == 0 ? NaN : total / cnt
end

# The ThêoBR factor as Stable32 computes it: SP1065's ladder (eq. 33) on the
# record, or on every `b`-th phase when it has more than 3000 points
# (Stable32's `MaxTheoBRPts`), plus 0.00055·(b − 1). The decimation and the
# added term are inferred from Stable32's output, which they reproduce to its
# five digits for b = 1 (1001 points) and b = 3 (8192 points).
function _stable32_theobr_ratio(x::Vector{Float64}, tau0::Float64)
    b = cld(length(x), 3000)
    return _theobr_ratio(x[1:b:end], b * tau0; exact=true, sp1065=true) + 0.00055 * (b - 1)
end

# ThêoH segment rule (SP1065 §5.2.16 eq. 34): requested grid points are
# τ-grid factors (target τ = m·τ0). Points at or below the crossover
# k = 20 % of the record span T = (N−1)·τ0 come from overlapping ADEV at
# that m; points above come from ThêoBR at the even Thêo1 factor whose
# effective τ = 0.75·m_theo·τ0 best matches m·τ0.
_theoh_is_adev(m::Int, N::Int) = m <= 0.2 * (N - 1)

# Even Thêo1 factor nearest to m/0.75 = 4m/3 (exact whenever 3 | m).
_theoh_theo_m(m::Int) = 2 * max(1, round(Int, m / 1.5))


"""
    _neff_counts(dev_sym::Symbol, N::Int, m_values::Vector{Int}) → Vector{Int}

Number of analysis windows the kernel for `dev_sym` averages at each
averaging factor `m`, the Stable32 "#" / allantools `ns` count. Mirrors the
loop bounds and small-`N` guards of the `_*_core` kernels above exactly;
entries are 0 precisely where the kernel returns NaN. The wrapper-only
deviations (`tdev`, `htdev`, `ttotdev`) inherit the wrapped estimator's
counts and never reach here.
"""
_neff_counts(dev_sym::Symbol, N::Int, m_values::Vector{Int}) =
    Int[_neff_count(dev_sym, N, m) for m in m_values]

function _neff_count(dev_sym::Symbol, N::Int, m::Int)
    if dev_sym === :adev || dev_sym === :pdev
        # _adev_core: i in 1:N−2m, NaN when N−2m < 2. _pdev_core shares the
        # bound: m = 1 falls back to _adev_core, m ≥ 2 sums M = N−2m windows
        # (plus an explicit NaN guard for m < 1).
        dev_sym === :pdev && m < 1 && return 0
        L = N - 2m
        return L >= 2 ? L : 0
    elseif dev_sym === :mdev
        # _mdev_core: i in 1:N−3m+1, NaN when < 2.
        Ne = N - 3m + 1
        return Ne >= 2 ? Ne : 0
    elseif dev_sym === :hdev
        # _hdev_core: i in 1:N−3m, NaN when < 2.
        L = N - 3m
        return L >= 2 ? L : 0
    elseif dev_sym === :mhdev
        # _mhdev_core: i in 1:N−4m+1, NaN when < 2.
        Ne = N - 4m + 1
        return Ne >= 2 ? Ne : 0
    elseif dev_sym === :totdev
        # _totdev_core: n in 2:N−1 for every m (the reflected extension
        # supplies the tails); NaN only when N ≤ 2.
        return N > 2 ? N - 2 : 0
    elseif dev_sym === :mtotdev
        # _mtotdev_core: nsubs = N−3m+1 subsequences, NaN when < 1.
        nsubs = N - 3m + 1
        return nsubs >= 1 ? nsubs : 0
    elseif dev_sym === :htotdev
        # _htotdev_core: m = 1 falls back to the HDEV sum over N−3
        # terms (NaN when ≤ 0); m ≥ 2 iterates n_iter = (N−1)−3m+1
        # subsequences of y = diff(x)/τ₀ (NaN when < 1).
        if m == 1
            L = N - 3
            return L > 0 ? L : 0
        end
        n_iter = (N - 1) - 3m + 1
        return n_iter >= 1 ? n_iter : 0
    elseif dev_sym === :mhtotdev
        # _mhtotdev_core: nsubs = N−4m+1 subsequences, NaN when < 1
        # (or m < 1). The public wrapper's global drift removal does not
        # change the record length, so the count is unaffected.
        nsubs = N - 4m + 1
        return (m >= 1 && nsubs >= 1) ? nsubs : 0
    elseif dev_sym === :mtie || dev_sym === :tierms
        # _mtie_core: N−m fully-populated windows of m+1 samples, NaN when ≤ 0.
        # _tierms_core: the same N−m spans of m·τ₀ (one squared difference
        # each), NaN when ≤ 0 (or m < 1).
        (dev_sym === :tierms && m < 1) && return 0
        L = N - m
        return L > 0 ? L : 0
    elseif dev_sym === :theo1
        # _theo1_dev_at: i in 1:N−m outer windows; NaN for odd m or m
        # outside [2, N−1]. (allantools reports the total inner-term count
        # (N−m)·m/2 instead; this counts outer analysis windows, matching the
        # deverr normalization √(N−m) both use.)
        (m < 2 || isodd(m) || m > N - 1) && return 0
        return N - m
    elseif dev_sym === :theoh
        # theoh: ADEV count below/at the 0.2·T crossover, Thêo1 count at the
        # mapped even factor above it. Mirrors the branch rule in the public
        # `theoh` wrapper exactly.
        if _theoh_is_adev(m, N)
            L = N - 2m
            return L >= 2 ? L : 0
        end
        mt = _theoh_theo_m(m)
        return (mt >= 2 && mt <= N - 1) ? N - mt : 0
    end
    throw(ArgumentError("_neff_count: unknown deviation :$dev_sym"))
end

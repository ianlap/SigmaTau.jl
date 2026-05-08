# core/total.jl — Core Total Stability Kernels

"""
    _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:legacy) → Vector{Float64}

Computes the Total Deviation (TOTDEV) for a set of averaging factors `m`.

`detrend` selects the boundary-handling recipe:
- `:howe` — no detrend, mean-flip endpoint reflection (Howe 1995, NIST SP1065 eqn 25)
- `:greenhall` — per-window half-mean slope removal + time-reverse extension (Greenhall 2003)
- `:linear` — per-window full LS detrend + time-reverse extension
- `:legacy` — current SigmaTau behavior: global LS detrend + mean-flip reflection
"""
function _totdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:legacy)
    detrend === :legacy && return _totdev_legacy(x, m_values, tau0)
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
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
    for (k, m) in enumerate(m_values)
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

"""
    _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall) → Vector{Float64}

Computes the Modified Total Deviation (MTOTDEV).

`detrend` selects the boundary-handling recipe (see `_totdev_core` docstring).
For MTOTDEV, `:legacy` is an alias for `:greenhall` (current implementation).
"""
function _mtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _mtotdev_greenhall(x, m_values, tau0)
    end
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end

function _mtotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_seg = 3 * max_m
    ext = Vector{Float64}(undef, 3 * max_seg)
    cs = Vector{Float64}(undef, 3 * max_seg + 1)

    for (k, m) in enumerate(m_values)
        nsubs = N - 3m + 1
        if nsubs < 1
            devs[k] = NaN
            continue
        end

        seg_len = 3m
        outer_sum = 0.0

        for n in 1:nsubs
            half_n = seg_len / 2.0
            if m == 1
                slope = (x[n+2] - x[n]) / (2.0 * tau0)
            else
                hi = floor(Int, half_n)
                s1 = 0.0
                @inbounds @simd for i in 1:hi
                    s1 += x[n-1+i]
                end
                s1 /= hi

                s2 = 0.0
                @inbounds @simd for i in hi+1:seg_len
                    s2 += x[n-1+i]
                end
                s2 /= (seg_len - hi)
                slope = (s2 - s1) / (half_n * tau0)
            end

            @inbounds for j in 1:seg_len
                val = x[n-1+j] - slope * tau0 * (j - 1)
                rev_val = x[n-1 + seg_len - j + 1] - slope * tau0 * (seg_len - j)

                ext[j] = rev_val
                ext[seg_len + j] = val
                ext[2seg_len + j] = rev_val
            end

            cs[1] = 0.0
            @inbounds for j in 1:3seg_len
                cs[j+1] = cs[j] + ext[j]
            end

            block_sum = 0.0
            @inbounds @simd for j in 0:(6m - 1)
                a1 = (cs[j+m+1]  - cs[j+1])
                a2 = (cs[j+2m+1] - cs[j+m+1])
                a3 = (cs[j+3m+1] - cs[j+2m+1])
                d2 = (a3 - 2.0*a2 + a1) / m
                block_sum += d2^2
            end
            outer_sum += block_sum / (6.0 * m)
        end

        devs[k] = sqrt(outer_sum / (2.0 * Float64(m)^2 * tau0^2 * nsubs))
    end

    return devs
end

"""
    _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall) → Vector{Float64}

Computes the Hadamard Total Deviation (HTOTDEV).

`detrend` selects the boundary-handling recipe (see `_totdev_core`). For HTOTDEV,
`:legacy` is an alias for `:greenhall`. The recipe operates on the frequency
series `y = diff(x) / tau0`.
"""
function _htotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:greenhall)
    if detrend === :greenhall || detrend === :legacy
        return _htotdev_greenhall(x, m_values, tau0)
    end
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end

function _htotdev_greenhall(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    y = Vector{Float64}(undef, N-1)
    @inbounds @simd for i in 1:N-1
        y[i] = (x[i+1] - x[i]) / tau0
    end
    Ny = length(y)

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_seg = 3 * max_m
    ext = Vector{Float64}(undef, 3 * max_seg)
    cs = Vector{Float64}(undef, 3 * max_seg + 1)

    for (k, m) in enumerate(m_values)
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
        dev_sum = 0.0

        for i in 0:(n_iter - 1)
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

            cs[1] = 0.0
            @inbounds for j in 1:3seg_len
                cs[j+1] = cs[j] + ext[j]
            end

            sq = 0.0
            @inbounds @simd for j in 0:(6m - 1)
                h1 = (cs[j+m+1]  - cs[j+1])
                h2 = (cs[j+2m+1] - cs[j+m+1])
                h3 = (cs[j+3m+1] - cs[j+2m+1])
                sq += ((h3 - 2.0*h2 + h1) / m)^2
            end
            dev_sum += sq / (6.0 * m)
        end

        devs[k] = sqrt(dev_sum / (6.0 * n_iter))
    end

    return devs
end

"""
    _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:linear) → Vector{Float64}

Computes the Modified Hadamard Total Deviation (MHTOTDEV).

`detrend` selects the boundary-handling recipe (see `_totdev_core`). MHTOTDEV
is novel to SigmaTau; `:legacy` aliases to `:linear` (current implementation).
The default is `:linear` in this phase; switches to `:greenhall` in Phase 4.
"""
function _mhtotdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64; detrend::Symbol=:linear)
    if detrend === :linear || detrend === :legacy
        return _mhtotdev_linear(x, m_values, tau0)
    end
    throw(ArgumentError("unknown detrend recipe: $detrend; valid: :howe, :greenhall, :linear, :legacy"))
end

function _mhtotdev_linear(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)
    N = length(x)
    devs = Vector{Float64}(undef, length(m_values))

    max_m = isempty(m_values) ? 0 : maximum(m_values)
    max_Lp = 3 * max_m + 1
    ext_len = 3 * max_Lp
    L3_max = ext_len - 3 * max_m
    ext = Vector{Float64}(undef, ext_len)
    d3_vec = Vector{Float64}(undef, L3_max)
    S = Vector{Float64}(undef, L3_max + 1)

    for (k, m) in enumerate(m_values)
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

        total_sum = 0.0
        for n in 1:nsubs
            Lp_float = Float64(Lp)
            sum_i = (Lp_float * (Lp_float + 1.0)) / 2.0
            sum_i2 = (Lp_float * (Lp_float + 1.0) * (2.0*Lp_float + 1.0)) / 6.0
            delta = Lp_float * sum_i2 - sum_i^2

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

            S[1] = 0.0
            @inbounds for j in 1:L3
                S[j+1] = S[j] + d3_vec[j]
            end

            n_avg = L3 + 1 - m
            if n_avg > 0
                block_var = 0.0
                @inbounds @simd for j in 1:n_avg
                    block_var += (S[j+m] - S[j])^2
                end
                block_var /= (n_avg * 6.0 * Float64(m)^2)
            else
                block_var = 0.0
            end

            total_sum += block_var
        end

        devs[k] = sqrt(total_sum / (nsubs * Float64(m)^2 * tau0^2))
    end

    return devs
end

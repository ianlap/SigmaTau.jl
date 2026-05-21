# test/legacy_kernels.jl — Extracted reference kernels from the legacy
# SigmaTau Julia codebase (gitignored under `legacy/julia/src/deviations/`).
# These are inlined verbatim and depend only on Base — no heavy externals —
# so the parity testset can run on CI where the legacy tree is absent.
#
# Each function returns the *variance* (legacy convention); take sqrt to
# compare against the new `_*_core` deviation outputs.
#
# Source: legacy/julia/src/{types.jl,deviations/{allan,hadamard,total}.jl}.
# If the legacy code drifts, regenerate this file; otherwise these are the
# canonical numerical reference for SP1065-equivalent point estimates.

module LegacyKernels

# ── detrend_linear! (legacy/types.jl) ─────────────────────────────────────────

function detrend_linear!(x::AbstractVector{T}) where T<:Real
    n = length(x)
    n < 2 && return x
    x_bar = (n + 1) / 2.0
    sum_num = 0.0
    sum_y   = 0.0
    for i in 1:n
        val = Float64(x[i])
        sum_num += (i - x_bar) * val
        sum_y   += val
    end
    ss_xx = n * (Float64(n)^2 - 1) / 12.0
    slope = sum_num / ss_xx
    y_bar = sum_y / n
    for i in 1:n
        x[i] -= T(y_bar + slope * (i - x_bar))
    end
    return x
end

# ── ADEV / MDEV (allan.jl) ────────────────────────────────────────────────────

function adev_var(x::AbstractVector{<:Real}, m::Int, tau0::Real)
    N = length(x); L = N - 2m
    L <= 0 && return NaN
    d2 = @view(x[1+2m:end]) .- 2 .* @view(x[1+m:end-m]) .+ @view(x[1:L])
    return sum(abs2, d2) / (L * 2.0 * Float64(m)^2 * tau0^2)
end

function mdev_var(x::AbstractVector{<:Real}, m::Int, tau0::Real, x_cs::AbstractVector{<:Real})
    N = length(x); Ne = N - 3m + 1
    Ne <= 0 && return NaN
    d = @view(x_cs[1+3m:Ne+3m]) .- 3 .* @view(x_cs[1+2m:Ne+2m]) .+
        3 .* @view(x_cs[1+m:Ne+m]) .- @view(x_cs[1:Ne])
    return sum(abs2, d) / (Ne * 2.0 * Float64(m)^4 * tau0^2)
end

# ── HDEV / MHDEV (hadamard.jl) ────────────────────────────────────────────────

function hdev_var(x::AbstractVector{<:Real}, m::Int, tau0::Real)
    N = length(x); L = N - 3m
    L <= 0 && return NaN
    d3 = @view(x[1+3m:end]) .- 3 .* @view(x[1+2m:end-m]) .+
         3 .* @view(x[1+m:end-2m]) .- @view(x[1:L])
    return sum(abs2, d3) / (L * 6.0 * Float64(m)^2 * tau0^2)
end

function mhdev_var(x::AbstractVector{<:Real}, m::Int, tau0::Real, x_cs::AbstractVector{<:Real})
    N = length(x); Ne = N - 4m + 1
    Ne <= 0 && return NaN
    d = @view(x_cs[1+4m:Ne+4m]) .- 4 .* @view(x_cs[1+3m:Ne+3m]) .+
        6 .* @view(x_cs[1+2m:Ne+2m]) .- 4 .* @view(x_cs[1+m:Ne+m]) .+
        @view(x_cs[1:Ne])
    return sum(abs2, d) / (Ne * 6.0 * Float64(m)^4 * tau0^2)
end

# ── TOTDEV (total.jl) ─────────────────────────────────────────────────────────

function totdev_var(x::AbstractVector{<:Real}, m::Int, tau0::Real)
    N = length(x)
    xd = copy(x); detrend_linear!(xd)

    x_star = Vector{Float64}(undef, 3N - 4)
    for i in 1:N-2; x_star[i]         = 2xd[1] - xd[i+1]; end
    for i in 1:N;   x_star[N-2+i]     = xd[i];             end
    for i in 1:N-2; x_star[2N-2+i]    = 2xd[N] - xd[N-i];  end

    off = N - 2
    D = 0.0; count = 0
    for i in 1:N
        lo = off + i; hi = off + i + 2m
        hi > length(x_star) && continue
        d2 = x_star[hi] - 2x_star[off + i + m] + x_star[lo]
        D += d2^2; count += 1
    end
    count == 0 && return NaN
    return D / (2 * (N - 2) * (m * tau0)^2)
end

# ── MTOTDEV (total.jl) ────────────────────────────────────────────────────────

function mtotdev_var(x::AbstractVector{<:Real}, m::Int, tau0::Real)
    N = length(x); nsubs = N - 3m + 1
    nsubs < 1 && return NaN

    seg_len = 3m
    seq     = Vector{Float64}(undef, seg_len)
    seq_det = Vector{Float64}(undef, seg_len)
    ext     = Vector{Float64}(undef, 3seg_len)
    cs      = Vector{Float64}(undef, 3seg_len + 1)

    outer_sum = 0.0
    for n in 1:nsubs
        copyto!(seq, 1, x, n, seg_len)

        half_n = seg_len / 2
        if m == 1
            slope = (seq[3] - seq[1]) / (2tau0)
        else
            hi = floor(Int, half_n)
            s1 = sum(@view(seq[1:hi])) / hi
            s2 = sum(@view(seq[hi+1:seg_len])) / (seg_len - hi)
            slope = (s2 - s1) / (half_n * tau0)
        end
        for j in 1:seg_len
            seq_det[j] = seq[j] - slope * tau0 * (j - 1)
        end

        for j in 1:seg_len
            ext[j]            = seq_det[seg_len - j + 1]
            ext[seg_len + j]  = seq_det[j]
            ext[2seg_len + j] = seq_det[seg_len - j + 1]
        end

        cs[1] = 0.0
        for j in 1:3seg_len
            cs[j+1] = cs[j] + ext[j]
        end

        block_sum = 0.0
        for j in 0:(6m - 1)
            a1 = (cs[j+m+1]  - cs[j+1])   / m
            a2 = (cs[j+2m+1] - cs[j+m+1]) / m
            a3 = (cs[j+3m+1] - cs[j+2m+1]) / m
            d2 = a3 - 2a2 + a1
            block_sum += d2^2
        end
        outer_sum += block_sum / (6m)
    end

    return outer_sum / (2 * (m * tau0)^2 * nsubs)
end

# ── HTOTDEV (total.jl) ────────────────────────────────────────────────────────

function htotdev_var(x::AbstractVector{<:Real}, m::Int, tau0::Real)
    N = length(x)
    if m == 1
        L = N - 3; L <= 0 && return NaN
        d3 = @view(x[4:end]) .- 3 .* @view(x[3:end-1]) .+
             3 .* @view(x[2:end-2]) .- @view(x[1:L])
        return sum(abs2, d3) / (L * 6 * tau0^2)
    end

    y  = diff(x) ./ tau0
    Ny = length(y); n_iter = Ny - 3m + 1
    n_iter < 1 && return NaN

    seg_len = 3m
    xs    = Vector{Float64}(undef, seg_len)
    x0    = Vector{Float64}(undef, seg_len)
    xstar = Vector{Float64}(undef, 3seg_len)
    cs    = Vector{Float64}(undef, 3seg_len + 1)

    dev_sum = 0.0
    for i in 0:(n_iter - 1)
        copyto!(xs, 1, y, i + 1, seg_len)

        hi       = floor(Int, seg_len / 2)
        lo_start = ceil(Int, seg_len / 2) + 1
        m1 = sum(@view(xs[1:hi])) / hi
        m2 = sum(@view(xs[lo_start:seg_len])) / (seg_len - lo_start + 1)
        slope = if isodd(seg_len)
            (m2 - m1) / (0.5(seg_len - 1) + 1)
        else
            (m2 - m1) / (0.5seg_len)
        end
        mid = floor(seg_len / 2)
        for j in 1:seg_len
            x0[j] = xs[j] - slope * (j - 1 - mid)
        end

        for j in 1:seg_len
            xstar[j]            = x0[seg_len - j + 1]
            xstar[seg_len + j]  = x0[j]
            xstar[2seg_len + j] = x0[seg_len - j + 1]
        end

        cs[1] = 0.0
        for j in 1:3seg_len
            cs[j+1] = cs[j] + xstar[j]
        end

        sq = 0.0
        for j in 0:(6m - 1)
            h1 = (cs[j+m+1]  - cs[j+1])    / m
            h2 = (cs[j+2m+1] - cs[j+m+1])  / m
            h3 = (cs[j+3m+1] - cs[j+2m+1]) / m
            sq += (h3 - 2h2 + h1)^2
        end
        dev_sum += sq / (6m)
    end

    return dev_sum / (6 * n_iter)
end

# ── MHTOTDEV (total.jl) ───────────────────────────────────────────────────────

function mhtotdev_var(x::AbstractVector{<:Real}, m::Int, tau0::Real)
    m >= 1 || throw(ArgumentError("averaging factor m must be >= 1"))
    N = length(x); nsubs = N - 4m + 1
    nsubs < 1 && return NaN

    Lp = 3m + 1; ext_len = 3Lp; L3 = ext_len - 3m
    pd     = Vector{Float64}(undef, Lp)
    ext    = Vector{Float64}(undef, ext_len)
    d3_vec = Vector{Float64}(undef, L3)
    S      = Vector{Float64}(undef, L3 + 1)

    total_sum = 0.0
    for n in 1:nsubs
        copyto!(pd, 1, x, n, Lp)
        detrend_linear!(pd)

        for j in 1:Lp
            ext[j]        = pd[Lp - j + 1]
            ext[Lp + j]   = pd[j]
            ext[2Lp + j]  = pd[Lp - j + 1]
        end

        for j in 1:L3
            d3_vec[j] = ext[j] - 3ext[j+m] + 3ext[j+2m] - ext[j+3m]
        end

        S[1] = 0.0
        for j in 1:L3
            S[j+1] = S[j] + d3_vec[j]
        end

        n_avg = L3 + 1 - m
        block_var = 0.0
        if n_avg > 0
            for j in 1:n_avg
                a = S[j+m] - S[j]
                block_var += a^2
            end
            block_var /= (n_avg * 6.0 * Float64(m)^2)
        end

        total_sum += block_var
    end

    return total_sum / (nsubs * (m * tau0)^2)
end

end # module

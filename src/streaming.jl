# Stability estimators updated one sample at a time.
#
# Dobrogowski and Kasznia (2007) keep, for each averaging factor, a running
# sum of squared second differences for ADEV (their Eqs. (6) to (9)) and,
# for TDEV, an inner sum of m second differences that slides by one sample,
# with a running sum of its squares (Eqs. (10) to (14)). A new sample costs
# O(1) per averaging factor. The same scheme is applied here to third
# differences for HDEV and MHDEV, whose inner sum slides by a fourth
# difference.
#
# After any number of samples, `snapshot(acc).dev` equals `_adev_core`,
# `_mdev_core`, `_hdev_core` or `_mhdev_core` on the samples so far, with
# the same terms and normalization (test/stab/streaming.jl).

# Samples behind the newest that each update reads, in units of m: the
# second difference of ADEV reaches x[t−2m], the third difference of HDEV
# and of the MDEV slide x[t−3m], the fourth difference of the MHDEV slide
# x[t−4m].
const _STREAMING_LAG = Dict{Symbol,Int}(:adev => 2, :mdev => 3, :hdev => 3, :mhdev => 4)

"""
    StreamingStability(dev, tau0, m_values)

Accumulator that updates ADEV, MDEV, HDEV or MHDEV as phase samples
arrive, without keeping the record. Feed it with `push!` or `append!` and
read the estimate at any moment with [`snapshot`](@ref).

After any number of samples the deviation equals that of [`adev`](@ref),
[`mdev`](@ref), [`hdev`](@ref) or [`mhdev`](@ref) on the samples so far.

The running sums are those of
[Dobrogowski and Kasznia (2007), Eqs. (6) to (14)](@cite dobrogowski-2007-realtime-adev),
who give them for ADEV and TDEV. The same scheme is applied here to the
third differences of HDEV and MHDEV. A sample costs a fixed number of
operations per averaging factor, and the accumulator holds the last
``k\\max(m) + 1`` samples, with ``k = 2`` for ADEV, 3 for MDEV and HDEV and
4 for MHDEV. The total and Thêo estimators read the whole record at every
averaging time and cannot be accumulated this way.

# Arguments

- `dev`: `:adev`, `:mdev`, `:hdev` or `:mhdev`.
- `tau0`: sample interval in seconds.
- `m_values`: averaging factors ``m``, with ``τ = mτ_0``.

!!! warning "Missing samples"
    Every sample must be a number. The estimators leave out the terms that
    a missing sample enters; the accumulator cannot, and after a `NaN`
    every estimate is `NaN`.

# Examples

```jldoctest
julia> using SigmaTau

julia> x = [0.0, 1.0, 3.0, 2.0, 5.0, 4.0, 7.0, 9.0];

julia> acc = StreamingStability(:adev, 1.0, [1, 2]);

julia> append!(acc, x);

julia> snapshot(acc).dev ≈ adev(PhaseData(x, 1.0), [1, 2]; ci=false).dev
true
```
"""
mutable struct StreamingStability
    const dev::Symbol
    const tau0::Float64
    const m_values::Vector{Int}
    "Ring buffer of the most recent `cap` samples; sample t lives at mod1(t, cap)."
    const buf::Vector{Float64}
    const cap::Int
    "Per-m running sum: Σd² (`:adev`/`:hdev`) or Σ Sᵢ(m)² (`:mdev`/`:mhdev`)."
    const sumsq::Vector{Float64}
    "Per-m sliding inner sum Sᵢ(m) (modified family only; unused otherwise)."
    const inner::Vector{Float64}
    "Total samples pushed so far."
    count::Int

    function StreamingStability(dev::Symbol, tau0::Real,
                                m_values::AbstractVector{<:Integer})
        haskey(_STREAMING_LAG, dev) || throw(ArgumentError(
            "StreamingStability: dev must be one of :adev, :mdev, :hdev, :mhdev " *
            "(totals and Thêo re-read the whole record and cannot stream), got :$dev"))
        tau0 > 0 || throw(ArgumentError(
            "StreamingStability: tau0 must be positive, got $tau0"))
        isempty(m_values) && throw(ArgumentError(
            "StreamingStability: m_values must be non-empty"))
        all(>=(1), m_values) || throw(ArgumentError(
            "StreamingStability: every m must be ≥ 1, got $(minimum(m_values))"))
        ms = Vector{Int}(m_values)
        cap = _STREAMING_LAG[dev] * maximum(ms) + 1
        return new(dev, Float64(tau0), ms, zeros(Float64, cap), cap,
                   zeros(Float64, length(ms)), zeros(Float64, length(ms)), 0)
    end
end

# Sample with absolute stream index `idx` (1-based since the first push).
# Valid only while idx is within the last `cap` samples.
@inline _ring_get(buf::Vector{Float64}, cap::Int, idx::Int) =
    @inbounds buf[mod1(idx, cap)]

"""
    push!(acc::StreamingStability, x) → acc

Add one phase sample `x`, in seconds, and update the running sum of every
averaging factor: the sum of squared second differences for `:adev`
([Dobrogowski and Kasznia (2007), Eq. (9)](@cite dobrogowski-2007-realtime-adev)),
the inner sum and the sum of its squares for `:mdev` (their Eqs. (11) to
(14)), and the same with third differences for `:hdev` and `:mhdev`.
"""
function Base.push!(acc::StreamingStability, x_next::Real)
    t = acc.count + 1
    buf = acc.buf
    cap = acc.cap
    @inbounds buf[mod1(t, cap)] = Float64(x_next)
    acc.count = t

    ms = acc.m_values
    sumsq = acc.sumsq
    dev = acc.dev
    if dev === :adev
        # Eq. (9): S_t(m) = S_{t−1}(m) + (x_t − 2x_{t−m} + x_{t−2m})².
        for (k, m) in enumerate(ms)
            t >= 2m + 1 || continue
            d = _ring_get(buf, cap, t) - 2.0 * _ring_get(buf, cap, t - m) +
                _ring_get(buf, cap, t - 2m)
            @inbounds sumsq[k] += d * d
        end
    elseif dev === :hdev
        # The same with the third difference.
        for (k, m) in enumerate(ms)
            t >= 3m + 1 || continue
            d = _ring_get(buf, cap, t) - 3.0 * _ring_get(buf, cap, t - m) +
                3.0 * _ring_get(buf, cap, t - 2m) - _ring_get(buf, cap, t - 3m)
            @inbounds sumsq[k] += d * d
        end
    elseif dev === :mdev
        inner = acc.inner
        for (k, m) in enumerate(ms)
            if t > 3m
                # Slide the inner sum (Eq. (14)): the oldest second difference
                # leaves and the newest enters, a single third difference.
                @inbounds inner[k] += _ring_get(buf, cap, t) -
                    3.0 * _ring_get(buf, cap, t - m) +
                    3.0 * _ring_get(buf, cap, t - 2m) -
                    _ring_get(buf, cap, t - 3m)
                @inbounds sumsq[k] += inner[k]^2     # Eq. (11)
            elseif t >= 2m + 1
                # Build the first inner sum (Eq. (13)), one second difference
                # at a time. It is complete at t = 3m.
                @inbounds inner[k] += _ring_get(buf, cap, t) -
                    2.0 * _ring_get(buf, cap, t - m) +
                    _ring_get(buf, cap, t - 2m)
                if t == 3m
                    @inbounds sumsq[k] += inner[k]^2 # Eq. (11), first term
                end
            end
        end
    else # :mhdev
        inner = acc.inner
        for (k, m) in enumerate(ms)
            if t > 4m
                # As for MDEV, with third differences: the inner sum slides by
                # the fourth difference x_t − 4x_{t−m} + 6x_{t−2m} − 4x_{t−3m}
                # + x_{t−4m}.
                @inbounds inner[k] += _ring_get(buf, cap, t) -
                    4.0 * _ring_get(buf, cap, t - m) +
                    6.0 * _ring_get(buf, cap, t - 2m) -
                    4.0 * _ring_get(buf, cap, t - 3m) +
                    _ring_get(buf, cap, t - 4m)
                @inbounds sumsq[k] += inner[k]^2
            elseif t >= 3m + 1
                # The first inner sum, of m third differences, is complete at
                # t = 4m.
                @inbounds inner[k] += _ring_get(buf, cap, t) -
                    3.0 * _ring_get(buf, cap, t - m) +
                    3.0 * _ring_get(buf, cap, t - 2m) -
                    _ring_get(buf, cap, t - 3m)
                if t == 4m
                    @inbounds sumsq[k] += inner[k]^2
                end
            end
        end
    end
    return acc
end

"""
    append!(acc::StreamingStability, xs) → acc

Add the phase samples `xs`, in seconds and in order. The result is that of
`push!` on each in turn.
"""
function Base.append!(acc::StreamingStability, xs::AbstractVector{<:Real})
    for v in xs
        push!(acc, v)
    end
    return acc
end

"""
    snapshot(acc::StreamingStability) → StabilityResult
    snapshot(live::LiveStability) → Vector{StabilityResult}

Deviation over all the samples received so far, as a
[`StabilityResult`](@ref); for a [`LiveStability`](@ref), one result per
deviation.

`dev` is `NaN` at an averaging factor with fewer than two terms, and `neff`
holds the number of terms. The result carries no noise type, EDF or
confidence interval: for those, run the estimator on the record.
"""
function snapshot(acc::StreamingStability)
    ms = acc.m_values
    t = acc.count
    tau0 = acc.tau0
    dev = acc.dev
    devs = Vector{Float64}(undef, length(ms))
    neff = Vector{Int}(undef, length(ms))
    for (k, m) in enumerate(ms)
        if dev === :adev
            L = t - 2m
            ok = L >= 2
            devs[k] = ok ? sqrt(acc.sumsq[k] / (2.0 * L * Float64(m)^2 * tau0^2)) : NaN
            neff[k] = ok ? L : 0
        elseif dev === :hdev
            L = t - 3m
            ok = L >= 2
            devs[k] = ok ? sqrt(acc.sumsq[k] / (6.0 * L * Float64(m)^2 * tau0^2)) : NaN
            neff[k] = ok ? L : 0
        elseif dev === :mdev
            Ne = t - 3m + 1
            ok = Ne >= 2
            devs[k] = ok ? sqrt(acc.sumsq[k] / (2.0 * Ne * Float64(m)^4 * tau0^2)) : NaN
            neff[k] = ok ? Ne : 0
        else # :mhdev
            Ne = t - 4m + 1
            ok = Ne >= 2
            devs[k] = ok ? sqrt(acc.sumsq[k] / (6.0 * Ne * Float64(m)^4 * tau0^2)) : NaN
            neff[k] = ok ? Ne : 0
        end
    end
    return StabilityResult(dev, ms .* tau0, devs, Symbol[], _CIBound[], Float64[], neff,
                           tau0, t, NaN)
end

"""
    nsamples(acc::StreamingStability) → Int

Number of phase samples the accumulator has received.
"""
nsamples(acc::StreamingStability) = acc.count

function Base.show(io::IO, acc::StreamingStability)
    print(io, "StreamingStability(:", acc.dev, ", ", length(acc.m_values),
          " m-values, τ₀=", _show4g(acc.tau0), " s, ", acc.count, " samples)")
end

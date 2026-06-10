# streaming.jl — real-time streaming stability accumulators.
#
# Implements the Dobrogowski–Kasznia real-time computation scheme (2007 IEEE
# Frequency Control Symposium): per averaging factor m, keep a running sum of
# squared overlapping differences (their eqs. 6–9 for ADEV), and for the
# modified family the cumulative inner-sum form (their eqs. 10–14 for TDEV),
# each updated in O(1) work per new sample per m. Generalized here to the
# third-difference (Hadamard) family: `:hdev` streams squared third
# differences, and `:mhdev` keeps a per-m sliding inner sum of third
# differences whose slide update is the binomial *fourth* difference
# (the third-difference analogue of their eq. 12).
#
# Correctness contract: at any sample count, `snapshot(acc).dev` reproduces
# the batch kernels `_adev_core` / `_mdev_core` / `_hdev_core` /
# `_mhdev_core` on the samples pushed so far — same windows, same
# normalization, same NaN guards. The streaming equivalence testset in
# test/stab/streaming.jl locks this in at rtol 1e-10.

# Deepest lag (in units of m) each kernel reads back from the newest sample:
# ADEV second difference reaches x_{t−2m}; MDEV's slide update (eq. 12) and
# HDEV's third difference reach x_{t−3m}; MHDEV's fourth-difference slide
# reaches x_{t−4m}.
const _STREAMING_LAG = Dict{Symbol,Int}(:adev => 2, :mdev => 3, :hdev => 3, :mhdev => 4)

"""
    StreamingStability(dev::Symbol, tau0, m_values)

Real-time streaming accumulator for the overlapping difference-family
deviations: feed phase samples one at a time with `push!` (or in chunks
with `append!`) and read off the running deviation estimate at any moment
with [`snapshot`](@ref) — without storing the record or recomputing
anything.

`dev` selects the estimator: `:adev`, `:mdev`, `:hdev`, or `:mhdev`. `tau0`
is the sample interval in seconds and `m_values` the averaging factors
(τ = m·τ₀), exactly as in the batch API. At every sample count the streamed
estimate equals the corresponding batch call on the samples seen so far
(`adev(PhaseData(x[1:k], tau0), m_values; ci=false)`, etc.).

Implements the real-time computation scheme of Dobrogowski & Kasznia
[dobrogowski-2007-realtime-adev](@cite). For `:adev`/`:hdev` each new
sample updates a per-m running sum of squared second/third differences
(their eqs. 6–9); for the modified family `:mdev`/`:mhdev` a per-m sliding
inner sum ``S_i(m)`` is maintained alongside the cumulative sum of its
squares (their eqs. 10–14, generalized to third differences for `:mhdev`).
Every update is O(1) per new sample per m — the per-sample cost is
O(`length(m_values)`) regardless of how long the stream runs.

Memory is bounded: only a ring buffer of the most recent samples is held,
sized to the selected kernel's deepest lag — `2·maximum(m_values) + 1`
samples for `:adev`, `3·maximum(m_values) + 1` for `:mdev`/`:hdev`, and
`4·maximum(m_values) + 1` for `:mhdev` (the deepest lag any kernel needs).

The total-family estimators (`totdev`, `mtotdev`, `ttotdev`, `htotdev`,
`mhtotdev`) and the Thêo family cannot stream: their extension/sampling
schemes re-read the *whole* record on every update, so they have no O(1)
per-sample form.

# Example

```julia
acc = StreamingStability(:adev, 1.0, [1, 10, 100])
for x in phase_source            # samples arrive one at a time
    push!(acc, x)
end
result = snapshot(acc)            # StabilityResult, same dev values as batch
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
            "StreamingStability: every m must be ≥ 1, got $(m_values)"))
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
    push!(acc::StreamingStability, x_next::Real) → acc

Feed one new phase sample (seconds) into the accumulator. Updates every
per-m running sum in O(1) per averaging factor, per the Dobrogowski–Kasznia
scheme [dobrogowski-2007-realtime-adev](@cite): the squared-difference sum
update of their eq. 9 for `:adev` (third-difference form for `:hdev`), the
inner-sum build/slide of their eqs. 12–13 plus the overall-sum update of
eq. 11 for `:mdev` (third-difference generalization for `:mhdev`).
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
        # Eq. 9: S_t(m) = S_{t−1}(m) + (x_t − 2x_{t−m} + x_{t−2m})².
        for (k, m) in enumerate(ms)
            t >= 2m + 1 || continue
            d = _ring_get(buf, cap, t) - 2.0 * _ring_get(buf, cap, t - m) +
                _ring_get(buf, cap, t - 2m)
            @inbounds sumsq[k] += d * d
        end
    elseif dev === :hdev
        # Eq. 9 generalized to the third difference.
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
                # Eq. 12 slide: drop the oldest second difference of the
                # window, pick up the newest — a single third difference.
                @inbounds inner[k] += _ring_get(buf, cap, t) -
                    3.0 * _ring_get(buf, cap, t - m) +
                    3.0 * _ring_get(buf, cap, t - 2m) -
                    _ring_get(buf, cap, t - 3m)
                @inbounds sumsq[k] += inner[k]^2     # eq. 11
            elseif t >= 2m + 1
                # Eq. 13 build phase: accumulate the first window's m second
                # differences one at a time; the window completes at t = 3m.
                @inbounds inner[k] += _ring_get(buf, cap, t) -
                    2.0 * _ring_get(buf, cap, t - m) +
                    _ring_get(buf, cap, t - 2m)
                if t == 3m
                    @inbounds sumsq[k] += inner[k]^2 # eq. 11, first window
                end
            end
        end
    else # :mhdev
        inner = acc.inner
        for (k, m) in enumerate(ms)
            if t > 4m
                # Third-difference generalization of eq. 12: the slide is the
                # binomial fourth difference x_t − 4x_{t−m} + 6x_{t−2m}
                # − 4x_{t−3m} + x_{t−4m}.
                @inbounds inner[k] += _ring_get(buf, cap, t) -
                    4.0 * _ring_get(buf, cap, t - m) +
                    6.0 * _ring_get(buf, cap, t - 2m) -
                    4.0 * _ring_get(buf, cap, t - 3m) +
                    _ring_get(buf, cap, t - 4m)
                @inbounds sumsq[k] += inner[k]^2
            elseif t >= 3m + 1
                # Eq. 13 build phase with third differences; the first window
                # completes at t = 4m.
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
    append!(acc::StreamingStability, xs::AbstractVector{<:Real}) → acc

Feed a chunk of phase samples (seconds), in order — exactly equivalent to
`push!`-ing them one at a time.
"""
function Base.append!(acc::StreamingStability, xs::AbstractVector{<:Real})
    for v in xs
        push!(acc, v)
    end
    return acc
end

"""
    snapshot(acc::StreamingStability) → StabilityResult

The deviation estimate over all samples streamed so far, as a plain
[`StabilityResult`](@ref). `dev` matches the batch kernel on the same record
at every m (NaN where the batch kernel has fewer than 2 analysis windows);
`neff` carries the current window counts (0 at NaN points). The
`ci`/`edf`/`noise_type` fields are empty — on-demand EDF/CI for streaming
accumulators is future work; for confidence intervals run the batch API on
the full record.
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
    return StabilityResult(dev, ms .* tau0, devs, Symbol[], _CIBound[], Float64[], neff)
end

"""
    nsamples(acc::StreamingStability) → Int

Number of phase samples streamed into the accumulator so far.
"""
nsamples(acc::StreamingStability) = acc.count

function Base.show(io::IO, acc::StreamingStability)
    print(io, "StreamingStability(:", acc.dev, ", ", length(acc.m_values),
          " m-values, τ₀=", _show4g(acc.tau0), " s, ", acc.count, " samples)")
end

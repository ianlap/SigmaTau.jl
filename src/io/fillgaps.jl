# io/fillgaps.jl — Howe & Schlossberger gap imputation that preserves the
# noise character (and AVAR shape) of the surrounding data.
#
# Reference: D. A. Howe and L. Schlossberger,
# "Characterizing Frequency Stability Measurements having Multiple Data Gaps",
# Proc. PTTI 2009.
#
# Algorithm: for each gap, take a number of samples equal to the gap length
# from one (or both) side(s), reflect-and-invert them, low-pass filter via FFT
# with exp(-k/N) shaping, then add a linear ramp matched to the endpoints of
# the surrounding (filtered) data. Iterates outward from the largest
# continuous run.

"""
    _howe_filter(x::Vector{Float64}) → Vector{Float64}

Howe's symmetric reflect-and-FFT-shape filter. `x` is padded as
`[reverse(x); x; reverse(x)]`, multiplied in the frequency domain by the
even-symmetric mask `exp(-|k|/N)`, then the central third is returned.
Used both for synthesising fill points and as a smoothing pass on existing
data so the linear ramp matches a low-passed endpoint instead of a noisy
single sample.
"""
function _howe_filter(x::AbstractVector{<:Real})
    n  = length(x)
    n == 0 && return Float64[]
    pad      = Vector{Float64}(undef, 3n)
    @inbounds for i in 1:n
        pad[i]       = x[n - i + 1]   # flip(x)
        pad[n + i]   = x[i]
        pad[2n + i]  = x[n - i + 1]   # flip(x)
    end
    N        = length(pad)
    F        = fft(pad)
    half     = fld(N, 2)
    filt     = Vector{Float64}(undef, N)
    @inbounds for k in 0:half-1
        filt[k + 1] = exp(-k / N)
    end
    if iseven(N)
        @inbounds for k in 0:half-1
            filt[N - k] = filt[k + 1]
        end
    else
        # odd: the centre bin gets filt[half] (k=half-1+1) repeated, then mirror.
        @inbounds filt[half + 1] = filt[half]
        @inbounds for k in 0:half-1
            filt[N - k] = filt[k + 1]
        end
    end
    Ffilt    = F .* filt
    ypad     = real.(ifft(Ffilt))
    lo       = round(Int, N / 3) + 1
    hi       = round(Int, 2N / 3)
    return ypad[lo:hi]
end

"""
    _filter_all(x::Vector{Float64}) → Vector{Float64}

Apply `_howe_filter` to every contiguous run of finite samples in `x`,
leaving the `NaN` slots alone. Mirrors the MATLAB `filter_all` helper.
"""
function _filter_all(x::AbstractVector{<:Real})
    y = fill(NaN, length(x))
    i = 1
    n = length(x)
    while i <= n
        if isnan(x[i])
            i += 1
            continue
        end
        j = i
        while j <= n && !isnan(x[j])
            j += 1
        end
        y[i:j-1] .= _howe_filter(@view x[i:j-1])
        i = j
    end
    return y
end

"""
    _identify_gaps(x::AbstractVector{<:Real}) → Vector{Tuple{Int,Int}}

Return a list of `(start, stop)` 1-based inclusive index ranges marking the
contiguous `NaN` runs in `x`. Empty if `x` has no gaps.
"""
function _identify_gaps(x::AbstractVector{<:Real})
    gaps = Tuple{Int,Int}[]
    n = length(x)
    i = 1
    while i <= n
        if isnan(x[i])
            j = i
            while j <= n && isnan(x[j])
                j += 1
            end
            push!(gaps, (i, j - 1))
            i = j
        else
            i += 1
        end
    end
    return gaps
end

"""
    _fill_singletons!(x::Vector{Float64}, gaps::Vector{Tuple{Int,Int}})

Fill any 1-sample gaps in `x` (mean of the two neighbours, or nearest
neighbour at the edges) and remove those entries from `gaps`. Returns the
mutated `gaps` vector.
"""
function _fill_singletons!(x::Vector{Float64}, gaps::Vector{Tuple{Int,Int}})
    keep_mask = trues(length(gaps))
    n = length(x)
    @inbounds for (k, g) in pairs(gaps)
        g[1] == g[2] || continue
        idx = g[1]
        if idx == 1
            x[idx] = x[2]
        elseif idx == n
            x[idx] = x[n - 1]
        else
            x[idx] = 0.5 * (x[idx - 1] + x[idx + 1])
        end
        keep_mask[k] = false
    end
    deleteat!(gaps, .!keep_mask)
    return gaps
end

"""
    _choose_initial_gap(gaps) → Int

1-based index into `gaps` of the gap immediately following the largest
continuous data run between consecutive gaps. Mirrors the MATLAB helper of
the same name. Returns `1` when there is only one gap.
"""
function _choose_initial_gap(gaps::Vector{Tuple{Int,Int}})
    length(gaps) <= 1 && return 1
    maxdiff = 0
    pick    = 1
    for k in 1:length(gaps) - 1
        d = gaps[k + 1][1] - gaps[k][2]
        if d > maxdiff
            maxdiff = d
            pick    = k + 1
        end
    end
    return pick
end

# --- per-gap point selection ----------------------------------------------
# Returns one of:
#   :left, view-into-x           — sizegap points to the left of the gap
#   :right, view-into-x          — sizegap points to the right
#   :both, (left_view, right_view) — half-sized chunks from each side
#   :none, nothing                — couldn't fill (caller advances/retreats)
#   :done, nothing                — globally unfillable, caller bails out
function _check_right(x, gaps, curgap, sizegap, gap_num, gap_total)
    rightlen = length(x) - curgap[2]
    if rightlen >= sizegap
        if gap_num == gap_total
            return :right, view(x, curgap[2] + 1 : curgap[2] + sizegap)
        else
            g_next = gaps[gap_num + 1]
            if g_next[1] - curgap[2] - 1 >= sizegap
                return :right, view(x, curgap[2] + 1 : curgap[2] + sizegap)
            end
        end
    end
    return _check_both_sides(x, gaps, curgap, sizegap, gap_num, gap_total)
end

function _check_left(x, gaps, curgap, sizegap, gap_num, gap_total)
    if gap_num > 1
        g_prev = gaps[gap_num - 1]
        if curgap[1] - g_prev[2] - 1 >= sizegap
            return :left, view(x, curgap[1] - sizegap : curgap[1] - 1)
        end
    elseif curgap[1] > sizegap
        return :left, view(x, curgap[1] - sizegap : curgap[1] - 1)
    end
    return _check_right(x, gaps, curgap, sizegap, gap_num, gap_total)
end

function _check_both_sides(x, gaps, curgap, sizegap, gap_num, gap_total)
    half = fld(sizegap, 2)
    half == 0 && return :none, nothing
    if gap_num == 1
        if gap_total == 1
            if curgap[1] > half && curgap[2] <= length(x) - half
                return :both,
                       (view(x, curgap[1] - half : curgap[1] - 1),
                        view(x, curgap[2] + 1   : curgap[2] + half))
            end
            return :done, nothing
        else
            g_next = gaps[gap_num + 1]
            if curgap[1] > half && g_next[1] - curgap[2] - 1 >= half
                return :both,
                       (view(x, curgap[1] - half : curgap[1] - 1),
                        view(x, curgap[2] + 1   : curgap[2] + half))
            end
            return :none, nothing
        end
    elseif gap_num == gap_total
        g_prev = gaps[gap_num - 1]
        if curgap[1] - g_prev[2] - 1 >= half && curgap[2] <= length(x) - half
            return :both,
                   (view(x, curgap[1] - half : curgap[1] - 1),
                    view(x, curgap[2] + 1   : curgap[2] + half))
        end
        return :none, nothing
    else
        g_prev = gaps[gap_num - 1]
        g_next = gaps[gap_num + 1]
        if curgap[1] - g_prev[2] - 1 >= half && g_next[1] - curgap[2] - 1 >= half
            return :both,
                   (view(x, curgap[1] - half : curgap[1] - 1),
                    view(x, curgap[2] + 1   : curgap[2] + half))
        end
        return :none, nothing
    end
end

# --- per-gap fill driver ---------------------------------------------------
# Returns (status, new_gap_num) where status is :ok, :advance, :done.
function _fill_one_gap!(x::Vector{Float64},
                        gaps::Vector{Tuple{Int,Int}},
                        gap_num::Int,
                        reverse::Bool)
    gap_total = length(gaps)
    curgap    = gaps[gap_num]
    sizegap   = curgap[2] - curgap[1] + 1

    side, pts = if reverse
        if curgap[2] > length(x) - sizegap
            _check_left(x, gaps, curgap, sizegap, gap_num, gap_total)
        else
            _check_right(x, gaps, curgap, sizegap, gap_num, gap_total)
        end
    else
        if curgap[1] <= sizegap
            _check_right(x, gaps, curgap, sizegap, gap_num, gap_total)
        else
            _check_left(x, gaps, curgap, sizegap, gap_num, gap_total)
        end
    end

    if side === :none
        return reverse ? (:advance, gap_num - 1) : (:advance, gap_num + 1)
    elseif side === :done
        return :done, gap_num
    end

    fillvec = if side === :both
        left_pts, right_pts = pts
        # reflect-and-invert each side
        left  = .-reverse_view(left_pts)
        right = .-reverse_view(right_pts)
        filt_left  = _howe_filter(left)
        filt_right = _howe_filter(right)
        shift      = filt_left[end] - filt_right[1]
        right    .+= shift
        if isodd(sizegap)
            midval = 0.5 * (left[end] + right[1])
            vcat(left, [midval], right)
        else
            vcat(left, right)
        end
    else
        # single-side: invert and reflect
        .-reverse_view(pts)
    end

    fillvec = _howe_filter(fillvec)

    # Endpoint ramp: anchor to the filtered surrounding data so the joins
    # don't kink off a single noisy boundary sample.
    filt_data = _filter_all(x)
    start_anchor = filt_data[curgap[1] - 1]
    end_anchor   = filt_data[curgap[2] + 1]
    ramp = range(start_anchor - fillvec[1], stop=end_anchor - fillvec[end], length=sizegap)

    @inbounds for (k, idx) in enumerate(curgap[1]:curgap[2])
        x[idx] = fillvec[k] + ramp[k]
    end

    deleteat!(gaps, gap_num)
    # Rightward sweep: the gap at gap_num+1 has just shifted into gap_num,
    # so keep the index. Leftward sweep: actively step left.
    new_gap_num = if reverse
        gap_num - 1
    else
        gap_num > length(gaps) ? length(gaps) : gap_num
    end
    return :ok, new_gap_num
end

reverse_view(v) = collect(Iterators.reverse(v))

"""
    _howe_fillgaps_core(x::Vector{Float64}) → (xfilled, filled_mask::BitVector)

In-place-equivalent Howe fill. Takes an equispaced phase/frequency vector
with `NaN` at missing samples; returns a finite vector and a bit-mask
flagging the indices that were imputed. Returns `(copy(x), falses(...))`
unchanged if `x` has no `NaN`.

Throws `ErrorException` if the gap topology cannot be filled (e.g. NaN
runs at both record ends with no sufficiently long internal run to mirror).
"""
function _howe_fillgaps_core(x::AbstractVector{<:Real})
    xv         = Vector{Float64}(x)
    n          = length(xv)
    nan_mask_0 = isnan.(xv)
    if !any(nan_mask_0)
        return xv, falses(n)
    end

    gaps = _identify_gaps(xv)
    _fill_singletons!(xv, gaps)

    gap_num = _choose_initial_gap(gaps)

    # Outer loop: rightward sweep, then leftward sweep, repeat until empty.
    while !isempty(gaps)
        # Rightward sweep
        while gap_num <= length(gaps) && gap_num >= 1
            status, gap_num = _fill_one_gap!(xv, gaps, gap_num, false)
            status === :done && error("fillgaps: gap topology cannot be filled (data record too short relative to gap span)")
            isempty(gaps) && break
        end
        isempty(gaps) && break

        # Leftward sweep
        gap_num = min(gap_num, length(gaps))
        while gap_num >= 1 && !isempty(gaps)
            status, gap_num = _fill_one_gap!(xv, gaps, gap_num, true)
            status === :done && error("fillgaps: gap topology cannot be filled (data record too short relative to gap span)")
            gap_num < 1 && break
        end
        isempty(gaps) && break

        # Reset to the gap right of the largest remaining run, otherwise the
        # outer loop would spin forever.
        gap_num = _choose_initial_gap(gaps)
    end

    filled_mask = nan_mask_0 .& .!isnan.(xv)
    return xv, filled_mask
end

"""
    _make_equispaced(t::AbstractVector, x::AbstractVector) → (tfilled, xfilled)

Snap an irregular `(t, x)` series to a uniform grid with spacing
`dt = minimum(diff(t))`, populating gaps with `NaN`. The grid spans
`t[1]:dt:t[end]`. Used by `read_phase` / `read_frequency` when `fillgaps=true`.
"""
function _make_equispaced(t::AbstractVector{<:Real}, x::AbstractVector{<:Real})
    length(t) == length(x) || throw(ArgumentError("_make_equispaced: t and x must have equal length"))
    length(t) >= 2 || throw(ArgumentError("_make_equispaced: need at least 2 samples"))
    dt = minimum(diff(t))
    dt > 0 || throw(ArgumentError("_make_equispaced: time vector is not strictly increasing"))
    tfilled = collect(t[1]:dt:t[end])
    xfilled = fill(NaN, length(tfilled))
    @inbounds for i in eachindex(t)
        idx = round(Int, (t[i] - t[1]) / dt) + 1
        xfilled[idx] = x[i]
    end
    return tfilled, xfilled
end

"""
    fillgaps(pd::PhaseData) → PhaseData
    fillgaps(fd::FrequencyData) → FrequencyData

Impute `NaN` samples using Howe & Schlossberger's reflect-and-FFT-filter
algorithm, which preserves the local noise character and the AVAR/MDEV
curves of the surrounding data. Returns a new `PhaseData` /
`FrequencyData` with the same `tau0`. The input is treated as already
equispaced — for raw `(t, x)` tables with irregular sampling, use
[`read_phase`](@ref) with `fillgaps=true`.

Single-sample gaps are filled with the mean of their neighbours; multi-sample
gaps are filled by reflecting `sizegap` points from the larger neighbouring
run (or `sizegap/2` from each side when neither side has enough), inverting,
low-pass filtering, and adding an endpoint-matched linear ramp.

Returns the input unchanged if there are no `NaN` samples.
"""
function fillgaps(pd::PhaseData)
    xfilled, _ = _howe_fillgaps_core(pd.x)
    return PhaseData(xfilled, pd.tau0)
end

function fillgaps(fd::FrequencyData)
    yfilled, _ = _howe_fillgaps_core(fd.y)
    return FrequencyData(yfilled, fd.tau0)
end

# Gap filling by reflection, after Howe and Schlossberger (2022): the samples
# beside a gap are reversed, inverted and placed in the gap, with a line added
# so that both ends join the record.
#
# For each gap, as many samples as are missing are taken from the run before
# it, or after it, or half from each side, according to which is long enough.
# The reflected and inverted segment is filtered (`_howe_filter`), and the
# line joins it to the filtered data on both sides. Gaps are filled outward
# from the longest run of data.
#
# Two steps differ from the paper and from the authors' script (Champagne,
# Howe and Schlossberger (2021)): there the reflected samples are inserted
# unfiltered, the filter serving only to find the ends the line joins, and
# the gain of the filter is exp(−16k/L) where this one has exp(−k/L).

"""
    _howe_filter(x) → Vector{Float64}

Low-pass filter of a segment. `x` is padded with its reflection on both
sides, `[reverse(x); x; reverse(x)]`, multiplied in the frequency domain by
the even mask `exp(-|k|/L)`, with `L` the padded length, and the middle
third is returned. The gain falls to `exp(-1/2)` at the Nyquist frequency.
Applied to the segment placed in a gap, and to the data beside it, whose
filtered ends the line joins.
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
        # odd: the center bin gets filt[half] (k=half-1+1) repeated, then mirror.
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

# The filtered value at `i`, the end of a run of samples, from up to `w`
# samples of that run reaching away from the gap (`dir` -1 leftward, +1 rightward).
function _edge_filtered(x::Vector{Float64}, i::Int, dir::Int, w::Int)
    j = i
    while abs(j - i) + 1 < w && 1 <= j + dir <= length(x) && !isnan(x[j + dir])
        j += dir
    end
    f = _howe_filter(view(x, min(i, j):max(i, j)))
    return dir < 0 ? f[end] : f[1]
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

Fill the gaps of one sample in `x` with the mean of the two neighbors, or
with the one neighbor at an end of the record, and remove them from `gaps`.
Returns `gaps`.
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
continuous data run between consecutive gaps. Returns `1` when there is
only one gap.
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
#   :left, view-into-x:            sizegap points to the left of the gap
#   :right, view-into-x:           sizegap points to the right
#   :both, (left_view, right_view): half-sized chunks from each side
#   :none, nothing:                 cannot fill (caller advances/retreats)
#   :done, nothing:                 globally unfillable, caller stops
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

    # Endpoint ramp: anchor to the filtered data beside the gap, at most
    # `sizegap` samples of it (the scale the fill is built at), so the joins
    # don't kink off a single noisy boundary sample and each gap costs its
    # own length, not the record's. A gap at either end of the record has
    # data on one side only, so it is shifted to join that side.
    start_off = curgap[1] > 1         ? _edge_filtered(x, curgap[1] - 1, -1, sizegap) - fillvec[1]   : NaN
    end_off   = curgap[2] < length(x) ? _edge_filtered(x, curgap[2] + 1, 1, sizegap) - fillvec[end] : NaN
    isnan(start_off) && (start_off = end_off)
    isnan(end_off)   && (end_off = start_off)
    ramp = range(start_off, stop=end_off, length=sizegap)

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
    _howe_fillgaps_core(x) → (xfilled, filled_mask::BitVector)

Fill the `NaN` samples of an equally spaced vector. Returns the filled
vector and a mask of the samples that were filled; a vector without `NaN`
comes back as a copy with an empty mask.

Throws an `ErrorException` if the gaps cannot be filled, as when both ends
of the record are missing and no run between them is long enough to
reflect.
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
        remaining = length(gaps)
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

        # A pass that filled nothing will never fill anything.
        length(gaps) == remaining &&
            error("fillgaps: gap topology cannot be filled (data record too short relative to gap span)")
        # Reset to the gap right of the largest remaining run, otherwise the
        # outer loop would spin forever.
        gap_num = _choose_initial_gap(gaps)
    end

    filled_mask = nan_mask_0 .& .!isnan.(xv)
    return xv, filled_mask
end

"""
    fillgaps(data) → PhaseData or FrequencyData

Copy of a record with its missing samples (`NaN`) filled, after the
reflection method of
[Howe and Schlossberger (2022)](@cite howe-2022-gap-imputation): the
samples before a gap, as many as are missing, are reversed in time,
inverted in sign and placed in the gap, with a line added so that both ends
join the record.

Where the run before a gap is too short, the samples after it are used, and
where both are, half from each side. A gap of one sample takes the mean of
its neighbors. Filling starts at the gap after the longest run of data and
proceeds to the right, then to the left, so a fill may rest on an earlier
one. An error is thrown if the gaps are too long for the data between them.
A record without missing samples is returned as a copy.

!!! warning "A filled record is not the record"
    Fill only for a function that needs a complete record, such as the
    spectra and the cross estimators. The stability estimators leave out
    the terms a missing sample enters (see [Timestamps and gaps](@ref)),
    which does not bias them.

    This implementation departs from the paper in two steps. The paper
    inserts the reflected samples as they are and filters only to find the
    ends that the line joins. Here the inserted samples are themselves
    low-pass filtered, and the filter is weaker than that of the authors'
    script ([Champagne et al. (2021)](@cite champagne-2021-gap-filling-script)).

    In simulation (4096 samples with one gap of 30 %), the Allan deviation
    of a filled phase record was 8 % low at ``τ_0`` and 1 to 5 % low at
    ``16τ_0`` on white PM to flicker FM, and 25 % high at ``τ_0`` on
    random-walk FM. That of a filled frequency record was 5 to 7 % low at
    ``τ_0`` on flicker and random-walk FM, and many times high at long
    ``τ`` on white PM, flicker PM and white FM: 5 to 50 times at
    ``256τ_0``. With the gap left in place it was within 1 % at ``τ_0``
    and within 3 % at ``256τ_0`` on every noise type.

# Examples

```jldoctest
julia> using SigmaTau

julia> x = PhaseData([0.0, 1.0, NaN, 3.0, 4.0], 1.0);

julia> print(fillgaps(x).x)
[0.0, 1.0, 2.0, 3.0, 4.0]
```
"""
function fillgaps(pd::PhaseData)
    xfilled, _ = _howe_fillgaps_core(pd.x)
    return PhaseData(xfilled, pd.tau0; source=pd.source, start=pd.start)
end

function fillgaps(fd::FrequencyData)
    yfilled, _ = _howe_fillgaps_core(fd.y)
    return FrequencyData(yfilled, fd.tau0; source=fd.source, start=fd.start)
end

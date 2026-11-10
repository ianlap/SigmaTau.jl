# Gap imputation by reflection, after D. A. Howe and N. Schlossberger,
# "Characterizing Frequency Stability Measurements Having Multiple Data
# Gaps", IEEE Trans. Ultrason. Ferroelectr. Freq. Control 69(2), 468–472
# (2022), doi:10.1109/TUFFC.2021.3137425, and the script that paper names as
# its algorithm (C. Champagne, D. A. Howe and N. Schlossberger, fillgaps.py,
# doi:10.5281/zenodo.5594587).
#
# Their method: a gap of n samples takes the n samples beside it (the run
# before it, the run after it when that one is too short, half from each
# when neither is long enough), reflected and inverted, and adds a straight
# line so that the low-pass-filtered fill meets the filtered data at both
# ends. Gaps are filled outward from the longest run, so a fill can serve
# the next gap.
#
# Where this file differs from the script, on purpose: the script's filter
# gain is one bin off in the mirrored half of the spectrum; its line stops
# short of the end offset; where half the fill comes from each side it
# filters the left half twice and adds a middle sample only for gaps of
# 4j + 1 samples; and it filters whole runs for the levels beside a gap,
# where a window serves.
#
# What the method leaves open is set from the record (see `fillgaps`): which
# series is filled, phase or frequency, and the span the filter averages
# over. The series' trend (a frequency offset and drift) is taken out first
# and put back after, since the filter pads a block with its mirror image,
# which turns a slope into a triangle wave and pulls the filtered endpoints
# off the data.

# The low-pass filter of the paper's eq. 3 as the script has it: `x` padded
# with its mirror image on both sides, its spectrum multiplied by
# exp(−span·|k|/L) over the L bins, the middle third returned. It averages
# over about `span` samples; the script's constant is 16.
function _howe_filter(x::AbstractVector{<:Real}, span::Int)
    n = length(x)
    n == 0 && return Float64[]
    L = 3n
    F = fft(Float64[x[i <= n ? n - i + 1 : i <= 2n ? i - n : L - i + 1] for i in 1:L])
    @inbounds for k in 0:L-1
        F[k + 1] *= exp(-span * min(k, L - k) / L)
    end
    return real.(ifft(F))[n+1:2n]
end

# The filtered value at `i`, the end of a run of samples, from up to `w`
# samples of that run reaching away from the gap (`dir` -1 leftward, +1 rightward).
function _edge_filtered(x::Vector{Float64}, i::Int, dir::Int, w::Int, span::Int)
    j = i
    while abs(j - i) + 1 < w && 1 <= j + dir <= length(x) && !isnan(x[j + dir])
        j += dir
    end
    f = _howe_filter(view(x, min(i, j):max(i, j)), span)
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
# Returns (status, new_gap_num) where status is :ok, :advance, :done. The
# filter averages over `span` samples. With `free` the level after the gap is
# unknown (the phase of a frequency record): the fill joins the data before
# it and the data after it are moved to meet the fill.
function _fill_one_gap!(x::Vector{Float64},
                        gaps::Vector{Tuple{Int,Int}},
                        gap_num::Int,
                        reverse::Bool, span::Int, free::Bool)
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

    # The reflected and inverted samples, inserted as they are.
    fillvec = if side === :both
        left_pts, right_pts = pts
        left  = .-reverse_view(left_pts)
        right = .-reverse_view(right_pts)
        right .+= _howe_filter(left, span)[end] - _howe_filter(right, span)[1]
        isodd(sizegap) ? vcat(left, 0.5 * (left[end] + right[1]), right) : vcat(left, right)
    else
        .-reverse_view(pts)
    end

    # The line that brings the filtered fill to the filtered data beside the
    # gap, read from a window of at least 64 samples. A gap at either end of
    # the record has data on one side only and is shifted to join that side.
    filtered = _howe_filter(fillvec, span)
    w = max(sizegap, 64, span)
    start_off = curgap[1] > 1         ? _edge_filtered(x, curgap[1] - 1, -1, w, span) - filtered[1]   : NaN
    end_off   = curgap[2] < length(x) ? _edge_filtered(x, curgap[2] + 1, 1, w, span) - filtered[end] : NaN
    if free && !isnan(start_off) && !isnan(end_off)
        @views x[curgap[2]+1:end] .+= start_off - end_off
        end_off = start_off
    end
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

# The reflection fill of every gap of `x` (NaN runs), with the filter span and
# `free` of `_fill_one_gap!`; gaps of one sample are filled before it (see
# `_fillgaps`). Throws when the layout cannot be filled, for instance gaps at
# both record ends with no run between them long enough to mirror.
function _howe_fill(x::Vector{Float64}, span::Int; free::Bool=false)
    xv = copy(x)
    gaps = _identify_gaps(xv)
    isempty(gaps) && return xv
    unfillable() = error("fillgaps: gap topology cannot be filled (data record too short relative to gap span)")
    gap_num = _choose_initial_gap(gaps)

    # Outer loop: rightward sweep, then leftward sweep, repeat until empty.
    while !isempty(gaps)
        remaining = length(gaps)
        while 1 <= gap_num <= length(gaps)
            status, gap_num = _fill_one_gap!(xv, gaps, gap_num, false, span, free)
            status === :done && unfillable()
        end
        isempty(gaps) && break

        gap_num = min(gap_num, length(gaps))
        while gap_num >= 1 && !isempty(gaps)
            status, gap_num = _fill_one_gap!(xv, gaps, gap_num, true, span, free)
            status === :done && unfillable()
        end
        isempty(gaps) && break

        # A pass that filled nothing will never fill anything.
        length(gaps) == remaining && unfillable()
        gap_num = _choose_initial_gap(gaps)
    end
    return xv
end

# The least-squares polynomial of `degree` 1 or 2 through the measured
# samples of `v`: the trend a fill is made around. A frequency offset is a
# line in phase and a constant in frequency; a frequency drift a parabola and
# a line.
function _measured_trend(v::Vector{Float64}, degree::Int)
    ok = isfinite.(v)
    count(ok) <= degree && return zeros(length(v))
    c = _index_polyfit(v, ok, degree)
    mid = (length(v) - 1) / 2
    return [c[1] + c[2] * (i - 1 - mid) + c[3] * (i - 1 - mid)^2 for i in eachindex(v)]
end

# The span, an octave of samples, to average the level of `v` over when it
# is carried across a gap, and whether `v` is redder than white at its
# shortest spans or beyond them (its running sum is then too smooth to fill;
# see `fillgaps`). Both are read from the two-sample deviation of `v` at the
# octaves m that at least 8m terms measure past the gaps. The error of one
# block mean as a prediction of the next is that deviation at the block
# length, so the span is the longest octave whose deviation is within 10 % of
# the lowest. The series is red when the deviation has turned up again by the
# last octave, or has nowhere fallen to 1/√2 of its first value, the fall of
# white noise in one octave.
function _level_span(v::Vector{Float64})
    n = count(isfinite, v)
    n < 4 && return 1, false
    ms = [2^k for k in 0:max(0, floor(Int, log2(n / 10)))]
    r = adev(FrequencyData(v, 1.0), ms; ci=false)
    ok = findall(i -> r.neff[i] >= 8 * ms[i], eachindex(ms))
    isempty(ok) && return 1, false
    d = view(r.dev, ok)
    lowest = minimum(d)
    k = findlast(<=(1.1 * lowest), d)
    return ms[ok[k]], k < length(d) || (length(d) > 1 && lowest > d[1] / sqrt(2))
end

# The missing phases of `x` from the filled first differences `y`
# (y[k] = x[k+1] − x[k]). Between two measured phases the differences are
# first corrected to sum to the measured change, by a parabola that is
# largest mid-gap and fades to the joins.
function _phase_from_steps!(x::Vector{Float64}, y::Vector{Float64})
    for (g1, g2) in _identify_gaps(x)
        if g1 == 1
            for i in g2:-1:1
                x[i] = x[i + 1] - y[i]
            end
            continue
        end
        if g2 < length(x)
            n = g2 - g1 + 2
            e = x[g2 + 1] - x[g1 - 1] - sum(view(y, g1-1:g2))
            for j in 1:n
                y[g1 - 2 + j] += e * (j - 0.5) * (n - j + 0.5) * 12 / (n * (2n^2 + 1))
            end
        end
        for i in g1:g2
            x[i] = x[i - 1] + y[i - 1]
        end
    end
    return x
end

# `v` (phase when `phase`, else frequency) with its gaps filled; the measured
# samples are returned as they came.
function _fillgaps(v::AbstractVector{<:Real}, phase::Bool, form::Symbol)
    form in (:auto, :phase, :frequency) ||
        throw(ArgumentError("fillgaps: form must be :auto, :phase or :frequency, got :$form"))
    out = Vector{Float64}(v)
    _fill_singletons!(out, _identify_gaps(out))
    gap = isnan.(out)
    any(gap) || return out
    y = phase ? diff(out) : out
    line = _measured_trend(y, 1)
    y0 = y .- line
    span, red = _level_span(y0)
    if form === :frequency || (form === :auto && red)
        yf = _howe_fill(y0, span) .+ line
        phase ? _phase_from_steps!(out, yf) : (out[gap] = yf[gap])
    elseif phase
        line = _measured_trend(out, 2)
        x0 = out .- line
        out[gap] = (_howe_fill(x0, first(_level_span(x0))) .+ line)[gap]
    else
        # Each run of frequency samples integrates to phase on its own level,
        # 0 at its start.
        x0 = fill(NaN, length(y0) + 1)
        for k in eachindex(y0)
            isnan(y0[k]) && continue
            isnan(x0[k]) && (x0[k] = 0.0)
            x0[k + 1] = x0[k] + y0[k]
        end
        out[gap] = (diff(_howe_fill(x0, first(_level_span(x0)); free=true)) .+ line)[gap]
    end
    return out
end

"""
    fillgaps(data::PhaseData;     form=:auto) → PhaseData
    fillgaps(data::FrequencyData; form=:auto) → FrequencyData

Fill the missing (`NaN`) samples of an equispaced record by the reflection
method of D. A. Howe and N. Schlossberger, IEEE Trans. UFFC 69(2), 468–472
(2022), doi:10.1109/TUFFC.2021.3137425, as their script has it (C. Champagne,
D. A. Howe and N. Schlossberger, `fillgaps.py`, doi:10.5281/zenodo.5594587),
so that an estimator which cannot skip gaps can run. The measured samples are
returned as they came, with the same `tau0`, `source` and `start`. The readers
put timestamped samples on their grid first, so `read_data(path;
fillgaps=true)` fills the readings a file is missing.

A gap of n samples takes the n samples beside it (the run before it, the run
after it when that one is too short, half from each when neither is long
enough), reflected and inverted, and a straight line is added so that the
low-pass-filtered fill meets the filtered data at both ends. A gap of one
sample takes the mean of its neighbours. Gaps are filled outward from the
longest run, so a fill can serve the next gap; a layout with no run long
enough to mirror throws an `ErrorException`.

The reflected samples carry the noise of the record only where the series
being filled lies between white and random walk. A smoother series shows the
added line as a step at each join (the phase of random-walk FM read 10 times
high at τ₀ after a 30 % fill), and in a rougher one the error of the level at
the joins is integrated across the gap (the frequency of white PM read 12
times high at 256 τ₀). Two things the paper leaves open are therefore set
from the record:

- `form`, the series that is filled. `:auto` reads the two-sample deviation
  of the frequency samples past the gaps, at the octaves of m samples that
  at least 8m terms measure. It fills the frequency when that deviation
  nowhere falls to 1/√2 of its value at τ₀, the fall of white FM in one
  octave, or turns up again by the last octave: flicker or random-walk FM
  is then present and the phase too smooth. It fills the phase otherwise.
  `:phase` and `:frequency` force the choice. A phase record filled as
  frequency has the filled differences corrected to sum to the measured
  phase change across each gap; a frequency record filled as phase leaves
  the phase step across each gap free.
- The span the filter averages over (16 samples in the script): the longest
  of those octaves at which the two-sample deviation of the series being
  filled is within 10 % of its lowest, so the level beside a gap is averaged
  for as long as averaging improves it.

The trend of the series, a least-squares parabola in phase or line in
frequency through the measured samples, is set aside while it is filled. The
filter mirrors each block at its ends, so a slope would pull the filtered
endpoints off the data: with a frequency offset of 100 σ_y(τ₀) a phase record
filled as scripted read 48 times high at τ₀.

On simulated records of 4096 samples (twelve mixes of the five power-law
noises, as phase and as frequency records, one to six gaps removing up to
30 %, with and without offset and drift, 100 records each) the Allan
deviation of the filled record stayed within 5 % of the original's from τ₀ to
256 τ₀, or 11 % at 256 τ₀ for a 30 % gap at the end of the record, which has
data on one side only; the modified Allan deviation stayed within 4 % to
16 τ₀. Three things a fill does not recover:

- MDEV at the longest τ where white PM dominates: 1.1 to 1.3 times high, and
  1.9 for a frequency record with six gaps. Averaging m phases resolves a
  level error that ADEV does not see.
- The phase step across a gap in a frequency record, which is not measured,
  so the errors of many gaps add. With gaps of 1 to 12 samples removing a
  quarter of a frequency record, white PM read 1.3 and 1.9 times high in ADEV
  at 64 and 256 τ₀ (4.8 and 19 in MDEV), flicker PM 1.5 at 256 τ₀ and white FM
  1.06; the same gaps in a phase record stayed within 8 % in ADEV.
- A periodic term that stands out of the noise, which neither form carries in
  every layout. A sine in the frequency (period 512 samples, amplitude
  3 σ_y(τ₀)) on white FM read 1.15 high in ADEV and 1.40 in MDEV at 256 τ₀
  across one 30 % gap, where the frequency is filled and reflection inverts
  the sine (1.00 with `form=:phase`), and 1.42 at τ₀ past gaps of 1 to 12
  samples, where the phase is filled (1.00 with `form=:frequency`; 1.05 at a
  third of the amplitude). Remove such a term before filling and put it back
  after.

The estimators that skip gaps (`adev`, `mdev`, `hdev`, `mhdev` and their time
forms, on the record with its `NaN`s) need no fill and have none of these
errors. The choice of form and span is this package's, checked by those
simulations; the paper prescribes neither.
"""
fillgaps(pd::PhaseData; form::Symbol=:auto) =
    PhaseData(_fillgaps(pd.x, true, form), pd.tau0; source=pd.source, start=pd.start)

fillgaps(fd::FrequencyData; form::Symbol=:auto) =
    FrequencyData(_fillgaps(fd.y, false, form), fd.tau0; source=fd.source, start=fd.start)

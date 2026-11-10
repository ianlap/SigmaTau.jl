# Gap-tolerant estimation for records with NaN samples.
#
# A missing phase sample is a missing reading on an otherwise intact time
# scale, so a phase-record term counts when every sample it reads was
# measured: ADEV, HDEV and TIE rms need their 3, 4 or 2 points, MDEV and
# MHDEV need the whole span they average over, and MTIE takes the
# peak-to-peak of the measured samples in each window. A missing frequency
# sample loses the phase step across it, so a frequency-record term counts
# only when no missing sample falls inside its span. This is the rule of
# allantools' `gradev`, extended to the other estimators.
#
# Noise identification runs on the record with its gaps bridged, and the EDF
# for each τ uses the record length that would give the number of terms
# actually summed.

struct _GapRecord
    x::Vector{Float64}     # phase, 0 at missing samples
    bad::Vector{Int}       # bad[k+1] = missing units among the first k
    edges::Bool            # units are frequency samples, not phase samples
    noise_x::Vector{Float64}  # gap-free phase for noise identification
    tau0::Float64
end

_has_gaps(data::PhaseData) = any(isnan, data.x)
_has_gaps(data::FrequencyData) = any(isnan, data.y)

function _missing_prefix(v::AbstractVector{<:Real})
    bad = Vector{Int}(undef, length(v) + 1)
    bad[1] = 0
    @inbounds for k in eachindex(v)
        bad[k+1] = bad[k] + !isfinite(v[k])
    end
    return bad
end

function _gap_record(data::PhaseData)
    x = Float64[isfinite(v) ? Float64(v) : 0.0 for v in data.x]
    noise_x = _gap_noise_record(_f64(data.x), false)
    return _GapRecord(x, _missing_prefix(data.x), false, noise_x, data.tau0)
end

function _gap_record(data::FrequencyData)
    y = Float64[isfinite(v) ? Float64(v) : 0.0 for v in data.y]
    x = _freq_to_phase(FrequencyData(y, data.tau0)).x
    noise_x = _gap_noise_record(_f64(data.y), true) .* data.tau0
    return _GapRecord(x, _missing_prefix(data.y), true, noise_x, data.tau0)
end

# A gap-free phase record for noise identification from `v` (phase, or
# frequency when `freq`): the measured span with interior gaps bridged by
# straight lines. Noise ID only needs the record's character, and this cannot
# fail on any gap layout.
function _gap_noise_record(v::Vector{Float64}, freq::Bool)
    ok = findall(isfinite, v)
    isempty(ok) && return Float64[]
    lo, hi = first(ok), last(ok)
    w = v[lo:hi]
    prev = 1
    for k in 2:length(w)
        isfinite(w[k]) || continue
        for j in prev+1:k-1
            w[j] = w[prev] + (w[k] - w[prev]) * (j - prev) / (k - prev)
        end
        prev = k
    end
    freq || return w
    x = Vector{Float64}(undef, length(w) + 1)
    x[1] = 0.0
    acc = 0.0
    @inbounds for k in eachindex(w)
        acc += w[k]
        x[k+1] = acc
    end
    return x
end

# No missing unit among the samples x[a..b] read by a term (for a frequency
# record, the steps between them).
@inline _span_ok(g::_GapRecord, a::Int, b::Int) =
    g.edges ? g.bad[b] == g.bad[a] : g.bad[b+1] == g.bad[a]

# x[k] was measured (always true for a frequency record, whose terms are
# checked by span).
@inline _point_ok(g::_GapRecord, k::Int) = g.edges || g.bad[k+1] == g.bad[k]

# Prefix sums of the phase with a least-squares line through the measured
# samples removed, for the averaging estimators. The line cancels in every
# term; removing it keeps the prefix sums small.
function _gap_prefix(g::_GapRecord)
    x = g.x
    N = length(x)
    n = 0; si = 0.0; sx = 0.0; sii = 0.0; six = 0.0
    @inbounds for i in 1:N
        _point_ok(g, i) || continue
        n += 1; si += i; sx += x[i]; sii += Float64(i)^2; six += i * x[i]
    end
    den = n * sii - si^2
    slope = den > 0 ? (n * six - si * sx) / den : 0.0
    icpt = n > 0 ? (sx - slope * si) / n : 0.0
    cs = Vector{Float64}(undef, N + 1)
    cs[1] = 0.0
    @inbounds for i in 1:N
        cs[i+1] = cs[i] + (_point_ok(g, i) ? x[i] - (icpt + slope * i) : 0.0)
    end
    return cs
end

# Sum of squared terms over the valid starts, and how many there were. `S`
# is the last offset the term reads; `pts` the offsets a pointwise term reads
# (empty for a term that averages over its whole span).
function _gap_sum(term, g::_GapRecord, S::Int, pts::NTuple{K,Int}) where {K}
    s = 0.0
    n = 0
    @inbounds for i in 1:length(g.x)-S
        ok = (K == 0 || g.edges) ? _span_ok(g, i, i + S) : all(o -> _point_ok(g, i + o), pts)
        ok || continue
        d = term(i)
        s += d * d
        n += 1
    end
    return s, n
end

# Deviation and term count per m for the σ_y estimators and TIE rms.
function _gap_core(est::Symbol, g::_GapRecord, m_values::Vector{Int})
    x = g.x
    cs = est === :mdev || est === :mhdev ? _gap_prefix(g) : Float64[]
    devs = Vector{Float64}(undef, length(m_values))
    ns = Vector{Int}(undef, length(m_values))
    for (k, m) in enumerate(m_values)
        if est === :adev
            s, n = _gap_sum(i -> x[i+2m] - 2x[i+m] + x[i], g, 2m, (0, m, 2m))
            v, nmin = s / (2.0 * n * Float64(m)^2 * g.tau0^2), 2
        elseif est === :hdev
            s, n = _gap_sum(i -> x[i+3m] - 3x[i+2m] + 3x[i+m] - x[i], g, 3m, (0, m, 2m, 3m))
            v, nmin = s / (6.0 * n * Float64(m)^2 * g.tau0^2), 2
        elseif est === :mdev
            s, n = _gap_sum(i -> cs[i+3m] - 3cs[i+2m] + 3cs[i+m] - cs[i], g, 3m - 1, ())
            v, nmin = s / (2.0 * n * Float64(m)^4 * g.tau0^2), 2
        elseif est === :mhdev
            s, n = _gap_sum(i -> cs[i+4m] - 4cs[i+3m] + 6cs[i+2m] - 4cs[i+m] + cs[i], g, 4m - 1, ())
            v, nmin = s / (6.0 * n * Float64(m)^4 * g.tau0^2), 2
        else # :tierms
            s, n = _gap_sum(i -> x[i+m] - x[i], g, m, (0, m))
            v, nmin = s / n, 1
        end
        ok = m >= 1 && n >= nmin
        devs[k] = ok ? sqrt(v) : NaN
        ns[k] = ok ? n : 0
    end
    return devs, ns
end

# MTIE: largest peak-to-peak over the measured samples of any window of
# m + 1 samples; for a frequency record, windows stay inside a gap-free run.
# Counts the windows holding at least two measured samples.
function _gap_mtie(g::_GapRecord, m_values::Vector{Int})
    N = length(g.x)
    devs = Vector{Float64}(undef, length(m_values))
    ns = Vector{Int}(undef, length(m_values))
    measured = [_point_ok(g, i) for i in 1:N]
    for (k, m) in enumerate(m_values)
        best, n = 0.0, 0
        lo = 1
        while lo <= N
            # A run of samples with no missing unit between them.
            hi = lo
            while hi < N && (g.edges ? g.bad[hi+1] == g.bad[hi] : true)
                hi += 1
            end
            b, c = _mtie_window_max(g.x, measured, lo, hi, m)
            best = max(best, b); n += c
            lo = hi + 1
        end
        devs[k] = n > 0 ? best : NaN
        ns[k] = n
    end
    return devs, ns
end

function _mtie_window_max(x::Vector{Float64}, measured::Vector{Bool}, lo::Int, hi::Int, m::Int)
    best, count = 0.0, 0
    hi - lo < m && return best, count
    # Monotone deques of indices, popped from the front by moving a head.
    maxq = Vector{Int}(undef, hi - lo + 1); minq = similar(maxq)
    h1 = h2 = 1; t1 = t2 = 0
    for j in lo:hi
        # Window [j − m, j]; drop what has left it.
        while h1 <= t1 && maxq[h1] < j - m; h1 += 1; end
        while h2 <= t2 && minq[h2] < j - m; h2 += 1; end
        if measured[j]
            while h1 <= t1 && x[maxq[t1]] <= x[j]; t1 -= 1; end
            maxq[t1 += 1] = j
            while h2 <= t2 && x[minq[t2]] >= x[j]; t2 -= 1; end
            minq[t2 += 1] = j
        end
        # Only full windows count; a partial one lies inside the first full one.
        j - lo >= m && h1 <= t1 && (best = max(best, x[maxq[h1]] - x[minq[h2]]))
    end
    # Windows [i, i + m] inside lo..hi with at least two measured samples.
    cm = cumsum([measured[j] for j in lo:hi])
    for i in lo:hi-m
        a = i - lo + 1; b = a + m
        (cm[b] - (a > 1 ? cm[a-1] : 0)) >= 2 && (count += 1)
    end
    return best, count
end

# Record length whose gap-free term count equals `n` terms of `est` at `m`.
function _gap_equivalent_n(est::Symbol, m::Int, n::Int)
    est === :adev  && return n + 2m
    est === :mdev  && return n + 3m - 1
    est === :hdev  && return n + 3m
    return n + 4m - 1   # :mhdev
end

"""
Gap-tolerant result for `est` (adev, mdev, tdev, hdev, mhdev, htdev, tierms
or mtie); see the note at the top of this file.
"""
function _gap_result(est::Symbol, g::_GapRecord, m_values::Vector{Int};
                     ci::Bool=true, remove_drift::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    _check_m(est, m_values)
    pd = PhaseData(g.x, g.tau0)
    taus = m_values .* g.tau0
    if est === :tdev || est === :htdev
        base = _gap_result(est === :tdev ? :mdev : :mhdev, g, m_values; ci=ci, remove_drift=remove_drift, confidence=confidence)
        factor = taus ./ sqrt(est === :tdev ? 3.0 : 10.0 / 3.0)
        ci || return _result(pd, confidence, est, taus, base.dev .* factor, Symbol[], _CIBound[], Float64[], base.neff)
        return _result(pd, confidence, est, taus, base.dev .* factor, base.noise_type,
                       _scale_ci(base.ci, factor), base.edf, base.neff)
    end
    devs, ns = est === :mtie ? _gap_mtie(g, m_values) : _gap_core(est, g, m_values)
    if !ci || est === :mtie || est === :tierms
        return _result(pd, confidence, est, taus, devs, Symbol[], _CIBound[], Float64[], ns)
    end
    noises = est === :adev || est === :mdev ? identify_noise(g.noise_x, m_values, dmin=0, dmax=2) :
                                              _hadamard_noise(g.noise_x, m_values, remove_drift)
    edfs = Vector{Float64}(undef, length(m_values))
    lower = similar(edfs); upper = similar(edfs)
    for k in eachindex(m_values)
        m = m_values[k]
        Neq = _gap_equivalent_n(est, m, ns[k])
        T = (Neq - 1) * g.tau0
        edfs[k] = only(calculate_edf(est, devs[k:k], noises[k:k], [m], taus[k:k], Neq, T))
        lo, hi = confidence_intervals(devs[k:k], edfs[k:k], confidence)
        lower[k] = only(lo); upper[k] = only(hi)
    end
    return _result(pd, confidence, est, taus, devs, noises, _zip_ci(lower, upper), edfs, ns)
end

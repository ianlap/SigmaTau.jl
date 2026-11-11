# The distribution of a record's samples and their autocorrelation: what
# Stable32's Statistics, Histogram and Autocorrelation functions show.

"""
    sample_statistics(data::PhaseData) → NamedTuple
    sample_statistics(data::FrequencyData) → NamedTuple

Basic statistics of a record's samples, the quantities of Stable32's
Statistics function
([Stable32 User Manual](@cite hamilton-2008-stable32-manual)):

- `n`: number of samples measured.
- `gaps`: number of `NaN` samples.
- `min`, `max`, `mean`, `median`: of the samples measured, phase in seconds
  or fractional frequency.
- `std`: their standard deviation, with ``n - 1`` in the denominator.
- `std_y`: that of the fractional frequency, the Std Dev that Stable32
  shows for either kind of record: of the samples for a frequency record,
  of ``(x_{i+1} - x_i)/τ_0`` where both were measured for a phase record.
- `adev`: the Allan deviation at ``τ_0``, from the terms whose samples were
  all measured, as [`adev`](@ref) takes them.
- `noise`: the dominant noise type at ``τ_0``, by [`identify_noise`](@ref)
  on the record with each gap bridged by a straight line.

The record needs three samples measured.

!!! note "Standard deviation of a phase record"
    `std` is that of the samples as given. Under any FM noise that of phase
    grows with the length of the record, so it describes this record only.
    `std_y` is the one to read against `adev`: for white FM the two agree,
    apart from the scatter of the estimates.

# Examples

```jldoctest
julia> using SigmaTau

julia> s = sample_statistics(FrequencyData([1.0, 3.0, 2.0, NaN, 6.0], 1.0));

julia> (s.n, s.gaps, s.min, s.max, s.mean, s.median)
(4, 1, 1.0, 6.0, 3.0, 2.5)

julia> s.adev ≈ sqrt(((3 - 1)^2 + (2 - 3)^2) / 4)
true
```
"""
sample_statistics(data::PhaseData) =
    _sample_statistics(_f64(data.x), data.tau0, false)
sample_statistics(data::FrequencyData) =
    _sample_statistics(_f64(data.y), data.tau0, true)

function _sample_statistics(v::Vector{Float64}, tau0::Float64, freq::Bool)
    w = filter(isfinite, v)
    n = length(w)
    n >= 3 || throw(ArgumentError("sample_statistics needs three samples measured, got $n"))
    # σ_y(τ₀) from the first differences of frequency or the second of phase,
    # each term only where every sample it reads was measured.
    s = 0.0; k = 0
    d = freq ? 1 : 2
    @inbounds for i in 1:length(v)-d
        t = freq ? v[i+1] - v[i] : (v[i+2] - 2v[i+1] + v[i]) / tau0
        isfinite(t) || continue
        s += t * t; k += 1
    end
    y = freq ? w : filter(isfinite, diff(v) ./ tau0)
    noise_x = _gap_noise_record(v, freq)
    return (n = n, gaps = length(v) - n, min = minimum(w), max = maximum(w),
            mean = mean(w), median = _middle(w), std = std(w),
            std_y = length(y) >= 2 ? std(y) : NaN,
            adev = k > 0 ? sqrt(s / (2k)) : NaN,
            noise = identify_noise(noise_x, [1])[1])
end

"""
    sample_histogram(data::PhaseData; bins) → NamedTuple
    sample_histogram(data::FrequencyData; bins) → NamedTuple

Histogram of a record's samples, the plot of Stable32's Histogram function
([Stable32 User Manual](@cite hamilton-2008-stable32-manual)), with what the
normal curve drawn over it needs. A distribution far from normal, two peaks
for instance, points to a problem with the record
([SP1065 §5.17](@cite riley-2008-sp1065)).

The bins have equal widths and span the samples measured. By default their
width follows Scott's rule,

```math
h = 3.49\\, s\\, n^{-1/3},
```

which minimizes the integrated mean squared error of the histogram of a
normal sample of ``n`` values with standard deviation ``s``
([Scott (1979)](@cite scott-1979-histograms)), with at most 1000 bins.

The normal curve with the mean and standard deviation of the samples gives
the expected count in a bin of width ``h`` centered on ``v``:

```math
c(v) = \\frac{n h}{s \\sqrt{2π}} \\exp\\left(-\\frac{(v - \\bar{v})^2}{2 s^2}\\right).
```

# Keywords

- `bins`: number of bins; Scott's rule by default.

# Returns

A named tuple:

- `edges`: the `bins + 1` bin edges, rising. A record whose samples are all
  equal has one bin, of zero width.
- `counts`: samples in each bin. A bin holds its lower edge, and the last
  one also its upper edge.
- `n`, `mean`, `std`: number of samples measured, their mean and standard
  deviation, as [`sample_statistics`](@ref) gives them.

# Examples

```jldoctest
julia> using SigmaTau

julia> h = sample_histogram(FrequencyData([0.0, 1.0, 1.5, 2.0, 4.0], 1.0); bins = 2);

julia> (h.edges, h.counts)
([0.0, 2.0, 4.0], [3, 2])
```
"""
sample_histogram(data::PhaseData; bins::Integer = _scott_bins(_f64(data.x))) =
    _sample_histogram(_f64(data.x), Int(bins))
sample_histogram(data::FrequencyData; bins::Integer = _scott_bins(_f64(data.y))) =
    _sample_histogram(_f64(data.y), Int(bins))

const _MAX_BINS = 1000

function _scott_bins(v::Vector{Float64})
    w = filter(isfinite, v)
    length(w) >= 2 || return 1
    lo, hi = extrema(w)
    h = 3.49 * std(w) * length(w)^(-1 / 3)
    h > 0 || return 1
    return clamp(ceil(Int, (hi - lo) / h), 1, _MAX_BINS)
end

function _sample_histogram(v::Vector{Float64}, bins::Int)
    bins >= 1 || throw(ArgumentError("bins must be at least 1, got $bins"))
    w = filter(isfinite, v)
    n = length(w)
    n >= 2 || throw(ArgumentError("sample_histogram needs two samples measured, got $n"))
    lo, hi = extrema(w)
    counts = zeros(Int, hi > lo ? bins : 1)
    nb = length(counts)
    edges = nb == 1 ? [lo, hi] : collect(range(lo, hi; length = nb + 1))
    width = (hi - lo) / nb
    @inbounds for x in w
        k = clamp(width > 0 ? floor(Int, (x - lo) / width) + 1 : 1, 1, nb)
        # Against the edges themselves, which round apart from the division.
        k < nb && x >= edges[k+1] && (k += 1)
        k > 1 && x < edges[k] && (k -= 1)
        counts[k] += 1
    end
    return (edges = edges, counts = counts, n = n, mean = mean(w), std = std(w))
end

"""
    acf(data::PhaseData; m=1, lags) → NamedTuple
    acf(data::FrequencyData; m=1, lags) → NamedTuple

Autocorrelation of a record at averaging factor ``m``, the plot of
Stable32's Autocorrelation function
([Stable32 User Manual](@cite hamilton-2008-stable32-manual)), estimated as

```math
r_k = \\frac{\\sum_{t=1}^{N-k} (z_t - \\bar{z})(z_{t+k} - \\bar{z})}
           {\\sum_{t=1}^{N} (z_t - \\bar{z})^2}, \\qquad k = 0, 1, …, K
```

([SP1065 §5.5.3, Eq. (47)](@cite riley-2008-sp1065)). The series ``z_t``
is the phase taken every ``m`` samples, or the frequency averaged over
``m`` samples without overlap: the series that the lag-1 noise
identification reads ([SP1065 §5.5.5](@cite riley-2008-sp1065)). White
frequency noise has ``r_k ≈ 0`` beyond lag 0, white phase noise read as
frequency has ``r_1 ≈ -1/2``, and the divergent noises stay near 1 over
many lags; a periodic component shows as an oscillation.

# Keywords

- `m=1`: averaging factor, ``τ = m τ_0``.
- `lags`: the largest lag ``K``, from 0 to ``N - 1``; by default
  ``N ÷ 2``, half the series, as Stable32 computes them.

# Returns

A named tuple:

- `lag`, `r`: the lags ``0, …, K`` and ``r_k``.
- `m`, `tau`, `n`: averaging factor, ``τ`` in seconds and ``N``, the
  number of points in the series.
- `alpha`: the lag-1 estimate of the power-law exponent ``α``, a real
  number that falls between two noise types for a mix of them
  ([Riley and Greenhall (2004)](@cite riley-2004-lag1-acf)). It is read as
  [`identify_noise`](@ref) reads it: on the phase without its quadratic
  drift, differenced until stationary. `NaN` when the phase holds fewer
  than ``30 m`` samples, where [`identify_noise`](@ref) stops using it.
- `noise`: the label [`identify_noise`](@ref) gives at ``m``.

A record with `NaN` gaps throws an error: fill them with
[`fillgaps`](@ref) first.

# Examples

```jldoctest
julia> using SigmaTau

julia> a = acf(FrequencyData([1.0, -1.0, 1.0, -1.0, 1.0, -1.0], 1.0); lags = 2);

julia> a.r ≈ [1.0, -5 / 6, 4 / 6]
true
```
"""
function acf(data::PhaseData; m::Integer = 1, lags::Integer = -1)
    _check_acf_m(m)
    x = _f64(data.x)
    return _acf(x, x[1:m:end], m, lags, data.tau0)
end

function acf(data::FrequencyData; m::Integer = 1, lags::Integer = -1)
    _check_acf_m(m)
    y = _f64(data.y)
    z = Float64[sum(@view y[(j - 1) * m + 1:j * m]) / m for j in 1:length(y) ÷ m]
    return _acf(_freq_to_phase(data).x, z, m, lags, data.tau0)
end

_check_acf_m(m::Integer) = m >= 1 || throw(ArgumentError("acf: m must be at least 1, got $m"))

function _acf(x::Vector{Float64}, z::Vector{Float64}, m::Integer, lags::Integer, tau0::Float64)
    _anynan(x) && throw(ArgumentError(
        "acf: the record has NaN gaps; fill them with fillgaps first"))
    N = length(z)
    N >= 3 || throw(ArgumentError("acf: m = $m leaves $N points; it needs 3"))
    K = lags == -1 ? N ÷ 2 : Int(lags)
    0 <= K <= N - 1 || throw(ArgumentError("acf: lags must lie between 0 and $(N - 1), got $lags"))
    zbar = mean(z)
    L = nextpow(2, N + K)   # padded so that no lag up to K wraps around
    pad = zeros(L)
    @inbounds for t in 1:N
        pad[t] = z[t] - zbar
    end
    c = irfft(abs2.(rfft(pad)), L)
    c[1] > 0 || throw(ArgumentError("acf: the series is constant"))
    alpha = length(x) ÷ m >= NEFF_RELIABLE ?
        _noise_id_lag1acf(_preprocess(x), Int(m), 0, 2; detrend=true)[1] : NaN
    return (lag = collect(0:K), r = c[1:K+1] ./ c[1], m = Int(m), tau = m * tau0, n = N,
            alpha = alpha, noise = identify_noise(x, [Int(m)])[1])
end

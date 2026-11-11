# Two-channel (cross) estimators: two simultaneous measurements of one
# device, each against its own reference or through its own channel. The
# cross-spectrum and the cross Allan variance keep what the channels share,
# the device, and average away what they do not, each reference and each
# channel's own noise, down to a floor set by the number of averages.

"""
    CrossSpectralResult

Result of [`cross_psd_y`](@ref) and [`cross_psd_x`](@ref).

# Fields

$(TYPEDFIELDS)
"""
struct CrossSpectralResult
    "Quantity estimated: `:Sy` (fractional frequency, 1/Hz) or `:Sx` (time, s²/Hz)."
    spectral_type::Symbol
    "Fourier frequencies ``f``, in Hz."
    freq::Vector{Float64}
    "Real part of the averaged cross-spectrum: the spectrum of the noise common to both channels."
    psd::Vector{Float64}
    """Rms magnitude of the residual that the noise the channels do not share
    leaves in the average, ``\\sqrt{S_{aa} S_{bb}/K}``. ``K`` is the number of
    segments reduced for the correlation of overlapping ones
    ([Welch (1967), p. 71](@cite welch-1967-psd)), 0.95 of `nseg` with the
    default half-overlapped Hann taper. Where the channels share nothing,
    `psd` scatters about zero with a standard deviation of `floor`/√2. A bin
    where `psd` is not well above the floor is not resolved."""
    floor::Vector{Float64}
    "Spectrum of channel a alone."
    psd_a::Vector{Float64}
    "Spectrum of channel b alone."
    psd_b::Vector{Float64}
    "Number of Welch segments averaged."
    nseg::Int
    "Welch segment length, in samples."
    nperseg::Int
    "Overlap of successive segments, in samples."
    noverlap::Int
    "Taper applied to each segment: `:hann`, `:hamming` or `:rectangular`."
    taper::Symbol
end

function Base.show(io::IO, r::CrossSpectralResult)
    print(io, "CrossSpectralResult(:", r.spectral_type, ", ", length(r.freq), " bins, ",
              r.nseg, " averages, nperseg=", r.nperseg, ")")
end

function _cross_pair(a::AbstractTimingData, b::AbstractTimingData; gaps::Bool=false)
    (a isa PhaseData) == (b isa PhaseData) ||
        throw(ArgumentError("cross: pass two PhaseData or two FrequencyData records"))
    gaps || (any(isnan, _samples(a)) || any(isnan, _samples(b))) &&
        throw(ArgumentError("cross: the records have gaps (NaN); fill them first with fillgaps"))
    a.tau0 == b.tau0 ||
        throw(ArgumentError("cross: the channels have different τ₀ ($(a.tau0) vs $(b.tau0))"))
    length(a) == length(b) ||
        throw(ArgumentError("cross: the channels differ in length ($(length(a)) vs $(length(b)))"))
end

function _cross_psd(type::Symbol, za::Vector{Float64}, zb::Vector{Float64}, tau0::Float64,
                    nperseg, noverlap, taper::Symbol)
    np, no = _welch_params(length(za), nperseg, noverlap)
    freq, Saa, Sbb, Sab, nseg = _welch(za, zb, 1 / tau0; nperseg=np, noverlap=no, window=taper)
    floor = sqrt.(Saa .* Sbb ./ _independent_averages(taper, np, np - no, nseg))
    psd = za === zb ? copy(Saa) : real.(Sab)       # one channel passed twice
    return CrossSpectralResult(type, freq, psd, floor, Saa, Sbb, nseg, np, no, taper)
end

"""
    cross_psd_y(a, b; nperseg=nothing, noverlap=nothing,
                taper=:hann) → CrossSpectralResult

Cross power spectral density of the fractional frequency of two
simultaneous records of one device, in 1/Hz.

The estimate is Welch's average of cross-periodograms
([Welch (1967), p. 73](@cite welch-1967-psd)): that of [`psd_y`](@ref) with
the product ``A_k^* B_k`` of the two channels' transforms in place of
``|Y_k|^2``. Its real part is reported
([Baudiquez et al. (2020)](@cite baudiquez-2020-cross-spectrum-statistics)).

Noise present in both channels survives the average. Noise in only one,
such as each reference or each channel's electronics, averages toward zero
and leaves the residual given in the `floor` field of the result. A device
measured against two independent references shows its own spectrum below
the noise of either reference, and one signal split into two channels
shows the noise of the source below the channels' own.

[Baudiquez et al. (2020)](@cite baudiquez-2020-cross-spectrum-statistics)
derive the distribution of the estimate, a variance-gamma distribution, and
from it an upper limit for a spectrum that is not resolved;
[Baudiquez et al. (2022)](@cite baudiquez-2022-cross-spectrum-multiple)
extend it to more than two instruments. Neither is implemented here.

# Arguments

- `a`, `b`: two `PhaseData` or two `FrequencyData` records with the same
  sample interval and length, and no missing samples.

The keywords are those of [`psd_y`](@ref).
"""
function cross_psd_y(a::AbstractTimingData, b::AbstractTimingData;
                     nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    _cross_pair(a, b)
    ya = _f64(_as_freq(a).y); yb = _f64(_as_freq(b).y)
    return _cross_psd(:Sy, ya, yb, a.tau0, nperseg, noverlap, taper)
end

"""
    cross_psd_x(a, b; nperseg=nothing, noverlap=nothing,
                taper=:hann) → CrossSpectralResult

Cross power spectral density of the time fluctuations of two simultaneous
records of one device, in s²/Hz: the estimate of [`cross_psd_y`](@ref), its
floor and the spectrum of each channel, each divided by the response of the
first difference as in [`psd_x`](@ref). Every field is `NaN` at ``f = 0``.

Arguments and keywords are those of [`cross_psd_y`](@ref), with `nperseg`
counted and its default chosen as in [`psd_x`](@ref).
"""
function cross_psd_x(a::AbstractTimingData, b::AbstractTimingData;
                     nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    np, no = _sx_welch_params(a, nperseg, noverlap)
    c = cross_psd_y(a, b; nperseg=np, noverlap=no, taper=taper)
    r = _difference_response(c.freq, a.tau0)
    return CrossSpectralResult(:Sx, c.freq, c.psd ./ r, c.floor ./ r, c.psd_a ./ r, c.psd_b ./ r,
                               c.nseg, np, no, taper)
end

# Equivalent number of independent averages of nseg segments of length L,
# `step` apart (Welch (1967)): nseg / (1 + 2 Σₖ (1 − k/nseg) ρ(k)²), with ρ(k)
# the taper's overlap correlation at a shift of k steps.
function _independent_averages(taper::Symbol, L::Int, step::Int, nseg::Int)
    w = _window(taper, L)
    U = sum(abs2, w)
    s = 0.0
    for k in 1:nseg-1
        k * step >= L && break
        ρ = sum(w[i] * w[i+k*step] for i in 1:L-k*step) / U
        s += (1 - k / nseg) * ρ^2
    end
    return nseg / (1 + 2s)
end

"""
    CrossStabilityResult

Result of [`cross_adev`](@ref).

# Fields

$(TYPEDFIELDS)
"""
struct CrossStabilityResult
    "Averaging times τ, in seconds."
    tau::Vector{Float64}
    "Cross Allan variance ``σ_{ab}^2(τ)``; negative where the common noise is below the scatter of the estimate."
    covar::Vector{Float64}
    "`√covar` where `covar > 0`, `NaN` elsewhere: the deviation common to both channels."
    dev::Vector{Float64}
    """Standard uncertainty of `covar`,
    ``\\sqrt{(σ_a^2 σ_b^2 + σ_{ab}^4)/ν}``: the variance of a mean of products
    of Gaussian second differences, as in the covariance model of
    [Vernotte and Lantz (2019)](@cite vernotte-2019-tch-gcov-uncertainty).
    ``ν`` is the geometric mean of the two channels' Allan-variance EDF for
    their identified noise types, a choice of this package. Where the
    channels share nothing this is the floor of the estimate; `covar` well
    above it is resolved."""
    uncertainty::Vector{Float64}
    "Allan deviation of channel a alone."
    dev_a::Vector{Float64}
    "Allan deviation of channel b alone."
    dev_b::Vector{Float64}
    "Number of overlapping second differences averaged at each τ."
    neff::Vector{Int}
    "Sample interval τ₀, in seconds."
    tau0::Float64
    "Number of phase samples per channel."
    n::Int
end

function Base.show(io::IO, r::CrossStabilityResult)
    print(io, "CrossStabilityResult(", length(r.tau), " τ, ", r.n, " samples per channel)")
end

"""
    cross_adev(a, b, m_values=Octave) → CrossStabilityResult

Cross Allan variance of two simultaneous records of one device: the
overlapping Allan variance with the square of each second difference
replaced by the product of the two channels' second differences,

```math
σ_{ab}^2(τ) = \\frac{1}{2τ^2 (N-2m)} \\sum_{i=1}^{N-2m}
\\left(a_{i+2m} - 2a_{i+m} + a_i\\right) \\left(b_{i+2m} - 2b_{i+m} + b_i\\right),
\\qquad τ = mτ_0 ,
```

with ``a_i`` and ``b_i`` the phase samples of the two channels. This is the
two-sample covariance that
[Vernotte et al. (2016), Eq. (7)](@cite vernotte-2016-gcov) call the
Groslambert covariance.

It estimates the Allan variance of the noise the channels share. Noise in
only one channel averages toward zero, to within the `uncertainty` given
at each ``τ``. With `b` equal to `a` the result is the overlapping Allan
variance. The estimate can be negative: `covar` keeps the sign, and `dev`
is `NaN` there.

# Arguments

- `a`, `b`: two `PhaseData` or two `FrequencyData` records with the same
  sample interval and length, and no missing samples.
- `m_values`: averaging factors, or a [`TauMode`](@ref) grid.
"""
function cross_adev(a::AbstractTimingData, b::AbstractTimingData, m_values::Vector{Int})
    _cross_pair(a, b)
    xa = _f64(_as_phase(a).x); xb = _f64(_as_phase(b).x)
    N = length(xa); τ₀ = Float64(a.tau0)
    all(m -> 1 <= m && 2m < N, m_values) ||
        throw(ArgumentError("cross_adev: every m must satisfy 1 ≤ m < N/2 (N = $N)"))
    K = length(m_values)
    covar = zeros(K); va = zeros(K); vb = zeros(K); neff = zeros(Int, K)
    for (j, m) in enumerate(m_values)
        n = N - 2m
        p, saa, sbb = _second_diff_products(xa, xb, m)
        c = 2 * (m * τ₀)^2
        covar[j] = sum(p) / (c * n); va[j] = saa / (c * n); vb[j] = sbb / (c * n)
        neff[j] = n
    end
    # Each channel's Allan-variance EDF for its identified noise type.
    ν = sqrt.(adev(PhaseData(xa, τ₀), m_values).edf .* adev(PhaseData(xb, τ₀), m_values).edf)
    unc = sqrt.((va .* vb .+ covar .^ 2) ./ ν)
    dev = [v > 0 ? sqrt(v) : NaN for v in covar]
    return CrossStabilityResult(m_values .* τ₀, covar, dev, unc, sqrt.(va), sqrt.(vb), neff, τ₀, N)
end

# Products Δ²x_a(i) Δ²x_b(i) of the overlapping second differences at lag m,
# and the two channels' sums of squares.
function _second_diff_products(xa::Vector{Float64}, xb::Vector{Float64}, m::Int)
    n = length(xa) - 2m
    p = Vector{Float64}(undef, n)
    saa = 0.0; sbb = 0.0
    @inbounds for i in 1:n
        da = xa[i+2m] - 2xa[i+m] + xa[i]
        db = xb[i+2m] - 2xb[i+m] + xb[i]
        p[i] = da * db; saa += da * da; sbb += db * db
    end
    return p, saa, sbb
end

cross_adev(a::AbstractTimingData, b::AbstractTimingData, taus::TauMode) =
    cross_adev(a, b, tau_values(taus, length(_as_phase(a).x), :adev))
cross_adev(a::AbstractTimingData, b::AbstractTimingData) =
    cross_adev(a, b, _default_m_values(length(_as_phase(a).x), :adev))

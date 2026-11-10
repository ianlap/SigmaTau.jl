# Dual-channel (cross) estimators: two simultaneous measurements of the same
# device, each against its own reference or through its own channel. Their
# cross-spectrum and cross Allan variance keep what the channels share (the
# device) and average away what they do not (each reference and each
# channel's own noise), down to a floor set by the number of averages.

"""
    CrossSpectralResult

Return type of [`cross_psd_y`](@ref) and [`cross_psd_x`](@ref).

$(TYPEDFIELDS)
"""
struct CrossSpectralResult
    "`:Sy` (fractional frequency, 1/Hz) or `:Sx` (phase, s²/Hz)."
    spectral_type::Symbol
    "One-sided frequency grid in Hz."
    freq::Vector{Float64}
    "Real part of the averaged cross-spectrum: the noise common to both channels."
    psd::Vector{Float64}
    """Floor of `psd`, `√(S_aa S_bb / K)`: the residual left by noise the
    channels do not share after `K` independent averages, `K` being `nseg`
    reduced for the correlation of overlapping segments (Welch 1967; ≈ `nseg`
    at the default 50 % Hann overlap). Bins where `|psd|` is not well above
    it are not resolved."""
    floor::Vector{Float64}
    "Spectrum of channel a alone."
    psd_a::Vector{Float64}
    "Spectrum of channel b alone."
    psd_b::Vector{Float64}
    "Number of Welch segments averaged."
    nseg::Int
    "Welch segment length (samples)."
    nperseg::Int
    "Welch segment overlap (samples)."
    noverlap::Int
    "Taper applied to each segment."
    taper::Symbol
end

function Base.show(io::IO, r::CrossSpectralResult)
    print(io, "CrossSpectralResult(:", r.spectral_type, ", ", length(r.freq), " bins, ",
              r.nseg, " averages, nperseg=", r.nperseg, ")")
end

function _cross_pair(a::AbstractTimingData, b::AbstractTimingData)
    (a isa PhaseData) == (b isa PhaseData) ||
        throw(ArgumentError("cross: pass two PhaseData or two FrequencyData records"))
    (any(isnan, _samples(a)) || any(isnan, _samples(b))) &&
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
    cross_psd_y(a, b; nperseg=nothing, noverlap=nothing, taper=:hann) → CrossSpectralResult

Cross power spectral density of the fractional frequency of two simultaneous
records of the same device, `a` and `b` (`PhaseData` or `FrequencyData`,
same τ₀ and length, both the same kind). Welch's method as in [`psd_y`](@ref), with the averaged
cross-periodogram `conj(A)·B` in place of `|A|²`.

Noise present in both channels survives the average; noise in only one
(each reference, each channel's electronics) averages toward zero, leaving a
residual of about `floor = √(S_aa S_bb / K)`, `K` the number of independent
averages; the full distribution of the estimate (a variance-gamma for two
channels) is derived in [baudiquez-2020-cross-spectrum-statistics](@cite) and,
for more instruments, [baudiquez-2022-cross-spectrum-multiple](@cite).
Measuring a device against
two independent references gives its own spectrum below either reference's
noise; splitting one signal into two channels gives the source's noise below
the channels' floor.
"""
function cross_psd_y(a::AbstractTimingData, b::AbstractTimingData;
                     nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    _cross_pair(a, b)
    ya = _f64(_as_freq(a).y); yb = _f64(_as_freq(b).y)
    return _cross_psd(:Sy, ya, yb, a.tau0, nperseg, noverlap, taper)
end

"""
    cross_psd_x(a, b; nperseg=nothing, noverlap=nothing, taper=:hann) → CrossSpectralResult

Cross power spectral density of the phase of two simultaneous records, in
s²/Hz: [`cross_psd_y`](@ref), its floor and the two channels' own spectra,
each divided by the response of the first difference as in [`psd_x`](@ref).
"""
function cross_psd_x(a::AbstractTimingData, b::AbstractTimingData;
                     nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    np, no = _sx_welch_params(a, nperseg, noverlap)
    c = cross_psd_y(a, b; nperseg=np, noverlap=no, taper=taper)
    r = _difference_response(c.freq, a.tau0)
    return CrossSpectralResult(:Sx, c.freq, c.psd ./ r, c.floor ./ r, c.psd_a ./ r, c.psd_b ./ r,
                               c.nseg, np, no, taper)
end

# Welch's equivalent number of independent averages of nseg segments of
# length L, `step` apart: nseg / (1 + 2 Σₖ (1 − k/nseg) ρ(k)²), with ρ(k) the
# window's overlap correlation at a shift of k steps.
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

Return type of [`cross_adev`](@ref).

$(TYPEDFIELDS)
"""
struct CrossStabilityResult
    "Analysis intervals τ in seconds."
    tau::Vector{Float64}
    "Cross Allan variance σ²_ab(τ); negative where the common noise is below the floor."
    covar::Vector{Float64}
    "`√covar` where `covar > 0`, else `NaN`: the deviation common to both channels."
    dev::Vector{Float64}
    """1σ statistical uncertainty of `covar`: `√((σ²_a σ²_b + covar²)/ν)`, the
    variance of a cross variance of Gaussian records with `ν` the geometric
    mean of the two channels' Allan-variance EDF for their identified noise
    types [vernotte-2019-tch-gcov-uncertainty](@cite). Where the channels share
    nothing this is the floor; `covar` well above it is resolved."""
    uncertainty::Vector{Float64}
    "Allan deviation of channel a alone."
    dev_a::Vector{Float64}
    "Allan deviation of channel b alone."
    dev_b::Vector{Float64}
    "Number of overlapping second differences averaged at each τ."
    neff::Vector{Int}
    "Sample interval τ₀ in seconds."
    tau0::Float64
    "Number of phase samples per channel."
    n::Int
end

function Base.show(io::IO, r::CrossStabilityResult)
    print(io, "CrossStabilityResult(", length(r.tau), " τ, ", r.n, " samples per channel)")
end

"""
    cross_adev(a, b, m_values=Octave) → CrossStabilityResult

Cross Allan variance of two simultaneous records of the same device: the
overlapping Allan variance with the square of each second difference
replaced by the product of the two channels' second differences,

    σ²_ab(τ) = Σᵢ Δ²x_a(i) Δ²x_b(i) / (2 τ² (N − 2m)),   τ = m τ₀.

It estimates the Allan variance of the noise the channels share; noise in
only one channel averages toward zero, to within the `uncertainty`
reported per τ.
With `b = a` it is exactly the overlapping Allan variance. Because the
estimate can go negative, `covar` keeps the sign and `dev` is `NaN` there.
`a`, `b` are two `PhaseData` or two `FrequencyData` with the same τ₀ and length;
`m_values` is a vector of averaging factors or a [`TauMode`](@ref) grid.
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

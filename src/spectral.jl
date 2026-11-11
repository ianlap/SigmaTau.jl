# Welch spectral estimates: the core, and the public `psd_y`, `psd_x`,
# `ssb_phase_noise` and `integrated_jitter`. The scaling is the one-sided
# density of `scipy.signal.welch` with `scaling="density"`, so the density
# summed over [0, fs/2] times the bin width approximates the variance of the
# input. Each segment loses its mean weighted by the taper, where scipy's
# `detrend="constant"` removes the plain mean: with the Hann and Hamming
# tapers the two estimates differ at f = 0 and in the first bin above it.

"""
    _window(kind::Symbol, n::Int) → Vector{Float64}

Taper of `n` points: `:hann`, `:hamming` or `:rectangular`. The cosine
tapers are in the periodic form, with `n` and not `n − 1` in the cosine's
denominator, as `scipy.signal.get_window` returns them by default.
"""
function _window(kind::Symbol, n::Int)
    if kind === :rectangular
        return ones(Float64, n)
    elseif kind === :hann
        return [0.5 - 0.5 * cospi(2 * k / n) for k in 0:n-1]
    elseif kind === :hamming
        return [0.54 - 0.46 * cospi(2 * k / n) for k in 0:n-1]
    else
        throw(ArgumentError("spectral: unknown taper $kind; " *
                            "expected :hann, :hamming, or :rectangular"))
    end
end

"""
    _welch_psd(z, fs; nperseg, noverlap, window) → (freq, psd)

One-sided Welch estimate of the power spectral density of the real series
`z` sampled at `fs` Hz, as stated in `psd_y`.

`z` is cut into segments of `nperseg` samples, `nperseg − noverlap` apart.
Each has its mean removed and is multiplied by the taper; its squared
transform is scaled by `1/(fs·Σwᵢ²)` and doubled at every bin but zero
frequency and, for even `nperseg`, the Nyquist frequency. The segments are
averaged.

The mean removed is the one weighted by the taper, `Σwᵢzᵢ/Σwᵢ`, which
leaves the tapered segment no component at zero frequency. With the Hann
and Hamming tapers the plain mean, scipy's `detrend="constant"`, gives
another value there and in the first bin, and the same value in every other
bin. The plain mean of a differenced series is the difference of its end
points over the segment length, and the taper carries it into the first
bin: white PM then reads `4/3 + nperseg/(12π²)` times its density there
under the Hann taper, against 13/12 times with the weighted mean.

Returns `freq = [0, fs/nperseg, …, fs/2]`, `nperseg ÷ 2 + 1` values, and the
density at each. Throws an `ArgumentError` for `nperseg < 2`,
`nperseg > length(z)` or `noverlap` outside `0:nperseg−1`.
"""
function _welch_psd(z::Vector{Float64}, fs::Float64;
                    nperseg::Int, noverlap::Int, window::Symbol)
    freq, psd, _, _, _ = _welch(z, z, fs; nperseg=nperseg, noverlap=noverlap, window=window)
    return freq, psd
end

"""
    _welch(a, b, fs; nperseg, noverlap, window) → (freq, Saa, Sbb, Sab, nseg)

Welch spectra and cross-spectrum of two real series of equal length, with
the segments, taper and scaling of `_welch_psd`: `Sab[k]` averages
`conj(A_k)·B_k` over the `nseg` segments. With `b === a` one transform per
segment is taken, `Sbb` is `Saa` and `Sab` is empty.
"""
function _welch(a::Vector{Float64}, b::Vector{Float64}, fs::Float64;
                nperseg::Int, noverlap::Int, window::Symbol)
    same = a === b
    N = length(a)
    length(b) == N ||
        throw(ArgumentError("spectral: the two signals differ in length ($N vs $(length(b)))"))
    nperseg >= 2 ||
        throw(ArgumentError("spectral: nperseg must be ≥ 2 (got $nperseg)"))
    nperseg <= N ||
        throw(ArgumentError("spectral: nperseg=$nperseg exceeds signal length $N"))
    0 <= noverlap < nperseg ||
        throw(ArgumentError("spectral: noverlap must be in [0, nperseg) " *
                            "(got $noverlap, nperseg=$nperseg)"))
    fs > 0 || throw(ArgumentError("spectral: fs must be > 0 (got $fs)"))

    win   = _window(window, nperseg)
    U     = sum(abs2, win)                  # window power
    step  = nperseg - noverlap
    nseg  = (N - noverlap) ÷ step           # number of full segments
    nseg >= 1 ||
        throw(ArgumentError("spectral: no full segment fits; " *
                            "reduce nperseg or noverlap"))

    nfreq = nperseg ÷ 2 + 1
    Saa   = zeros(Float64, nfreq)
    Sbb   = same ? Saa : zeros(Float64, nfreq)
    Sab   = zeros(ComplexF64, same ? 0 : nfreq)
    seg   = Vector{Float64}(undef, nperseg)

    wsum = sum(win)
    function taper!(z, off)
        μ = 0.0
        @inbounds for i in 1:nperseg
            μ += win[i] * z[off + i]
        end
        μ /= wsum
        @inbounds for i in 1:nperseg
            seg[i] = (z[off + i] - μ) * win[i]
        end
        return rfft(seg)
    end

    for s in 0:nseg-1
        A = taper!(a, s * step)
        if same
            @inbounds for k in 1:nfreq
                Saa[k] += abs2(A[k])
            end
        else
            B = taper!(b, s * step)
            @inbounds for k in 1:nfreq
                Saa[k] += abs2(A[k])
                Sbb[k] += abs2(B[k])
                Sab[k] += conj(A[k]) * B[k]
            end
        end
    end

    # Average, scale to a density, and fold to one-sided: double all but DC
    # and (for even nperseg) Nyquist.
    last_doubled = iseven(nperseg) ? nfreq - 1 : nfreq
    @inbounds for k in 1:nfreq
        c = (2 <= k <= last_doubled ? 2.0 : 1.0) / (fs * U * nseg)
        Saa[k] *= c
        same && continue
        Sbb[k] *= c
        Sab[k] *= c
    end

    freq = [(k - 1) * fs / nperseg for k in 1:nfreq]
    return freq, Saa, Sbb, Sab, nseg
end

"""
    _default_nperseg(N::Int) → Int

Default Welch segment length: the largest power of two up to `N ÷ 4`, at
least 16 and at most `N`. From 64 samples on, seven or more half-overlapped
segments are then averaged.
"""
function _default_nperseg(N::Int)
    np = 16
    while 2np <= N ÷ 4
        np *= 2
    end
    return min(np, N)
end

# Public estimators: each wraps the Welch core in a `SpectralResult`. The
# quantities are defined in IEEE Std 1139-2022, Annex A.2.

# Segment length and overlap, with the defaults filled in: `_default_nperseg`
# and half the segment.
function _welch_params(n::Int, nperseg, noverlap)
    np = nperseg === nothing ? _default_nperseg(n) : Int(nperseg)
    no = noverlap === nothing ? np ÷ 2 : Int(noverlap)
    return np, no
end

"""
    psd_y(data; nperseg=nothing, noverlap=nothing, taper=:hann) → SpectralResult

One-sided power spectral density of the fractional frequency, ``S_y(f)`` in
1/Hz ([IEEE Std 1139-2022, Eq. (A.8)](@cite ieee1139-2022-definitions)),
estimated by Welch's method ([Welch (1967)](@cite welch-1967-psd)).

The record is cut into ``K`` overlapping segments of ``L`` samples. Each
segment has its mean removed and is multiplied by the taper ``w``, and the
squared magnitudes of the segments' discrete Fourier transforms ``Y^{(j)}``
are averaged:

```math
\\hat S_y(f_k) = \\frac{2τ_0}{K \\sum_{i=1}^{L} w_i^2} \\sum_{j=1}^{K} \\left| Y^{(j)}_k \\right|^2,
\\qquad f_k = \\frac{k}{Lτ_0}, \\quad k = 0, \\ldots, \\lfloor L/2 \\rfloor .
```

The factor 2 folds the negative frequencies into the one-sided density. It
is left out at ``f = 0`` and, for even ``L``, at the Nyquist frequency
``1/(2τ_0)``. With this scaling the density summed over the bins, times the
bin width ``1/(Lτ_0)``, approximates the variance of ``y``.

The mean removed from a segment is the one weighted by the taper,
``\\sum_i w_i y_i / \\sum_i w_i``. The estimate is then zero at ``f = 0``,
and a frequency offset leaves no trace in it. Segments and scaling are
those of `scipy.signal.welch` with `scaling="density"`, whose
`detrend="constant"` removes the plain mean instead. With the Hann and
Hamming tapers the two estimates differ at ``f = 0`` and in the first bin,
``1/(Lτ_0)``, and are equal from the second bin on, where the validation
suite compares them.

# Arguments

- `data`: a `FrequencyData` record, or a `PhaseData` record, which is
  differenced to frequency first. A missing sample makes the whole spectrum
  `NaN`: fill gaps first with [`fillgaps`](@ref).

# Keywords

- `nperseg`: segment length ``L``, in samples. The default is the largest
  power of two up to a quarter of the record, at least 16 and no longer
  than the record. A longer segment resolves lower frequencies; a shorter
  one averages more segments and scatters less.
- `noverlap`: samples shared by successive segments; half a segment by
  default.
- `taper=:hann`: `:hann`, `:hamming` or `:rectangular`.

!!! note "Lowest bins, tapers and drift"
    The taper spreads each frequency over its neighbors, and the first bin
    shares the main lobe of the taper with the mean that was removed. In
    simulation with the Hann taper it reads between 0.56 times the spectrum,
    on flicker FM, and 1.08 times, on white PM, at any segment length, and
    should not be read. From the second bin on, the expected value of the
    Hann estimate is the spectrum itself on flicker PM and white FM. It is
    8 % high at the second bin and 4 % at the third on white PM, 11 % and
    4 % on flicker FM, and on random-walk FM 48 %, 15 % and 7 % high at the
    second to fourth bins and within 2 % from the eighth.

    The side lobes of the Hamming and rectangular tapers fall slowly, and
    they carry the top of a PM spectrum down into the lowest bins, the more
    so the longer the segment. With ``L = 1024`` on white PM the Hamming
    estimate is 28 % high at the second bin and 7 % at the fourth, and the
    rectangular one is high by 20 % or more up to the 16th bin. The
    rectangular taper is also a factor of 2 high at every bin on
    random-walk FM.

    A linear frequency drift is not removed, and it raises the lowest
    bins. Take it out first with [`detrend`](@ref): `method=:quadratic` for
    a phase record, `:linear` for a frequency record.

# Examples

A sinusoid of unit amplitude has the variance 1/2:

```jldoctest
julia> using SigmaTau

julia> y = FrequencyData([sinpi(k / 4) for k in 0:1023], 1.0);

julia> s = psd_y(y);

julia> (s.nperseg, s.noverlap, length(s.freq))
(256, 128, 129)

julia> sum(s.psd) * (s.freq[2] - s.freq[1]) ≈ 0.5
true

julia> s.freq[argmax(s.psd)]
0.125
```
"""
function psd_y(data::FrequencyData; nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    y  = _f64(data.y)
    np, no = _welch_params(length(y), nperseg, noverlap)
    freq, psd = _welch_psd(y, 1 / data.tau0; nperseg=np, noverlap=no, window=taper)
    return SpectralResult(:Sy, freq, psd, :per_Hz, np, no, taper)
end

psd_y(data::PhaseData; kwargs...) = psd_y(FrequencyData(data); kwargs...)

# Response of the first difference y[k] = (x[k+1] − x[k])/τ₀ at each frequency,
# |1 − e^{−2πifτ₀}|²/τ₀² = 4 sin²(πfτ₀)/τ₀², which is (2πf)² at low f. S_x is
# S_y over it, and has no value at f = 0.
_difference_response(freq::Vector{Float64}, tau0::Float64) =
    [f > 0 ? 4 * sinpi(f * tau0)^2 / tau0^2 : NaN for f in freq]

# Welch parameters of a phase spectrum: the default segment length comes from
# the number of phase samples, and fits the one-fewer frequency samples.
function _sx_welch_params(data::AbstractTimingData, nperseg, noverlap)
    n = length(data) + (data isa FrequencyData)
    np = nperseg === nothing ? min(_default_nperseg(n), n - 1) : Int(nperseg)
    return np, noverlap === nothing ? np ÷ 2 : Int(noverlap)
end

"""
    psd_x(data; nperseg=nothing, noverlap=nothing, taper=:hann) → SpectralResult

One-sided power spectral density of the time fluctuations, ``S_x(f)`` in
s²/Hz ([IEEE Std 1139-2022, Eq. (A.14)](@cite ieee1139-2022-definitions)):
the estimate of ``S_y(f)`` by [`psd_y`](@ref), divided by the response of
the first difference that takes phase to frequency.

```math
\\hat S_x(f_k) = \\frac{τ_0^2}{4 \\sin^2(π f_k τ_0)}\\, \\hat S_y(f_k),
\\qquad k = 1, \\ldots, \\lfloor L/2 \\rfloor .
```

The two densities are related by ``S_y(f) = (2πf)^2 S_x(f)``
([IEEE Std 1139-2022, Eq. (A.16)](@cite ieee1139-2022-definitions)), and
the response of the first difference over ``τ_0`` is ``(2πf)^2`` well below
the Nyquist frequency. The value at ``f = 0`` is `NaN`.

The spectrum is not estimated on the phase itself, because the phase
spectrum of FM noise is steep, ``f^{-2}`` for white FM to ``f^{-4}`` for
random-walk FM, and what a taper leaks from the lowest frequencies then
exceeds the spectrum. In simulation a Welch estimate on the phase, with the
Hann taper, reads 3.5 times high at the second bin on flicker FM and more
than 100 times on random-walk FM, and a frequency offset swamps every low
bin. The spectrum of the frequency is less steep by ``f^2``, and a
frequency offset is only its mean. At each bin the estimate is high or low
by the factor that [`psd_y`](@ref) is, and the note there on the lowest
bins, the tapers and drift applies.

# Arguments

- `data`: a `PhaseData` or `FrequencyData` record. A missing sample makes
  the whole spectrum `NaN`, as in `psd_y`.

# Keywords

Those of [`psd_y`](@ref), with one difference. A segment of `nperseg`
frequency samples spans `nperseg + 1` phase samples, and the default is the
largest power of two up to a quarter of the number of phase samples, at
least 16. `psd_y` counts the frequency samples, one fewer, so for a phase
record of ``2^j ≥ 128`` samples its default segment is half as long.
"""
function psd_x(data::AbstractTimingData; nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    np, no = _sx_welch_params(data, nperseg, noverlap)
    sy = psd_y(data; nperseg=np, noverlap=no, taper=taper)
    return SpectralResult(:Sx, sy.freq, sy.psd ./ _difference_response(sy.freq, data.tau0),
                          :s2_per_Hz, np, no, taper)
end

"""
    ssb_phase_noise(data; f_carrier, nperseg=nothing, noverlap=nothing,
                    taper=:hann) → SpectralResult

Phase noise ``ℒ(f)`` in dBc/Hz of a signal at the carrier frequency
`f_carrier`, in Hz, from its phase or frequency record.

The time spectrum of [`psd_x`](@ref) gives the phase spectrum in rad²/Hz,
``S_φ(f) = (2πν_0)^2 S_x(f)`` with ``ν_0`` the carrier frequency, and the
phase noise is half of it, ``ℒ(f) = S_φ(f)/2``, reported as
``10 \\log_{10} ℒ(f)``
([IEEE Std 1139-2022, Eqs. (A.13) and (A.15)](@cite ieee1139-2022-definitions)).
So defined, ``ℒ(f)`` is the single-sideband noise power relative to the
carrier only for small phase modulation (clause 4 of the standard).

The result runs over the offset frequencies from ``1/(Lτ_0)`` to the
Nyquist frequency ``1/(2τ_0)``; the bin at ``f = 0`` is dropped. The other
keywords are those of [`psd_x`](@ref), and the note of [`psd_y`](@ref) on
the lowest bins, the tapers and drift applies: the first point of the
result is its first bin.
"""
function ssb_phase_noise(data::AbstractTimingData; f_carrier::Real, nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    f_carrier > 0 ||
        throw(ArgumentError("ssb_phase_noise: f_carrier must be > 0 Hz (got $f_carrier)"))
    sx = psd_x(data; nperseg=nperseg, noverlap=noverlap, taper=taper)
    # Drop the bin at f = 0; S_x → S_φ → ℒ(f) in dBc/Hz.
    freq = sx.freq[2:end]
    sphi = (2π * f_carrier)^2 .* sx.psd[2:end]
    script_l = 10 .* log10.(0.5 .* sphi)
    return SpectralResult(:L, freq, script_l, :dBc_per_Hz, sx.nperseg, sx.noverlap, taper)
end

"""
    integrated_jitter(s::SpectralResult; f_carrier, band=nothing) → NamedTuple
    integrated_jitter(f, S; quantity=:L, f_carrier, band=nothing) → NamedTuple

Phase jitter, time jitter and residual FM of a signal at the carrier
frequency `f_carrier`, in Hz: its phase spectrum ``S_φ(f)`` integrated over
the band ``f_1`` to ``f_2`` of offset frequencies. The result has the fields

- `phase`: ``Δφ = \\left(\\int_{f_1}^{f_2} S_φ(f)\\,df\\right)^{1/2}``, in
  radians rms
  ([IEEE Std 1139-2022, Eq. (A.17)](@cite ieee1139-2022-definitions)).
- `time`: ``Δφ/(2πν_0)``, in seconds rms, the same integral of ``S_x(f)``.
- `residual_fm`: ``Δf = \\left(\\int_{f_1}^{f_2} f^2 S_φ(f)\\,df\\right)^{1/2}``,
  in hertz rms.

Between neighboring points the spectrum is taken as a power law,
``S_φ(f) = K f^x``, a straight line on a log-log plot, and each segment is
integrated in closed form
([SP1065 §5.2.19, Eqs. (38) to (43)](@cite riley-2008-sp1065)).

# Arguments

- `s`: a result of [`ssb_phase_noise`](@ref), [`psd_x`](@ref) or
  [`psd_y`](@ref).
- `f`, `S`: spectrum points, such as the ``ℒ(f)`` table of a phase-noise
  analyzer or a data sheet. `f` is in Hz and increases.

# Keywords

- `f_carrier`: the carrier frequency ``ν_0``, in Hz.
- `band=nothing`: the limits `(f₁, f₂)`, in Hz, inside the range of the
  spectrum. The whole spectrum is integrated by default.
- `quantity=:L`: what `S` holds: `:L`, ``ℒ(f)`` in dBc/Hz, with
  ``S_φ = 2ℒ``; `:Sx`, ``S_x(f)`` in s²/Hz; or `:Sy`, ``S_y(f)`` in 1/Hz,
  converted with ``(2πf)^2``. A [`psd_y`](@ref) result is converted with
  the response of the first difference instead (see [`psd_x`](@ref)).

Points at ``f = 0`` and points without a positive spectrum are left out.

This is jitter as the integral of the phase noise over a band.
[Howe and Tasset (2003)](@cite howe-2003-clock-jitter) treat jitter as a
function of the delay ``τ``, the rms of ``x(t+τ) - x(t)``, which
[`tierms`](@ref) estimates from a phase record.

# Examples

A phase noise that falls by 10 dB per decade, ``x = -1``, on a 10 MHz
carrier; the results are in rad, ps and Hz:

```jldoctest
julia> using SigmaTau

julia> j = integrated_jitter([1e3, 1e4, 1e5], [-100.0, -110.0, -120.0];
                             f_carrier=10e6);

julia> r4(v) = round(v; sigdigits=4);

julia> r4(j.phase), r4(j.time * 1e12), r4(j.residual_fm)
(0.0009597, 15.27, 31.62)
```
"""
function integrated_jitter(s::SpectralResult; f_carrier::Real, band=nothing)
    τ₀ = 1 / (s.nperseg * (s.freq[2] - s.freq[1]))
    return _jitter(s.freq, s.psd, s.spectral_type, f_carrier, band, τ₀)
end

integrated_jitter(f::AbstractVector{<:Real}, S::AbstractVector{<:Real};
                  quantity::Symbol=:L, f_carrier::Real, band=nothing) =
    _jitter(Float64.(f), Float64.(S), quantity, f_carrier, band, NaN)

function _jitter(f::Vector{Float64}, S::Vector{Float64}, quantity::Symbol,
                 f_carrier::Real, band, τ₀::Float64)
    length(f) == length(S) || throw(DimensionMismatch(
        "integrated_jitter: $(length(f)) frequencies but $(length(S)) values"))
    f_carrier > 0 || throw(ArgumentError("integrated_jitter: f_carrier must be > 0 Hz (got $f_carrier)"))
    ω₀² = (2π * f_carrier)^2
    Sφ = quantity === :L  ? 2 .* exp10.(S ./ 10) :
         quantity === :Sx ? ω₀² .* S :
         quantity === :Sy ? ω₀² .* S ./ _powerlaw_to_sy(:Sx, f, ones(length(f)), nothing, τ₀) :
         throw(ArgumentError("integrated_jitter: quantity must be :L, :Sx or :Sy, got :$quantity"))
    keep = findall(k -> f[k] > 0 && isfinite(Sφ[k]) && Sφ[k] > 0, eachindex(f))
    fk, Sk = f[keep], Sφ[keep]
    length(fk) >= 2 && issorted(fk; lt=<=) ||
        throw(ArgumentError("integrated_jitter: need at least 2 increasing frequencies with a positive spectrum"))
    lo, hi = band === nothing ? (fk[1], fk[end]) : Float64.(band)
    fk[1] <= lo < hi <= fk[end] || throw(ArgumentError(
        "integrated_jitter: band ($lo, $hi) Hz must lie inside the spectrum, $(fk[1]) to $(fk[end]) Hz"))
    φ² = 0.0; ν² = 0.0
    for k in 1:length(fk)-1
        a, b = max(fk[k], lo), min(fk[k+1], hi)
        a < b || continue
        x = log(Sk[k+1] / Sk[k]) / log(fk[k+1] / fk[k])
        φ² += _powerlaw_segment(Sk[k] * fk[k], x, a / fk[k], b / fk[k])
        ν² += _powerlaw_segment(Sk[k] * fk[k]^3, x + 2, a / fk[k], b / fk[k])
    end
    return (phase = √φ², time = √φ² / (2π * f_carrier), residual_fm = √ν²)
end

# ∫ S₁ (f/f₁)^x df from a·f₁ to b·f₁, with c = S₁·f₁ (SP1065 §5.2.19,
# Eqs. (41) and (43)).
_powerlaw_segment(c::Float64, x::Float64, a::Float64, b::Float64) =
    abs(x + 1) < 1e-9 ? c * log(b / a) : c * (b^(x + 1) - a^(x + 1)) / (x + 1)

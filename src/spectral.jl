# Welch PSD core and the public `psd_y` / `psd_x` / `ssb_phase_noise` estimators.
# One-sided "density" normalization (as scipy.signal.welch with
# `scaling="density", return_onesided=True`), so the integral of the PSD over
# [0, fs/2] recovers the variance of the per-segment mean-removed input.

"""
    _window(kind::Symbol, n::Int) → Vector{Float64}

Length-`n` analysis window. Supports `:hann`, `:hamming`, and `:rectangular`
(the no-op boxcar). Periodic ("DFT-even") convention, with `n` rather than
`n-1` in the cosine denominator, which is the correct form for spectral
estimation (avoids the half-sample bias of the symmetric window).
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

One-sided Welch PSD estimate of the real signal `z` sampled at `fs` Hz.

Splits `z` into overlapping segments of length `nperseg` stepped by
`nperseg - noverlap`, removes each segment's mean (constant detrend), applies
the chosen `window`, periodograms each segment with a real FFT, and averages.
The normalization is the "density" convention: each periodogram is scaled by
`1 / (fs · Σ wᵢ²)` and every bin except DC (and the Nyquist bin when `nperseg`
is even) is doubled to fold negative-frequency power into the one-sided
spectrum. With this scaling `Σ psd · Δf ≈ var(z)` over the analysis band.

Returns the frequency grid `freq = [0, fs/nperseg, …, fs/2]` (length
`nperseg ÷ 2 + 1`) and the averaged PSD.

Throws `ArgumentError` on `nperseg < 2`, `nperseg > length(z)`, or
`noverlap ∉ [0, nperseg)`.
"""
function _welch_psd(z::Vector{Float64}, fs::Float64;
                    nperseg::Int, noverlap::Int, window::Symbol)
    freq, psd, _, _, _ = _welch(z, z, fs; nperseg=nperseg, noverlap=noverlap, window=window)
    return freq, psd
end

"""
    _welch(a, b, fs; nperseg, noverlap, window) → (freq, Saa, Sbb, Sab, nseg)

Welch auto- and cross-spectra of two equal-length real signals, with the
segmenting, detrend, window and one-sided density scaling of
[`_welch_psd`](@ref): `Sab[k]` averages `conj(A_k)·B_k` over the `nseg`
segments. When `b === a` only one FFT per segment is taken, `Sbb` is
`Saa` and `Sab` is empty.
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

    function taper!(z, off)
        μ = 0.0
        @inbounds for i in 1:nperseg
            μ += z[off + i]
        end
        μ /= nperseg
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

Default Welch segment length when the caller does not supply one: the largest
power of two up to `N ÷ 4`, at least 16 and at most `N`. The spectrum then
reaches within a factor of four of the record's lowest frequency and still
averages seven half-overlapped segments. Tune upward for finer low-frequency
resolution, downward for more segments (lower variance).
"""
function _default_nperseg(N::Int)
    np = 16
    while 2np <= N ÷ 4
        np *= 2
    end
    return min(np, N)
end

# Public API: `psd_y` / `psd_x` / `ssb_phase_noise` take `PhaseData` / `FrequencyData` and return a
# `SpectralResult` wrapping the `_welch_psd` core.
# Reference: IEEE Std 1139-2022 §3.3–3.5; NIST SP1065 (Riley & Howe 2008).

# Resolve Welch parameters from explicit kwargs, filling defaults from the
# sequence length. `noverlap` defaults to 50 % of the (resolved) segment.
function _welch_params(n::Int, nperseg, noverlap)
    np = nperseg === nothing ? _default_nperseg(n) : Int(nperseg)
    no = noverlap === nothing ? np ÷ 2 : Int(noverlap)
    return np, no
end

"""
    psd_y(data::FrequencyData; nperseg=nothing, noverlap=nothing, taper=:hann) → SpectralResult
    psd_y(data::PhaseData;     nperseg=nothing, noverlap=nothing, taper=:hann) → SpectralResult

One-sided power spectral density of the fractional frequency `y(t)`,
`S_y(f)` in units of `1/Hz` (IEEE Std 1139-2022 §3.4), estimated by Welch's
method at sample rate `fs = 1/τ₀`.

`PhaseData` input is differenced to fractional frequency
(`y[k] = (x[k+1] − x[k])/τ₀`) first. `nperseg` defaults to the largest
power of two up to `N ÷ 4` (at least 16), `noverlap` to half of it; raise `nperseg` for finer
low-frequency resolution, lower it for more averaging. `taper` is `:hann`,
`:hamming`, or `:rectangular`.

The estimate is variance-preserving: `Σ S_y(f)·Δf ≈ var(y)` over `[0, fs/2]`.
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
    psd_x(data::PhaseData;     nperseg=nothing, noverlap=nothing, taper=:hann) → SpectralResult
    psd_x(data::FrequencyData; nperseg=nothing, noverlap=nothing, taper=:hann) → SpectralResult

One-sided power spectral density of the phase (time) residual `x(t)`,
`S_x(f)` in units of `s²/Hz` (IEEE Std 1139-2022 §3.3), at sample rate
`fs = 1/τ₀`.

It is the Welch estimate of `S_y(f)` (see [`psd_y`](@ref), same keyword
arguments) divided by the response of the first difference that takes phase to
frequency, `4 sin²(πfτ₀)/τ₀²`, the sampled form of `S_x(f) = S_y(f)/(2πf)²`.
The phase of FM noise has a spectrum as steep as f⁻² to f⁻⁴, which leaks
through any taper into the lowest bins (a Welch estimate on the phase itself
reads 3 times high at the second bin on flicker FM and about 100 times on
random-walk FM, and is swamped at every low bin by a frequency offset); its
first difference is two powers of f flatter and a frequency offset is only its
mean. The value at `f = 0` is `NaN`, and the first bin above it, which shares
the taper's main lobe with the removed mean, is the least reliable: on white
PM it reads about 10 times high. A segment holds `nperseg` frequency samples;
the default is the largest power of two up to a quarter of the phase samples.
"""
function psd_x(data::AbstractTimingData; nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    np, no = _sx_welch_params(data, nperseg, noverlap)
    sy = psd_y(data; nperseg=np, noverlap=no, taper=taper)
    return SpectralResult(:Sx, sy.freq, sy.psd ./ _difference_response(sy.freq, data.tau0),
                          :s2_per_Hz, np, no, taper)
end

"""
    ssb_phase_noise(data::PhaseData;     f_carrier, nperseg=nothing, noverlap=nothing, taper=:hann) → SpectralResult
    ssb_phase_noise(data::FrequencyData; f_carrier, nperseg=nothing, noverlap=nothing, taper=:hann) → SpectralResult

Single-sideband phase noise `ℒ(f)` in `dBc/Hz` (IEEE Std 1139-2022 §3.5).

Converts the seconds-valued phase residual to phase in radians at the carrier
`f_carrier` (Hz, required), `S_φ(f) = (2π·f_carrier)²·S_x(f)`, and reports
`ℒ(f) = (1/2)·S_φ(f)` on a decibel scale, `10·log₁₀(ℒ_linear)`. Useful for
direct comparison against oscillator-vendor phase-noise plots.

The DC bin (`f = 0`) is dropped, since `ℒ(0)` is not meaningful and the dB
conversion would be `-Inf`. Keyword arguments otherwise match [`psd_y`](@ref).
"""
function ssb_phase_noise(data::AbstractTimingData; f_carrier::Real, nperseg=nothing, noverlap=nothing, taper::Symbol=:hann)
    f_carrier > 0 ||
        throw(ArgumentError("ssb_phase_noise: f_carrier must be > 0 Hz (got $f_carrier)"))
    sx = psd_x(data; nperseg=nperseg, noverlap=noverlap, taper=taper)
    # Drop DC; convert S_x → S_φ → ℒ(f) → dBc/Hz.
    freq = sx.freq[2:end]
    sphi = (2π * f_carrier)^2 .* sx.psd[2:end]
    script_l = 10 .* log10.(0.5 .* sphi)
    return SpectralResult(:L, freq, script_l, :dBc_per_Hz, sx.nperseg, sx.noverlap, taper)
end

"""
    integrated_jitter(s::SpectralResult; f_carrier, band=nothing) → (; phase, time, residual_fm)
    integrated_jitter(f, S; quantity=:L, f_carrier, band=nothing) → (; phase, time, residual_fm)

Integrated phase jitter and residual FM over the offset band `band = (f₁, f₂)`
in Hz (default: the whole spectrum), NIST SP1065 §5.2.19
[riley-2008-sp1065](@cite) and IEEE Std 1139-1999 eq. (A.16)
[ieee1139-1999-definitions](@cite). Between neighbouring points the phase
spectrum `S_φ(f)` is taken as the power law `K·f^x` through them, a straight line on a
log-log plot, and integrated exactly:

- `phase`: `Δφ = √∫ S_φ(f) df`, rms radians;
- `time`: `Δφ / (2π·f_carrier)`, rms seconds (`S_x = S_φ/(2πν₀)²`, IEEE 1139
  eq. A.14) [howe-2003-clock-jitter](@cite);
- `residual_fm`: `Δf = √∫ f²·S_φ(f) df`, rms hertz.

`s` is a [`ssb_phase_noise`](@ref), [`psd_x`](@ref) or [`psd_y`](@ref)
result. The second form takes spectrum points, such as a phase-noise
analyzer's ℒ(f) in dBc/Hz (`quantity=:L`, `S_φ = 2·ℒ`, IEEE 1139 eq. A.12),
or `:Sx` (s²/Hz) or `:Sy` (1/Hz); `f` must increase. `f_carrier` is the carrier in Hz.

```jldoctest
julia> f = [1e3, 1e4, 1e5];  # ℒ falling 10 dB a decade: x = -1

julia> j = integrated_jitter(f, [-100.0, -110.0, -120.0]; f_carrier=10e6);

julia> round(j.phase; sigdigits=4), round(j.time * 1e12; sigdigits=4), round(j.residual_fm; sigdigits=4)
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

# ∫ S₁ (f/f₁)^x df from a·f₁ to b·f₁, with c = S₁·f₁ (SP1065 eqs 41 and 43).
_powerlaw_segment(c::Float64, x::Float64, a::Float64, b::Float64) =
    abs(x + 1) < 1e-9 ? c * log(b / a) : c * (b^(x + 1) - a^(x + 1)) / (x + 1)

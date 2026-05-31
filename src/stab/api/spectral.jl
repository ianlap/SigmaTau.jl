# stab/api/spectral.jl — public spectral-density entry points.
#
# `Sy` / `Sx` / `L` take `PhaseData` / `FrequencyData` and return a
# `SpectralResult`, wrapping the `_welch_psd` core in `stab/spectral.jl`.
# Reference: IEEE Std 1139-2022 §3.3–3.5; NIST SP1065 (Riley & Howe 2008).

# Resolve Welch parameters from explicit kwargs, filling defaults from the
# sequence length. `noverlap` defaults to 50 % of the (resolved) segment.
function _welch_params(n::Int, nperseg, noverlap)
    np = nperseg === nothing ? _default_nperseg(n) : Int(nperseg)
    no = noverlap === nothing ? np ÷ 2 : Int(noverlap)
    return np, no
end

"""
    Sy(data::FrequencyData; nperseg=nothing, noverlap=nothing, window=:hann) → SpectralResult
    Sy(data::PhaseData;     nperseg=nothing, noverlap=nothing, window=:hann) → SpectralResult

One-sided power spectral density of the fractional frequency `y(t)`,
`S_y(f)` in units of `1/Hz` (IEEE Std 1139-2022 §3.4), estimated by Welch's
method at sample rate `fs = 1/τ₀`.

`PhaseData` input is differenced to fractional frequency
(`y[k] = (x[k+1] − x[k])/τ₀`) first. `nperseg` / `noverlap` default to the
scipy convention (`min(N, 256)` with 50 % overlap); raise `nperseg` for finer
low-frequency resolution, lower it for more averaging. `window` is `:hann`,
`:hamming`, or `:rectangular`.

The estimate is variance-preserving: `Σ S_y(f)·Δf ≈ var(y)` over `[0, fs/2]`.
"""
function Sy(data::FrequencyData; nperseg=nothing, noverlap=nothing, window::Symbol=:hann)
    y  = _f64(data.y)
    np, no = _welch_params(length(y), nperseg, noverlap)
    freq, psd = _welch_psd(y, 1 / data.tau0; nperseg=np, noverlap=no, window=window)
    return SpectralResult(:Sy, freq, psd, :per_Hz, np, no, window)
end

Sy(data::PhaseData; kwargs...) = Sy(_phase_to_freq(data); kwargs...)

"""
    Sx(data::PhaseData;     nperseg=nothing, noverlap=nothing, window=:hann) → SpectralResult
    Sx(data::FrequencyData; nperseg=nothing, noverlap=nothing, window=:hann) → SpectralResult

One-sided power spectral density of the phase (time) residual `x(t)`,
`S_x(f)` in units of `s²/Hz` (IEEE Std 1139-2022 §3.3), estimated by Welch's
method at sample rate `fs = 1/τ₀`.

`FrequencyData` input is integrated to phase (`x[k] = τ₀·Σⱼ y[j]`) first. The
estimators are related by `S_x(f) = S_y(f) / (2πf)²`; the two agree under that
identity at low frequency, where the discrete integrator matches the continuous
one. Keyword arguments match [`Sy`](@ref).
"""
function Sx(data::PhaseData; nperseg=nothing, noverlap=nothing, window::Symbol=:hann)
    x  = _f64(data.x)
    np, no = _welch_params(length(x), nperseg, noverlap)
    freq, psd = _welch_psd(x, 1 / data.tau0; nperseg=np, noverlap=no, window=window)
    return SpectralResult(:Sx, freq, psd, :s2_per_Hz, np, no, window)
end

Sx(data::FrequencyData; kwargs...) = Sx(_freq_to_phase(data); kwargs...)

"""
    L(data::PhaseData;     f_carrier, nperseg=nothing, noverlap=nothing, window=:hann) → SpectralResult
    L(data::FrequencyData; f_carrier, nperseg=nothing, noverlap=nothing, window=:hann) → SpectralResult

Single-sideband phase noise `ℒ(f)` in `dBc/Hz` (IEEE Std 1139-2022 §3.5).

Converts the seconds-valued phase residual to phase in radians at the carrier
`f_carrier` (Hz, required), `S_φ(f) = (2π·f_carrier)²·S_x(f)`, and reports
`ℒ(f) = (1/2)·S_φ(f)` on a decibel scale, `10·log₁₀(ℒ_linear)`. Useful for
direct comparison against oscillator-vendor phase-noise plots.

The DC bin (`f = 0`) is dropped, since `ℒ(0)` is not meaningful and the dB
conversion would be `-Inf`. Keyword arguments otherwise match [`Sy`](@ref).
"""
function L(data::PhaseData; f_carrier::Real, nperseg=nothing, noverlap=nothing, window::Symbol=:hann)
    f_carrier > 0 ||
        throw(ArgumentError("L: f_carrier must be > 0 Hz (got $f_carrier)"))
    sx = Sx(data; nperseg=nperseg, noverlap=noverlap, window=window)
    # Drop DC; convert S_x → S_φ → ℒ(f) → dBc/Hz.
    freq = sx.freq[2:end]
    sphi = (2π * f_carrier)^2 .* sx.psd[2:end]
    script_l = 10 .* log10.(0.5 .* sphi)
    return SpectralResult(:L, freq, script_l, :dBc_per_Hz, sx.nperseg, sx.noverlap, window)
end

L(data::FrequencyData; kwargs...) = L(_freq_to_phase(data); kwargs...)

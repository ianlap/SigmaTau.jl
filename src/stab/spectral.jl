# stab/spectral.jl — Welch power-spectral-density core.
#
# Pure `Vector{Float64}` → arrays kernel shared by the public `Sy` / `Sx` / `L`
# API in `stab/api/spectral.jl`. One-sided "density" normalization (matching
# scipy.signal.welch with `scaling="density", return_onesided=True`) so that
# the integral of the returned PSD over [0, fs/2] recovers the variance of the
# (per-segment mean-removed) input. Uses the real FFT from FFTW (imported once
# in `SigmaTau.jl`).

"""
    _window(kind::Symbol, n::Int) → Vector{Float64}

Length-`n` analysis window. Supports `:hann`, `:hamming`, and `:rectangular`
(the no-op boxcar). Periodic ("DFT-even") convention — the `n` in the cosine
denominator is `n`, not `n-1` — which is the correct form for spectral
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
        throw(ArgumentError("spectral: unknown window $kind; " *
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
    N = length(z)
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
    scale = 1.0 / (fs * U)
    psd   = zeros(Float64, nfreq)
    seg   = Vector{Float64}(undef, nperseg)

    for s in 0:nseg-1
        off = s * step
        @inbounds begin
            μ = 0.0
            for i in 1:nperseg
                μ += z[off + i]
            end
            μ /= nperseg
            for i in 1:nperseg
                seg[i] = (z[off + i] - μ) * win[i]
            end
        end
        Z = rfft(seg)
        @inbounds for k in 1:nfreq
            psd[k] += abs2(Z[k]) * scale
        end
    end
    psd ./= nseg

    # Fold to one-sided: double all but DC and (for even nperseg) Nyquist.
    last_doubled = iseven(nperseg) ? nfreq - 1 : nfreq
    @inbounds for k in 2:last_doubled
        psd[k] *= 2
    end

    freq = [(k - 1) * fs / nperseg for k in 1:nfreq]
    return freq, psd
end

"""
    _default_nperseg(N::Int) → Int

Default Welch segment length when the caller does not supply one: `min(N, 256)`,
the scipy default. Tune upward for finer low-frequency resolution, downward for
more segments (lower variance).
"""
_default_nperseg(N::Int) = min(N, 256)

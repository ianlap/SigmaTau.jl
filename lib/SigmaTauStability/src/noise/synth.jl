# noise/synth.jl — Power-law phase-noise synthesizer (non-exported helper).
#
# Generates synthetic phase data with a chosen fractional-frequency power-law
# slope α ∈ {-2, -1, 0, 1, 2} via a frequency-domain `f^(α/2)` shaping of
# white Gaussian noise, integrated to phase. Used by the test suite to
# exercise stability deviations across all five SP1065 noise types.
#
# Requires an `AbstractFFTs` backend (e.g. `FFTW`) to be loaded by the caller.
# `SigmaTauStability` declares `AbstractFFTs` in its [deps] but no concrete
# backend; pulling in `using FFTW` at the call site supplies the methods.

using AbstractFFTs

"""
    _gen_powerlaw_phase(alpha, N; tau0=1.0) → Vector{Float64}

Synthesize an N-sample phase residual vector whose fractional-frequency
sequence has power spectral density ∝ `f^alpha`. Mapping:

| α  | Noise type | Phase spectrum |
|---:|:-----------|:---------------|
|  2 | WPM        | flat           |
|  1 | FLPM       | 1/f            |
|  0 | WHFM       | 1/f²           |
| -1 | FLFM       | 1/f³           |
| -2 | RWFM       | 1/f⁴           |

For deterministic output, seed the global RNG (`Random.seed!`) before calling.
The DC component is zeroed before the inverse transform so the output has
zero mean. Calls into `AbstractFFTs.fft`/`ifft` — caller must have an FFT
backend loaded (e.g. `using FFTW`).
"""
function _gen_powerlaw_phase(alpha::Real, N::Int; tau0::Real = 1.0)
    # White Gaussian noise → fractional-frequency series shaped to f^alpha.
    w = randn(N)
    W = fft(w)

    # Symmetric FFT-bin frequency magnitude. Bin 0 (DC) gets a placeholder
    # so the f^(alpha/2) term is finite; we zero the DC bin afterwards.
    f = abs.(AbstractFFTs.fftfreq(N, 1.0))
    f[1] = 1.0

    Y    = W .* f .^ (alpha / 2)
    Y[1] = 0.0
    y    = real.(ifft(Y))

    # Phase: x[k] = τ₀ · Σⱼ₌₁ᵏ y[j]
    return cumsum(y) .* tau0
end

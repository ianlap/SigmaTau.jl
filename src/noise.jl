# noise.jl — power-law noise identification and synthesis.
#
# Merged from stab/noise/{lag1, synth, gen}.jl (lag-1 ACF / B1 / R(n) noise-type
# ID, the internal spectral shaper, and the public noise_gen generator).


# ──────────────────────────────────────────────────────────────────────
# ── stab/noise/lag1.jl ──────────────────────────────────────────────────

# noise/lag1.jl — Lag-1 Autocorrelation for noise identification
using Statistics
using StaticArrays

const NEFF_RELIABLE = 30

"""
    identify_noise(x::Vector{Float64}, m_values::Vector{Int}; dmin::Int=0, dmax::Int=2, detrend::Bool=true) → Vector{Symbol}

Identifies dominant power-law noise type using lag-1 autocorrelation for
phase data with B1-ratio fallback.

The full record always passes through a 5σ outlier filter before the
per-m loop; that step is unconditional and unrelated to the polynomial
detrend below.

## `detrend`

When `true` (default), each decimated-by-m subseries is quadratically
detrended before the lag-1 ACF / B1 ratio is computed. Allantools'
`autocorr_noise_id` does the same thing (`detrend(x, deg=2)` in
[`ci.py`](https://github.com/aewallin/allantools/blob/master/allantools/ci.py))
to remove a frequency offset (linear phase ramp) and drift (quadratic
phase term) from the slow content of each averaging window.

When `false`, the per-m polynomial fit is skipped. This matches
Stable32's noise-ID convention; use it when you need α values that line
up with a Stable32-generated reference fixture point-for-point.
"""
function identify_noise(x::Vector{Float64}, m_values::Vector{Int}; dmin::Int=0, dmax::Int=2, detrend::Bool=true)
    x_clean = _preprocess(x)
    N = length(x_clean)
    noises = Vector{Symbol}(undef, length(m_values))
    last_reliable = :unknown

    for (k, m) in enumerate(m_values)
        N_eff = N ÷ m
        alpha = NaN

        try
            if N_eff >= NEFF_RELIABLE
                alpha, _, _, _ = _noise_id_lag1acf(x_clean, m, dmin, dmax; detrend=detrend)
            else
                alpha, _, _ = _noise_id_b1rn(x_clean, m; detrend=detrend)
            end
        catch
            alpha = NaN
        end

        if !isnan(alpha)
            # Map alpha to symbol
            if round(alpha) == 2
                noises[k] = :WHPM
            elseif round(alpha) == 1
                noises[k] = :FLPM
            elseif round(alpha) == 0
                noises[k] = :WHFM
            elseif round(alpha) == -1
                noises[k] = :FLFM
            elseif round(alpha) == -2
                noises[k] = :RWFM
            else
                noises[k] = :unknown
            end
            last_reliable = noises[k]
        else
            noises[k] = last_reliable
        end
    end

    return noises
end

function _preprocess(x::Vector{Float64})
    # 5σ outlier filter only — no polynomial detrending. The pre-record
    # linear detrend that lived here previously was distinct to SigmaTau
    # (Stable32 and allantools both leave the full record untouched and
    # detrend per-m instead). Removing it lets our default α track the
    # external references; the per-m polynomial detrend still runs in
    # `_noise_id_lag1acf` / `_noise_id_b1rn` (toggleable via the
    # `detrend` kwarg on `identify_noise`).
    x_std = std(x)
    x_std < eps() && return x  # degenerate input — nothing to filter against.
    x_mean = mean(x)
    thresh = 5.0 * x_std
    # Two-pass count-then-fill avoids the N-length boolean mask + temporary
    # `z` array the broadcast form allocated.
    nkeep = 0
    @inbounds for v in x
        nkeep += (abs(v - x_mean) < thresh) ? 1 : 0
    end
    nkeep == length(x) && return x
    out = Vector{Float64}(undef, nkeep)
    j = 1
    @inbounds for v in x
        if abs(v - x_mean) < thresh
            out[j] = v
            j += 1
        end
    end
    return out
end

function _detrend_quadratic(x::Vector{Float64})
    N = length(x)
    N_float = Float64(N)

    # One pass for the three weighted-sum moments. Comprehensions in the
    # previous form allocated a generator state per `sum(... for i in 1:N)`.
    X1 = 0.0; X2 = 0.0; X3 = 0.0
    @inbounds @simd for i in 1:N
        v = x[i]; fi = Float64(i)
        X1 += v
        X2 += fi * v
        X3 += fi * fi * v
    end

    S1 = N_float
    S2 = (N_float * (N_float + 1.0)) / 2.0
    S3 = (N_float * (N_float + 1.0) * (2.0*N_float + 1.0)) / 6.0
    S4 = (N_float^2 * (N_float + 1.0)^2) / 4.0
    S5 = (N_float * (N_float + 1.0) * (2.0*N_float + 1.0) * (3.0*N_float^2 + 3.0*N_float - 1.0)) / 30.0

    # 3×3 LS solve via StaticArrays — heap-free (vs `M = [S1 S2 S3; …]`
    # which allocated the matrix, the V vector, and the LU work).
    M = @SMatrix [S1 S2 S3; S2 S3 S4; S3 S4 S5]
    V = @SVector [X1, X2, X3]
    C = M \ V
    a, b, c = C[1], C[2], C[3]

    out = Vector{Float64}(undef, N)
    @inbounds @simd for i in 1:N
        fi = Float64(i)
        out[i] = x[i] - (a + b*fi + c*fi*fi)
    end
    return out
end

function _noise_id_lag1acf(x::Vector{Float64}, m::Int, dmin::Int = 0, dmax::Int = 2; detrend::Bool=true)
    x_dec = m > 1 ? x[1:m:end] : x
    x_det = detrend ? _detrend_quadratic(x_dec) : copy(x_dec)

    d = 0
    while true
        r1  = _lag1_acf(x_det)
        rho = r1 / (1 + r1)

        if d >= dmin && (rho < 0.25 || d >= dmax)
            p       = -2 * (rho + d)
            alpha   = p + 2
            isnan(alpha) && return (NaN, 0, d, rho)
            return (alpha, round(Int, alpha), d, rho)
        end

        x_det = diff(x_det)
        d += 1
        length(x_det) >= 5 || throw(ArgumentError("Data too short after differencing"))
    end
end

function _lag1_acf(x::Vector{Float64})
    # Three single-pass loops over `x` rather than allocating a centred
    # copy and a broadcasted product. Equivalent up to FP reduction order.
    N = length(x)
    N < 2 && return NaN
    s = 0.0
    @inbounds @simd for v in x; s += v; end
    mu = s / N
    ssx = 0.0
    raw = 0.0
    @inbounds @simd for v in x
        d = v - mu
        ssx += d*d
        raw += v*v
    end
    # Scale-invariant degeneracy guard: bail out only when the centred
    # signal really has no power (or float-zero relative to the
    # raw-magnitude scale). The previous form `ssx < eps(Float64) * N`
    # mixed dimensions — it tripped on any data with std below
    # ~√eps ≈ 1.5e-8 regardless of N, so phase records in seconds
    # (typically 1e-9 .. 1e-12) produced spurious NaN classifications
    # even though there was 12+ orders of dynamic range left in the
    # Float64 representation.
    raw > 0 && ssx <= eps(Float64) * raw && return NaN
    ssx == 0.0 && return NaN
    num = 0.0
    @inbounds @simd for i in 1:N-1
        num += (x[i] - mu) * (x[i+1] - mu)
    end
    return num / ssx
end

function _noise_id_b1rn(x::Vector{Float64}, m::Int; detrend::Bool=true)
    x_dec = x[1:m:end]
    if detrend
        x_dec = _detrend_quadratic(x_dec)
    end

    avar_val = _simple_avar(x_dec, 1) / Float64(m)^2
    # Howe 2005 / Barnes 1971: the N in B1_theory(N, r=1, μ) refers to the
    # input data run, not the post-averaging sample count. Using the latter
    # collapses the dynamic range of the theoretical B1 values at small N_eff
    # and misclassifies long-τ points. The original frequency-sample count
    # `length(x) - 1` keeps the boundaries well-separated regardless of AF.
    N_avar   = length(x) - 1

    dx = diff(x)
    Nd = (length(dx) ÷ m) * m
    Nd < m && return (NaN, -2, NaN)

    y_blocks = reshape(dx[1:Nd], m, :)
    y_avg = vec(mean(y_blocks, dims=1))
    # Bessel-corrected (N-1 divisor) — Howe/Barnes B1 theory assumes the
    # unbiased population-variance estimator. The N-divisor variant
    # systematically understates B1 at small N_eff and pushes long-τ
    # classifications into the PM region. Matches Stable32 exactly.
    var_class = var(y_avg)

    (isnan(avar_val) || avar_val <= 0) && return (NaN, -2, NaN)

    B1_obs = var_class / avar_val

    mu_list    = [1, 0, -1, -2]
    alpha_list = [-2, -1, 0, 2]
    b1_vals    = [_b1_theory(N_avar, mu) for mu in mu_list]

    mu_best   = mu_list[end]
    alpha_int = alpha_list[end]

    for i in 1:(length(mu_list) - 1)
        boundary = sqrt(b1_vals[i] * b1_vals[i+1])
        if B1_obs > boundary
            mu_best   = mu_list[i]
            alpha_int = alpha_list[i]
            break
        end
    end

    if mu_best == -2
        adev_val = sqrt(avar_val)
        mdev_val = _simple_mdev(x, m, 1.0)
        if !isnan(mdev_val) && adev_val > 0
            Rn_obs = (mdev_val / adev_val)^2
            R_hi = _rn_theory(m, 0)
            R_lo = _rn_theory(m, -1)
            alpha_int = Rn_obs > sqrt(R_hi * R_lo) ? 1 : 2
        end
    end

    return (alpha_int, mu_best, B1_obs)
end

function _b1_theory(N::Int, mu::Int)
    if mu == 2;  return N * (N + 1) / 6
    elseif mu == 1;  return N / 2
    elseif mu == 0;  return N * log(N) / (2 * (N - 1) * log(2))
    elseif mu == -1; return 1.0
    elseif mu == -2; return (N^2 - 1) / (1.5 * N * (N - 1))
    else;            return (N * (1 - N^mu)) / (2 * (N - 1) * (1 - 2^mu))
    end
end

function _rn_theory(af::Int, b::Int)
    if b == 0
        return 1.0 / af
    elseif b == -1
        avar = (1.038 + 3 * log(2π * 0.5 * af)) / (4π^2)
        mvar = 3 * log(256 / 27) / (8π^2)
        return mvar / avar
    else
        return 1.0
    end
end

function _simple_avar(x::Vector{Float64}, m::Int)
    N = length(x)
    L = N - 2m
    L <= 0 && return NaN
    d2 = @view(x[1+2m:N]) .- 2 .* @view(x[1+m:N-m]) .+ @view(x[1:L])
    return sum(abs2, d2) / (2.0 * Float64(m)^2 * L)
end

function _simple_mdev(x::Vector{Float64}, m::Int, tau0::Float64)
    N  = length(x)
    Ne = N - 3m + 1
    Ne <= 0 && return NaN
    # One length-(N+1) buffer for the prefix sum, then fold the three
    # running-sum differences into a single accumulator. Previously three
    # broadcast subtractions and `cumsum(pushfirst!(copy(x), 0.0))` allocated
    # six length-N(ish) intermediates per call.
    cs = Vector{Float64}(undef, N + 1)
    cs[1] = 0.0
    acc = 0.0
    @inbounds for i in 1:N
        acc += x[i]
        cs[i+1] = acc
    end
    sum_sq = 0.0
    @inbounds for i in 1:Ne
        s1 = cs[i+m]  - cs[i]
        s2 = cs[i+2m] - cs[i+m]
        s3 = cs[i+3m] - cs[i+2m]
        d  = (s3 - 2.0*s2 + s1) / m
        sum_sq += d * d
    end
    return sqrt(sum_sq / (2.0 * Float64(m)^2 * tau0^2 * Ne))
end


# ──────────────────────────────────────────────────────────────────────
# ── stab/noise/synth.jl ─────────────────────────────────────────────────

# noise/synth.jl — Power-law phase-noise synthesizer (non-exported helper).
#
# Generates synthetic phase data with a chosen fractional-frequency power-law
# slope α ∈ {-2, -1, 0, 1, 2} via a frequency-domain `f^(α/2)` shaping of
# white Gaussian noise, integrated to phase. Used by the test suite to
# exercise stability deviations across all five SP1065 noise types.
#
# Requires an `AbstractFFTs` backend (e.g. `FFTW`) to be loaded by the caller.
# `SigmaTau` declares `AbstractFFTs` in its [deps] but no concrete backend;
# pulling in `using FFTW` at the call site supplies the methods.

using AbstractFFTs

"""
    _gen_powerlaw_y(alpha, N; rng=nothing) → Vector{Float64}

Internal primitive. Synthesize an N-sample fractional-frequency vector `y`
with power spectral density ∝ `f^alpha`. The DC component is zeroed so the
output has zero mean; the absolute amplitude is whatever the f^(α/2) shaper
happens to produce on unit-variance white noise. Callers that need a
calibrated level should rescale (e.g. via [`noise_gen`](@ref)).

Pass an `AbstractRNG` as `rng` to draw from a specific stream (e.g. for a
Monte Carlo with independent per-realization seeds); the default `nothing`
draws from the global RNG, so existing seeded callers are unaffected.
"""
function _gen_powerlaw_y(alpha::Real, N::Int; rng=nothing)
    # White Gaussian noise → fractional-frequency series shaped to f^alpha.
    w = rng === nothing ? randn(N) : randn(rng, N)
    W = fft(w)

    # Symmetric FFT-bin frequency magnitude. Bin 0 (DC) gets a placeholder
    # so the f^(alpha/2) term is finite; we zero the DC bin afterwards.
    f = abs.(AbstractFFTs.fftfreq(N, 1.0))
    f[1] = 1.0

    Y    = W .* f .^ (alpha / 2)
    Y[1] = 0.0
    return real.(ifft(Y))
end

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

For deterministic output, seed the global RNG (`Random.seed!`) before calling,
or pass an `AbstractRNG` as `rng`. The DC component is zeroed before the inverse
transform so the output has zero mean. Calls into `AbstractFFTs.fft`/`ifft` —
caller must have an FFT backend loaded (e.g. `using FFTW`).
"""
function _gen_powerlaw_phase(alpha::Real, N::Int; tau0::Real = 1.0, rng=nothing)
    # Phase: x[k] = τ₀ · Σⱼ₌₁ᵏ y[j]
    return cumsum(_gen_powerlaw_y(alpha, N; rng=rng)) .* tau0
end


# ──────────────────────────────────────────────────────────────────────
# ── stab/noise/gen.jl ───────────────────────────────────────────────────

# noise/gen.jl — Public power-law clock-noise generator.
#
# Wraps the internal `_gen_powerlaw_y` spectral shaper with a magnitude-aware
# API. Users specify the noise mixture by α exponent and either a target
# Allan deviation σ_y(τ=τ₀) per noise type or the corresponding fractional-
# frequency PSD coefficient h_α. Contributions for different α are drawn
# independently and summed.

# σ_y(τ=τ₀) from h_α via SP1065 Table 3 with the Nyquist bandwidth
# f_h = 1/(2τ₀) for the bandwidth-dependent WPM (α=2) and FPM (α=1) cases.
function _h_to_sigma1(alpha::Int, h::Float64, tau0::Float64)
    if alpha == 2
        return sqrt(3 * h / (8 * π^2 * tau0^3))
    elseif alpha == 1
        # 2π·f_h·τ = π at τ=τ₀, f_h=1/(2τ₀); SP1065 Table 3 numerator
        # collapses to 1.038 + 3·ln(π).
        coef = 1.038 + 3 * log(π)
        return sqrt(coef * h / (4 * π^2 * tau0^2))
    elseif alpha == 0
        return sqrt(h / (2 * tau0))
    elseif alpha == -1
        return sqrt(2 * log(2) * h)
    elseif alpha == -2
        return sqrt(2 * π^2 * h * tau0 / 3)
    elseif alpha == -4
        throw(ArgumentError(
            "noise_gen: σ_y(τ₀) is undefined for α = -4 — the Allan-variance " *
            "integral ∫S_y(f)|H(f)|²df diverges at low frequency for α ≤ -3, so " *
            "there is no sigma1 ↔ h₋₄ conversion to apply. Specify this component " *
            "by its PSD coefficient instead: h = Dict(-4 => h₋₄)."))
    elseif alpha <= -3
        throw(ArgumentError(
            "noise_gen: α = $alpha is not synthesizable — σ_y(τ₀) is undefined " *
            "(the Allan-variance integral diverges at low frequency for α ≤ -3) " *
            "and there is no generator for it either. Supported: α ∈ {-2, -1, 0, " *
            "1, 2} via `sigma1` or `h`, and α = -4 via `h` only."))
    else
        throw(ArgumentError(
            "noise_gen: unsupported α = $alpha; expected ∈ {-2, -1, 0, 1, 2} " *
            "(α = -4 is supported through `h` only)"))
    end
end

# α = -4 (RRFM, random-run FM) generator — h-coefficient path only.
#
# The f^(α/2) spectral shaper cannot be used here: it zeroes the DC bin of a
# periodic record, which caps the lowest Fourier component exactly where an f⁻⁴
# law puts nearly all of its power. Double discrete integration of white noise
# realizes the law down to the record's own fundamental instead.
#
# Amplitude derivation. Let w[k] be white Gaussian with per-sample variance q,
# sampled at τ₀, and y = cumsum(cumsum(w)). One discrete integrator has transfer
# function H₁(z) = 1/(1 − z⁻¹), so on the unit circle z = exp(i2πfτ₀)
#
#     |H₁(f)|² = 1 / |1 − e^{-i2πfτ₀}|² = 1 / (4 sin²(π f τ₀)),
#
# and the cascade of two gives |H₂(f)|² = 1 / (16 sin⁴(π f τ₀)). White w has the
# one-sided PSD S_w(f) = 2 q τ₀ on 0 ≤ f ≤ 1/(2τ₀) (its integral over that band
# is q, the sample variance), so
#
#     S_y(f) = 2 q τ₀ / (16 sin⁴(π f τ₀)) = q τ₀ / (8 sin⁴(π f τ₀)).
#
# As f → 0, sin(π f τ₀) → π f τ₀ and therefore
#
#     S_y(f) → [q / (8 π⁴ τ₀³)] · f⁻⁴,
#
# so matching S_y(f) = h₋₄ f⁻⁴ fixes the white driver variance at
#
#     q = 8 π⁴ τ₀³ h₋₄.
#
# The same algebra at α = -2 gives h₋₂ = q/(2π²τ₀); as there, the discrete
# construction reproduces the continuous power law for f ≪ 1/(2τ₀) — i.e. at
# τ ≫ τ₀, which is the only regime where an f⁻⁴ specification is meaningful.
# The record mean is removed: every stability estimator is a difference
# operator, so the mean of y is unobservable, and dropping it keeps the
# synthesized series centred.
function _gen_rrfm_y(N::Int, tau0::Float64, h::Float64; rng=nothing)
    q = 8 * π^4 * tau0^3 * h
    w = rng === nothing ? randn(N) : randn(rng, N)
    y = cumsum(cumsum(w)) .* sqrt(q)
    return y .- mean(y)
end

# Empirical σ_y(τ=τ₀) of a fractional-frequency vector via overlapping ADEV
# at m=1. Used to rescale a raw spectral-shaper realization to a known level.
function _measure_sigma1(y::Vector{Float64}, tau0::Float64)
    x = cumsum(y) .* tau0
    return _adev_core(x, [1], tau0)[1]
end

function _noise_gen_y(N::Int, tau0::Float64,
                     sigma1::AbstractDict, h::AbstractDict; rng=nothing)
    !isempty(sigma1) && !isempty(h) &&
        throw(ArgumentError("noise_gen: pass either `sigma1` or `h`, not both"))
    isempty(sigma1) && isempty(h) &&
        throw(ArgumentError("noise_gen: noise mixture is empty; pass `sigma1` or `h`"))
    N >= 4 || throw(ArgumentError("noise_gen: N must be ≥ 4 (got $N)"))
    tau0 > 0 || throw(ArgumentError("noise_gen: tau0 must be > 0 (got $tau0)"))

    # `targets`: α → σ_y(τ₀), realized by the spectral shaper and rescaled to hit
    # the level exactly. `direct`: α → h_α, realized at absolute scale because no
    # finite σ_y(τ₀) exists to rescale against (α = -4 only).
    targets = Dict{Int, Float64}()
    direct  = Dict{Int, Float64}()
    if !isempty(sigma1)
        for (α, s) in sigma1
            α isa Integer ||
                throw(ArgumentError("noise_gen: α keys must be Int, got $(typeof(α))"))
            # α ≤ -3 has no σ_y(τ₀) at all; `_h_to_sigma1` owns that explanation.
            if α <= -3
                _h_to_sigma1(Int(α), 1.0, tau0)
            end
            -2 <= α <= 2 ||
                throw(ArgumentError("noise_gen: α must be ∈ {-2,…,2}, got $α"))
            s >= 0 ||
                throw(ArgumentError("noise_gen: σ values must be ≥ 0, got $s for α=$α"))
            targets[Int(α)] = Float64(s)
        end
    else
        for (α, hα) in h
            α isa Integer ||
                throw(ArgumentError("noise_gen: α keys must be Int, got $(typeof(α))"))
            hα >= 0 ||
                throw(ArgumentError("noise_gen: h values must be ≥ 0, got $hα for α=$α"))
            if α == -4
                direct[-4] = Float64(hα)
            elseif -2 <= α <= 2
                targets[Int(α)] = _h_to_sigma1(Int(α), Float64(hα), tau0)
            else
                throw(ArgumentError(
                    "noise_gen: α must be ∈ {-2,…,2} or -4, got $α"))
            end
        end
    end

    y_total = zeros(Float64, N)
    for (α, σ_target) in targets
        σ_target == 0 && continue
        y_raw = _gen_powerlaw_y(α, N; rng=rng)
        σ_raw = _measure_sigma1(y_raw, tau0)
        σ_raw > 0 || continue
        y_total .+= y_raw .* (σ_target / σ_raw)
    end
    for (_, hα) in direct
        hα == 0 && continue
        y_total .+= _gen_rrfm_y(N, tau0, hα; rng=rng)
    end
    return y_total
end

"""
    noise_gen(::Type{PhaseData},     N, tau0; sigma1=Dict(), h=Dict()) → PhaseData
    noise_gen(::Type{FrequencyData}, N, tau0; sigma1=Dict(), h=Dict()) → FrequencyData

Synthesize a length-`N` clock record whose fractional-frequency power spectrum
is a sum of user-specified power laws `S_y(f) = h_α f^α`. Returns either a
`PhaseData` (integrated to phase with sample interval `tau0`) or a
`FrequencyData` (raw `y` sequence) depending on the first argument.

# α index

| α  | Noise type                | Phase PSD slope |
|---:|:--------------------------|:----------------|
|  2 | WPM — white phase         | flat            |
|  1 | FPM — flicker phase       | 1/f             |
|  0 | WFM — white frequency     | 1/f²            |
| -1 | FFM — flicker frequency   | 1/f³            |
| -2 | RWFM — random-walk freq.  | 1/f⁴            |
| -4 | RRFM — random-run freq.   | 1/f⁶            |

# Specifying the mixture

Pass *either* `sigma1` *or* `h` (not both):

- `sigma1[α] = σ` — target Allan deviation σ_y(τ=τ₀) of the α-component.
  Convenient when working from a clock spec sheet that quotes
  σ_y(1 s) directly. Available for α ∈ {-2,…,2} only.
- `h[α] = h_α` — fractional-frequency PSD coefficient in
  `S_y(f) = h_α f^α`. Unambiguous everywhere.

# α = -4 (RRFM)

RRFM is reachable through `h` only. The Allan variance of an α ≤ -3 process
diverges, so `sigma1[-4]` has no finite value to target and raises an
`ArgumentError` pointing at `h`. The component is synthesized by double
discrete integration of white noise with driver variance `q = 8π⁴τ₀³h₋₄`,
which realizes `S_y(f) → h₋₄ f⁻⁴` down to the record's own fundamental (the
f^(α/2) spectral shaper zeroes the DC bin and cannot). RRFM is invisible to
the Allan family: use a Hadamard-family estimator (`hdev`, `mhdev`, `htotdev`,
`mhtotdev`) — its third-difference kernel converges for α > -5.

For WPM (α=2) and FPM (α=1) the SP1065 σ_y(τ) ↔ h_α relation depends on a
measurement bandwidth `f_h`. This function uses the Nyquist convention
`f_h = 1/(2τ₀)` — the highest frequency representable on a record sampled
at `τ₀`.

Each per-α component is rescaled to hit its requested amplitude exactly for
the drawn realization, so the empirical σ_y(τ=τ₀) of an isolated component
equals the requested value to numerical precision; the natural slope of the
power-law carries through to larger τ. Components from different α are
statistically independent.

Seed the global RNG (`Random.seed!`) before calling for reproducible output,
or pass an `AbstractRNG` as `rng` to draw from a specific stream (e.g. for a
Monte Carlo with independent per-realization seeds).

# Examples

```julia
using SigmaTau, Random
Random.seed!(0)

# Composite: WFM at σ_y(1 s) = 1e-12 plus RWFM at σ_y(1 s) = 1e-14
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))

# Same record but as fractional-frequency
y = noise_gen(FrequencyData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))

# Specify by h_α (here pure WFM with h_0 = 2e-24 → σ_y(1 s) ≈ 1e-12)
noise_gen(PhaseData, 8192, 1.0; h = Dict(0 => 2e-24))

# WFM plus a random-run FM floor (α = -4 is an h-only component)
noise_gen(PhaseData, 8192, 1.0; h = Dict(0 => 2e-24, -4 => 1e-40))
```
"""
function noise_gen(::Type{PhaseData}, N::Int, tau0::Real;
                   sigma1::AbstractDict = Dict{Int,Float64}(),
                   h::AbstractDict      = Dict{Int,Float64}(),
                   rng = nothing)
    τ₀ = Float64(tau0)
    y  = _noise_gen_y(N, τ₀, sigma1, h; rng=rng)
    return PhaseData(cumsum(y) .* τ₀, τ₀)
end

function noise_gen(::Type{FrequencyData}, N::Int, tau0::Real;
                   sigma1::AbstractDict = Dict{Int,Float64}(),
                   h::AbstractDict      = Dict{Int,Float64}(),
                   rng = nothing)
    τ₀ = Float64(tau0)
    y  = _noise_gen_y(N, τ₀, sigma1, h; rng=rng)
    return FrequencyData(y, τ₀)
end


# ──────────────────────────────────────────────────────────────────────
# ── NNLS power-law decomposition ────────────────────────────────────────

# noise/nnls.jl — nonnegative least-squares decomposition of a measured
# variance profile σ̂²(τ) into per-α power-law components.
#
# Where `identify_noise` answers "which single α dominates at this τ", this
# answers "how much of each α is present across the whole τ range" — one fit to
# the entire profile rather than an independent classification per point. The
# nonnegativity constraint is what makes it well posed: variance contributions
# cannot be negative, and without the constraint the near-collinear power-law
# columns trade large positive and negative amplitudes.
#
# The deviation API this dispatches into lives in deviations.jl, which is
# included after this file; that is fine, the call is resolved at run time.

"""
    _nnls_exact(A, b) → (a, rss)

Exact nonnegative least squares — minimize ‖Aa − b‖² subject to `a ≥ 0` — by
exhaustive active-set enumeration. Every one of the `2ⁿ − 1` nonempty supports
is solved as an unconstrained least-squares problem restricted to that support
and kept only if the solution is feasible; `a = 0` seeds the search.

This is exact, not heuristic: the optimum restricted to its own support *is*
the unconstrained minimizer over that subspace, so the optimum appears among
the candidates, and every candidate is a feasible point of the constrained
problem. Cost is exponential in the number of columns, which is why the caller
is capped at 20 — the intended `n` here is 5–7 power laws.

The restricted solves go through `\\` (QR) rather than explicit normal
equations: power-law columns are near-collinear and squaring the condition
number would lose the smallest amplitudes.
"""
function _nnls_exact(A::AbstractMatrix{<:Real}, b::AbstractVector{<:Real})
    n = size(A, 2)
    n <= 20 || throw(ArgumentError(
        "_nnls_exact: $n columns exceeds the exhaustive-enumeration cap of 20"))
    size(A, 1) == length(b) || throw(DimensionMismatch(
        "_nnls_exact: A has $(size(A, 1)) rows but b has $(length(b)) entries"))
    Af = Matrix{Float64}(A)
    bf = Vector{Float64}(b)
    best_a = zeros(n)
    best_r = sum(abs2, bf)             # empty support — always feasible
    for mask in 1:((1 << n) - 1)
        S = [j for j in 1:n if (mask >> (j - 1)) & 1 == 1]
        AS = Af[:, S]
        aS = try
            AS \ bf
        catch
            continue                   # rank-deficient support — skip
        end
        (all(isfinite, aS) && all(>=(0.0), aS)) || continue
        r = sum(abs2, AS * aS .- bf)
        if r < best_r
            best_r = r
            fill!(best_a, 0.0)
            best_a[S] = aS
        end
    end
    return best_a, best_r
end

# Estimator family: `:modified` for the phase-averaged (modified) kernels,
# `:unmodified` for the plain second/third-difference kernels.
function _nnls_family(deviation::Symbol)
    (deviation === :mdev || deviation === :mhdev || deviation === :mhtotdev) &&
        return :modified
    deviation === :hdev && return :unmodified
    throw(ArgumentError(
        "noise_id_nnls: unsupported deviation $deviation; " *
        "expected :mdev, :mhdev, :hdev or :mhtotdev"))
end

# Difference order: 2 for the Allan family (converges for α > -3), 3 for the
# Hadamard family (converges for α > -5).
_nnls_diff_order(deviation::Symbol) = deviation === :mdev ? 2 : 3

"""
    _nnls_mu(family, alpha) → Float64

Response exponent μ of the *variance*: σ²_D(τ) ∝ τ^μ for a pure α power law.

- Modified family (MDEV / MHDEV / MHTOTDEV): the phase pre-average separates
  the two PM noises, so a single relation covers the whole range,
  `μ = -α - 1` — α=2 → -3, α=1 → -2, α=0 → -1, α=-1 → 0, α=-2 → 1, α=-3 → 2,
  α=-4 → 3.
- Unmodified family (ADEV / HDEV): `μ = -α - 1` for α ≤ 0, but WPM and FPM are
  *degenerate* at μ = -2. AVAR/HVAR of white PM is exactly τ⁻² (with a
  bandwidth factor) and of flicker PM is τ⁻² times a slowly varying
  `ln(2π f_h τ)`, which no power-law column can resolve. A basis for an
  unmodified estimator therefore carries at most one of α ∈ {1, 2}.
"""
function _nnls_mu(family::Symbol, alpha::Int)
    family === :modified && return -alpha - 1.0
    return alpha >= 1 ? -2.0 : -alpha - 1.0
end

# Convergent, non-degenerate default basis per estimator.
function _nnls_default_alphas(deviation::Symbol)
    # MHTOTDEV's per-window detrend is quadratic (matched to the third
    # difference's annihilation order), so its kernel converges over the full
    # Hadamard range α > -5 like MHDEV — verified on pure RRFM (variance slope
    # 3.00, level N-stable). It therefore carries the full modified-family
    # basis. (The pre-quadratic kernel had a second-difference low-frequency
    # response and excluded α ≤ -3 here.)
    family = _nnls_family(deviation)
    if _nnls_diff_order(deviation) == 2
        # Allan family: the second-difference kernel diverges for α ≤ -3.
        return family === :modified ? [2, 1, 0, -1, -2] : [2, 0, -1, -2]
    end
    # Hadamard family: third differences converge down to α = -4.
    return family === :modified ? [2, 1, 0, -1, -2, -3, -4] : [2, 0, -1, -2, -3, -4]
end

# Smallest averaging factor kept by the default grid — see `noise_id_nnls`.
const NNLS_M_MIN = 4

function _nnls_symbol(alpha::Int)
    alpha == 2  && return :WHPM
    alpha == 1  && return :FLPM
    alpha == 0  && return :WHFM
    alpha == -1 && return :FLFM
    alpha == -2 && return :RWFM
    alpha == -3 && return :FWFM       # flicker-walk FM
    alpha == -4 && return :RRFM       # random-run FM
    return :unknown
end

# Deviation dispatch: the RAW profile, uniform kwargs. `ci=false` skips the
# lag-1/B1 classifier and `correct_bias=false` skips MHTOTDEV's bias divisor —
# both of which are functions of `identify_noise`, i.e. of exactly the
# classification this routine is meant to replace. Fitting a bias-corrected
# profile would make the decomposition depend on the answer it is computing.
function _nnls_deviation(deviation::Symbol, data::PhaseData, m_values::Vector{Int})
    deviation === :mdev  && return mdev(data, m_values; ci=false)
    deviation === :mhdev && return mhdev(data, m_values; ci=false)
    deviation === :hdev  && return hdev(data, m_values; ci=false)
    deviation === :mhtotdev &&
        return mhtotdev(data, m_values; ci=false, correct_bias=false)
    throw(ArgumentError("noise_id_nnls: unsupported deviation $deviation"))
end

# The shipped EDF tables stop at α = -2 (`_alpha_from_noise` maps anything else
# to α = 0), so the α ≤ -3 rows borrow the RWFM model — the closest available.
# This only sets the χ² weights; it does not enter the fitted amplitudes.
_nnls_edf_symbol(alpha::Int) = alpha <= -3 ? :RWFM : _nnls_symbol(alpha)

"""
    noise_id_nnls(data::PhaseData [, m_values]; deviation=:mhdev, alphas=…, weights=:edf)
    noise_id_nnls(data::FrequencyData [, m_values]; …)

Decompose a record's measured variance profile into nonnegative power-law
components. One weighted nonnegative least-squares fit over the whole τ range
returns the amplitude of every α in the basis, plus the dominant noise type at
each τ.

Returns a `NamedTuple`:

| field | meaning |
|:------|:--------|
| `deviation` | estimator used |
| `tau` | τ grid actually fitted (s) |
| `var` | measured σ̂²(τ) |
| `fit` | fitted Σᵢ aᵢ (τ/τ₀)^μᵢ |
| `alphas` | basis exponents |
| `noise_types` | `:WHPM`, `:FLPM`, `:WHFM`, `:FLFM`, `:RWFM`, `:FWFM`, `:RRFM` per α |
| `amplitudes` | aᵢ ≥ 0 — the α-component's *variance at τ = τ₀* |
| `contributions` | `length(tau) × length(alphas)` matrix aᵢ (τ/τ₀)^μᵢ |
| `dominant` | dominant noise type at each τ |
| `residual` | weighted sum of squared residuals |
| `weights` | weight vector used |

!!! warning "Experimental return shape"
    The set and naming of these fields is experimental and may change in a
    minor release; `amplitudes`, `alphas` and `dominant` are the stable core.

The fit runs on the **raw** deviation profile: `ci=false` and, for MHTOTDEV,
`correct_bias=false`. Both the CI machinery and the total-family bias divisor
are functions of `identify_noise`, so fitting a corrected profile would make
the decomposition depend on the classification it is computing.

# Basis

The columns are analytic power laws `Kᵢ(τ) = (τ/τ₀)^μ(αᵢ)`, with μ the
*variance* response exponent — `μ = -α-1` for the modified family
(MDEV/MHDEV/MHTOTDEV) and, for the unmodified family (HDEV), `μ = -α-1` for
α ≤ 0 with WPM and FPM both collapsing onto μ = -2. The defaults therefore
carry α ∈ {2, 1, 0, -1, -2} for `:mdev`, {2, 0, -1, -2, -3, -4} for `:hdev`,
and {2, 1, 0, -1, -2, -3, -4} for `:mhdev` and `:mhtotdev`. The convergence
condition per estimator is:

| estimator | converges for | reason |
|:----------|:--------------|:-------|
| `:mdev` | α > -3 | second-difference kernel |
| `:mhdev`, `:hdev`, `:mhtotdev` | α > -5 | third-difference kernel (MHTOTDEV's per-window quadratic detrend matches the kernel's annihilation order, preserving the range through the boundary extension) |

An unmodified basis additionally cannot carry both PM columns — they are the
same column. Passing a basis that violates either rule is an error.

# Weights

- `:edf` (default) — `wₖ = νₖ / (2 σ̂⁴ₖ)`, the inverse of the χ² sampling
  variance of a variance estimator with νₖ degrees of freedom. The per-τ noise
  type that the EDF model needs comes from a first `:relative`-weighted pass of
  this same fit, not from [`identify_noise`](@ref) — the routine is meant to
  replace that classifier, so it does not lean on it. Points where the EDF model
  returns a non-finite or non-positive ν fall back to the `:relative` weight,
  and α ≤ -3 rows borrow the RWFM EDF (the shipped tables stop at α = -2).
- `:relative` — `wₖ = 1 / σ̂⁴ₖ`, i.e. equal *fractional* weight per τ.
- any `AbstractVector` of per-τ weights.

# Caveats

The analytic basis assumes each component is in its asymptotic power-law
regime over the whole fitted range. Two things break that:

- **Total estimators.** MHTOTDEV's reflected/extended windows bend the response
  away from τ^μ as τ → T — a bend, not a different power: the per-window
  quadratic detrend matches the third difference's annihilation order, so the
  kernel converges over the full α > -5 range (a pure RRFM record gives the
  MHDEV-family τ³ with an N-stable level). The bend still distorts the last
  octaves. The `tools/nnls_noise_id.jl` study uses an *empirical* basis
  (ensemble-mean σ̂²(τ) of synthesized unit noises) for exactly this reason;
  prefer that route when the long-τ octaves matter.
- **Flicker PM under an unmodified estimator.** The `ln(2π f_h τ)` factor is
  absorbed into the α=2 column, so an FPM record fitted with `:hdev` reports
  white PM.
- **Small τ.** Every estimator is still in transition toward its asymptotic
  power law at m = 1–2 (the modified family runs 2–3× high there, settling by
  m ≈ 8), so the default grid starts at m = 4. Fitting the dropped octaves
  lets a pure WFM record buy off its m = 1 excess with a spurious WPM column.

Amplitudes are variances at τ = τ₀; `sqrt.(result.amplitudes)` is the
corresponding σ_y(τ₀) per component.

# Examples

```julia
using SigmaTau, Random
Random.seed!(0)

p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))
r = noise_id_nnls(p)
collect(zip(r.noise_types, r.amplitudes))   # amplitude per α
r.dominant                                  # dominant type per τ
```

See also [`identify_noise`](@ref) for the per-τ lag-1 ACF classifier and
[`noise_gen`](@ref) for the matching synthesizer.
"""
function noise_id_nnls(data::PhaseData, m_values::Vector{Int} = Int[];
                       deviation::Symbol = :mhdev,
                       alphas::AbstractVector{<:Integer} = _nnls_default_alphas(deviation),
                       weights::Union{Symbol, AbstractVector{<:Real}} = :edf)
    x  = _f64(data.x)
    τ₀ = data.tau0
    # Default grid: the estimator's own octave grid with the first two octaves
    # dropped. Every estimator's response to a power law is still in transition
    # at m = 1–2 (the modified family runs 2–3× above its asymptotic τ^μ there,
    # converging by m ≈ 8), and those are exactly the points the EDF weighting
    # trusts most — leaving them in lets a pure WFM record buy off its m = 1
    # excess with a spurious WPM column. Pass `m_values` explicitly to override.
    ms = isempty(m_values) ? filter(>=(NNLS_M_MIN), _default_m_values(length(x), deviation)) :
                             m_values
    isempty(ms) && throw(ArgumentError(
        "noise_id_nnls: record of $(length(x)) samples is too short for the default " *
        "m ≥ $NNLS_M_MIN grid; pass `m_values` explicitly"))

    family = _nnls_family(deviation)
    αs     = Int.(collect(alphas))
    isempty(αs) && throw(ArgumentError("noise_id_nnls: `alphas` is empty"))
    αmin   = _nnls_diff_order(deviation) == 2 ? -2 : -4
    for α in αs
        α <= 2 || throw(ArgumentError("noise_id_nnls: α = $α exceeds the α = 2 ceiling"))
        α >= αmin || throw(ArgumentError(
            "noise_id_nnls: $deviation does not converge at α = $α — its " *
            "$(_nnls_diff_order(deviation))-difference kernel requires α ≥ $αmin. " *
            "Use a Hadamard-family estimator (:hdev, :mhdev, :mhtotdev) for α < -2."))
    end
    μs = [_nnls_mu(family, α) for α in αs]
    if length(unique(μs)) != length(μs)
        throw(ArgumentError(
            "noise_id_nnls: the requested basis is degenerate for $deviation — " *
            "α = $αs map to response exponents μ = $μs. The unmodified family " *
            "cannot separate WPM (α=2) from FPM (α=1); drop one of them."))
    end

    res  = _nnls_deviation(deviation, PhaseData(x, τ₀), ms)
    keep = findall(k -> isfinite(res.dev[k]) && res.dev[k] > 0, eachindex(res.dev))
    length(keep) >= length(αs) || throw(ArgumentError(
        "noise_id_nnls: only $(length(keep)) usable τ points for a " *
        "$(length(αs))-column basis; shorten `alphas` or lengthen the record"))

    taus = res.tau[keep]
    v    = res.dev[keep] .^ 2

    # Weighted, column-equilibrated design matrix. The power-law columns span
    # many orders of magnitude over a wide τ grid; equilibrating to unit 2-norm
    # keeps the restricted QR solves well scaled, and the amplitudes are
    # unscaled afterwards.
    K = [(taus[k] / τ₀) ^ μs[j] for k in eachindex(taus), j in eachindex(μs)]
    function solve(w)
        sw = sqrt.(w)
        A  = sw .* K
        b  = sw .* v
        d  = [sqrt(sum(abs2, @view A[:, j])) for j in eachindex(μs)]
        all(dj -> isfinite(dj) && dj > 0, d) ||
            throw(ArgumentError("noise_id_nnls: degenerate basis column (norms $d)"))
        â, rss = _nnls_exact(A ./ d', b)
        return â ./ d, rss
    end

    w = if weights isa AbstractVector
        length(weights) == length(keep) || throw(DimensionMismatch(
            "noise_id_nnls: `weights` has $(length(weights)) entries, expected $(length(keep))"))
        Float64.(collect(weights))
    elseif weights === :relative
        1.0 ./ v .^ 2
    elseif weights === :edf
        # Self-consistent χ² weights: a first pass with equal fractional weight
        # supplies the per-τ noise type, the EDF model turns that into ν, and the
        # second pass uses wₖ = νₖ/(2 σ̂⁴ₖ). Deriving the types from the fit rather
        # than from `identify_noise` keeps the whole routine independent of the
        # classifier it is meant to replace.
        a0, _ = solve(1.0 ./ v .^ 2)
        syms  = [_nnls_edf_symbol(αs[argmax(vec(K[k, :]) .* a0)]) for k in eachindex(taus)]
        ν     = calculate_edf(deviation, sqrt.(v), syms, ms[keep], taus,
                              length(x), (length(x) - 1) * τ₀)
        [(isfinite(ν[k]) && ν[k] > 0 ? ν[k] / 2 : 1.0) / v[k] ^ 2 for k in eachindex(v)]
    else
        throw(ArgumentError(
            "noise_id_nnls: unknown weights $weights; expected :edf, :relative or a vector"))
    end
    all(wk -> isfinite(wk) && wk > 0, w) ||
        throw(ArgumentError("noise_id_nnls: weights must be finite and positive"))

    a, rss = solve(w)

    contributions = K .* a'
    dominant = [_nnls_symbol(αs[argmax(@view contributions[k, :])])
                for k in eachindex(taus)]

    return (; deviation, tau = taus, var = v, fit = vec(sum(contributions; dims=2)),
            alphas = αs, noise_types = _nnls_symbol.(αs), amplitudes = a,
            contributions, dominant, residual = rss, weights = w)
end

function noise_id_nnls(data::FrequencyData, m_values::Vector{Int} = Int[]; kwargs...)
    return noise_id_nnls(_freq_to_phase(data), m_values; kwargs...)
end

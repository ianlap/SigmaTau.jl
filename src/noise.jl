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
    else
        throw(ArgumentError(
            "noise_gen: unsupported α = $alpha; expected ∈ {-2, -1, 0, 1, 2}"))
    end
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

    targets = Dict{Int, Float64}()
    if !isempty(sigma1)
        for (α, s) in sigma1
            α isa Integer ||
                throw(ArgumentError("noise_gen: α keys must be Int, got $(typeof(α))"))
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
            -2 <= α <= 2 ||
                throw(ArgumentError("noise_gen: α must be ∈ {-2,…,2}, got $α"))
            hα >= 0 ||
                throw(ArgumentError("noise_gen: h values must be ≥ 0, got $hα for α=$α"))
            targets[Int(α)] = _h_to_sigma1(Int(α), Float64(hα), tau0)
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

# Specifying the mixture

Pass *either* `sigma1` *or* `h` (not both):

- `sigma1[α] = σ` — target Allan deviation σ_y(τ=τ₀) of the α-component.
  Convenient when working from a clock spec sheet that quotes
  σ_y(1 s) directly.
- `h[α] = h_α` — fractional-frequency PSD coefficient in
  `S_y(f) = h_α f^α`. Unambiguous everywhere.

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

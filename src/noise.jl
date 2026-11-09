# Power-law noise identification (lag-1 ACF with B1 / R(n) fallback),
# calibrated power-law noise synthesis and the NNLS variance decomposition.

const NEFF_RELIABLE = 30

"""
    identify_noise(x::Vector{Float64}, m_values::Vector{Int}; dmin::Int=0, dmax::Int=2, remove_drift::Bool=true, b1_length::Symbol=:record) → Vector{Symbol}
    identify_noise(data::PhaseData, m_values; kwargs...) → Vector{Symbol}
    identify_noise(data::FrequencyData, m_values; kwargs...) → Vector{Symbol}

The record forms take the same keywords and accept any integer vector or
[`TauMode`](@ref) for `m_values`; a `FrequencyData` record is integrated to
phase first. Use `dmax=2` for the Allan family and `dmax=3` for the Hadamard
family, as the deviation functions do.

Identify the dominant power-law noise type at each averaging factor by the
lag-1 autocorrelation of the differenced phase, with a B1/R(n) fallback when
fewer than `NEFF_RELIABLE` effective samples remain.

Every entry is one of `:WHPM`, `:FLPM`, `:WHFM`, `:FLFM`, `:RWFM`; the
function never returns `:unknown`. Estimates outside α ∈ [-2, 2] are
clamped to the nearest end, an averaging factor whose estimate fails
takes the type of the successful one nearest to it in log(m), and a
record with no usable estimate at all (e.g. constant) is reported as
`:WHFM`.

The full record always passes through a 5σ outlier filter before the
per-m loop; that step is unconditional and unrelated to the polynomial
detrend below.

## `remove_drift`

When `true` (default), each decimated-by-m subseries is quadratically
detrended before the lag-1 ACF / B1 ratio is computed. Allantools'
`autocorr_noise_id` does the same thing (`detrend(x, deg=2)` in
[`ci.py`](https://github.com/aewallin/allantools/blob/master/allantools/ci.py))
to remove a frequency offset (linear phase ramp) and drift (quadratic
phase term) from the slow content of each averaging window.

When `false`, the per-m polynomial fit is skipped. This matches
Stable32's noise-ID convention; use it when you need α values that line
up with a Stable32-generated reference fixture point-for-point. (Before
0.6 this keyword was called `detrend`; `detrend=false` still works, with a
deprecation warning.)

## `b1_length`

Sets the record length `N` at which the B1 fallback (used when fewer than
`NEFF_RELIABLE` averaged points remain) evaluates the theoretical
B1(N, μ) values that separate the noise types.

- `:record` (default): the full record's frequency-sample count, as
  Stable32 does. α values line up with Stable32 point for point.
- `:averaged`: the number of m-averaged frequency points the observed B1
  is computed from. This is the length the B1 theory describes, and it
  identifies random-walk FM far more reliably at long τ (in Monte Carlo
  runs with N = 1024 and m = 40–128, 55–75 % correct instead of 1–4 %),
  at the cost of somewhat weaker white-FM identification and less
  agreement with Stable32.
"""
function identify_noise(x::Vector{Float64}, m_values::Vector{Int}; dmin::Int=0, dmax::Int=2,
                        remove_drift::Bool=true, detrend::Bool=true, b1_length::Symbol=:record)
    b1_length in (:record, :averaged) ||
        throw(ArgumentError("b1_length must be :record or :averaged, got :$b1_length"))
    detrend || Base.depwarn("identify_noise: `detrend` is deprecated, use `remove_drift`", :identify_noise)
    detrend = remove_drift && detrend
    x_clean = _preprocess(x)
    N = length(x_clean)
    alphas = Vector{Int}(undef, length(m_values))
    ok = falses(length(m_values))

    for (k, m) in enumerate(m_values)
        N_eff = N ÷ m
        alpha = NaN

        try
            if N_eff >= NEFF_RELIABLE
                alpha, _, _, _ = _noise_id_lag1acf(x_clean, m, dmin, dmax; detrend=detrend)
            else
                alpha, _, _ = _noise_id_b1rn(x_clean, m; detrend=detrend, b1_length=b1_length)
            end
        catch
            alpha = NaN
        end

        isnan(alpha) && continue
        # The lag-1 estimate is continuous and noisy: a white-PM record can
        # come out at α ≈ 2.6, and at d = dmax a ρ near its 0.5 ceiling gives
        # α ≈ -3 (and FWFM/RRFM-dominated data, or an infinite ratio, can go
        # further). None of those has a symbol here, so clamp to the
        # supported range and take the nearest power law instead of
        # reporting :unknown.
        alphas[k] = round(Int, clamp(alpha, -2.0, 2.0))
        ok[k] = true
    end

    return _fill_noise_symbols(alphas, ok, m_values)
end

# Record forms: phase taken from the record (frequency integrated), any integer
# vector (or a TauMode, defined in deviations.jl) accepted for the averaging factors.
identify_noise(data::PhaseData, m_values::AbstractVector{<:Integer}; kwargs...) =
    identify_noise(Vector{Float64}(data.x), collect(Int, m_values); kwargs...)
identify_noise(data::FrequencyData, m_values::AbstractVector{<:Integer}; kwargs...) =
    identify_noise(_freq_to_phase(data), m_values; kwargs...)

const _NOISE_SYMBOLS = (:RWFM, :FLFM, :WHFM, :FLPM, :WHPM)  # α = -2 … 2

# Map integer α to noise symbols. An averaging factor whose estimate failed
# (too few points after differencing, a degenerate decimated series) takes
# the type of the successful averaging factor closest to it in log(m), the
# natural distance on octave or decade τ grids; a tie goes to the shorter m
# because it rests on more data. With no successful estimate at all
# (constant or very short record), fall back to :WHFM, the same α = 0
# default `_alpha_from_noise` uses for EDF.
function _fill_noise_symbols(alphas::Vector{Int}, ok::AbstractVector{Bool},
                             m_values::AbstractVector{<:Integer})
    n = length(alphas)
    noises = Vector{Symbol}(undef, n)
    any(ok) || return fill!(noises, :WHFM)
    for k in 1:n
        j = k
        if !ok[k]
            best = Inf
            for i in 1:n
                ok[i] || continue
                dist = abs(log(m_values[i] / m_values[k]))
                if dist < best || (dist == best && m_values[i] < m_values[j])
                    best = dist
                    j = i
                end
            end
        end
        noises[k] = _NOISE_SYMBOLS[alphas[j] + 3]
    end
    return noises
end

function _preprocess(x::Vector{Float64})
    # 5σ outlier filter only; Stable32 and allantools detrend per m, as the
    # identifiers below do.
    x_std = std(x)
    x_std < eps() && return x  # degenerate input; nothing to filter against
    x_mean = mean(x)
    thresh = 5.0 * x_std
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

function _noise_id_lag1acf(x::Vector{Float64}, m::Int, dmin::Int = 0, dmax::Int = 2; detrend::Bool=true)
    x_dec = m > 1 ? x[1:m:end] : x
    x_det = detrend ? _remove_quadratic(x_dec) : copy(x_dec)

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
    # Degenerate record: no power relative to the raw magnitude (phase in
    # seconds can be 1e-12, so no absolute threshold).
    raw > 0 && ssx <= eps(Float64) * raw && return NaN
    ssx == 0.0 && return NaN
    num = 0.0
    @inbounds @simd for i in 1:N-1
        num += (x[i] - mu) * (x[i+1] - mu)
    end
    return num / ssx
end

function _noise_id_b1rn(x::Vector{Float64}, m::Int; detrend::Bool=true, b1_length::Symbol=:record)
    x_dec = x[1:m:end]
    detrend && (x_dec = _remove_quadratic(x_dec))

    avar_val = _adev_core(x_dec, [1], 1.0)[1]^2 / Float64(m)^2
    dx = diff(x)
    Nd = (length(dx) ÷ m) * m
    Nd < m && return (NaN, -2, NaN)

    y_blocks = reshape(dx[1:Nd], m, :)
    y_avg = vec(mean(y_blocks, dims=1))
    # Unbiased variance (N − 1), as B1 theory and Stable32 assume.
    var_class = var(y_avg)

    # N in B1_theory(N, μ). `:record` uses the full frequency-sample count
    # `length(x) - 1`, as Stable32 does; `:averaged` uses the number of
    # averaged points B1_obs is computed from, which is what the theory
    # describes (see the `b1_length` docs of `identify_noise`).
    N_avar = b1_length === :averaged ? length(y_avg) : length(x) - 1
    N_avar < 2 && return (NaN, -2, NaN)

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
        mdev_val = _mdev_core(x, [m], 1.0)[1]
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

# ── Power-law noise synthesis ─────────────────────────────────────────────

# Fractional frequency with S_y ∝ |2 sin(πf)|^α, the discrete-time f^α: white
# Gaussian noise through the Kasdin–Walter filter (1 − z⁻¹)^(α/2), exact at
# the integer orders (α = 2 is the first difference of white noise, α = 0
# white, α = −2 a cumulative sum). Zero-mean and uncalibrated; `noise_gen`
# rescales. `rng = nothing` draws from the global RNG.
function _gen_powerlaw_y(alpha::Real, N::Int; rng=nothing)
    w = rng === nothing ? randn(N) : randn(rng, N)
    h = Vector{Float64}(undef, N)
    h[1] = 1.0
    @inbounds for k in 1:N-1
        h[k+1] = h[k] * (k - 1 - alpha / 2) / k
    end
    # Linear (not circular) convolution through zero-padded FFTs.
    M = 2N
    W = fft(vcat(w, zeros(M - N)))
    H = fft(vcat(h, zeros(M - N)))
    y = real.(ifft(W .* H))[1:N]
    return y .- (sum(y) / N)
end

# Phase x = τ₀ · cumsum(y) of a `_gen_powerlaw_y` draw.
function _gen_powerlaw_phase(alpha::Real, N::Int; tau0::Real = 1.0, rng=nothing)
    return cumsum(_gen_powerlaw_y(alpha, N; rng=rng)) .* tau0
end


# σ_y(τ₀) of a fractional-frequency vector, to rescale a raw draw to a level.
function _measure_sigma1(y::Vector{Float64}, tau0::Float64)
    x = cumsum(y) .* tau0
    return _adev_core(x, [1], tau0)[1]
end

# Power-law exponent α of a noise_gen key: an integer α or a noise-type name.
const _NOISE_ALPHA = Dict(:WHPM => 2, :FLPM => 1, :WHFM => 0, :FLFM => -1,
                          :RWFM => -2, :FWFM => -3, :RRFM => -4)
_noise_alpha(α::Integer) = Int(α)
_noise_alpha(s::Symbol) = haskey(_NOISE_ALPHA, s) ? _NOISE_ALPHA[s] :
    throw(ArgumentError("noise_gen: unknown noise type :$s; use one of " * join(keys(_NOISE_ALPHA), ", ")))
_noise_alpha(k) = throw(ArgumentError("noise_gen: keys must be an Int α or a noise-type Symbol, got $(typeof(k))"))

function _noise_gen_y(N::Int, tau0::Float64,
                     sigma1::AbstractDict, h::AbstractDict; rng=nothing)
    !isempty(sigma1) && !isempty(h) &&
        throw(ArgumentError("noise_gen: pass either `sigma1` or `h`, not both"))
    isempty(sigma1) && isempty(h) &&
        throw(ArgumentError("noise_gen: noise mixture is empty; pass `sigma1` or `h`"))
    N >= 4 || throw(ArgumentError("noise_gen: N must be ≥ 4 (got $N)"))
    tau0 > 0 || throw(ArgumentError("noise_gen: tau0 must be > 0 (got $tau0)"))

    check(α, v, what) = begin
        -4 <= α <= 2 ||
            throw(ArgumentError("noise_gen: α must be ∈ {-4,…,2}, got $α"))
        v >= 0 ||
            throw(ArgumentError("noise_gen: $what values must be ≥ 0, got $v for α=$α"))
    end

    y_total = zeros(Float64, N)
    if !isempty(sigma1)
        # Empirical calibration: rescale the realization so its overlapping
        # ADEV at τ₀ equals the request exactly.
        for (key, s) in sigma1
            α = _noise_alpha(key)
            check(α, s, "σ")
            s == 0 && continue
            y_raw = _gen_powerlaw_y(α, N; rng=rng)
            σ_raw = _measure_sigma1(y_raw, tau0)
            σ_raw > 0 || continue
            y_total .+= y_raw .* (Float64(s) / σ_raw)
        end
    else
        # Analytic calibration of the PSD coefficient: the Kasdin–Walter
        # filter driven by unit white noise sampled at τ₀ has one-sided PSD
        # 2τ₀·|2 sin(πfτ₀)|^α = 2τ₀·(2πτ₀f)^α at low frequency, so scaling by
        # √(h_α / (2τ₀(2πτ₀)^α)) gives S_y(f) = h_α f^α there.
        for (key, hα) in h
            α = _noise_alpha(key)
            check(α, hα, "h")
            hα == 0 && continue
            y_raw = _gen_powerlaw_y(α, N; rng=rng)
            y_total .+= y_raw .* sqrt(Float64(hα) / (2 * tau0 * (2π * tau0)^α))
        end
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
|  2 | WHPM, white phase          | flat            |
|  1 | FLPM, flicker phase        | 1/f             |
|  0 | WHFM, white frequency      | 1/f²            |
| -1 | FLFM, flicker frequency    | 1/f³            |
| -2 | RWFM, random-walk frequency | 1/f⁴           |
| -3 | FWFM, flicker-walk frequency | 1/f⁵          |
| -4 | RRFM, random-run frequency  | 1/f⁶           |

Keys are the integer α or the noise-type name (`Dict(:WHFM => 1e-12)`).

# Specifying the mixture

Pass *either* `sigma1` *or* `h` (not both):

- `sigma1[α] = σ`: target Allan deviation σ_y(τ=τ₀) of the α-component.
  Convenient when working from a clock spec sheet that quotes
  σ_y(1 s) directly.
- `h[α] = h_α`: fractional-frequency PSD coefficient in
  `S_y(f) = h_α f^α`. Unambiguous everywhere.

`sigma1` components are rescaled so the empirical σ_y(τ=τ₀) of the drawn
realization equals the request exactly. `h` components are scaled
analytically: the Kasdin–Walter filter driven by unit white noise has
one-sided PSD `2τ₀ |2 sin(πfτ₀)|^α`, which is `2τ₀ (2πτ₀ f)^α` at low
frequency, so the component is multiplied by `√(h_α / (2τ₀ (2πτ₀)^α))` and
its spectrum is

    S_y(f) = h_α · (sin(π f τ₀) / (π τ₀))^α,

which is `h_α f^α` up to the sine-versus-linear roll-off near Nyquist. For
α ≤ −1 the process is not stationary and this is the spectrum of its
differences: the d-th difference of `y` has spectrum
`(2 sin(π f τ₀))^(2d) S_y(f)`. The test suite checks the spectrum, through
differences where needed, for every α from +2 to −4. The empirical σ_y(τ₀) of an `h` component therefore fluctuates
with the realization (≈ √(2/N) relative), and for the coloured noises it is
not the continuous-time SP1065 value, because a sampled discrete process is
not a band-limited continuous one. Components from different α are
statistically independent.

Seed the global RNG (`Random.seed!`) before calling for reproducible output,
or pass an `AbstractRNG` as `rng` to draw from a specific stream (e.g. for a
Monte Carlo with independent per-realization seeds).

# Examples

```julia
using SigmaTau, Random
Random.seed!(0)

# Composite: WHFM at σ_y(1 s) = 1e-12 plus RWFM at σ_y(1 s) = 1e-14
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))

# Same record but as fractional-frequency
y = noise_gen(FrequencyData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))

# Specify by h_α (here pure WHFM with h_0 = 2e-24 → σ_y(1 s) ≈ 1e-12)
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
    # MHTOTDEV's per-window detrend is quadratic, matched to the third
    # difference, so it converges over the full Hadamard range like MHDEV.
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

# Weighted NNLS of v ≈ K·a. The power-law columns span many orders of
# magnitude over a wide τ grid; equilibrating them to unit 2-norm keeps the
# restricted QR solves well scaled, and the amplitudes are unscaled afterwards.
function _nnls_weighted(K::Matrix{Float64}, v::Vector{Float64}, w::Vector{Float64})
    sw = sqrt.(w)
    A  = sw .* K
    d  = [sqrt(sum(abs2, @view A[:, j])) for j in axes(A, 2)]
    all(dj -> isfinite(dj) && dj > 0, d) ||
        throw(ArgumentError("degenerate power-law basis column"))
    â, rss = _nnls_exact(A ./ d', sw .* v)
    return â ./ d, rss
end

# Self-consistent χ² weights wₖ = νₖ/(2 σ̂⁴ₖ): a first pass with equal
# fractional weight supplies the per-τ noise type, and the EDF model turns
# that into ν, with `Ns[k]` the record length at the k-th τ. Deriving the types from the fit rather than from
# `identify_noise` keeps the fit independent of that classifier. Points where
# the EDF model gives no finite positive ν keep the fractional weight.
function _nnls_edf_weights(deviation::Symbol, K::Matrix{Float64}, v::Vector{Float64},
                           αs::Vector{Int}, ms::Vector{Int}, taus::Vector{Float64},
                           Ns::Vector{Int}, τ₀::Float64)
    a0, _ = _nnls_weighted(K, v, 1.0 ./ v .^ 2)
    map(eachindex(v)) do k
        sym = _nnls_edf_symbol(αs[argmax(vec(K[k, :]) .* a0)])
        ν = only(calculate_edf(deviation, [sqrt(v[k])], [sym], [ms[k]], [taus[k]], Ns[k], (Ns[k] - 1) * τ₀))
        (isfinite(ν) && ν > 0 ? ν / 2 : 1.0) / v[k] ^ 2
    end
end

"""
    noise_id_nnls(data [, m_values]; deviation=:mhdev, alphas=…, weights=:edf) → NamedTuple

Decompose a record's measured variance profile into nonnegative power-law
components: one weighted nonnegative least-squares fit of
σ̂²(τ) ≈ Σᵢ aᵢ (τ/τ₀)^μ(αᵢ) over the whole τ range, giving the amplitude of
every α in the basis and the dominant noise type at each τ. The fit runs on
the raw deviation (`ci=false`, no bias correction), so it does not lean on
[`identify_noise`](@ref), which it is meant to replace at run time. `data`
is a `PhaseData` or `FrequencyData`.

Returns a `NamedTuple` with `deviation`, `tau` (the τ fitted), `var` (σ̂²),
`fit` (the model at each τ), `alphas`, `noise_types`, `amplitudes` (aᵢ ≥ 0,
each component's variance at τ = τ₀), `contributions` (a τ × α matrix),
`dominant` (per τ), `residual` and `weights`.

`deviation` is `:mdev`, `:hdev`, `:mhdev` (default) or `:mhtotdev`. `alphas`
defaults to every α the estimator converges for: 2…−2 for `:mdev`, 2…−4 for
the Hadamard family, which also tolerates drift; `:hdev` cannot tell white
from flicker PM and takes one of them. The default grid is the estimator's
octave grid from m = 4, since every estimator is still in transition toward
its power law at m = 1–2. `weights` is `:edf` (χ² weights νₖ/(2σ̂⁴ₖ), the
noise type per τ taken from a first `:relative` pass), `:relative` (equal
fractional weight) or a vector.

The basis assumes each component is in its asymptotic regime over the
fitted range; the total estimators bend away from τ^μ as τ → T, so prefer
`:mhdev` when the long-τ octaves matter. For calibrated h_α see
[`powerlaw_fit`](@ref), which fits each component's exact response.

# Examples

```julia
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 1e-14))
r = noise_id_nnls(p)
collect(zip(r.noise_types, r.amplitudes))   # amplitude per α
r.dominant                                  # dominant type per τ
```
"""
function noise_id_nnls(data::PhaseData, m_values::Vector{Int} = Int[];
                       deviation::Symbol = :mhdev,
                       alphas::AbstractVector{<:Integer} = _nnls_default_alphas(deviation),
                       weights::Union{Symbol, AbstractVector{<:Real}} = :edf)
    x  = _f64(data.x)
    τ₀ = data.tau0
    # The estimator's octave grid from m = NNLS_M_MIN; see the docstring.
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

    K = [(taus[k] / τ₀) ^ μs[j] for k in eachindex(taus), j in eachindex(μs)]

    w = if weights isa AbstractVector
        length(weights) == length(keep) || throw(DimensionMismatch(
            "noise_id_nnls: `weights` has $(length(weights)) entries, expected $(length(keep))"))
        Float64.(collect(weights))
    elseif weights === :relative
        1.0 ./ v .^ 2
    elseif weights === :edf
        _nnls_edf_weights(deviation, K, v, αs, ms[keep], taus, fill(length(x), length(keep)), τ₀)
    else
        throw(ArgumentError(
            "noise_id_nnls: unknown weights $weights; expected :edf, :relative or a vector"))
    end
    all(wk -> isfinite(wk) && wk > 0, w) ||
        throw(ArgumentError("noise_id_nnls: weights must be finite and positive"))

    a, rss = _nnls_weighted(K, v, w)

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

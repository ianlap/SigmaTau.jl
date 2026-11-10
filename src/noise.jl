# Power-law noise identification (lag-1 ACF with B1 / R(n) fallback),
# calibrated power-law noise synthesis and the NNLS variance decomposition.

# Fewest decimated points the lag-1 test is run on: Riley and Greenhall
# (2004) find it reliable from about 30, and Stable32 switches to B1 below.
const NEFF_RELIABLE = 30
const B1_MIN_POINTS = 6   # fewest averaged points the averaged-length B1 test is run on

"""
    identify_noise(x, m_values; dmin=0, dmax=2, remove_drift=true,
                   b1_length=:averaged) → Vector{Symbol}
    identify_noise(data, m_values; kwargs...) → Vector{Symbol}

Dominant power-law noise type of a record at each averaging factor: `:WHPM`,
`:FLPM`, `:WHFM`, `:FLFM` or `:RWFM`, and with `dmax=3` also `:FWFM`
(flicker-walk FM, ``α = -3``) and `:RRFM` (random-run FM, ``α = -4``). The
estimators use it for the noise type that sets their EDF and bias
correction.

# Method

For each ``m`` the phase is decimated by ``m``, and the lag-1
autocorrelation ``r_1`` of the result gives ``δ = r_1/(1 + r_1)``. The
series is differenced until ``δ < 0.25`` or `dmax` differences have been
taken; after ``d`` differences, ``α = 2 - 2(δ + d)``
([Riley and Greenhall (2004)](@cite riley-2004-lag1-acf);
[SP1065 §5.5.6](@cite riley-2008-sp1065)). This package adds four steps:

- For ``m ≥ 2`` the limits between the FM labels are moved. A decimated
  record is a sampled continuous process, on which the estimate reads
  ``α = 0``, −1.446 and −2.4 for white, flicker and random-walk FM, and
  after a third difference −3.0, −3.93 and −4.57 for random-walk,
  flicker-walk and random-run FM. The limits sit halfway between those
  readings. Rounding to the nearest integer, as the paper does and as is
  done here at ``m = 1``, would split flicker FM almost evenly between its
  own label and random-walk FM.
- A white-PM result at ``m ≥ 4`` is confirmed with ``R(n)``, the ratio of
  the modified to the plain Allan variance
  ([SP1065 §5.5.2](@cite riley-2008-sp1065)), because decimation whitens
  flicker PM.
- Where fewer than 30 decimated points remain, the lag-1 test gives way to
  the ``B_1`` ratio of the standard variance to the Allan variance of the
  ``m``-averaged frequency, with ``R(n)`` to separate the PM types
  ([SP1065 §5.5.2](@cite riley-2008-sp1065)), as in
  [Stable32](@cite hamilton-2008-stable32-manual). That test stops at
  random-walk FM.
- From ``m = 16`` a PM result of either test stands only if ``R(n)`` lies
  nearer its flicker-PM value than 1/2, its value on white FM (it is 27/40
  on flicker FM and 33/40 on random-walk FM). One that fails is reported
  as white FM. ``R(n)`` uses every sample where the tests read one in
  ``m``, and with few decimated points they take FM for PM often: white FM
  in one record in four at 32 points. A PM label there makes the interval
  several times too narrow.

Samples farther than five standard deviations from the record mean are
dropped before the analysis.

The function always returns a label. An estimate outside the labeled range
counts as the nearest end of it, an averaging factor whose estimate fails
takes the label of the nearest one in ``\\log m`` that succeeded, and a
record with no usable estimate, such as a constant one, is reported as
`:WHFM`.

# Arguments

- `x`: phase samples, or `data`, a `PhaseData` or `FrequencyData` record
  (frequency is integrated to phase).
- `m_values`: averaging factors. With a record, any integer vector or a
  [`TauMode`](@ref).

# Keywords

- `dmin=0`, `dmax=2`: the fewest and the most differences taken. Use
  `dmax=2` for the Allan family and `dmax=3` for the Hadamard family, as
  the estimators do.
- `remove_drift=true`: remove a least-squares quadratic, that is a frequency
  offset and a linear frequency drift, from each decimated series before
  the test, as allantools does. `false` skips it.
- `b1_length=:averaged`: the record length ``N`` at which the expected
  ``B_1(N, μ)`` of each noise type is evaluated in the test for fewer than
  30 points. `:averaged` is the number of ``m``-averaged frequency values
  the ratio is computed from, which is the ``N`` of the definition
  ([SP1065 §5.7](@cite riley-2008-sp1065)). `:record` is the number of
  frequency samples in the whole record, which follows Stable32's labels
  more closely; at that length the limit for random-walk FM lies far above
  anything a few averaged values can reach.

!!! note "Labels from fewer than 30 points"
    In simulation (``N = 1024``, 6 to 24 averaged values) `:averaged` finds
    random-walk FM in 60 to 74 % of records and flicker FM in 33 to 79 %,
    where `:record` finds them in 1 to 2 % and 23 to 46 %. White FM pays
    for it: 62 to 96 % against 86 to 100 %. With ten values or fewer,
    flicker FM and random-walk FM get their own label in a third to two
    thirds of records, so a label at the longest averaging times is weak
    evidence.

    Below 6 averaged values `:averaged` reports the label found at the
    longest averaging factor that still has 6, because with 2 to 4 values
    the ratio is close to a coin toss. The labels still cost white FM
    something there: its bias-corrected total variance reads 4 % high with
    four averaged values left and 18 % high with two, where the uncorrected
    one is within 1 %. The value 6 comes from simulation, not from the
    literature.

!!! compat "SigmaTau 0.6"
    `remove_drift` was called `detrend` before version 0.6. `detrend=false`
    still works and gives a deprecation warning.
"""
function identify_noise(x::Vector{Float64}, m_values::Vector{Int}; dmin::Int=0, dmax::Int=2,
                        remove_drift::Bool=true, detrend::Bool=true, b1_length::Symbol=:averaged)
    b1_length in (:record, :averaged) ||
        throw(ArgumentError("b1_length must be :record or :averaged, got :$b1_length"))
    detrend || Base.depwarn("identify_noise: `detrend` is deprecated, use `remove_drift`", :identify_noise)
    detrend = remove_drift && detrend
    x_clean = _preprocess(x)
    N = length(x_clean)
    alphas = Vector{Int}(undef, length(m_values))
    ok = falses(length(m_values))
    # Drift-removed record, shared by m = 1 and R(n), and its ADEV and MDEV at
    # every m for R(n); each made on first use.
    x_rn = nothing
    rn_adev = rn_mdev = nothing

    for (k, m) in enumerate(m_values)
        N_eff = N ÷ m
        alpha = NaN

        try
            lag1 = N_eff >= NEFF_RELIABLE
            # Below B1_MIN_POINTS the label is the one at the longest m that has them.
            mr = lag1 || b1_length !== :averaged ? m : clamp(N ÷ B1_MIN_POINTS, 1, m)
            if lag1 && m == 1
                x_rn === nothing && (x_rn = detrend ? _noise_detrend(x_clean) : x_clean)
                alpha = _noise_id_lag1acf(x_rn, 1, dmin, dmax; detrend=false)[1]
            elseif lag1
                alpha = _noise_id_lag1acf(x_clean, m, dmin, dmax; detrend=detrend)[1]
            else
                alpha = _noise_id_b1rn(x_clean, mr; detrend=detrend, b1_length=b1_length)[1]
            end
            # A PM result is checked with R(n) = MVAR/AVAR (SP1065 §5.5.2, step 5),
            # which uses every sample where the tests above use one in m. From
            # m = 16 it must lie nearer its flicker-PM value than the 1/2 of white
            # FM, or the noise is not PM. And decimation whitens flicker PM, so a
            # white-PM result of the lag-1 test is confirmed with it from m = 4.
            if alpha >= 0.5 && mr >= 4
                x_rn === nothing && (x_rn = detrend ? _noise_detrend(x_clean) : x_clean)
                if rn_mdev === nothing
                    rn_adev = _adev_core(x_rn, m_values, 1.0)
                    rn_mdev = _mdev_core(x_rn, m_values, 1.0)
                end
                md, ad = mr == m ? (rn_mdev[k], rn_adev[k]) :
                                   (_mdev_core(x_rn, [mr], 1.0)[1], _adev_core(x_rn, [mr], 1.0)[1])
                if mr >= 16 && _rn_fm(md, ad, mr)
                    alpha = 0.0
                elseif lag1 && alpha >= 1.5 && _rn_flicker(md, ad, mr)
                    alpha = 1.0
                end
            end
        catch
            alpha = NaN
        end

        isnan(alpha) && continue
        alphas[k] = _lag1_label(alpha, m, dmax >= 3 ? -4 : -2)
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

# Decimating the phase by m leaves a sampled continuous process, on which the
# lag-1 estimate does not settle on the integers. After two differences flicker
# FM and random-walk FM read α = −1.446 and −2.4: their second differences
# correlate at (9 ln 3 − 16 ln 2)/(8 ln 2) and 1/4. After three, random-walk,
# flicker-walk and random-run FM read −3.0, −3.927 and −4.565, from
# −1/3, (128 ln 2 − 81 ln 3)/(27 ln 3 − 32 ln 2) and 13/33. All follow from the
# generalized autocovariances of Greenhall and Riley (2003).
_lag1_reading(r1, d) = 2 - 2 * (d + r1 / (1 + r1))
const _LAG1_FLFM = _lag1_reading((9log(3) - 16log(2)) / (8log(2)), 2)
const _LAG1_RWFM = _lag1_reading(1 / 4, 2)
const _LAG1_RWFM3 = _lag1_reading(-1 / 3, 3)
const _LAG1_FWFM = _lag1_reading((128log(2) - 81log(3)) / (27log(3) - 32log(2)), 3)
const _LAG1_RRFM = _lag1_reading(13 / 33, 3)

# Integer α of a lag-1 estimate, no lower than `amin` (−2, or −4 when a third
# difference is allowed). The estimate is continuous and noisy (white PM can
# come out at 2.6), so it is clamped to the labeled range. At m = 1 it rounds
# to the nearest integer (Riley and Greenhall (2004)); for m ≥ 2 the limits on
# the FM side sit halfway between the sampled-process readings above, at
# −0.72, −1.92, −3.46 and −4.25.
function _lag1_label(alpha::Real, m::Int, amin::Int)
    m == 1 && return round(Int, clamp(alpha, amin, 2))
    amin <= -4 && alpha < (_LAG1_FWFM + _LAG1_RRFM) / 2 && return -4
    amin <= -3 && alpha < (_LAG1_RWFM3 + _LAG1_FWFM) / 2 && return -3
    alpha < (_LAG1_FLFM + _LAG1_RWFM) / 2 && return -2
    alpha < _LAG1_FLFM / 2 && return -1
    return round(Int, clamp(alpha, 0.0, 2.0))
end

const _NOISE_SYMBOLS = (:RRFM, :FWFM, :RWFM, :FLFM, :WHFM, :FLPM, :WHPM)  # α = -4 … 2

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
        noises[k] = _NOISE_SYMBOLS[alphas[j] + 5]
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

# Drift removal for the noise identifier: the least-squares quadratic in the
# centred index t in closed form (the odd moments of t vanish), with one
# correction pass on the residual. It agrees with `_remove_quadratic` to
# rounding at about a seventh of the cost; a record with gaps or under three
# points takes that general fit.
function _noise_detrend(x::Vector{Float64})
    N = length(x)
    (N < 3 || !all(isfinite, x)) && return _remove_quadratic(x)
    mu = (N - 1) / 2
    S2 = S4 = 0.0
    @inbounds @simd for i in 1:N
        t2 = ((i - 1) - mu)^2
        S2 += t2; S4 += t2 * t2
    end
    det = N * S4 - S2^2
    r = copy(x)
    for _ in 1:2
        a = b = c = 0.0
        @inbounds @simd for i in 1:N
            t = (i - 1) - mu
            a += r[i]; b += t * r[i]; c += t * t * r[i]
        end
        c0, c1, c2 = (S4 * a - S2 * c) / det, b / S2, (N * c - S2 * a) / det
        @inbounds @simd for i in 1:N
            t = (i - 1) - mu
            r[i] -= c0 + (c1 + c2 * t) * t
        end
    end
    return r
end

function _noise_id_lag1acf(x::Vector{Float64}, m::Int, dmin::Int = 0, dmax::Int = 2; detrend::Bool=true)
    x_dec = m > 1 ? x[1:m:end] : x
    x_det = detrend ? _noise_detrend(x_dec) : copy(x_dec)

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

function _noise_id_b1rn(x::Vector{Float64}, m::Int; detrend::Bool=true, b1_length::Symbol=:averaged)
    x_dec = x[1:m:end]
    detrend && (x_dec = _noise_detrend(x_dec))

    avar_val = _adev_core(x_dec, [1], 1.0)[1]^2 / Float64(m)^2
    dx = diff(x)
    Nd = (length(dx) ÷ m) * m
    Nd < m && return (NaN, -2, NaN)

    y_blocks = reshape(dx[1:Nd], m, :)
    y_avg = vec(mean(y_blocks, dims=1))
    # Unbiased variance (N − 1), as B1 is defined.
    var_class = var(y_avg)

    # N in B1_theory(N, μ): the number of averaged points B1_obs is computed
    # from, which is what the theory describes, or with `:record` the full
    # frequency-sample count (see the `b1_length` docs of `identify_noise`).
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

    mu_best == -2 && (alpha_int = _rn_flicker(_mdev_core(x, [m], 1.0)[1], sqrt(avar_val), m) ? 1 : 2)

    return (alpha_int, mu_best, B1_obs)
end

# R(n) = MVAR/AVAR at averaging factor m separates the PM noises, which ADEV
# and the decimated phase cannot: true when it lies nearer its flicker-PM
# value than 1/m, its white-PM value.
function _rn_flicker(mdev_val::Float64, adev_val::Float64, m::Int)
    (isnan(mdev_val) || !(adev_val > 0)) && return false
    return (mdev_val / adev_val)^2 > sqrt(_rn_theory(m, 0) * _rn_theory(m, -1))
end

# And PM from FM: true when R(n) lies nearer 1/2, its white-FM value, than its
# flicker-PM value.
_rn_fm(mdev_val::Float64, adev_val::Float64, m::Int) =
    adev_val > 0 && (mdev_val / adev_val)^2 > sqrt(_rn_theory(m, -1) / 2)

# Barnes's B1(N, μ) without dead time: the expected ratio of the N-sample
# standard variance to the Allan variance when σ_y²(τ) ∝ τ^μ (SP1065 §5.5.2
# and §5.7, Eq. (48)).
function _b1_theory(N::Int, mu::Int)
    if mu == 2;  return N * (N + 1) / 6
    elseif mu == 1;  return N / 2
    elseif mu == 0;  return N * log(N) / (2 * (N - 1) * log(2))
    elseif mu == -1; return 1.0
    elseif mu == -2; return (N^2 - 1) / (1.5 * N * (N - 1))
    else;            return (N * (1 - N^mu)) / (2 * (N - 1) * (1 - 2^mu))
    end
end

# Expected R(n) = MVAR/AVAR (SP1065 §5.10, Eq. (51)) at averaging factor `af`:
# 1/af for white PM (b = 0); for flicker PM (b = −1) the ratio of the two
# conversion formulas with the bandwidth at the Nyquist frequency of τ₀.
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
# Gaussian noise through the fractional differencing filter (1 − z⁻¹)^(α/2)
# of Kasdin and Walter (1992), Eqs. (36) and (42), exact at the integer
# orders (α = 2 is the first difference of white noise, α = 0 white, α = −2 a
# cumulative sum). Zero-mean and not yet scaled; `noise_gen` sets the level.
# `rng = nothing` draws from the global RNG.
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
        # Analytic calibration of the PSD coefficient: the filter driven by
        # unit white noise sampled at τ₀ has one-sided PSD
        # 2τ₀·|2 sin(πfτ₀)|^α (Kasdin and Walter (1992), Eq. (37)), which is
        # 2τ₀·(2πτ₀f)^α at low frequency, so scaling by
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
    noise_gen(PhaseData, N, tau0; sigma1=Dict(), h=Dict(), rng=nothing) → PhaseData
    noise_gen(FrequencyData, N, tau0; sigma1=Dict(), h=Dict(),
              rng=nothing) → FrequencyData

Simulated clock record of `N` samples at sample interval `tau0` whose
fractional-frequency spectrum is a sum of power laws,
``S_y(f) = \\sum_α h_α f^α``. The first argument selects the kind of record
returned.

Each component is white Gaussian noise passed through the fractional
differencing filter ``(1 - z^{-1})^{α/2}`` of
[Kasdin and Walter (1992), Eqs. (36) and (42)](@cite kasdin-1992-discrete-simulation),
which is exact at the integer orders: a first difference for ``α = 2``, a
running sum for ``α = -2``. The components are independent.

| ``α`` | Noise type      | Label   | Slope of ``S_x(f)`` |
|------:|:----------------|:--------|:--------------------|
|     2 | white PM        | `:WHPM` | flat                |
|     1 | flicker PM      | `:FLPM` | ``f^{-1}``          |
|     0 | white FM        | `:WHFM` | ``f^{-2}``          |
|    −1 | flicker FM      | `:FLFM` | ``f^{-3}``          |
|    −2 | random-walk FM  | `:RWFM` | ``f^{-4}``          |
|    −3 | flicker-walk FM | `:FWFM` | ``f^{-5}``          |
|    −4 | random-run FM   | `:RRFM` | ``f^{-6}``          |

# Keywords

Give `sigma1` or `h`, not both. The keys are the integer ``α`` or the
label, as in `Dict(:WHFM => 1e-12)`.

- `sigma1`: the Allan deviation ``σ_y(τ_0)`` of each component, for when a
  specification quotes a deviation at one averaging time. Each component is
  rescaled so that its overlapping Allan deviation at ``τ_0`` equals the
  request exactly, in every realization.
- `h`: the coefficient ``h_α`` of each component. The filter driven by unit
  white noise has the one-sided spectrum ``2τ_0\\,|2\\sin(πfτ_0)|^α`` (their
  Eq. (37)), so the component is scaled by
  ``\\sqrt{h_α/(2τ_0(2πτ_0)^α)}`` and has
  ``S_y(f) = h_α\\,(\\sin(πfτ_0)/(πτ_0))^α``, which is ``h_α f^α`` well
  below the Nyquist frequency.
- `rng`: a random number generator, such as `Xoshiro(1)`, for a
  reproducible or an independent stream. The global generator is used by
  default.

For ``α ≤ -1`` the process is not stationary, and the spectrum above is
that of its differences: the ``d``-th difference of ``y`` has the spectrum
``(2\\sin(πfτ_0))^{2d} S_y(f)``.

The record is a sampled process, not a band-limited continuous one. The
Allan variance of an `h` component follows the usual conversion
([SP1065 §7.1, Eq. (67)](@cite riley-2008-sp1065), with the bandwidth
``f_h = 1/(2τ_0)``) at every ``τ`` for white PM and white FM, and only at
long ``τ`` for the other types. [`powerlaw_sigma`](@ref) gives the exact
deviation of such a record and states the difference.

# Examples

```jldoctest
julia> using SigmaTau, Random

julia> p = noise_gen(PhaseData, 4096, 1.0; sigma1=Dict(:WHFM => 1e-12),
                     rng=Xoshiro(1));

julia> adev(p, [1]; ci=false).dev[1] ≈ 1e-12
true
```

White FM with random-walk FM, by deviation and by coefficient:

```julia
p = noise_gen(PhaseData, 8192, 1.0; sigma1=Dict(0 => 1e-12, -2 => 1e-14))
y = noise_gen(FrequencyData, 8192, 1.0; h=Dict(0 => 2e-24, -2 => 1e-28))
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

# Nonnegative least-squares decomposition of a measured variance profile
# σ̂²(τ) into power-law components.
#
# `identify_noise` finds the one α that dominates at each τ. This fit gives
# the amount of each α over the whole τ range, from one fit to the entire
# profile. The constraint makes it well posed: a variance contribution cannot
# be negative, and without the constraint the nearly collinear power-law
# columns trade large positive and negative amplitudes.
#
# The estimators this calls are defined in deviations.jl, included after
# this file; the calls are resolved at run time.

"""
    _nnls_exact(A, b) → (a, rss)

Nonnegative least squares, minimizing ‖Aa − b‖² subject to `a ≥ 0`, by
enumerating the active sets. Each of the `2ⁿ − 1` nonempty supports is
solved as an unconstrained least-squares problem on that support and kept
if its solution is feasible; `a = 0` seeds the search.

The result is exact: the constrained optimum is the unconstrained minimizer
on its own support, so it is among the candidates, and every candidate is
feasible. The cost is exponential in the number of columns, capped at 20;
the intended use has five to seven power laws.

The restricted problems are solved with `\\` (QR), not the normal
equations: power-law columns are nearly collinear, and squaring the
condition number would lose the smallest amplitudes.
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
    best_r = sum(abs2, bf)             # the empty support is always feasible
    for mask in 1:((1 << n) - 1)
        S = [j for j in 1:n if (mask >> (j - 1)) & 1 == 1]
        AS = Af[:, S]
        aS = try
            AS \ bf
        catch
            continue                   # rank-deficient support
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

# Estimator family: `:modified` for the kernels that average the phase,
# `:unmodified` for the plain second- and third-difference kernels.
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

Exponent μ of the variance response, σ²(τ) ∝ τ^μ, to the power law `alpha`.

- Modified family (MDEV, MHDEV, MHTOTDEV): `μ = −α − 1` over the whole
  range, from −3 for white PM to 3 for random-run FM, because the phase
  average separates the two PM types.
- Unmodified family (ADEV, HDEV): `μ = −α − 1` for `α ≤ 0`. White and
  flicker PM both give μ = −2: the variance of white PM is proportional to
  τ⁻², and that of flicker PM to τ⁻² times a slowly varying `ln(2π f_h τ)`,
  which a power-law column cannot resolve. A basis for an unmodified
  estimator therefore holds at most one of α = 1 and α = 2.
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

# Smallest averaging factor of the default grid (see `noise_id_nnls`).
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

# The raw profile: `ci=false` skips the lag-1/B1 classifier and
# `correct_bias=false` MHTOTDEV's bias divisor. Both depend on
# `identify_noise`, and fitting a profile that depends on a noise label would
# make the decomposition depend on the answer it is computing.
function _nnls_deviation(deviation::Symbol, data::PhaseData, m_values::Vector{Int})
    deviation === :mdev  && return mdev(data, m_values; ci=false)
    deviation === :mhdev && return mhdev(data, m_values; ci=false)
    deviation === :hdev  && return hdev(data, m_values; ci=false)
    deviation === :mhtotdev &&
        return mhtotdev(data, m_values; ci=false, correct_bias=false)
    throw(ArgumentError("noise_id_nnls: unsupported deviation $deviation"))
end

# The χ² weights of the α ≤ −3 rows use the random-walk FM EDF, the reddest
# type MDEV's model covers. This sets the weights only; it does not enter the
# fitted amplitudes.
_nnls_edf_symbol(alpha::Int) = alpha <= -3 ? :RWFM : _nnls_symbol(alpha)

# Weighted NNLS of v ≈ K·a. The power-law columns span many orders of
# magnitude over a wide τ grid; equilibrating them to unit 2-norm keeps the
# restricted QR solves well scaled, and the amplitudes are unscaled afterward.
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
# fractional weight supplies the noise type at each τ, and the EDF model
# turns that into ν, with `Ns[k]` the record length at the k-th τ. Taking the
# types from the fit, not from `identify_noise`, keeps the fit independent of
# that test. A point where the EDF model gives no finite positive ν keeps the
# fractional weight.
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
    noise_id_nnls(data, m_values; deviation=:mhdev, alphas, weights=:edf) → NamedTuple

Amount of each power-law noise type in a record, from one nonnegative
least-squares fit of the measured variance to
``σ^2(τ) ≈ \\sum_i a_i\\,(τ/τ_0)^{μ(α_i)}`` over the whole range of ``τ``.
Where [`identify_noise`](@ref) names the dominant type at each ``τ``, this
fit gives an amplitude for every type in the basis. It runs on the
deviation without bias correction or noise identification, so it does not
depend on that test.

For coefficients ``h_α`` use [`powerlaw_fit`](@ref), which fits the exact
response of each component. The basis here is the asymptotic power law
``τ^μ`` of each type, and the total estimators bend away from it as ``τ``
approaches the record length, so prefer `:mhdev` when the longest averaging
times matter.

# Arguments

- `data`: a `PhaseData` or `FrequencyData` record.
- `m_values`: averaging factors. The default is the octave grid of the
  estimator from ``m = 4``, because at ``m = 1`` and 2 the estimators have
  not yet reached their power laws.

# Keywords

- `deviation=:mhdev`: the estimator fitted, one of `:mdev`, `:hdev`,
  `:mhdev` and `:mhtotdev`.
- `alphas`: the exponents in the basis. The default is every type the
  estimator converges for: 2 to −2 for `:mdev` and 2 to −4 for the Hadamard
  family, which is also insensitive to drift. `:hdev` cannot tell white PM
  from flicker PM and takes one of them, white PM by default.
- `weights=:edf`: `:edf` weights each point by ``ν/(2σ^4)``, with the EDF
  ``ν`` for the type that dominates a first fit there; `:relative` gives
  every point the same fractional weight; a vector gives the weights
  directly.

# Returns

A named tuple with `deviation`, `tau` (the averaging times fitted), `var`
(the measured variance), `fit` (the model at each ``τ``), `alphas`,
`noise_types`, `amplitudes` (``a_i ≥ 0``, the variance of each component at
``τ_0``), `contributions` (a matrix, ``τ`` by type), `dominant` (the
largest contribution at each ``τ``), `residual` and `weights`.

# Examples

```julia
p = noise_gen(PhaseData, 8192, 1.0; sigma1=Dict(0 => 1e-12, -2 => 1e-14))
r = noise_id_nnls(p)
collect(zip(r.noise_types, r.amplitudes))   # amplitude of each type
r.dominant                                  # dominant type at each τ
```
"""
function noise_id_nnls(data::PhaseData, m_values::AbstractVector{<:Integer} = Int[];
                       deviation::Symbol = :mhdev,
                       alphas::AbstractVector{<:Integer} = _nnls_default_alphas(deviation),
                       weights::Union{Symbol, AbstractVector{<:Real}} = :edf)
    x  = _f64(data.x)
    τ₀ = data.tau0
    # The estimator's octave grid from m = NNLS_M_MIN; see the docstring.
    ms = isempty(m_values) ? filter(>=(NNLS_M_MIN), _default_m_values(length(x), deviation)) :
                             collect(Int, m_values)
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

function noise_id_nnls(data::FrequencyData, m_values::AbstractVector{<:Integer} = Int[]; kwargs...)
    return noise_id_nnls(_freq_to_phase(data), m_values; kwargs...)
end

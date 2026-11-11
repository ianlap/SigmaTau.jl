# Barnes's bias functions B1, B2 and B3 for frequency averages taken with dead
# time between them (NBS TN 375, 1969; NIST TN 1318, 1990), and the dead-time
# correction of the Allan deviation built on them (SP1065 §5.15).
#
# With F(A) = 2A^(μ+2) − (A+1)^(μ+2) − |A−1|^(μ+2) of TN 1318 Eq. (10), the
# expected squared difference of two averages over τ that start A·τ apart is
# proportional to τ^μ·(2 + F(A)), and each bias function is a ratio of sums of
# such terms. Because 2A² − (A+1)² − (A−1)² = −2,
#
#   (2 + F(A))/μ = 2A²·e(A) − (A+1)²·e(A+1) − (A−1)²·e(|A−1|),  e(t) = (t^μ − 1)/μ,
#
# and e(t) → ln t at μ = 0. The functions below work with this quotient: the
# factor 1/μ cancels in every ratio, the 0/0 that the published forms have at
# μ = 0 does not arise, and nothing is lost to cancellation near it.

# e(t) = (t^μ − 1)/μ.
_dt_e(t::Float64, μ::Float64) = μ == 0 ? log(t) : expm1(μ * log(t)) / μ

# t²·e(t). At t = 0 the term is 0 for every μ, 0^(μ+2) = 0 at μ = −2 included
# (TN 375 p. 6, which gives B2(1, −2) = 1 and 2/3 for any other r > 0).
_dt_term(t::Float64, μ::Float64) = t == 0 ? 0.0 : t^2 * _dt_e(t, μ)

# (2 + F(A))/μ. A lag within rounding of 1 is the adjacent one, where white PM
# (μ = −2) is discontinuous. From A = 32 the three terms cancel to less than a
# thousandth of their size, and the function is summed instead as its series
# in 1/A², −2e(A) − A^μ·[3 + μ + 2·Σ c_k A^(−2(k−1))] over k = 2 … 5 with
# c_k = C(μ+2, 2k)/μ, which is exact to rounding there.
function _dt_g(A::Float64, μ::Float64)
    if A >= 32
        u = 1 / A^2
        c = (μ + 2) * (μ + 1) * (μ - 1) / 24 * u
        s = c * (1 + (μ - 2) * (μ - 3) / 30 * u * (1 + (μ - 4) * (μ - 5) / 56 * u * (1 + (μ - 6) * (μ - 7) / 90 * u)))
        return -2 * _dt_e(A, μ) - (3 + μ + 2s) * A^μ
    end
    d = abs(A - 1)
    return 2 * _dt_term(A, μ) - _dt_term(A + 1, μ) - _dt_term(d < 8eps() ? 0.0 : d, μ)
end

# Σ (M − |d|)·g((kM + d)r) over d = −(M−1) … M−1, less 2·Σ (M − n)·g(nr) over
# n = 1 … M−1: M² times the expected squared difference of two averages of M
# measurements that start kM measurement periods apart (TN 1318 Appendix,
# where k = 1 gives the numerator of Eq. (13)).
function _dt_distributed(k::Int, M::Int, r::Float64, μ::Float64)
    s = M * _dt_g(k * M * r, μ)
    for n in 1:M-1
        s += (M - n) * (_dt_g((k * M + n) * r, μ) + _dt_g((k * M - n) * r, μ) - 2 * _dt_g(n * r, μ))
    end
    return s
end

function _dt_args(f::Symbol, r::Real, μ::Real; zero_ok::Bool=false)
    -2 <= μ <= 2 || throw(ArgumentError("$f: μ must lie in −2 … 2, got $μ"))
    (isfinite(r) && (r > 0 || (zero_ok && r == 0))) || throw(ArgumentError(
        "$f: the ratio r = T/τ must be finite and " * (zero_ok ? "≥ 0" : "positive") * ", got $r"))
    return Float64(r), Float64(μ)
end

"""
    _noise_mu(noise::Symbol) → Float64

Exponent ``μ`` of ``σ_y^2(τ) ∝ τ^μ`` for a noise type the Allan variance
converges for: −2 for white and flicker PM, −1, 0 and 1 for white, flicker
and random-walk FM (TN 375, Fig. 1).
"""
function _noise_mu(noise::Symbol)
    α = get(_NOISE_ALPHA, noise, -9)
    -2 <= α <= 2 || throw(ArgumentError(
        "the dead-time bias functions cover :WHPM, :FLPM, :WHFM, :FLFM and :RWFM; got :$noise"))
    return Float64(max(-α - 1, -2))
end

"""
    bias_b1(N, r, μ) → Float64

Barnes's bias function ``B_1``: the expected ``N``-sample variance of
frequency averages over the expected two-sample variance of the same
averages,

```math
B_1(N, r, μ) = \\frac{⟨σ^2(N, T, τ)⟩}{⟨σ^2(2, T, τ)⟩}
             = \\frac{1 + \\sum_{n=1}^{N-1} \\frac{N-n}{N(N-1)} F(nr)}{1 + \\frac{1}{2}F(r)} ,
\\qquad F(A) = 2A^{μ+2} - (A+1)^{μ+2} - |A-1|^{μ+2} ,
```

for averages over ``τ`` that start every ``T = rτ``, on noise whose Allan
variance goes as ``τ^μ``
([Barnes (1969), Eq. (8)](@cite barnes-1969-bias-functions);
[Barnes and Allan (1990), Eqs. (10) and (11)](@cite barnes-1990-dead-time)).
``σ^2(N, T, τ)`` is the sample variance with divisor ``N - 1``. At
``r = 1`` there is no dead time and ``B_1 = N(1 - N^μ)/[2(N-1)(1 - 2^μ)]``,
the ratio of the standard variance to the Allan variance that
[`identify_noise`](@ref) tests
([SP1065 §5.7, Eq. (48)](@cite riley-2008-sp1065)).

# Arguments
- `N`: the number of averages, an integer ``≥ 2``, or `Inf` for the limit of
  a long record: ``1/[1 + F(r)/2]`` for ``μ < 0``, and infinite for
  ``μ ≥ 0``, where the standard variance does not converge.
- `r`: the ratio ``T/τ > 0`` of the measurement period to the averaging
  time; ``r > 1`` is dead time, ``r < 1`` overlapping averages.
- `μ`: the exponent of ``σ_y^2(τ) ∝ τ^μ``, any real number in ``-2 ≤ μ ≤ 2``,
  or a noise type: `:WHPM` and `:FLPM` (``μ = -2``), `:WHFM` (−1), `:FLFM`
  (0) or `:RWFM` (1). In terms of ``S_y(f) ∝ f^α``, ``μ = -α - 1`` for
  ``-3 < α ≤ 1``.

The expressions are evaluated in a form that is finite at ``μ = 0``, where
the published ones read 0/0.

!!! note "White and flicker PM"
    At ``μ = -2`` the functions are those of white PM sampled without a
    bandwidth limit, with ``0^0 = 0`` as in the tables. Flicker PM has the
    same ``μ`` and a variance that grows with the logarithm of the
    measurement bandwidth, which the functions leave out.

# Examples

```jldoctest
julia> using SigmaTau

julia> bias_b1(4, 1, :RWFM)          # N/2 without dead time
2.0

julia> round(bias_b1(8, 2, 0.2); sigdigits=4)
1.629
```
"""
function bias_b1(N::Real, r::Real, μ::Real)
    r, μ = _dt_args(:bias_b1, r, μ)
    g1 = _dt_g(r, μ)
    if N == Inf
        return μ < 0 ? 2 / (μ * g1) : Inf
    end
    (isinteger(N) && N >= 2) || throw(ArgumentError("bias_b1: N must be an integer ≥ 2 or Inf, got $N"))
    n_avg = Int(N)
    s = 0.0
    for n in 1:n_avg-1
        s += (n_avg - n) * _dt_g(n * r, μ)
    end
    return 2s / (n_avg * (n_avg - 1) * g1)
end

"""
    bias_b2(r, μ) → Float64

Barnes's bias function ``B_2``: the expected two-sample variance of
frequency averages with dead time over the Allan variance at the same
averaging time,

```math
B_2(r, μ) = \\frac{⟨σ^2(2, T, τ)⟩}{⟨σ^2(2, τ, τ)⟩}
          = \\frac{1 + \\frac{1}{2}F(r)}{2(1 - 2^μ)} ,
```

with ``r = T/τ`` and ``F`` as in [`bias_b1`](@ref)
([Barnes (1969), Eq. (10)](@cite barnes-1969-bias-functions);
[Barnes and Allan (1990), Eq. (12)](@cite barnes-1990-dead-time);
[SP1065 §5.8, Eq. (49)](@cite riley-2008-sp1065)). Dividing a two-sample
variance measured with the dead time ``T - τ`` by ``B_2`` gives the Allan
variance at ``τ``.

``B_2(1, μ) = 1``. For ``r ≥ 1`` the function is ``r^2`` at ``μ = 2``,
``(3r - 1)/2`` for random-walk FM, 1 for white FM and 2/3 for white PM at
any ``r > 1``. `r` may be any ratio ``≥ 0``, and `μ` is an exponent or a
noise type, as in [`bias_b1`](@ref).

# Examples

```jldoctest
julia> using SigmaTau

julia> bias_b2(36, :RWFM)            # 100 s measured once an hour
53.5

julia> round(sqrt(bias_b2(36, :FLFM)); digits=3)
1.915
```
"""
function bias_b2(r::Real, μ::Real)
    r, μ = _dt_args(:bias_b2, r, μ; zero_ok=true)
    return _dt_g(r, μ) / _dt_g(1.0, μ)
end

"""
    bias_b3(N, M, r, μ) → Float64

Barnes and Allan's bias function ``B_3``: the expected ``N``-sample variance
of averages of ``M`` measurements, with the dead time spread among them,
over that of single measurements ``M`` times as long with the dead time
gathered at the end,

```math
B_3(N, M, r, μ) = \\frac{⟨σ^2(N, M, T, τ)⟩}{⟨σ^2(N, T, τ)⟩} ,
\\qquad τ = Mτ_0 , \\quad T = MT_0 , \\quad r = T_0/τ_0 ,
```

where the measurements average over ``τ_0`` and start every ``T_0``
([Barnes and Allan (1990), Eq. (7)](@cite barnes-1990-dead-time);
[SP1065 §5.9, Eq. (50)](@cite riley-2008-sp1065)). For ``N = 2``,

```math
B_3(2, M, r, μ) = \\frac{2M + M F(Mr) - \\sum_{n=1}^{M-1} (M-n)
    \\left[2F(nr) - F((M+n)r) - F((M-n)r)\\right]}{M^{μ+2}\\,[F(r) + 2]} ,
```

with ``F`` as in [`bias_b1`](@ref)
([Barnes and Allan (1990), Eq. (13)](@cite barnes-1990-dead-time)).

The product ``B_2(r, μ)\\,B_3(2, M, r, μ)`` is the ratio of the two-sample
variance of ``M``-averages of dead-time data to the Allan variance at
``τ = Mτ_0``; [`adev`](@ref) divides by it under `dead_time`.

``B_3 = 1`` at ``M = 1``, at ``r = 1``, at ``μ = 2`` and, for ``r ≥ 1``, on
white FM; on white PM with dead time ``B_3(2, M, r, -2) = M``. `N` and `M`
are integers, ``N ≥ 2`` and ``M ≥ 1``; `r` and `μ` are as in
[`bias_b1`](@ref).

!!! note "More than two samples"
    The reference gives the expression and its tables for ``N = 2``. For
    ``N > 2`` the function evaluates the definition by the same steps: the
    ``N``-sample variance is the mean squared difference of every pair of
    averages, and the pair ``k`` averages apart contributes the two-sample
    term with its lags moved to ``(kM ± n)r``.

# Examples

```jldoctest
julia> using SigmaTau

julia> bias_b3(2, 8, 36, :WHPM)
8.0

julia> round(bias_b3(2, 4, 2, :FLFM); sigdigits=4)
0.7052
```
"""
function bias_b3(N::Integer, M::Integer, r::Real, μ::Real)
    N >= 2 || throw(ArgumentError("bias_b3: N must be ≥ 2, got $N"))
    M >= 1 || throw(ArgumentError("bias_b3: M must be ≥ 1, got $M"))
    r, μ = _dt_args(:bias_b3, r, μ)
    num = den = 0.0
    for k in 1:Int(N)-1
        num += (N - k) * _dt_distributed(k, Int(M), r, μ)
        den += (N - k) * _dt_g(k * r, μ)
    end
    return num / (Float64(M)^(μ + 2) * den)
end

bias_b1(N::Real, r::Real, noise::Symbol) = bias_b1(N, r, _noise_mu(noise))
bias_b2(r::Real, noise::Symbol) = bias_b2(r, _noise_mu(noise))
bias_b3(N::Integer, M::Integer, r::Real, noise::Symbol) = bias_b3(N, M, r, _noise_mu(noise))

# ── The Allan deviation of averages taken with dead time ──────────────────

# B2(r, μ)·B3(2, m, r, μ): the two-sample variance of m-averages of frequency
# measurements with dead time between them, over the Allan variance at
# τ = m·τ₀ (TN 1318 §5; SP1065 §5.15).
_dead_time_bias(m::Int, r::Float64, μ::Float64) = bias_b2(r, μ) * bias_b3(2, m, r, μ)

# ADEV of a frequency record whose samples average over τ₀ and start every
# τ₀ + dead_time. The record summed to phase with step τ₀ gives the two-sample
# deviation of the m-averages, which `_dead_time_scaled` turns into σ_y(m·τ₀).
# The noise type cannot be identified on such a record, whose σ(τ) slope is
# not that of the noise, so it has to be stated.
function _dead_time_adev(data::FrequencyData, m_values::Vector{Int}, dead_time::Float64;
                         ci::Bool=true, overlapping::Bool=true, noise::_NoiseChoice=:auto,
                         confidence::Float64=DEFAULT_CONFIDENCE)
    noise === :auto && throw(ArgumentError(
        "adev: dead_time needs the noise type stated with noise=; it cannot be identified on a record with dead time"))
    _has_gaps(data) && throw(ArgumentError("adev: dead_time needs a record without NaN gaps"))
    res = adev(_freq_to_phase(data), m_values; ci=ci, overlapping=overlapping, noise=noise, confidence=confidence)
    return _dead_time_scaled(res, m_values, dead_time, noise)
end

# `res`, the ADEV of dead-time readings at `m_values`, over √(B2·B3) for the
# stated noise type at each τ. The interval is scaled with the deviation.
function _dead_time_scaled(res::StabilityResult, m_values::Vector{Int}, dead_time::Float64, noise::_NoiseChoice)
    r = 1 + dead_time / res.tau0
    (isfinite(r) && r > 0) || throw(ArgumentError(
        "adev: dead_time must be finite and above −tau0 = $(-res.tau0), got $dead_time"))
    noises = _stated_noises(noise, :adev, res.tau, res.tau0)
    scale = Float64[1 / sqrt(_dead_time_bias(m, r, _noise_mu(n))) for (m, n) in zip(m_values, noises)]
    return StabilityResult(:adev, res.tau, res.dev .* scale, res.noise_type, _scale_ci(res.ci, scale),
                           res.edf, res.neff, res.tau0, res.n, res.confidence)
end

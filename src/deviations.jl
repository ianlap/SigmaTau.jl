# Public deviation API: PhaseData / FrequencyData in, StabilityResult out.
# Wraps the kernels in kernels.jl with noise identification, EDF and
# confidence intervals from noise.jl and edf.jl. The methods here take phase
# data and a `Vector{Int}` of averaging factors; registry.jl adds the other
# call forms (frequency data, a `TauMode` grid or none, any integer vector).

# Every estimator builds its result through here so the record metadata
# (τ₀, sample count, confidence level) rides along.
_result(data::PhaseData, confidence::Real, dt::Symbol, tau, dev, noise, ci, edf, neff) =
    StabilityResult(dt, tau, dev, noise, ci, edf, neff, data.tau0, length(data.x),
                    isempty(ci) ? NaN : Float64(confidence))

# Averaging factors count from 1; the kernels index with them unchecked.
_check_m(est::Symbol, m_values::Vector{Int}) = all(>=(1), m_values) ||
    throw(ArgumentError("$est: averaging factors must be ≥ 1, got $(minimum(m_values))"))

# Noise types for the Hadamard family. With `remove_drift` they are identified
# on the record less its frequency drift (a least-squares quadratic in phase):
# these statistics ignore drift, so the noise type that sets their EDF and
# bias correction should too (SP1065 §5.5.6 removes deterministic components
# before noise identification). The B1 test is not used on that record: its
# limits are for a record with the drift left in, where the ratio is about
# twice as large, and it stops at random-walk FM.
_hadamard_noise(x::Vector{Float64}, m_values::Vector{Int}, remove_drift::Bool) =
    remove_drift ? identify_noise(_noise_detrend(x), m_values, dmin=0, dmax=3, b1=false) :
                   identify_noise(x, m_values, dmin=0, dmax=3)

# What `noise=` takes: `:auto`, one noise type, one per τ, or a fitted model.
const _NoiseChoice = Union{Symbol,AbstractVector{Symbol},PowerLawFit}

# Reddest α the kernel of an estimator converges for.
_kernel_amin(kernel::Symbol) = kernel === :adev || kernel === :mdev ? -2 : -4

# A stated noise type in `kernel`; one redder than the kernel covers counts as
# its reddest, as in identification.
function _stated_noise(n::Symbol, kernel::Symbol)
    n in _NOISE_SYMBOLS || throw(ArgumentError(
        "noise must be :auto, a noise type (:WHPM, :FLPM, :WHFM, :FLFM, :RWFM, :FWFM, :RRFM), " *
        "a vector of them or a PowerLawFit; got :$n"))
    return _NOISE_SYMBOLS[max(_NOISE_ALPHA[n], _kernel_amin(kernel)) + 5]
end

# The noise type at each τ from a `noise` other than `:auto`. `kernel` (:adev,
# :mdev, :hdev or :mhdev) is the estimator's kind.
_stated_noises(noise::Symbol, kernel::Symbol, taus, tau0::Float64) =
    fill(_stated_noise(noise, kernel), length(taus))

function _stated_noises(noise::AbstractVector{Symbol}, kernel::Symbol, taus, tau0::Float64)
    length(noise) == length(taus) || throw(ArgumentError(
        "noise has $(length(noise)) entries for $(length(taus)) averaging factors"))
    return Symbol[_stated_noise(n, kernel) for n in noise]
end

# The component of a power-law model that gives the most variance at each τ.
function _stated_noises(fit::PowerLawFit, kernel::Symbol, taus, tau0::Float64)
    amin = _kernel_amin(kernel)
    noises = Vector{Symbol}(undef, length(taus))
    for (i, τ) in enumerate(taus)
        m = max(1, round(Int, τ / tau0))
        best, vbest = 0, -Inf
        for (α, h) in zip(fit.alphas, fit.h)
            # A term the Allan kernels diverge for is weighed in the Hadamard kernel of the same kind.
            k = α >= amin ? kernel : kernel === :adev ? :hdev : :mhdev
            v = h * _powerlaw_var(k, α, m, tau0)
            v > vbest && ((best, vbest) = (α, v))
        end
        noises[i] = _NOISE_SYMBOLS[max(best, amin) + 5]
    end
    return noises
end

# EDF of the non-overlapped ADEV (d = 2) and HDEV (d = 3): the algorithm of
# Greenhall and Riley (2003) with stride factor S = 1 in place of the
# overlapped S = m (their §3.1).
_nonoverlapped_edf(d::Int, noises::Vector{Symbol}, m_values::Vector{Int}, N::Int) =
    Float64[_calc_edf_core(_alpha_from_noise(n), d, m, m, 1, N) for (n, m) in zip(noises, m_values)]

# Terms the non-overlapped kernel sums: i = 1, 1 + m, … ≤ N − d·m, 0 below two.
_nonoverlapped_counts(d::Int, N::Int, m_values::Vector{Int}) =
    Int[(n = N - d * m >= 1 ? (N - d * m - 1) ÷ m + 1 : 0; n >= 2 ? n : 0) for m in m_values]

# The phase samples, for an estimator that cannot skip gaps.
function _gapfree(est::Symbol, data::PhaseData, m_values::Vector{Int})
    _check_m(est, m_values)
    _has_gaps(data) && throw(ArgumentError(
        "$est: the record has NaN gaps, which this estimator cannot skip; fill them with fillgaps first"))
    return _f64(data.x)
end

"""
    adev(data, m_values; ci=true, overlapping=true, noise=:auto,
         confidence=DEFAULT_CONFIDENCE) → StabilityResult

Overlapping Allan deviation ``σ_y(τ)`` (ADEV), the standard time-domain
measure of frequency stability.

For ``N`` phase samples ``x_i`` at sample interval ``τ_0`` and averaging time
``τ = mτ_0``,

```math
σ_y^2(τ) = \\frac{1}{2τ^2(N-2m)} \\sum_{i=1}^{N-2m} \\left(x_{i+2m} - 2x_{i+m} + x_i\\right)^2 .
```

This is the estimator of
[IEEE Std 1139-2022, Eq. (A.23)](@cite ieee1139-2022-definitions) and
[SP1065 §5.2.4, Eq. (11)](@cite riley-2008-sp1065) for the two-sample
variance of [Allan (1966)](@cite allan-1966-statistics). The standard
recommends it over the non-overlapped form, which has the same expectation
and fewer degrees of freedom.

# Arguments
- `data`: a [`PhaseData`](@ref) or [`FrequencyData`](@ref) record. Frequency
  data is integrated to phase first.
- `m_values`: the averaging factors ``m ≥ 1``, as an integer vector or range;
  a [`TauMode`](@ref) grid such as `Decade`; or nothing, for octave spacing
  up to the largest ``m`` the estimator supports. A factor that leaves fewer
  than two terms in the sum gives `NaN`.

# Keywords
- `ci=true`: identify the noise type at each ``τ`` and report the
  equivalent degrees of freedom (EDF) and a χ² confidence interval. The EDF
  is that of the full algorithm of
  [Greenhall and Riley (2003), §3.6](@cite greenhall-2003-edf-stability).
  `false` computes the deviation alone.
- `overlapping=true`: `false` selects the non-overlapped estimator of
  [SP1065 §5.2.2, Eq. (7)](@cite riley-2008-sp1065), one second difference
  every ``m`` samples (Stable32's "Normal" Allan deviation, allantools'
  `adev`). Its EDF comes from the same algorithm with stride 1, and it needs
  a record without gaps.
- `noise=:auto`: the power-law noise type that sets the EDF. `:auto`
  identifies it at each ``τ`` with [`identify_noise`](@ref). A symbol
  (`:WHPM`, `:FLPM`, `:WHFM`, `:FLFM`, `:RWFM`, `:FWFM` or `:RRFM`) states one
  type for every ``τ``, a vector states one per averaging factor, and a
  [`PowerLawFit`](@ref) gives its strongest component at each ``τ``. A type
  redder than the estimator converges for counts as the reddest it does.
  State the noise when the identified types contradict what is known of the
  clock: the test reads each ``τ`` on its own, and a slow periodic or
  wandering term can mislead it where few averaged points remain.
- `confidence=DEFAULT_CONFIDENCE`: confidence level of the interval.

Missing samples (`NaN`) are skipped: a term counts when every sample it
reads was measured, and the EDF follows the number of terms summed (see
[Timestamps and gaps](@ref)).

The other estimators take `data` and `m_values` in the same forms, and
those that report an interval share `ci`, `noise` and `confidence`.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData(collect(1.0:10.0), 1.0);

julia> r = adev(p, [1, 2]; ci=false);

julia> round.(r.dev; sigdigits=4)
2-element Vector{Float64}:
 0.0
 0.0
```
"""
function adev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, overlapping::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    _has_gaps(data) && return _gap_result(:adev, _gap_record(data), m_values; ci=ci, overlapping=overlapping, noise=noise, confidence=confidence)
    _check_m(:adev, m_values)
    x = _f64(data.x)
    raw_devs = _adev_core(x, m_values, data.tau0; overlapping=overlapping)
    taus = m_values .* data.tau0
    neff = overlapping ? _neff_counts(:adev, length(x), m_values) : _nonoverlapped_counts(2, length(x), m_values)

    if !ci
        return _result(data, confidence, :adev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = noise === :auto ? identify_noise(x, m_values, dmin=0, dmax=2) :
                               _stated_noises(noise, :adev, taus, data.tau0)
    edfs = overlapping ? calculate_edf(:adev, raw_devs, noises, m_values, taus, length(x), (length(x) - 1) * data.tau0) :
                         _nonoverlapped_edf(2, noises, m_values, length(x))
    lower, upper = confidence_intervals(raw_devs, edfs, confidence)

    return _result(data, confidence, :adev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    mdev(data, m_values; ci=true, noise=:auto,
         confidence=DEFAULT_CONFIDENCE) → StabilityResult

Modified Allan deviation (MDEV): the second difference of phase averaged
over ``m`` adjacent samples.

```math
\\mathrm{Mod}\\,σ_y^2(τ) = \\frac{1}{2m^4τ_0^2(N-3m+1)} \\sum_{j=1}^{N-3m+1}
    \\left[\\sum_{i=j}^{j+m-1} \\left(x_{i+2m} - 2x_{i+m} + x_i\\right)\\right]^2
```

This is the estimator of
[IEEE Std 1139-2022, Eq. (A.27)](@cite ieee1139-2022-definitions) and
[SP1065 §5.2.5, Eq. (14)](@cite riley-2008-sp1065), introduced by
[Allan and Barnes (1981)](@cite allan-1981-modified). The averaging narrows
the measurement bandwidth as ``τ`` grows, so white PM falls as ``τ^{-3/2}``
and flicker PM as ``τ^{-1}``, two noise types that [`adev`](@ref) cannot
tell apart.

Arguments and keywords are those of [`adev`](@ref), less `overlapping`.
Missing samples are skipped in the same way.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData(collect(1.0:10.0), 1.0);

julia> r = mdev(p, [1, 2]; ci=false);

julia> round.(r.dev; sigdigits=4)
2-element Vector{Float64}:
 0.0
 0.0
```
"""
function mdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    _has_gaps(data) && return _gap_result(:mdev, _gap_record(data), m_values; ci=ci, noise=noise, confidence=confidence)
    _check_m(:mdev, m_values)
    x = _f64(data.x)
    raw_devs = _mdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:mdev, length(x), m_values)

    if !ci
        return _result(data, confidence, :mdev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = noise === :auto ? identify_noise(x, m_values, dmin=0, dmax=2) :
                               _stated_noises(noise, :mdev, taus, data.tau0)
    edfs = calculate_edf(:mdev, raw_devs, noises, m_values, taus, length(x), (length(x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, confidence)

    return _result(data, confidence, :mdev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    tdev(data, m_values; ci=true, noise=:auto,
         confidence=DEFAULT_CONFIDENCE) → StabilityResult

Time deviation ``σ_x(τ)`` (TDEV), in seconds: the time-stability form of
[`mdev`](@ref).

```math
σ_x(τ) = \\frac{τ}{\\sqrt{3}}\\,\\mathrm{Mod}\\,σ_y(τ)
```

([IEEE Std 1139-2022, Eq. (A.29)](@cite ieee1139-2022-definitions);
[SP1065 §5.2.6, Eq. (15)](@cite riley-2008-sp1065)). With this factor, TDEV
on white PM is the standard deviation of the ``m``-sample mean of the phase.
It is the usual measure for time transfer and network timing. For a clock
with frequency drift use [`htdev`](@ref).

Arguments and keywords are those of [`mdev`](@ref). The interval is MDEV's
scaled by the same factor, with the same EDF. See
[Choosing an estimator](@ref) for when to use it.
"""
function tdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    res = mdev(data, m_values; ci=ci, noise=noise, confidence=confidence)
    factor = res.tau ./ sqrt(3.0)

    if !ci
        return _result(data, confidence, :tdev, res.tau, res.dev .* factor, Symbol[], _CIBound[], Float64[], res.neff)
    end

    return _result(data, confidence, :tdev, res.tau, res.dev .* factor, res.noise_type,
                           _scale_ci(res.ci, factor), res.edf, res.neff)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData;
# a record with gaps goes to the gap-tolerant path (src/gaps.jl).


"""
    hdev(data, m_values; ci=true, remove_drift=true, overlapping=true,
         noise=:auto, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Overlapping Hadamard deviation (HDEV): a three-sample deviation built on the
third difference of phase.

```math
\\mathrm{H}σ_y^2(τ) = \\frac{1}{6τ^2(N-3m)} \\sum_{i=1}^{N-3m} \\left(x_{i+3m} - 3x_{i+2m} + 3x_{i+m} - x_i\\right)^2
```

This is the estimator of
[SP1065 §5.2.9, Eq. (20)](@cite riley-2008-sp1065);
[IEEE Std 1139-2022, Eq. (A.31)](@cite ieee1139-2022-definitions) gives the
same variance in its non-overlapped form. The third difference cancels
linear frequency drift, which [`adev`](@ref) reads as a term proportional to
``τ``, and it converges for noise as red as random-run FM (``α = -4``). Like
`adev` it cannot tell white PM from flicker PM; [`mhdev`](@ref) can.

Arguments and the keywords `ci`, `noise` and `confidence` are those of
[`adev`](@ref). Missing samples are skipped in the same way. See
[Choosing an estimator](@ref) for when to use it.

# Keywords
- `remove_drift=true`: identify the noise type on the record less its
  least-squares frequency drift (a quadratic in phase), so that a drift the
  statistic ignores cannot change the EDF. The deviation is the same either
  way. On that record the ``B_1`` ratio of [`identify_noise`](@ref) is not
  used: where fewer than 30 averaged points remain, an FM type is the one
  found at the longest ``τ`` that has 30, and a change from one FM type to
  another beyond that ``τ`` has to be stated with `noise`. `false`
  identifies the noise on the record as given, as Stable32 does.
- `overlapping=true`: `false` selects the non-overlapped estimator of
  [SP1065 §5.2.8, Eq. (18)](@cite riley-2008-sp1065), one third difference
  every ``m`` samples (Stable32's "Normal" Hadamard deviation, allantools'
  `hdev`), under the conditions given for [`adev`](@ref).
"""
function hdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, remove_drift::Bool=true, overlapping::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    _has_gaps(data) && return _gap_result(:hdev, _gap_record(data), m_values; ci=ci, remove_drift=remove_drift, overlapping=overlapping, noise=noise, confidence=confidence)
    _check_m(:hdev, m_values)
    x = _f64(data.x)
    raw_devs = _hdev_core(x, m_values, data.tau0; overlapping=overlapping)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = overlapping ? _neff_counts(:hdev, length(x), m_values) : _nonoverlapped_counts(3, length(x), m_values)

    if !ci
        return _result(data, confidence, :hdev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = noise === :auto ? _hadamard_noise(x, m_values, remove_drift) :
                               _stated_noises(noise, :hdev, taus, data.tau0)
    edfs = overlapping ? calculate_edf(:hdev, raw_devs, noises, m_values, taus, length(x), T) :
                         _nonoverlapped_edf(3, noises, m_values, length(x))
    lower, upper = confidence_intervals(raw_devs, edfs, confidence)

    return _result(data, confidence, :hdev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    mhdev(data, m_values; ci=true, remove_drift=true, noise=:auto,
          confidence=DEFAULT_CONFIDENCE) → StabilityResult

Modified Hadamard deviation (MHDEV): the third difference of phase averaged
over ``m`` adjacent samples.

```math
\\mathrm{Mod}\\,σ_H^2(τ) = \\frac{1}{6m^4τ_0^2(N-4m+1)} \\sum_{j=1}^{N-4m+1}
    \\left[\\sum_{i=j}^{j+m-1} \\left(x_{i+3m} - 3x_{i+2m} + 3x_{i+m} - x_i\\right)\\right]^2
```

This is the estimator of
[SP1065 §5.2.10, Eq. (23)](@cite riley-2008-sp1065), the third-difference
case of the modified Hadamard variance of
[Bregni and Jmoda (2006), Eq. (6)](@cite bregni-2006-mhvar). It ignores
linear frequency drift as [`hdev`](@ref) does and separates white PM from
flicker PM as [`mdev`](@ref) does. The EDF is that of the
[Greenhall and Riley (2003)](@cite greenhall-2003-edf-stability) algorithm
for a modified third-difference variance.

Arguments and keywords are those of [`hdev`](@ref), less `overlapping`.
Missing samples are skipped in the same way. See
[Choosing an estimator](@ref) for when to use it.
"""
function mhdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, remove_drift::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    _has_gaps(data) && return _gap_result(:mhdev, _gap_record(data), m_values; ci=ci, remove_drift=remove_drift, noise=noise, confidence=confidence)
    _check_m(:mhdev, m_values)
    x = _f64(data.x)
    raw_devs = _mhdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:mhdev, length(x), m_values)

    if !ci
        return _result(data, confidence, :mhdev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = noise === :auto ? _hadamard_noise(x, m_values, remove_drift) :
                               _stated_noises(noise, :mhdev, taus, data.tau0)
    edfs = calculate_edf(:mhdev, raw_devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, confidence)

    return _result(data, confidence, :mhdev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    htdev(data, m_values; ci=true, remove_drift=true, noise=:auto,
          confidence=DEFAULT_CONFIDENCE) → StabilityResult

Hadamard time deviation (HTDEV), in seconds: the time-stability form of
[`mhdev`](@ref), as [`tdev`](@ref) is of [`mdev`](@ref).

```math
σ_{x,H}(τ) = \\frac{τ}{\\sqrt{10/3}}\\,\\mathrm{Mod}\\,σ_H(τ)
```

The statistic is defined in this package:
[SP1065](@cite riley-2008-sp1065) and
[IEEE Std 1139-2022](@cite ieee1139-2022-definitions) define no Hadamard
time deviation. The factor is chosen as TDEV's is. On white PM the result is
the standard deviation of the ``m``-sample mean of the phase, so `htdev` and
`tdev` agree there, and `htdev` keeps the drift rejection of MHDEV. Use it in
place of `tdev` when the clock drifts.

Arguments and keywords are those of [`mhdev`](@ref). The interval is
MHDEV's scaled by the same factor, with the same EDF. See
[Choosing an estimator](@ref) for when to use it.
"""
function htdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, remove_drift::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    res = mhdev(data, m_values; ci=ci, remove_drift=remove_drift, noise=noise, confidence=confidence)
    factor = res.tau ./ sqrt(10.0 / 3.0)

    if !ci
        return _result(data, confidence, :htdev, res.tau, res.dev .* factor, Symbol[], _CIBound[], Float64[], res.neff)
    end

    return _result(data, confidence, :htdev, res.tau, res.dev .* factor, res.noise_type,
                           _scale_ci(res.ci, factor), res.edf, res.neff)
end


# Deviation-scale unbias divisor √B, with NaN for B ≤ 0: :totvar's
# B = 1 − a·τ/T goes negative for τ > T/a, which gives a NaN row instead of
# a DomainError for the whole call.
_unbias_divisor(B::AbstractVector{<:Real}) = [b > 0 ? sqrt(b) : NaN for b in B]

# What `correct_bias` asks of one estimator: :exact (`true`), :none (`false`),
# :published, or for `:stable32` the mode that the estimator passes: the one of
# those three that Stable32 applies to it, or :stable32 where Stable32 has a
# correction of its own (Thêo1).
function _bias_mode(correct_bias::Union{Bool,Symbol}, stable32::Symbol)
    correct_bias === true && return :exact
    correct_bias === false && return :none
    correct_bias === :published && return :published
    correct_bias === :stable32 && return stable32
    throw(ArgumentError("correct_bias must be true, false, :published or :stable32; got :$correct_bias"))
end

_unbias(raw, mode::Symbol, noises, var_type::Symbol, m_values, N::Int) = mode == :none ? raw :
    raw ./ _unbias_divisor(bias_correction(noises, var_type, m_values, N; published=mode == :published))

"""
    totdev(data, m_values; ci=true, correct_bias=true, noise=:auto,
           confidence=DEFAULT_CONFIDENCE) → StabilityResult

Total deviation (TOTDEV): an estimator of the Allan deviation with more
degrees of freedom at long averaging times.

The record is extended at each end by reflection about the end sample,
``x^*_{1-j} = 2x_1 - x_{1+j}`` and ``x^*_{N+j} = 2x_N - x_{N-j}``, and a
second difference is taken at every interior sample:

```math
σ_{\\mathrm{total}}^2(τ) = \\frac{1}{2τ^2(N-2)} \\sum_{i=2}^{N-1} \\left(x^*_{i-m} - 2x^*_i + x^*_{i+m}\\right)^2
```

([Greenhall et al. (1999), Eqs. (2) and (3)](@cite greenhall-1999-totvar);
[SP1065 §5.2.11, Eq. (25)](@cite riley-2008-sp1065);
[IEEE Std 1139-2022, Eqs. (A.32) and (A.33)](@cite ieee1139-2022-definitions)).

Arguments and the keywords `ci`, `noise` and `confidence` are those of
[`adev`](@ref). The record must have no gaps; fill them with
[`fillgaps`](@ref) first. Each extension holds ``N - 2`` samples, so the
deviation exists for ``m ≤ N - 1``. A larger averaging factor gives `NaN`
with `neff` 0, and the grids stop at ``m = \\lfloor (N-1)/2 \\rfloor``.

# Keywords
- `correct_bias=true`: divide the deviation by ``\\sqrt{B}``, where ``B`` is
  the expected ratio of the total variance to the Allan variance for the
  noise type at that ``τ`` (see [`bias_correction`](@ref)). The total
  variance reads low on flicker FM and random-walk FM, ``B = 1 - aτ/T`` for
  a record of length ``T``
  ([Greenhall et al. (1999), Eq. (25)](@cite greenhall-1999-totvar)), and
  high on white and flicker PM
  (see [Total estimators on PM noise](@ref)). `:published` and `:stable32`
  correct FM noise only, as SP1065 and Stable32 do. `false` returns the
  uncorrected value, as allantools does.

The EDF on FM noise is ``b\\,T/τ - c``
([Greenhall et al. (1999), Eq. (26)](@cite greenhall-1999-totvar), with the
coefficients as rounded in SP1065 Table 7). On PM noise it follows the model
of [Total estimators on PM noise](@ref). At ``m = 1`` no term reaches the
extension, the estimate equals ADEV, and ADEV's EDF is used.
"""
function totdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Union{Bool,Symbol}=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _gapfree(:totdev, data, m_values)
    raw_devs = _totdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:totdev, length(x), m_values)

    # Noise IDs are needed for either path: bias correction and CIs both read α.
    mode = _bias_mode(correct_bias, :published)
    need_noise = mode != :none || ci
    noises = !need_noise     ? Symbol[] :
             noise === :auto ? identify_noise(x, m_values, dmin=0, dmax=2) :
                               _stated_noises(noise, :adev, taus, data.tau0)

    devs = _unbias(raw_devs, mode, noises, :totvar, m_values, length(x))

    if !ci
        return _result(data, confidence, :totdev, taus, devs, noises, _CIBound[], Float64[], neff)
    end

    edfs = calculate_edf(:totdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, confidence)

    return _result(data, confidence, :totdev, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    mtotdev(data, m_values; ci=true, correct_bias=true, noise=:auto,
            confidence=DEFAULT_CONFIDENCE) → StabilityResult

Modified total deviation (MTOTDEV): the total form of [`mdev`](@ref).

Each of the ``N - 3m + 1`` subsequences of ``3m`` phase samples is treated
on its own. A linear trend is removed, with its slope taken from the means
of the two halves of the subsequence. The residual is extended at both ends
by even reflection to ``9m`` samples. The squared second differences of
``m``-sample phase means are then averaged over the ``6m`` positions of the
extension. The variance is the mean over the subsequences, divided by
``2τ^2``
([Howe and Vernotte (1999), Eqs. (5) and (7)](@cite howe-1999-modtotvar);
[Howe and Peppler (2001), Eqs. (1) to (5)](@cite howe-2001-total-estimators);
[SP1065 §5.2.12, Eq. (27)](@cite riley-2008-sp1065)).

Arguments and the keywords `ci`, `noise` and `confidence` are those of
[`adev`](@ref). The record must have no gaps.

# Keywords
- `correct_bias=true`: divide the deviation by ``\\sqrt{B}``, where ``B`` is
  the expected ratio of the modified total variance to the modified Allan
  variance (see [`bias_correction`](@ref)). The estimator reads low:
  [Howe and Peppler (2001), Table 2](@cite howe-2001-total-estimators) give
  ``B = 0.995``, 0.851, 0.771, 0.717 and 0.679 for white PM, flicker PM,
  white FM, flicker FM and random-walk FM at ``τ ≥ 16τ_0``. Below that the
  exact value at each ``m`` is used. At ``m = 1`` the estimate is
  MDEV/``\\sqrt{2}`` for any record and ``B = 1/2``. `:published` applies
  the table at every ``m > 1``. `false` or `:stable32` returns the
  uncorrected value, which is what Stable32 and allantools report.

The EDF is MDEV's below ``m = 16``, as in Stable32, and
``(T/τ)/(b_0 + b_1τ/T)`` with the coefficients of Howe and Peppler's Table 2
from there. SP1065 Table 8 lists a different approximation, which is not
used.
"""
function mtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Union{Bool,Symbol}=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _gapfree(:mtotdev, data, m_values)
    raw_devs = _mtotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:mtotdev, length(x), m_values)

    mode = _bias_mode(correct_bias, :none)
    need_noise = mode != :none || ci
    noises = !need_noise     ? Symbol[] :
             noise === :auto ? identify_noise(x, m_values, dmin=0, dmax=2) :
                               _stated_noises(noise, :mdev, taus, data.tau0)

    # At m = 1 the detrended three-point subsequence is (0, z, 0) with z half
    # the second difference, and its reflection averages to half of MVAR.
    devs = _unbias(raw_devs, mode, noises, :mtot, m_values, length(x))

    if !ci
        return _result(data, confidence, :mtotdev, taus, devs, noises, _CIBound[], Float64[], neff)
    end

    edfs = calculate_edf(:mtotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, confidence)

    return _result(data, confidence, :mtotdev, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    htotdev(data, m_values; ci=true, correct_bias=true, remove_drift=true,
            noise=:auto, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Hadamard total deviation (HTOTDEV): the total form of [`hdev`](@ref).

The computation is that of [`mtotdev`](@ref) applied to the ``N - 1``
fractional-frequency values ``y_i = (x_{i+1} - x_i)/τ_0``. Each of the
``N - 3m`` subsequences of ``3m`` values has a linear trend removed, with
its slope taken from the means of its two halves, and is extended at both
ends by even reflection to ``9m`` values. The squared second differences of
``m``-sample frequency means are averaged over the ``6m`` positions of the
extension. The variance is the mean over the subsequences, divided by 6
([Howe et al. (2005), Eqs. (3) to (5)](@cite howe-2005-tothvar);
[Howe and Peppler (2001)](@cite howe-2001-total-estimators);
[SP1065 §5.2.14, Eq. (29)](@cite riley-2008-sp1065)). At ``m = 1`` the
estimate is the overlapping Hadamard deviation.

Arguments and the keywords `ci`, `noise` and `confidence` are those of
[`adev`](@ref). The record must have no gaps.

# Keywords
- `correct_bias=true`: divide the deviation by ``\\sqrt{B}``, where ``B`` is
  the expected ratio of the Hadamard total variance to the Hadamard variance
  (see [`bias_correction`](@ref)). On FM noise the estimator reads low,
  ``B = 0.995``, 0.851, 0.771, 0.717 and 0.679 for white FM, flicker FM,
  random-walk FM, flicker-walk FM and random-run FM
  ([Howe et al. (2005), Table I](@cite howe-2005-tothvar)). Those values
  hold from ``m = 16``; below that the exact value at each ``m`` is used. On
  white and flicker PM it reads high, by up to 38 % in variance, and ``B``
  is its exact expectation (see [Total estimators on PM noise](@ref)).
  ``B = 1`` at ``m = 1``. `:published` or `:stable32` applies the published
  FM factor at every ``m > 1`` and none on PM, as Stable32 does. `false`
  returns the uncorrected value, as allantools does.
- `remove_drift=true`: as in [`hdev`](@ref). Here the noise type also
  selects ``B``, so the keyword can change the deviation where the two
  identifications differ.

The EDF on FM noise is ``(T/τ)/(b_0 + b_1τ/T)``
([Howe et al. (2005), Eq. (7)](@cite howe-2005-tothvar)). On PM noise it
follows the model of [Total estimators on PM noise](@ref), and at ``m = 1``
it is HDEV's.
"""
function htotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Union{Bool,Symbol}=true, remove_drift::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _gapfree(:htotdev, data, m_values)
    raw_devs = _htotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:htotdev, length(x), m_values)

    mode = _bias_mode(correct_bias, :published)
    need_noise = mode != :none || ci
    noises = !need_noise     ? Symbol[] :
             noise === :auto ? _hadamard_noise(x, m_values, remove_drift) :
                               _stated_noises(noise, :hdev, taus, data.tau0)

    devs = _unbias(raw_devs, mode, noises, :htot, m_values, length(x))

    if !ci
        return _result(data, confidence, :htotdev, taus, devs, noises, _CIBound[], Float64[], neff)
    end

    edfs = calculate_edf(:htotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, confidence)

    return _result(data, confidence, :htotdev, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    ttotdev(data, m_values; ci=true, correct_bias=true, noise=:auto,
            confidence=DEFAULT_CONFIDENCE) → StabilityResult

Time total deviation (TTOTDEV), in seconds: the time-stability form of
[`mtotdev`](@ref), as [`tdev`](@ref) is of [`mdev`](@ref).

```math
σ_{x,\\mathrm{total}}(τ) = \\frac{τ}{\\sqrt{3}}\\,\\mathrm{Mod}\\,σ_{\\mathrm{total}}(τ)
```

([Howe and Peppler (2001)](@cite howe-2001-total-estimators);
[SP1065 §5.2.13](@cite riley-2008-sp1065), whose Eq. (28) prints ``τ^3/3``
where the variance takes ``τ^2/3``).

Arguments and keywords are those of [`mtotdev`](@ref). The interval is
MTOTDEV's scaled by the same factor, with the same EDF.
"""
function ttotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Union{Bool,Symbol}=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    res = mtotdev(data, m_values; ci=ci, correct_bias=correct_bias, noise=noise, confidence=confidence)
    factor = res.tau ./ sqrt(3.0)

    if !ci
        return _result(data, confidence, :ttotdev, res.tau, res.dev .* factor, Symbol[], _CIBound[], Float64[], res.neff)
    end

    return _result(data, confidence, :ttotdev, res.tau, res.dev .* factor, res.noise_type,
                           _scale_ci(res.ci, factor), res.edf, res.neff)
end

"""
    mhtotdev(data, m_values; ci=true, correct_bias=true, remove_drift=true,
             noise=:auto, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Modified Hadamard total deviation (MHTOTDEV): the total form of
[`mhdev`](@ref).

Each of the ``N - 4m + 1`` windows of ``4m`` phase samples has its
least-squares quadratic removed and is extended at both ends by even
reflection to ``12m`` samples. The squared third differences of
``m``-sample phase means are averaged over the ``8m + 1`` positions of the
extension. The variance is the mean over the windows, divided by ``6τ^2``.
Removing a quadratic from each window makes the estimate insensitive to
frequency offset and linear frequency drift, and it converges for every
noise type MHDEV does, down to random-run FM (``α = -4``).

The construction is defined in this package.
[Vernotte and Howe (2000)](@cite vernotte-2000-total-structure) outline a
total form of the same third-difference variance, with a linear trend
removed from each subsequence and a quadratic from the whole record, and
give no bias or EDF for it. Bias and EDF for the estimator here are derived
in [MHTOTDEV bias and EDF](@ref).

Arguments and the keywords `ci`, `noise` and `confidence` are those of
[`adev`](@ref). The record must have no gaps.

# Keywords
- `correct_bias=true`: divide the deviation by ``\\sqrt{B}``, where ``B`` is
  the expected ratio of this variance to the modified Hadamard variance
  (see [`bias_correction`](@ref)): 0.78 to 0.90 from ``m = 16``, the exact
  value at each smaller ``m``, and 124/225 at ``m = 1``. `false` returns the
  uncorrected value. `:published` and `:stable32` act as `true`, since no
  reference defines the estimator.
- `remove_drift=true`: remove the least-squares frequency drift of the whole
  record before the analysis, so that the noise type is identified on the
  record without it, as in [`hdev`](@ref). The uncorrected deviation does
  not depend on it; the corrected one changes only where the identified
  noise type, which selects ``B``, does.

The EDF is a model fitted to the exact EDF of the statistic, and MHDEV's at
``m = 1``.
"""
function mhtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Union{Bool,Symbol}=true, remove_drift::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _gapfree(:mhtotdev, data, m_values)
    if remove_drift
        x = _remove_quadratic(x)
    end
    raw_devs = _mhtotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:mhtotdev, length(x), m_values)

    mode = _bias_mode(correct_bias, :exact)
    need_noise = mode != :none || ci
    noises = !need_noise     ? Symbol[] :
             noise === :auto ? _hadamard_noise(x, m_values, remove_drift) :
                               _stated_noises(noise, :mhdev, taus, data.tau0)

    # At m = 1 the detrended four-point window is c·(−1, 3, −3, 1), and the
    # nine third differences of its reflection average to 124/225 of MHVAR.
    devs = _unbias(raw_devs, mode, noises, :mhtot, m_values, length(x))

    if !ci
        return _result(data, confidence, :mhtotdev, taus, devs, noises, _CIBound[], Float64[], neff)
    end

    edfs = calculate_edf(:mhtotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, confidence)

    return _result(data, confidence, :mhtotdev, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end


"""
    mtie(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Maximum time interval error (MTIE), in seconds: the largest peak-to-peak
phase excursion within any observation interval ``τ = mτ_0``, a window of
``m + 1`` samples.

```math
\\mathrm{MTIE}(mτ_0) = \\max_{1 \\le k \\le N-m} \\left(\\max_{k \\le i \\le k+m} x_i - \\min_{k \\le i \\le k+m} x_i\\right)
```

([ITU-T G.810, clause 4.5.15 and Appendix II.5](@cite ITU-G810);
[SP1065 §5.2.17, Eq. (35)](@cite riley-2008-sp1065)).

MTIE is a peak statistic, not a variance: one outlier or phase step sets it,
and it has no χ² interval. The fields `noise_type`, `ci` and `edf` of the
result are empty. The keywords `ci` and `confidence` are accepted, so that
every estimator can be called alike, and ignored. With missing samples, each
window gives the peak-to-peak of the samples measured in it. Compare the
result with a network limit using [`mask_check`](@ref).

`data` and `m_values` are as in [`adev`](@ref).

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData([0.0, 1.0, 0.5, 2.0, 1.5], 1.0);

julia> r = mtie(p, [1, 3]);

julia> r.dev
2-element Vector{Float64}:
 1.5
 2.0
```
"""
function mtie(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    _has_gaps(data) && return _gap_result(:mtie, _gap_record(data), m_values; ci=ci, confidence=confidence)
    _check_m(:mtie, m_values)
    raw_devs = _mtie_core(_f64(data.x), m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:mtie, length(data.x), m_values)
    return _result(data, confidence, :mtie, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
end


"""
    tierms(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Root-mean-square time interval error (TIE rms), in seconds: the rms phase
change over an interval ``τ = mτ_0``.

```math
\\mathrm{TIE}_{\\mathrm{rms}}(mτ_0) = \\sqrt{\\frac{1}{N-m} \\sum_{i=1}^{N-m} \\left(x_{i+m} - x_i\\right)^2}
```

([SP1065 §5.2.18, Eq. (37)](@cite riley-2008-sp1065);
[ITU-T G.810, Appendix II.4](@cite ITU-G810)). Stable32 and allantools
(`tierms`) compute the same estimator.

Unlike [`tdev`](@ref), TIE rms includes any frequency offset in the record,
which contributes the offset times ``τ``, and it does not converge for
flicker FM and random-walk FM (G.810, Appendix II.4). Like [`mtie`](@ref)
it has no χ² interval: `noise_type`, `ci` and `edf` are empty, and the
keywords `ci` and `confidence` are accepted and ignored. A term with a
missing sample is skipped.

`data` and `m_values` are as in [`adev`](@ref).

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData([0.0, 1.0, 0.5, 2.0, 1.5], 1.0);

julia> r = tierms(p, [1, 3]);

julia> round.(r.dev; sigdigits=4)
2-element Vector{Float64}:
 0.9682
 1.458
```
"""
function tierms(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    _has_gaps(data) && return _gap_result(:tierms, _gap_record(data), m_values; ci=ci, confidence=confidence)
    _check_m(:tierms, m_values)
    raw_devs = _tierms_core(_f64(data.x), m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:tierms, length(data.x), m_values)
    return _result(data, confidence, :tierms, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
end


"""
    pdev(data, m_values; ci=true, noise=:auto,
         confidence=DEFAULT_CONFIDENCE) → StabilityResult

Parabolic deviation (PDEV): the two-sample deviation of frequencies
estimated by linear regression on the phase, as an Ω counter measures them
([IEEE Std 1139-2022, A.3.4](@cite ieee1139-2022-definitions)).

For ``m ≥ 2``,

```math
σ_P^2(τ) = \\frac{72}{m^4τ^2(N-2m)} \\sum_{i=1}^{N-2m} \\left[\\sum_{k=0}^{m-1}
    \\left(\\frac{m-1}{2} - k\\right)\\left(x_{i+k} - x_{i+k+m}\\right)\\right]^2
```

([Vernotte et al. (2016), Eq. (15)](@cite vernotte-2016-pvar)), and at
``m = 1`` it is the Allan deviation by definition (their Eq. (14)). Like
[`mdev`](@ref), PDEV falls as ``τ^{-3/2}`` on white PM and ``τ^{-1}`` on
flicker PM. It is the more efficient of the two on white PM, because linear
regression is the best estimate of a frequency in white phase noise.

The EDF follows
[Vernotte et al. (2021)](@cite vernotte-2021-pvar-noninteger):
``ν ≈ 35/(A\\,m/M - 12\\,(m/M)^2)`` with ``M = N - 2m`` up to ``m ≈ N/4``
(their Eqs. (16) and (18)), and their interpolation in ``\\ln m`` down to
``ν = 1`` near ``m = N/2`` beyond (Eqs. (22) to (24)). The authors give its
accuracy as ±10 % except at ``m = 1`` and 2. At ``m = 1`` the EDF of
[`adev`](@ref) is used.

Arguments and keywords are those of [`adev`](@ref), less `overlapping`. The
record must have no gaps.

# Examples

```jldoctest
julia> using SigmaTau

julia> p = PhaseData(collect(1.0:10.0), 1.0);

julia> r = pdev(p, [1, 2]; ci=false);

julia> round.(r.dev; sigdigits=4)
2-element Vector{Float64}:
 0.0
 0.0
```
"""
function pdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _gapfree(:pdev, data, m_values)
    raw_devs = _pdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:pdev, length(x), m_values)

    if !ci
        return _result(data, confidence, :pdev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = noise === :auto ? identify_noise(x, m_values, dmin=0, dmax=2) :
                               _stated_noises(noise, :mdev, taus, data.tau0)
    edfs = calculate_edf(:pdev, raw_devs, noises, m_values, taus, length(x), (length(x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, confidence)

    return _result(data, confidence, :pdev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end


"""
    theo1(data, m_values; ci=true, correct_bias=true, noise=:auto,
          confidence=DEFAULT_CONFIDENCE) → StabilityResult

Thêo1 deviation, bias-removed by default (ThêoBR): an estimator of the Allan
deviation that reaches averaging times up to 75 % of the record length.

For an even averaging factor ``m``,

```math
\\begin{aligned}
\\mathrm{Thêo1}(m, τ_0, N) = {} & \\frac{1}{0.75\\,(N-m)(mτ_0)^2} \\sum_{i=1}^{N-m} \\sum_{δ=0}^{m/2-1} \\frac{1}{m/2-δ} \\\\
& \\times \\left[\\left(x_i - x_{i-δ+m/2}\\right) + \\left(x_{i+m} - x_{i+δ+m/2}\\right)\\right]^2
\\end{aligned}
```

at the effective averaging time ``τ = 0.75\\,mτ_0``
([Howe and Peppler (2003), Eq. (1)](@cite howe-2003-theo1);
[Howe and Tasset (2004), Eq. (1)](@cite howe-2004-theo1);
[IEEE Std 1139-2022, Eq. (A.35)](@cite ieee1139-2022-definitions)).
[SP1065 §5.2.15, Eq. (30)](@cite riley-2008-sp1065) prints the sum without
the factor 0.75, which the white-FM bias of 1 in its Table 2 requires.

# Arguments
- `data`: as in [`adev`](@ref). The record must have no gaps.
- `m_values`: even Thêo1 factors with ``2 ≤ m ≤ N - 1``; any other value
  gives `NaN`. The published definition takes ``m ≥ 10``. A
  [`TauMode`](@ref) grid is rounded up to even factors. The result reports
  ``τ = 0.75\\,mτ_0``, as Stable32 does; allantools reports the same
  deviations at ``mτ_0``.

# Keywords
- `correct_bias=true`: apply the ThêoBR correction. The Thêo1 variance is
  multiplied by the mean ratio of the Allan variance at ``m = 9 + 3i`` to
  Thêo1 at ``m = 12 + 4i``, two factors with the same averaging time, over
  ``i = 0, …, n`` with ``n = \\lfloor N/6 - 3 \\rfloor``
  ([SP1065 §5.2.16, Eq. (33)](@cite riley-2008-sp1065);
  [IEEE Std 1139-2022, Eq. (A.37)](@cite ieee1139-2022-definitions)). The
  ratio is measured on the record, so no noise type is identified. One ratio
  removes the bias where it changes little with ``m``, as on FM noise of one
  type. On white and flicker PM the bias changes with ``m`` throughout
  ([Howe (2006), Table 3](@cite howe-2006-theoh)), and the corrected
  variance was 1.04 to 1.25 times the Allan variance at 0.2 to 0.7 of the
  record length in simulation (``N = 2048``). It needs ``N ≥ 19``. `false`
  returns Thêo1 uncorrected.
- `correct_bias=:published`: the same correction with
  ``n = \\lfloor N/30 - 3 \\rfloor``, as
  [Howe (2006), Eq. (6)](@cite howe-2006-theoh) prints it, which stops at a
  tenth of the record. Where the noise type changes within the record
  beyond that point, the ratio is that of the short averaging times at which
  it is measured: with PM giving way to FM the corrected variance at long
  ``τ`` was 0.18 to 0.77 times the Allan variance in simulation.
- `correct_bias=:stable32`: the correction of a Stable32 Thêo1 or ThêoH
  run with automatic noise identification. The
  [Stable32 User Manual](@cite hamilton-2008-stable32-manual) gives it as
  the default ratio, computed from at most 3000 phase samples (its setting
  `MaxTheoBRPts`), and names
  [Taylor and Howe (2010)](@cite taylor-2010-fast-theobr) as the method.
  For a longer record this package takes every ``b``-th sample,
  ``b = \\lceil N/3000 \\rceil``, and adds ``0.00055\\,(b - 1)`` to the
  ratio. The term is the linear correction that Taylor and Howe fitted on
  simulated FM noise for a record averaged in blocks of ``b`` samples
  (their Eq. (4)). The deviations agree with Stable32's to the five digits
  that it prints, at each of the 181 averaging factors of a Thêo1 run on
  its 1001-point record `PHASE.DAT` (``b = 1``) and at the two Thêo1
  points of a ThêoH run on an 8192-point record (``b = 3``), where block
  averages give deviations 6.5 % higher. That Stable32 takes every
  ``b``-th sample, and the value of ``b``, are inferred from that
  agreement, and no other ``b`` has been compared.

  With a noise type entered by hand, Stable32 applies the bias that its
  manual tabulates for the type, which is that of
  [SP1065 Table 2](@cite riley-2008-sp1065). It is 1 for white FM, so
  `correct_bias=false` matches a Stable32 run made for white FM. This
  package has no option for the tabulated bias of the other noise types.
- `ci`, `noise`, `confidence`: as in [`adev`](@ref). The EDF is that of the
  empirical formulas of [SP1065 Table 3](@cite riley-2008-sp1065), with or
  without the correction.

!!! note "Length of the ThêoBR average"
    SP1065 Eq. (33), the Stable32 manual and IEEE Std 1139-2022 Eq. (A.37)
    print ``n = \\lfloor N/6 - 3 \\rfloor``, which runs the Allan variance
    terms out to half the record. Howe (2006),
    [Howe and Tasset (2004), Eq. (4)](@cite howe-2004-theo1),
    [McGee and Howe (2007)](@cite mcgee-2007-theoh) and
    [Taylor and Howe (2010), Eq. (1)](@cite taylor-2010-fast-theobr) all give
    ``n = \\lfloor 0.1N/3 - 3 \\rfloor``, which keeps them within the first
    tenth. This package takes the longer average, and the shorter one under
    `correct_bias=:published`. In simulation (``N = 2048``) the two fail in
    different places. Where white or flicker PM gave way to FM within the
    record, the shorter average read 0.18 to 0.77 of the Allan variance at
    long ``τ`` and the longer 0.56 to 1.12. On a record of one FM noise type
    the last terms of the longer average hold one or two Allan variance
    samples each: its ratio scattered 2.5 times more (white FM, 5 to 95 %:
    0.58 to 1.41 against 0.83 to 1.16) and averaged 7 to 9 % low, at every
    record length tried, so ThêoBR at short ``τ`` is less certain than ADEV
    there.

The double sum costs ``O(Nm)`` per averaging factor. On long records the
ThêoBR ratio computes Thêo1 at 48 of its ``n + 1`` factors and interpolates
the rest, which stays within 0.7 % of the full average; the longer average
still makes the default call about 2.5 times slower than `:published` on
very long records. The `:stable32` ratio computes every factor.
"""
function theo1(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Union{Bool,Symbol}=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _gapfree(:theo1, data, m_values)
    N = length(x)
    raw_devs = _theo1_core(x, m_values, data.tau0)
    taus = 0.75 .* m_values .* data.tau0
    T = (N - 1) * data.tau0
    neff = _neff_counts(:theo1, N, m_values)

    mode = _bias_mode(correct_bias, :stable32)
    devs = mode == :none ? raw_devs :
           mode == :stable32 ? raw_devs .* sqrt(_stable32_theobr_ratio(x, data.tau0)) :
                               raw_devs .* sqrt(_theobr_ratio(x, data.tau0; howe=mode == :published))

    if !ci
        return _result(data, confidence, :theo1, taus, devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = noise === :auto ? identify_noise(x, m_values, dmin=0, dmax=2) :
                               _stated_noises(noise, :adev, taus, data.tau0)
    edfs = calculate_edf(:theo1, devs, noises, m_values, taus, N, T)
    lower, upper = confidence_intervals(devs, edfs, confidence)

    return _result(data, confidence, :theo1, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    theoh(data, m_values; ci=true, noise=:auto,
          confidence=DEFAULT_CONFIDENCE) → StabilityResult

ThêoH: the overlapping Allan deviation at short averaging times joined to
ThêoBR, the bias-removed Thêo1, at long ones, so that one curve runs from
``τ_0`` to 75 % of the record length
([Howe (2006), Eq. (7)](@cite howe-2006-theoh);
[SP1065 §5.2.16, Eq. (34)](@cite riley-2008-sp1065)).

Each requested factor ``m``, any positive integer, is a target averaging
time ``mτ_0`` and is computed according to its place in the record of length
``T = (N-1)τ_0``:

- ``mτ_0 ≤ 0.2\\,T``: [`adev`](@ref) at ``m``.
- ``mτ_0 > 0.2\\,T``: ThêoBR (see [`theo1`](@ref)) at the even Thêo1 factor
  nearest ``4m/3``, reported at its own effective averaging time, which can
  differ from ``mτ_0`` by one rounding step.

The crossover at 20 % of the record is that of SP1065. Howe (2006) places it
at 10 %. The interval uses the EDF of ADEV below the crossover and of Thêo1
above it. The bias removal is always applied, with the default ratio of
`theo1`, and needs ``N ≥ 19``. For uncorrected Thêo1, or for the ratio of
a Stable32 ThêoH run, call `theo1` with `correct_bias=false` or
`correct_bias=:stable32`.

`data` and the keywords are as in [`adev`](@ref), less `overlapping`. The
record must have no gaps.
"""
function theoh(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, noise::_NoiseChoice=:auto, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _gapfree(:theoh, data, m_values)
    N = length(x)
    tau0 = data.tau0
    T = (N - 1) * tau0
    n_pts = length(m_values)

    # Split the requested τ grid at the SP1065 Eq. (34) crossover (20 % of T).
    adev_idx = Int[]
    theo_idx = Int[]
    for (j, m) in enumerate(m_values)
        push!(_theoh_is_adev(m, N) ? adev_idx : theo_idx, j)
    end
    adev_ms = m_values[adev_idx]
    theo_ms = Int[_theoh_theo_m(m) for m in m_values[theo_idx]]

    taus = Float64[m * tau0 for m in m_values]
    devs = fill(NaN, n_pts)

    adev_devs = _adev_core(x, adev_ms, tau0)
    devs[adev_idx] = adev_devs

    # Skip the ThêoBR ratio (itself a sampled Thêo1 ladder) when no point
    # falls on the Thêo1 segment.
    theo_devs = isempty(theo_idx) ? Float64[] :
                _theo1_core(x, theo_ms, tau0) .* sqrt(_theobr_ratio(x, tau0))
    devs[theo_idx] = theo_devs
    taus[theo_idx] = 0.75 .* theo_ms .* tau0

    neff = _neff_counts(:theoh, N, m_values)

    if !ci
        return _result(data, confidence, :theoh, taus, devs, Symbol[], _CIBound[], Float64[], neff)
    end

    # Per-segment noise ID and EDF: ADEV EDF below the crossover, Thêo1
    # EDF above, merged back into request order.
    noises = noise === :auto ? fill(:WHFM, n_pts) : _stated_noises(noise, :adev, taus, tau0)
    edfs = fill(NaN, n_pts)
    if !isempty(adev_idx)
        noise === :auto && (noises[adev_idx] = identify_noise(x, adev_ms, dmin=0, dmax=2))
        edfs[adev_idx] = calculate_edf(:adev, adev_devs, noises[adev_idx], adev_ms,
                                       taus[adev_idx], N, T)
    end
    if !isempty(theo_idx)
        noise === :auto && (noises[theo_idx] = identify_noise(x, theo_ms, dmin=0, dmax=2))
        edfs[theo_idx] = calculate_edf(:theo1, theo_devs, noises[theo_idx], theo_ms,
                                       taus[theo_idx], N, T)
    end
    lower, upper = confidence_intervals(devs, edfs, confidence)

    return _result(data, confidence, :theoh, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end


# Dynamic (time-resolved) deviations: the overlapping kernel (`_adev_core`,
# `_hdev_core`) evaluated at each position of a sliding window, the estimator
# of Galleani and Tavella (2009), Eq. (34). One window buffer is allocated
# and each window is copied into it, so the cost is one static-kernel call on
# `window` samples per window.

# Shared driver. `core` is one of the `_*_core` kernels; `min_window` is the
# smallest window that supports m = 1 for that kernel (window − span·1 ≥ 2).
function _dynamic_dev(core::F, dev_sym::Symbol, data::PhaseData,
                      m_values::Vector{Int};
                      window::Int, step::Int, min_window::Int) where {F}
    x = _gapfree(dev_sym, data, m_values)
    N = length(x)
    window <= N || throw(ArgumentError(
        "$dev_sym: window=$window exceeds the record length N=$N"))
    window >= min_window || throw(ArgumentError(
        "$dev_sym: window=$window is too short to support any averaging factor (need ≥ $min_window samples)"))
    step >= 1 || throw(ArgumentError("$dev_sym: step must be ≥ 1, got $step"))

    starts = 1:step:(N - window + 1)
    taus = m_values .* data.tau0
    # Window covering samples s … s+window−1 is centered at coordinate time
    # t_c = (s − 1 + (window − 1)/2)·τ₀, with sample x[1] at t = 0.
    centers = Float64[(s - 1 + (window - 1) / 2) * data.tau0 for s in starts]

    dev = Matrix{Float64}(undef, length(starts), length(m_values))
    buf = Vector{Float64}(undef, window)
    for (i, s) in enumerate(starts)
        copyto!(buf, 1, x, s, window)
        dev[i, :] = core(buf, m_values, data.tau0)
    end

    return DynamicStabilityResult(dev_sym, centers, taus, dev, window, data.tau0)
end

"""
    dadev(data, m_values=Octave; window, step=window÷2) → DynamicStabilityResult

Dynamic Allan deviation ``σ_y(t, τ)`` (DADEV): the overlapping Allan
deviation of a sliding window of the record, as a function of the position
of the window.

For a window of ``N_w`` samples centered at time ``t``,

```math
σ_y^2(t, τ) = \\frac{1}{2τ^2(N_w-2m)} \\sum_{i} \\left(x_{i+2m} - 2x_{i+m} + x_i\\right)^2 ,
```

summed over the ``N_w - 2m`` second differences that fit in the window
([Galleani and Tavella (2009), Eq. (34)](@cite galleani-2009-davar);
[IEEE Std 1139-2022, Eq. (A.40)](@cite ieee1139-2022-definitions)). The map
shows changes of stability with time, such as a step in noise level or a
burst, that a single ``σ_y(τ)`` curve averages away.
[McKelvy et al. (2026)](@cite mckelvy-2026-telemetry-stability) take it as
the target of models that estimate the stability of hydrogen masers from
their telemetry.

# Arguments
- `data`: a [`PhaseData`](@ref) or [`FrequencyData`](@ref) record without
  gaps.
- `m_values`: averaging factors shared by every window, or a
  [`TauMode`](@ref) grid, which is resolved against the window length and
  not the record length. A factor the window cannot support
  (``N_w - 2m < 2``) gives a `NaN` column.

# Keywords
- `window`: the window length ``N_w`` in samples. The default is the largest
  power of two up to ``N/8``, and at least 16, which gives about fifteen
  windows.
- `step=window÷2`: the number of samples between successive windows. The
  cost is one Allan deviation of `window` samples per window, so halving
  `step` doubles the run time.

The result holds one row per window and one column per ``τ``. No confidence
interval is reported. [Galleani (2011)](@cite galleani-2011-davar-confidence)
derives confidence surfaces for the dynamic Allan variance, which are not
implemented here. Each row is an ordinary Allan deviation of `window`
samples, so [`adev`](@ref) on that slice of the record gives its interval.
"""
function dadev(data::PhaseData, m_values::Vector{Int};
               window::Int = _default_window(length(data.x)), step::Int = window ÷ 2)
    # min window: _adev_core needs N_w − 2m ≥ 2 ⇒ window ≥ 4 at m = 1.
    return _dynamic_dev(_adev_core, :dadev, data, m_values;
                        window=window, step=step, min_window=4)
end

"""
    dhdev(data, m_values=Octave; window, step=window÷2) → DynamicStabilityResult

Dynamic Hadamard deviation (DHDEV): the overlapping Hadamard deviation
([`hdev`](@ref)) of a sliding window, built as [`dadev`](@ref) is
([SP1065 §5.2.20](@cite riley-2008-sp1065)). The third difference cancels
linear frequency drift within each window, so the map of a drifting clock
shows its noise and not its drift.

Arguments and keywords are those of [`dadev`](@ref). A factor the window
cannot support (``N_w - 3m < 2``) gives a `NaN` column. No confidence
interval is reported.
"""
function dhdev(data::PhaseData, m_values::Vector{Int};
               window::Int = _default_window(length(data.x)), step::Int = window ÷ 2)
    # min window: _hdev_core needs N_w − 3m ≥ 2 ⇒ window ≥ 5 at m = 1.
    return _dynamic_dev(_hdev_core, :dhdev, data, m_values;
                        window=window, step=step, min_window=5)
end


# TauMode grid selector: resolve against the *window* length (every window
# computes the same m grid; the record length is irrelevant to the clamp).
# With no grid, the octave grid; with no window, `_default_window`.
dadev(data::AbstractTimingData, taus::TauMode = Octave; window::Int = _default_window(length(data)), kwargs...) =
    dadev(data, tau_values(taus, window, :adev); window=window, kwargs...)
dhdev(data::AbstractTimingData, taus::TauMode = Octave; window::Int = _default_window(length(data)), kwargs...) =
    dhdev(data, tau_values(taus, window, :hdev); window=window, kwargs...)

# Default dynamic window: the largest power of two up to N ÷ 8 (at least 16,
# at most N), so the map has about fifteen half-overlapped windows.
function _default_window(N::Int)
    w = 16
    while 2w <= N ÷ 8
        w *= 2
    end
    return min(w, N)
end

identify_noise(data::AbstractTimingData, taus::TauMode; kwargs...) =
    identify_noise(data, tau_values(taus, length(data), :adev); kwargs...)

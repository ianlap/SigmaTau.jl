# api/total.jl — User wrappers for Total stability calculations

# Map variance-scale bias factor B to its deviation-scale unbias divisor √B,
# treating non-positive B as NaN. Only :totvar can yield B ≤ 0 in practice
# (its B = 1 − a·τ/T formula goes negative if a caller passes an oversized
# τ — e.g., τ > T/a for FFM); the :mtot and :htot tables are bounded
# positive. Without this guard a single out-of-range τ would throw
# DomainError from sqrt and abort the whole deviation call instead of
# producing a NaN row for just that τ.
_unbias_divisor(B::Vector{Float64}) = [b > 0 ? sqrt(b) : NaN for b in B]

"""
    totdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Total Deviation for the given PhaseData, using the canonical
Howe 1995 / NIST SP1065 eqn 25 endpoint mean-flip extension (see `_totdev_core`).

`correct_bias=true` (default) applies the SP1065 noise-type-dependent
unbias correction `σ_unbiased = σ_raw / √B`, where `B = E[TOTVAR]/AVAR`
is the variance-scale bias factor from [`bias_correction`](@ref). TOTVAR
is biased low for FFM and RWFM as τ approaches T, so the correction
increases σ at long τ. Set `correct_bias=false` to return the raw
kernel value (Stable32 actually applies the correction by default,
contrary to older notes — verify against the build you compare with).
"""
function totdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _totdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0

    # Noise IDs needed for either path — bias correction reads α, CIs read α.
    need_noise = correct_bias || ci
    noises = need_noise ? identify_noise(x, m_values, dmin=0, dmax=2) : Symbol[]

    devs = correct_bias ? raw_devs ./ _unbias_divisor(bias_correction(noises, :totvar, taus, T)) : raw_devs

    if !ci
        return StabilityResult(:totdev, taus, devs, noises, Float64[], Float64[], Float64[])
    end

    edfs = calculate_edf(:totdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(x), confidence)

    return StabilityResult(:totdev, taus, devs, noises, lower, upper, edfs)
end

"""
    mtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Modified Total Deviation for the given PhaseData, using the
canonical Greenhall 2003 per-window time-reverse extension (see `_mtotdev_core`).

`correct_bias=true` (default) applies the SP1065 Table 11 unbias
correction `σ_unbiased = σ_raw / √B`, where `B ∈ {1.06, 1.17, 1.27,
1.30, 1.31}` for α ∈ {2, 1, 0, -1, -2}. MTOT is biased high, so the
correction drops σ by roughly 3 – 13 %. Pass `correct_bias=false` to
return the raw kernel — matches Stable32 and allantools, which do not
apply this correction.
"""
function mtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _mtotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0

    need_noise = correct_bias || ci
    noises = need_noise ? identify_noise(x, m_values, dmin=0, dmax=2) : Symbol[]

    devs = correct_bias ? raw_devs ./ _unbias_divisor(bias_correction(noises, :mtot, taus, T)) : raw_devs

    if !ci
        return StabilityResult(:mtotdev, taus, devs, noises, Float64[], Float64[], Float64[])
    end

    edfs = calculate_edf(:mtotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(x), confidence)

    return StabilityResult(:mtotdev, taus, devs, noises, lower, upper, edfs)
end

"""
    htotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Hadamard Total Deviation for the given PhaseData, using the
canonical Greenhall 2003 per-window time-reverse extension (see `_htotdev_core`).

`correct_bias=true` (default) applies the FCS 2001 (Howe & Tasset)
Table I unbias correction `σ_unbiased = σ_raw / √B`, where
`B = 1 + a ∈ {0.995, 0.851, 0.771, 0.717, 0.679}` for
α ∈ {0, -1, -2, -3, -4}. HTOT is biased *low* (a < 0), so the
correction raises σ — substantially for divergent FM. PM noises
(α > 0) get B = 1 (no published model). Pass
`correct_bias=false` to return the raw kernel value (Stable32 actually
applies the correction by default — see TOTDEV docstring for the same
caveat).
"""
function htotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _htotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0

    need_noise = correct_bias || ci
    noises = need_noise ? identify_noise(x, m_values, dmin=0, dmax=3) : Symbol[]

    devs = correct_bias ? raw_devs ./ _unbias_divisor(bias_correction(noises, :htot, taus, T)) : raw_devs

    if !ci
        return StabilityResult(:htotdev, taus, devs, noises, Float64[], Float64[], Float64[])
    end

    edfs = calculate_edf(:htotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(x), confidence)

    return StabilityResult(:htotdev, taus, devs, noises, lower, upper, edfs)
end

"""
    ttotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Time-Total Deviation. Wraps [`mtotdev`](@ref) and rescales by `τ/√3`,
analogous to [`tdev`](@ref) wrapping [`mdev`](@ref). TTOTDEV has units
of seconds (it is a `σ_x` quantity, not `σ_y`) and gives a
time-deviation summary of long-τ stability with MTOTDEV's
per-subsegment extended window — useful for telecom / time-transfer
analyses on records too short for ordinary TDEV at the τ of interest.

The `correct_bias`, `ci`, and `confidence` kwargs flow through to
the underlying `mtotdev` call unchanged. Confidence-interval bounds
inherit MTOTDEV's χ²/Gaussian limits scaled by the same `τ / √3` factor;
the EDF column is reused as-is (a time rescaling does not change the
degrees of freedom).
"""
function ttotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    res = mtotdev(data, m_values; ci=ci, correct_bias=correct_bias, confidence=confidence)
    factor = res.tau ./ sqrt(3.0)

    if !ci
        return StabilityResult(:ttotdev, res.tau, res.dev .* factor, Symbol[], Float64[], Float64[], Float64[])
    end

    return StabilityResult(:ttotdev, res.tau, res.dev .* factor, res.noise_type,
                           res.ci_lower .* factor, res.ci_upper .* factor, res.edf)
end

"""
    mhtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Modified Hadamard Total Deviation, using the Greenhall methodology SigmaTau
adopts for this novel estimator (see `_mhtotdev_core`).

`correct_bias=true` (default) applies the Monte-Carlo-measured unbias
correction `σ_unbiased = σ_raw / √B`, where `B = b0 + b1·(τ/T)` is the
τ/T-linear bias from [`bias_correction`](@ref)`(…, :mhtot, …)`. MHTOTDEV is
≈ unbiased for white/flicker noise (B ≈ 1) and reads progressively high for
redder FM (RWFM B ≈ 1.9 at small τ, where the correction lowers σ); it is
*not* exactly unbiased, contrary to the earlier assumption. Pass
`correct_bias=false` for the raw kernel value. EDF uses the measured
`_coeff_mhtot` fit. See the "MHTOTDEV bias and EDF" theory page for the
measurement methodology.
"""
function mhtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _mhtotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0

    need_noise = correct_bias || ci
    noises = need_noise ? identify_noise(x, m_values, dmin=0, dmax=3) : Symbol[]

    devs = correct_bias ? raw_devs ./ _unbias_divisor(bias_correction(noises, :mhtot, taus, T)) : raw_devs

    if !ci
        return StabilityResult(:mhtotdev, taus, devs, noises, Float64[], Float64[], Float64[])
    end

    edfs = calculate_edf(:mhtotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(x), confidence)

    return StabilityResult(:mhtotdev, taus, devs, noises, lower, upper, edfs)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
totdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)   = totdev(_freq_to_phase(data),   m_values; kwargs...)
mtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = mtotdev(_freq_to_phase(data),  m_values; kwargs...)
ttotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = ttotdev(_freq_to_phase(data),  m_values; kwargs...)
htotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = htotdev(_freq_to_phase(data),  m_values; kwargs...)
mhtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mhtotdev(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`.
totdev(data::PhaseData,       taus::TauMode; kwargs...) = totdev(data,   tau_values(taus, length(data.x), :totdev);   kwargs...)
totdev(data::FrequencyData,   taus::TauMode; kwargs...) = totdev(data,   tau_values(taus, length(data.y), :totdev);   kwargs...)
mtotdev(data::PhaseData,      taus::TauMode; kwargs...) = mtotdev(data,  tau_values(taus, length(data.x), :mtotdev);  kwargs...)
mtotdev(data::FrequencyData,  taus::TauMode; kwargs...) = mtotdev(data,  tau_values(taus, length(data.y), :mtotdev);  kwargs...)
ttotdev(data::PhaseData,      taus::TauMode; kwargs...) = ttotdev(data,  tau_values(taus, length(data.x), :ttotdev);  kwargs...)
ttotdev(data::FrequencyData,  taus::TauMode; kwargs...) = ttotdev(data,  tau_values(taus, length(data.y), :ttotdev);  kwargs...)
htotdev(data::PhaseData,      taus::TauMode; kwargs...) = htotdev(data,  tau_values(taus, length(data.x), :htotdev);  kwargs...)
htotdev(data::FrequencyData,  taus::TauMode; kwargs...) = htotdev(data,  tau_values(taus, length(data.y), :htotdev);  kwargs...)
mhtotdev(data::PhaseData,     taus::TauMode; kwargs...) = mhtotdev(data, tau_values(taus, length(data.x), :mhtotdev); kwargs...)
mhtotdev(data::FrequencyData, taus::TauMode; kwargs...) = mhtotdev(data, tau_values(taus, length(data.y), :mhtotdev); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to each kernel's
# algorithmic m-max (see `_default_m_values`). All kwargs pass through.
totdev(data::PhaseData;       kwargs...) = totdev(data,   _default_m_values(length(data.x), :totdev);   kwargs...)
totdev(data::FrequencyData;   kwargs...) = totdev(data,   _default_m_values(length(data.y), :totdev);   kwargs...)
mtotdev(data::PhaseData;      kwargs...) = mtotdev(data,  _default_m_values(length(data.x), :mtotdev);  kwargs...)
mtotdev(data::FrequencyData;  kwargs...) = mtotdev(data,  _default_m_values(length(data.y), :mtotdev);  kwargs...)
ttotdev(data::PhaseData;      kwargs...) = ttotdev(data,  _default_m_values(length(data.x), :ttotdev);  kwargs...)
ttotdev(data::FrequencyData;  kwargs...) = ttotdev(data,  _default_m_values(length(data.y), :ttotdev);  kwargs...)
htotdev(data::PhaseData;      kwargs...) = htotdev(data,  _default_m_values(length(data.x), :htotdev);  kwargs...)
htotdev(data::FrequencyData;  kwargs...) = htotdev(data,  _default_m_values(length(data.y), :htotdev);  kwargs...)
mhtotdev(data::PhaseData;     kwargs...) = mhtotdev(data, _default_m_values(length(data.x), :mhtotdev); kwargs...)
mhtotdev(data::FrequencyData; kwargs...) = mhtotdev(data, _default_m_values(length(data.y), :mhtotdev); kwargs...)

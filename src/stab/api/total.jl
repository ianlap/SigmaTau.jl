# api/total.jl — User wrappers for Total stability calculations

"""
    totdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:howe, calc_ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Total Deviation for the given PhaseData. See `_totdev_core` for
the meaning of `detrend`.

`correct_bias=true` (default) applies the SP1065 noise-type-dependent
bias factor `B(α)` to the raw kernel output. Set `correct_bias=false`
to return the raw kernel value, which matches Stable32 and allantools'
default outputs (neither of those tools applies the bias correction).
"""
function totdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:howe, calc_ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    raw_devs = _totdev_core(data.x, m_values, data.tau0; detrend=detrend)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    # Noise IDs needed for either path — bias correction reads α, CIs read α.
    need_noise = correct_bias || calc_ci
    noises = need_noise ? identify_noise(data.x, m_values, dmin=0, dmax=2) : Symbol[]

    devs = correct_bias ? raw_devs ./ bias_correction(noises, :totvar, taus, T) : raw_devs

    if !calc_ci
        return StabilityResult(:totdev, taus, devs, noises, Float64[], Float64[], Float64[])
    end

    edfs = calculate_edf(:totdev, devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:totdev, taus, devs, noises, lower, upper, edfs)
end

"""
    mtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Modified Total Deviation for the given PhaseData. See `_mtotdev_core` for
the meaning of `detrend`.

`correct_bias=true` (default) applies the SP1065 bias factor `B(α)` to
the raw kernel output (B≈1.27 for white FM). Pass `correct_bias=false`
to return the raw kernel — matches Stable32 and allantools, which do
not apply this correction. The resulting centerline drops by ~6–30%
depending on τ.
"""
function mtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    raw_devs = _mtotdev_core(data.x, m_values, data.tau0; detrend=detrend)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    need_noise = correct_bias || calc_ci
    noises = need_noise ? identify_noise(data.x, m_values, dmin=0, dmax=2) : Symbol[]

    devs = correct_bias ? raw_devs ./ bias_correction(noises, :mtot, taus, T) : raw_devs

    if !calc_ci
        return StabilityResult(:mtotdev, taus, devs, noises, Float64[], Float64[], Float64[])
    end

    edfs = calculate_edf(:mtotdev, devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:mtotdev, taus, devs, noises, lower, upper, edfs)
end

"""
    htotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Hadamard Total Deviation for the given PhaseData. See `_htotdev_core` for
the meaning of `detrend`.

`correct_bias=true` (default) applies the FCS 2001 bias factor (B≈1.005
for white FM, larger for divergent noises). Pass `correct_bias=false`
to match Stable32 and allantools, which report the raw kernel value.
"""
function htotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    raw_devs = _htotdev_core(data.x, m_values, data.tau0; detrend=detrend)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    need_noise = correct_bias || calc_ci
    noises = need_noise ? identify_noise(data.x, m_values, dmin=0, dmax=3) : Symbol[]

    devs = correct_bias ? raw_devs ./ bias_correction(noises, :htot, taus, T) : raw_devs

    if !calc_ci
        return StabilityResult(:htotdev, taus, devs, noises, Float64[], Float64[], Float64[])
    end

    edfs = calculate_edf(:htotdev, devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:htotdev, taus, devs, noises, lower, upper, edfs)
end

"""
    mhtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Modified Hadamard Total Deviation. See `_mhtotdev_core` for
the meaning of `detrend`.

No bias correction is applied for any value of `correct_bias` — the
kwarg is accepted for API symmetry with the other total-family
functions but is a documented no-op. FCS 2001 and NIST SP1065 publish
no bias-correction model for MHTOTDEV; the estimator is treated as
unbiased (B = 1) by policy, matching Stable32 and AllanLab.
`bias_correction(:mhtot, …)` returns ones for the same reason. EDF
uses the empirical SP1065 fit coefficients (`_coeff_mhtot`).
"""
function mhtotdev(data::PhaseData, m_values::Vector{Int}; detrend::Symbol=:greenhall, calc_ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    _ = correct_bias  # accepted for API symmetry; B = 1 by FCS 2001 / SP1065 policy. See docstring.
    raw_devs = _mhtotdev_core(data.x, m_values, data.tau0; detrend=detrend)
    taus = m_values .* data.tau0
    T = (length(data.x) - 1) * data.tau0

    if !calc_ci
        return StabilityResult(:mhtotdev, taus, raw_devs, Symbol[], Float64[], Float64[], Float64[])
    end

    noises = identify_noise(data.x, m_values, dmin=0, dmax=3)
    edfs = calculate_edf(:mhtotdev, raw_devs, noises, m_values, taus, length(data.x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(data.x), confidence)

    return StabilityResult(:mhtotdev, taus, raw_devs, noises, lower, upper, edfs)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
totdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)   = totdev(_freq_to_phase(data),   m_values; kwargs...)
mtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = mtotdev(_freq_to_phase(data),  m_values; kwargs...)
htotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = htotdev(_freq_to_phase(data),  m_values; kwargs...)
mhtotdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mhtotdev(_freq_to_phase(data), m_values; kwargs...)

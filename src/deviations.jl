# Public deviation API: PhaseData / FrequencyData in, StabilityResult out.
# Wraps the kernels in kernels.jl with noise identification, EDF and
# confidence intervals from noise.jl and edf.jl.

"""
    adev(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Overlapping Allan deviation σ_y(τ) of a phase record, per IEEE 1139-2022 §C.2.
The EDF for the χ² confidence interval uses the truncated-sum algorithm of
[Greenhall & Riley 2003](@cite greenhall-2003-edf-stability).

`m_values` selects the analysis-interval factors (τ = m·τ₀). When
`ci=true`, the result populates per-τ noise type, χ²-based confidence
bounds, and equivalent degrees of freedom.

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
function adev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _adev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:adev, length(x), m_values)

    if !ci
        return StabilityResult(:adev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:adev, raw_devs, noises, m_values, taus, length(x), (length(x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(x), confidence)

    return StabilityResult(:adev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    mdev(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Modified Allan deviation Mod σ_y(τ) of a phase record, per IEEE 1139-2022 §C.3.
The phase-averaged second difference separates white PM from flicker PM, which
`adev` cannot.

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
function mdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _mdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:mdev, length(x), m_values)

    if !ci
        return StabilityResult(:mdev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:mdev, raw_devs, noises, m_values, taus, length(x), (length(x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(x), confidence)

    return StabilityResult(:mdev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    tdev(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Time deviation σ_x(τ) = (τ/√3) · Mod σ_y(τ), per IEEE 1139-2022 §C.4. A σ_x
quantity in seconds. Confidence bounds are MDEV's scaled by the same factor.
"""
function tdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    res = mdev(data, m_values; ci=ci, confidence=confidence)
    factor = res.tau ./ sqrt(3.0)

    if !ci
        return StabilityResult(:tdev, res.tau, res.dev .* factor, Symbol[], _CIBound[], Float64[], res.neff)
    end

    return StabilityResult(:tdev, res.tau, res.dev .* factor, res.noise_type,
                           _scale_ci(res.ci, factor), res.edf, res.neff)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
adev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = adev(_freq_to_phase(data), m_values; kwargs...)
mdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mdev(_freq_to_phase(data), m_values; kwargs...)
tdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = tdev(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`.
adev(data::PhaseData,     taus::TauMode; kwargs...) = adev(data, tau_values(taus, length(data.x), :adev); kwargs...)
adev(data::FrequencyData, taus::TauMode; kwargs...) = adev(data, tau_values(taus, length(data.y), :adev); kwargs...)
mdev(data::PhaseData,     taus::TauMode; kwargs...) = mdev(data, tau_values(taus, length(data.x), :mdev); kwargs...)
mdev(data::FrequencyData, taus::TauMode; kwargs...) = mdev(data, tau_values(taus, length(data.y), :mdev); kwargs...)
tdev(data::PhaseData,     taus::TauMode; kwargs...) = tdev(data, tau_values(taus, length(data.x), :tdev); kwargs...)
tdev(data::FrequencyData, taus::TauMode; kwargs...) = tdev(data, tau_values(taus, length(data.y), :tdev); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to the kernel's algorithmic
# m-max (see `_default_m_values`). All kwargs pass through unchanged.
adev(data::PhaseData;     kwargs...) = adev(data, _default_m_values(length(data.x), :adev); kwargs...)
adev(data::FrequencyData; kwargs...) = adev(data, _default_m_values(length(data.y), :adev); kwargs...)
mdev(data::PhaseData;     kwargs...) = mdev(data, _default_m_values(length(data.x), :mdev); kwargs...)
mdev(data::FrequencyData; kwargs...) = mdev(data, _default_m_values(length(data.y), :mdev); kwargs...)
tdev(data::PhaseData;     kwargs...) = tdev(data, _default_m_values(length(data.x), :tdev); kwargs...)
tdev(data::FrequencyData; kwargs...) = tdev(data, _default_m_values(length(data.y), :tdev); kwargs...)


"""
    hdev(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Overlapping Hadamard deviation of a phase record, per NIST SP1065 §5.2.6. The
third difference removes linear frequency drift.
"""
function hdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _hdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:hdev, length(x), m_values)

    if !ci
        return StabilityResult(:hdev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=3)
    edfs = calculate_edf(:hdev, raw_devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(x), confidence)

    return StabilityResult(:hdev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    mhdev(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Modified Hadamard deviation of a phase record: the phase-averaged third
difference, drift-immune like `hdev` and able to separate white PM from
flicker PM like `mdev`.
"""
function mhdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _mhdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:mhdev, length(x), m_values)

    if !ci
        return StabilityResult(:mhdev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=3)
    edfs = calculate_edf(:mhdev, raw_devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(x), confidence)

    return StabilityResult(:mhdev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    htdev(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Hadamard time deviation σ_x,HT(τ) = τ/√(10/3) · σ_y,MH(τ), the σ_x form of
`mhdev` in seconds, as TDEV is of MDEV. The construction is original to this
package; SP1065, IEEE 1139-2022 and NBS TN-1337 do not define it.
"""
function htdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    res = mhdev(data, m_values; ci=ci, confidence=confidence)
    factor = res.tau ./ sqrt(10.0 / 3.0)

    if !ci
        return StabilityResult(:htdev, res.tau, res.dev .* factor, Symbol[], _CIBound[], Float64[], res.neff)
    end

    return StabilityResult(:htdev, res.tau, res.dev .* factor, res.noise_type,
                           _scale_ci(res.ci, factor), res.edf, res.neff)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
hdev(data::FrequencyData, m_values::Vector{Int}; kwargs...)  = hdev(_freq_to_phase(data),  m_values; kwargs...)
mhdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mhdev(_freq_to_phase(data), m_values; kwargs...)
htdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = htdev(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`.
hdev(data::PhaseData,      taus::TauMode; kwargs...) = hdev(data,  tau_values(taus, length(data.x), :hdev);  kwargs...)
hdev(data::FrequencyData,  taus::TauMode; kwargs...) = hdev(data,  tau_values(taus, length(data.y), :hdev);  kwargs...)
mhdev(data::PhaseData,     taus::TauMode; kwargs...) = mhdev(data, tau_values(taus, length(data.x), :mhdev); kwargs...)
mhdev(data::FrequencyData, taus::TauMode; kwargs...) = mhdev(data, tau_values(taus, length(data.y), :mhdev); kwargs...)
htdev(data::PhaseData,     taus::TauMode; kwargs...) = htdev(data, tau_values(taus, length(data.x), :htdev); kwargs...)
htdev(data::FrequencyData, taus::TauMode; kwargs...) = htdev(data, tau_values(taus, length(data.y), :htdev); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to each kernel's
# algorithmic m-max (see `_default_m_values`). All kwargs pass through.
hdev(data::PhaseData;      kwargs...) = hdev(data,  _default_m_values(length(data.x), :hdev);  kwargs...)
hdev(data::FrequencyData;  kwargs...) = hdev(data,  _default_m_values(length(data.y), :hdev);  kwargs...)
mhdev(data::PhaseData;     kwargs...) = mhdev(data, _default_m_values(length(data.x), :mhdev); kwargs...)
mhdev(data::FrequencyData; kwargs...) = mhdev(data, _default_m_values(length(data.y), :mhdev); kwargs...)
htdev(data::PhaseData;     kwargs...) = htdev(data, _default_m_values(length(data.x), :htdev); kwargs...)
htdev(data::FrequencyData; kwargs...) = htdev(data, _default_m_values(length(data.y), :htdev); kwargs...)


# Map variance-scale bias factor B to its deviation-scale unbias divisor √B,
# treating non-positive B as NaN. Only :totvar can yield B ≤ 0 in practice
# (its B = 1 − a·τ/T formula goes negative if a caller passes an oversized
# τ, e.g. τ > T/a for FLFM); the :mtot and :htot tables are bounded
# positive. Without this guard a single out-of-range τ would throw
# DomainError from sqrt and abort the whole deviation call instead of
# producing a NaN row for just that τ.
_unbias_divisor(B::AbstractVector{<:Real}) = [b > 0 ? sqrt(b) : NaN for b in B]

"""
    totdev(data, m_values; ci=true, correct_bias=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Total deviation of a phase record, per NIST SP1065 eq. 25 (Howe 1995), with
the whole-record mean-flip extension.

`correct_bias=true` (default) applies the SP1065 noise-type-dependent
unbias correction `σ_unbiased = σ_raw / √B`, where `B = E[TOTVAR]/AVAR`
is the variance-scale bias factor from [`bias_correction`](@ref). TOTVAR
is biased low for FLFM and RWFM as τ approaches T, so the correction
increases σ at long τ. Stable32 applies the same correction by default; `correct_bias=false`
returns the raw kernel value, which is what allantools reports.
"""
function totdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _totdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:totdev, length(x), m_values)

    # Noise IDs are needed for either path: bias correction and CIs both read α.
    need_noise = correct_bias || ci
    noises = need_noise ? identify_noise(x, m_values, dmin=0, dmax=2) : Symbol[]

    devs = correct_bias ? raw_devs ./ _unbias_divisor(bias_correction(noises, :totvar, taus, T)) : raw_devs

    if !ci
        return StabilityResult(:totdev, taus, devs, noises, _CIBound[], Float64[], neff)
    end

    edfs = calculate_edf(:totdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(x), confidence)

    return StabilityResult(:totdev, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    mtotdev(data, m_values; ci=true, correct_bias=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Modified total deviation of a phase record, per NIST SP1065 §5.2.11, with the
per-subsegment detrend and time-reverse extension of Howe and Vernotte.

`correct_bias=true` (default) applies the SP1065 Table 11 correction
`σ_unbiased = σ_raw / √B`, where `B ∈ {1.06, 1.17, 1.27, 1.30, 1.31}` for
α ∈ {2, 1, 0, -1, -2}. MTOT is biased high, so the correction lowers σ by 3
to 13 %. `correct_bias=false` returns the raw kernel, which is what Stable32
and allantools report.
"""
function mtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _mtotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:mtotdev, length(x), m_values)

    need_noise = correct_bias || ci
    noises = need_noise ? identify_noise(x, m_values, dmin=0, dmax=2) : Symbol[]

    devs = correct_bias ? raw_devs ./ _unbias_divisor(bias_correction(noises, :mtot, taus, T)) : raw_devs

    if !ci
        return StabilityResult(:mtotdev, taus, devs, noises, _CIBound[], Float64[], neff)
    end

    edfs = calculate_edf(:mtotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(x), confidence)

    return StabilityResult(:mtotdev, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    htotdev(data, m_values; ci=true, correct_bias=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Hadamard total deviation of a phase record, per NIST SP1065 §5.2.13, with the
per-subsegment detrend and time-reverse extension of Howe et al. 2000.

`correct_bias=true` (default) applies the FCS 2001 (Howe & Tasset)
Table I unbias correction `σ_unbiased = σ_raw / √B`, where
`B = 1 + a ∈ {0.995, 0.851, 0.771, 0.717, 0.679}` for
α ∈ {0, -1, -2, -3, -4}. HTOT is biased *low* (a < 0), so the
correction raises σ, substantially for divergent FM. PM noises (α > 0) get
B = 1 (no published model). Stable32 applies the same correction by default; `correct_bias=false`
returns the raw kernel value, which is what allantools reports.
"""
function htotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _htotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:htotdev, length(x), m_values)

    need_noise = correct_bias || ci
    noises = need_noise ? identify_noise(x, m_values, dmin=0, dmax=3) : Symbol[]

    devs = correct_bias ? raw_devs ./ _unbias_divisor(bias_correction(noises, :htot, taus, T)) : raw_devs

    if !ci
        return StabilityResult(:htotdev, taus, devs, noises, _CIBound[], Float64[], neff)
    end

    edfs = calculate_edf(:htotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(x), confidence)

    return StabilityResult(:htotdev, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    ttotdev(data, m_values; ci=true, correct_bias=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Time total deviation (τ/√3) · MTOTDEV, per NIST SP1065 §5.2.12: the σ_x form of
[`mtotdev`](@ref) in seconds, as [`tdev`](@ref) is of [`mdev`](@ref). Keyword
arguments pass through to `mtotdev`; confidence bounds are scaled by the same
factor and the EDF is unchanged.
"""
function ttotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    res = mtotdev(data, m_values; ci=ci, correct_bias=correct_bias, confidence=confidence)
    factor = res.tau ./ sqrt(3.0)

    if !ci
        return StabilityResult(:ttotdev, res.tau, res.dev .* factor, Symbol[], _CIBound[], Float64[], res.neff)
    end

    return StabilityResult(:ttotdev, res.tau, res.dev .* factor, res.noise_type,
                           _scale_ci(res.ci, factor), res.edf, res.neff)
end

"""
    mhtotdev(data, m_values; ci=true, correct_bias=true, remove_drift=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Modified Hadamard total deviation of a phase record: the total (boundary-
extended) form of [`mhdev`](@ref), using the per-subsegment detrend and
time-reverse extension of the other total estimators. Original to this
package.

`remove_drift=true` (default) removes the record's least-squares frequency
drift (the quadratic phase term) before the analysis, and the per-τ noise
identification runs on the drift-removed record. Unlike MHDEV, the total
estimator is not drift-immune: residual drift re-enters through the
reflected window boundaries and reads high at long τ. Removing a
deterministic drift discards no statistical information and is the practice
SP1065 prescribes before any total estimator. Pass `remove_drift=false` for
data that is already detrended.

`correct_bias=true` (default) applies the Monte-Carlo-measured unbias
correction `σ_unbiased = σ_raw / √B`, where `B = b0 + b1·(τ/T)` is the
τ/T-linear bias from [`bias_correction`](@ref)`(…, :mhtot, …)`; the
coefficients were measured with the default drift removal in the loop. Pass
`correct_bias=false` for the raw kernel value. The EDF uses the same Monte
Carlo fit; see the "MHTOTDEV bias and EDF" page for the method.
"""
function mhtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, remove_drift::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    if remove_drift
        x = _remove_quadratic(x)
    end
    raw_devs = _mhtotdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    T = (length(x) - 1) * data.tau0
    neff = _neff_counts(:mhtotdev, length(x), m_values)

    need_noise = correct_bias || ci
    noises = need_noise ? identify_noise(x, m_values, dmin=0, dmax=3) : Symbol[]

    devs = correct_bias ? raw_devs ./ _unbias_divisor(bias_correction(noises, :mhtot, taus, T)) : raw_devs

    if !ci
        return StabilityResult(:mhtotdev, taus, devs, noises, _CIBound[], Float64[], neff)
    end

    edfs = calculate_edf(:mhtotdev, devs, noises, m_values, taus, length(x), T)
    lower, upper = confidence_intervals(devs, edfs, noises, length(x), confidence)

    return StabilityResult(:mhtotdev, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
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


"""
    mtie(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Maximum time interval error (MTIE) of a phase record, per ITU-T G.810.
For each averaging factor `m`, returns the largest peak-to-peak phase
excursion observed in any sliding window of `m+1` samples (spanning
`τ = m·τ₀`).

MTIE is a σ_x quantity (units of seconds), reported as a single
deterministic envelope rather than a statistical σ; there is no
published EDF / χ² confidence model for it, so `noise_type`, `ci`,
and `edf` are returned empty even when `ci=true` (`neff` is still
populated). `ci` and `confidence` are accepted for API uniformity with
the other deviations but are no-ops here.

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
    raw_devs = _mtie_core(_f64(data.x), m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:mtie, length(data.x), m_values)
    return StabilityResult(:mtie, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
end

mtie(data::FrequencyData, m_values::Vector{Int}; kwargs...) = mtie(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`.
mtie(data::PhaseData,     taus::TauMode; kwargs...) = mtie(data, tau_values(taus, length(data.x), :mtie); kwargs...)
mtie(data::FrequencyData, taus::TauMode; kwargs...) = mtie(data, tau_values(taus, length(data.y), :mtie); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to MTIE's algorithmic
# m-max (`N − 1`, see `_default_m_values`).
mtie(data::PhaseData;     kwargs...) = mtie(data, _default_m_values(length(data.x), :mtie); kwargs...)
mtie(data::FrequencyData; kwargs...) = mtie(data, _default_m_values(length(data.y), :mtie); kwargs...)


"""
    tierms(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

RMS time interval error (TIE rms) of a phase record, per NIST SP1065
§5.2.18 eq. 37 [riley-2008-sp1065](@cite). For each averaging factor `m`,
returns the root-mean-square phase change over all `N − m` spans of
`τ = m·τ₀`:

```math
\\mathrm{TIE\\,rms}(m\\tau_0) = \\sqrt{\\frac{1}{N-m}
    \\sum_{i=1}^{N-m} \\bigl(x_{i+m} - x_i\\bigr)^2}
```

This is the estimator Stable32 and allantools (`tierms`) compute. TIE rms
is a σ_x quantity (units of seconds) used by the telecommunications
industry alongside MTIE; for a record with no frequency offset it
approaches the standard deviation of the fractional-frequency
fluctuations multiplied by τ, similar in behavior to TDEV.

Like `mtie`, TIE rms has no published EDF / χ² confidence model, so
`noise_type`, `ci`, and `edf` are returned empty even when `ci=true`
(`neff` carries the `N − m` span count). `ci` and `confidence` are
accepted for API uniformity with the other deviations but are no-ops
here.

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
    raw_devs = _tierms_core(_f64(data.x), m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:tierms, length(data.x), m_values)
    return StabilityResult(:tierms, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
end

tierms(data::FrequencyData, m_values::Vector{Int}; kwargs...) = tierms(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`.
tierms(data::PhaseData,     taus::TauMode; kwargs...) = tierms(data, tau_values(taus, length(data.x), :tierms); kwargs...)
tierms(data::FrequencyData, taus::TauMode; kwargs...) = tierms(data, tau_values(taus, length(data.y), :tierms); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to TIE rms' algorithmic
# m-max (`N − 1`, see `_default_m_values`).
tierms(data::PhaseData;     kwargs...) = tierms(data, _default_m_values(length(data.x), :tierms); kwargs...)
tierms(data::FrequencyData; kwargs...) = tierms(data, _default_m_values(length(data.y), :tierms); kwargs...)


"""
    pdev(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Parabolic deviation σ_PDEV(τ) of a phase record, per Vernotte–Lenczner–
Bourgeois–Rubiola (IEEE T-UFFC 63(4), 2016) and Vernotte–Chen–Rubiola 2020. PDEV
is built from a least-squares parabolic fit and is the recommended estimator when
evaluating the uncertainty of an ω-averaged frequency measurement. At `τ = τ₀`
(i.e. `m = 1`) PDEV reduces to overlapping ADEV.

When `ci=true`, the result populates per-τ noise type, EDF, and χ²-based
confidence bounds using the PVAR EDF model of Vernotte–Chen–Rubiola 2020
(arXiv:2005.13631), `ν ≈ 35 / (A(α)·(m/M) − 12·(m/M)²)` with `M = N − 2m`. At
`m = 1` the EDF is ADEV's (the PVAR(τ₀) ≡ AVAR(τ₀) identity); the large-τ region
(`m ≳ N/4`) uses the paper's semi-log interpolation down to `ν = 1` at `m = N/2`.
The approximation is accurate to ≈ ±10 % except at `m = 2`, where it reads high.
PVAR is an unbiased wavelet variance, so no bias correction is applied.

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
function pdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    raw_devs = _pdev_core(x, m_values, data.tau0)
    taus = m_values .* data.tau0
    neff = _neff_counts(:pdev, length(x), m_values)

    if !ci
        return StabilityResult(:pdev, taus, raw_devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:pdev, raw_devs, noises, m_values, taus, length(x), (length(x) - 1) * data.tau0)
    lower, upper = confidence_intervals(raw_devs, edfs, noises, length(x), confidence)

    return StabilityResult(:pdev, taus, raw_devs, noises, _zip_ci(lower, upper), edfs, neff)
end

pdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = pdev(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`.
pdev(data::PhaseData,     taus::TauMode; kwargs...) = pdev(data, tau_values(taus, length(data.x), :pdev); kwargs...)
pdev(data::FrequencyData, taus::TauMode; kwargs...) = pdev(data, tau_values(taus, length(data.y), :pdev); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to PDEV's algorithmic
# m-max (`(N − 2) ÷ 2`, see `_default_m_values`).
pdev(data::PhaseData;     kwargs...) = pdev(data, _default_m_values(length(data.x), :pdev); kwargs...)
pdev(data::FrequencyData; kwargs...) = pdev(data, _default_m_values(length(data.y), :pdev); kwargs...)


"""
    theo1(data, m_values; ci=true, correct_bias=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

Thêo1 deviation per NIST SP1065 (Riley 2008) §5.2.15 eq. 30, a
two-sample variance that mimics the Allan variance while extending the
usable averaging-time range to 75 % of the record span (see
`_theo1_core` for the kernel).

Conventions (SP1065 / Stable32):

- `m_values` are Thêo1 averaging factors and must be **even** with
  `2 ≤ m ≤ N − 1`; other entries produce NaN rows. SP1065 states the
  statistic for `m ≥ 10`. Use `tau_values(mode, N, :theo1)` for an
  even-rounded grid.
- The reported `r.tau` is the **effective** averaging time
  `τ = 0.75·m·τ0` (SP1065 §5.2.15), as plotted by Stable32. allantools'
  `theo1` computes the same deviation values (eq. 30 with the
  Howe & Peppler 0.75 normalization, which SP1065's printed eq. 30
  omits, although its own Table 2 requires it) but reports them at `τ = m·τ0`;
  SigmaTau resolves that definitional offset in favor of
  SP1065/Stable32, so `r.dev` matches allantools at the same `m` while
  `r.tau` is allantools' τ scaled by 0.75.
- `correct_bias=true` (default) applies the ThêoBR automatic bias
  removal of SP1065 §5.2.16 eq. 33: the Thêo1 variance is multiplied by
  the average of the `AVAR(m=9+3i)/Theo1(m=12+4i)` ratios over
  `i = 0…⌊N/6 − 3⌋` (matched-τ ladders spanning the long-τ ADEV overlap
  region), removing the noise-dependent Thêo1 bias without explicit
  noise identification. ThêoBR needs `N ≥ 19` samples; shorter records
  yield NaN rows (pass `correct_bias=false` for raw Thêo1, which
  Stable32 reports as its "Theo1" run).

The EDF for the χ²-based CI uses the empirical per-noise-type Thêo1
formulas of SP1065 §5.2.15 Table 3 (with `r = 0.75m`), for both the raw
and bias-removed forms. The kernel's double sum is O(N·m) per τ, far
costlier than the O(N) Allan kernels at large m.
"""
function theo1(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    N = length(x)
    raw_devs = _theo1_core(x, m_values, data.tau0)
    taus = 0.75 .* m_values .* data.tau0
    T = (N - 1) * data.tau0
    neff = _neff_counts(:theo1, N, m_values)

    devs = correct_bias ? raw_devs .* sqrt(_theobr_ratio(x, data.tau0)) : raw_devs

    if !ci
        return StabilityResult(:theo1, taus, devs, Symbol[], _CIBound[], Float64[], neff)
    end

    noises = identify_noise(x, m_values, dmin=0, dmax=2)
    edfs = calculate_edf(:theo1, devs, noises, m_values, taus, N, T)
    lower, upper = confidence_intervals(devs, edfs, noises, N, confidence)

    return StabilityResult(:theo1, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

"""
    theoh(data, m_values; ci=true, confidence=DEFAULT_CONFIDENCE) → StabilityResult

ThêoH ("hybrid-ThêoBR") per NIST SP1065 (Riley 2008) §5.2.16 eq. 34: a
composite σ_y(τ) curve that is overlapping ADEV at short averaging
times and the bias-removed Thêo1 (ThêoBR, eq. 33) beyond, extending a
single continuous stability plot from `τ0` out to 75 % of the record
span.

`m_values` are **τ-grid factors** (target `τ = m·τ0`, any positive
integers; unlike [`theo1`](@ref) they need not be even). Each point is
assigned per eq. 34's crossover `k = 20 %` of the record span
`T = (N−1)·τ0`:

- `m·τ0 ≤ k`: overlapping ADEV at factor `m`; reported at `τ = m·τ0`.
- `m·τ0 > k`: ThêoBR at the even Thêo1 factor nearest `4m/3`
  (exact whenever `m` is divisible by 3); reported at its true
  effective `τ = 0.75·m_theo·τ0`.

`r.tau` therefore always holds the realized τ of each point, which on
the Thêo1 segment can differ from the requested `m·τ0` by up to one
even-factor rounding step. The ThêoBR scaling needs `N ≥ 19`; shorter
records yield NaN on the Thêo1 segment. CIs use the ADEV
(Greenhall–Riley) EDF on the ADEV segment and the SP1065 Table 3 Thêo1
EDF on the Thêo1 segment. ThêoBR/ThêoH always carry the eq. 33 bias
removal; there is no `correct_bias` switch (raw Thêo1 is available via
`theo1(...; correct_bias=false)`).
"""
function theoh(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)
    x = _f64(data.x)
    N = length(x)
    tau0 = data.tau0
    T = (N - 1) * tau0
    n_pts = length(m_values)

    # Split the requested τ grid at the eq. 34 crossover (20 % of T).
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

    # Skip the ThêoBR ratio (itself a Thêo1 ladder, O(N²)) when no point
    # falls on the Thêo1 segment.
    theo_devs = isempty(theo_idx) ? Float64[] :
                _theo1_core(x, theo_ms, tau0) .* sqrt(_theobr_ratio(x, tau0))
    devs[theo_idx] = theo_devs
    taus[theo_idx] = 0.75 .* theo_ms .* tau0

    neff = _neff_counts(:theoh, N, m_values)

    if !ci
        return StabilityResult(:theoh, taus, devs, Symbol[], _CIBound[], Float64[], neff)
    end

    # Per-segment noise ID and EDF: ADEV EDF below the crossover, Thêo1
    # EDF above, merged back into request order.
    noises = fill(:WHFM, n_pts)
    edfs = fill(NaN, n_pts)
    if !isempty(adev_idx)
        noises_a = identify_noise(x, adev_ms, dmin=0, dmax=2)
        noises[adev_idx] = noises_a
        edfs[adev_idx] = calculate_edf(:adev, adev_devs, noises_a, adev_ms,
                                       taus[adev_idx], N, T)
    end
    if !isempty(theo_idx)
        noises_t = identify_noise(x, theo_ms, dmin=0, dmax=2)
        noises[theo_idx] = noises_t
        edfs[theo_idx] = calculate_edf(:theo1, theo_devs, noises_t, theo_ms,
                                       taus[theo_idx], N, T)
    end
    lower, upper = confidence_intervals(devs, edfs, noises, N, confidence)

    return StabilityResult(:theoh, taus, devs, noises, _zip_ci(lower, upper), edfs, neff)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
theo1(data::FrequencyData, m_values::Vector{Int}; kwargs...) = theo1(_freq_to_phase(data), m_values; kwargs...)
theoh(data::FrequencyData, m_values::Vector{Int}; kwargs...) = theoh(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve to the explicit m_values form via `tau_values`
# (the :theo1 grid is rounded to even averaging factors there).
theo1(data::PhaseData,     taus::TauMode; kwargs...) = theo1(data, tau_values(taus, length(data.x), :theo1); kwargs...)
theo1(data::FrequencyData, taus::TauMode; kwargs...) = theo1(data, tau_values(taus, length(data.y), :theo1); kwargs...)
theoh(data::PhaseData,     taus::TauMode; kwargs...) = theoh(data, tau_values(taus, length(data.x), :theoh); kwargs...)
theoh(data::FrequencyData, taus::TauMode; kwargs...) = theoh(data, tau_values(taus, length(data.y), :theoh); kwargs...)

# Zero-arg convenience: octave-spaced m_values up to each kernel's
# algorithmic m-max (see `_default_m_values`). All kwargs pass through.
theo1(data::PhaseData;     kwargs...) = theo1(data, _default_m_values(length(data.x), :theo1); kwargs...)
theo1(data::FrequencyData; kwargs...) = theo1(data, _default_m_values(length(data.y), :theo1); kwargs...)
theoh(data::PhaseData;     kwargs...) = theoh(data, _default_m_values(length(data.x), :theoh); kwargs...)
theoh(data::FrequencyData; kwargs...) = theoh(data, _default_m_values(length(data.y), :theoh); kwargs...)


# Dynamic (time-resolved) deviations.
# Both slide an analysis window of `window` phase samples across the record
# and evaluate an existing overlapping kernel (`_adev_core` / `_hdev_core`)
# inside each window, producing the σ_y(t, τ) map of Galleani & Tavella 2009
# as applied by McKelvy et al. 2025. No new kernel math is introduced: one
# window buffer is allocated, each window is copied into it, and the core is
# called per window. Simple and provably consistent with the static
# estimators; the cost is (number of windows) × (one static-kernel call on
# `window` samples), so halving `step` doubles the run time.

# Shared driver. `core` is one of the `_*_core` kernels; `min_window` is the
# smallest window that supports m = 1 for that kernel (window − span·1 ≥ 2).
function _dynamic_dev(core::F, dev_sym::Symbol, data::PhaseData,
                      m_values::Vector{Int};
                      window::Int, step::Int, min_window::Int) where {F}
    x = _f64(data.x)
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
    dadev(data, m_values; window, step=window÷2) → DynamicStabilityResult

Dynamic Allan deviation σ_y(t, τ), a time-resolved stability map
(Galleani & Tavella 2009; McKelvy et al. 2025
[mckelvy-2025-telemetrystability](@cite)). An analysis window `W(t_c)` of
`window` phase samples slides across the record in increments of `step`
samples, and the overlapping Allan deviation is computed inside each
window:

```math
\\sigma_y^2(t_c, \\tau) = \\frac{1}{2(N_w - 2m)(m\\tau_0)^2}
    \\sum_{i \\in W(t_c)} (x_{i+2m} - 2x_{i+m} + x_i)^2 .
```

`m_values` selects the averaging factors (τ = m·τ₀) shared by every
window; entries the window cannot support (`window − 2m < 2`) yield `NaN`
columns. Pass a [`TauMode`](@ref) instead to get a grid clamped to the
*window* length (`tau_values(mode, window, :adev)`), not the record length.
The result's rows are windows (centered at the times in `t`), columns are τ.

No confidence intervals are produced: the literature gives no EDF model for
the time-resolved map (each window's value alone is an ordinary ADEV with
`window` samples; use [`adev`](@ref) on a slice when a CI is needed).
Each window is evaluated by a full `_adev_core` call on a copied buffer, so
runtime scales with the window count: halving `step` doubles the cost.
"""
function dadev(data::PhaseData, m_values::Vector{Int};
               window::Int, step::Int = window ÷ 2)
    # min window: _adev_core needs N_w − 2m ≥ 2 ⇒ window ≥ 4 at m = 1.
    return _dynamic_dev(_adev_core, :dadev, data, m_values;
                        window=window, step=step, min_window=4)
end

"""
    dhdev(data, m_values; window, step=window÷2) → DynamicStabilityResult

Dynamic Hadamard deviation: the same sliding-window construction as
[`dadev`](@ref) with the overlapping third-difference (Hadamard) kernel
inside each window, so the map rejects linear frequency drift *within each
window* (a drifting clock shows its noise floor, not its drift ramp).
Entries the window cannot support (`window − 3m < 2`) yield `NaN` columns.
Pass a [`TauMode`](@ref) for a grid clamped to the window length
(`tau_values(mode, window, :hdev)`).

Like `dadev`, no confidence intervals are produced; there is no published
EDF model for the time-resolved map.
"""
function dhdev(data::PhaseData, m_values::Vector{Int};
               window::Int, step::Int = window ÷ 2)
    # min window: _hdev_core needs N_w − 3m ≥ 2 ⇒ window ≥ 5 at m = 1.
    return _dynamic_dev(_hdev_core, :dhdev, data, m_values;
                        window=window, step=step, min_window=5)
end

# FrequencyData entry points: convert via _freq_to_phase, dispatch to PhaseData.
dadev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = dadev(_freq_to_phase(data), m_values; kwargs...)
dhdev(data::FrequencyData, m_values::Vector{Int}; kwargs...) = dhdev(_freq_to_phase(data), m_values; kwargs...)

# TauMode grid selector: resolve against the *window* length (every window
# computes the same m grid; the record length is irrelevant to the clamp).
dadev(data::PhaseData,     taus::TauMode; window::Int, kwargs...) = dadev(data, tau_values(taus, window, :adev); window=window, kwargs...)
dadev(data::FrequencyData, taus::TauMode; window::Int, kwargs...) = dadev(data, tau_values(taus, window, :adev); window=window, kwargs...)
dhdev(data::PhaseData,     taus::TauMode; window::Int, kwargs...) = dhdev(data, tau_values(taus, window, :hdev); window=window, kwargs...)
dhdev(data::FrequencyData, taus::TauMode; window::Int, kwargs...) = dhdev(data, tau_values(taus, window, :hdev); window=window, kwargs...)

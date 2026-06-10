# deviations.jl — public deviation API (PhaseData/FrequencyData ->
# StabilityResult). Merged from stab/api/{allan, hadamard, total, mtie,
# pdev}.jl. Wraps the kernels in kernels.jl and the stats in edf.jl / noise.jl.


# ──────────────────────────────────────────────────────────────────────
# ── stab/api/allan.jl ───────────────────────────────────────────────────

# api/allan.jl — User wrappers for stability calculations

"""
$(SIGNATURES)

Overlapping Allan deviation σ_y(τ) for a phase record, per IEEE 1139-2022 §C.2.
EDF for the χ²-based CI uses the truncated-sum algorithm of [Greenhall & Riley 2003](@cite greenhall-2003-edf-stability).

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
$(SIGNATURES)

Modified Allan deviation Mod σ_y(τ) for a phase record, per IEEE 1139-2022 §C.3.

Uses a phase-averaged second difference; better than `adev` at separating
white-PM from flicker-PM noise.

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
    tdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Time Deviation (TDEV) for the given PhaseData.

TDEV has units of seconds (it is a σ_x quantity, not σ_y), defined as
`σ_x(τ) = (τ/√3) · σ_y,MDEV(τ)`. Confidence-interval bounds inherit
MDEV's χ²/Gaussian limits scaled by the same `τ/√3` factor.
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


# ──────────────────────────────────────────────────────────────────────
# ── stab/api/hadamard.jl ────────────────────────────────────────────────

# api/hadamard.jl — User wrappers for Hadamard stability calculations

"""
    hdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Overlapping Hadamard Deviation for the given PhaseData.
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
    mhdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Modified Hadamard Deviation for the given PhaseData.
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
    htdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Computes the Hadamard time deviation (HTDEV) for the given `PhaseData`
— a third-difference time deviation, related to MHDEV by
`σ_x,HT(τ) = τ / √(10/3) · σ_y,MH(τ)`. HTDEV has units of seconds
(it is a σ_x quantity); HTDEV is to MHDEV what TDEV is to MDEV. The
construction is original to this package; SP1065, IEEE 1139-2022,
and NBS-TN-1337 do not define it.
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


# ──────────────────────────────────────────────────────────────────────
# ── stab/api/total.jl ───────────────────────────────────────────────────

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
    neff = _neff_counts(:totdev, length(x), m_values)

    # Noise IDs needed for either path — bias correction reads α, CIs read α.
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
        return StabilityResult(:ttotdev, res.tau, res.dev .* factor, Symbol[], _CIBound[], Float64[], res.neff)
    end

    return StabilityResult(:ttotdev, res.tau, res.dev .* factor, res.noise_type,
                           _scale_ci(res.ci, factor), res.edf, res.neff)
end

"""
    mhtotdev(data::PhaseData, m_values::Vector{Int}; ci::Bool=true, correct_bias::Bool=true, remove_drift::Bool=true, confidence::Float64=DEFAULT_CONFIDENCE)

Modified Hadamard Total Deviation, using the Greenhall methodology SigmaTau
adopts for this novel estimator (see `_mhtotdev_core`).

`remove_drift=true` (default) removes the record's least-squares frequency
drift (the quadratic phase term) before the analysis, and the per-τ noise
identification runs on the drift-removed record. Unlike MHDEV, the total
estimator is not drift-immune: its per-window detrend removes only the local
frequency offset, so residual deterministic drift re-enters through the
reflected window boundaries and reads several times high at long τ. Removing
the global least-squares drift is exact for a deterministic drift (the fit
is a linear projection) and discards no statistical information; it is the
same practice SP1065 prescribes before TOTVAR / MTOTVAR. Per-window drift
removal was evaluated and rejected — both the phase-domain and
frequency-domain variants measurably damage the statistic at the PM noise
types (see the "MHTOTDEV bias and EDF" theory page). Pass
`remove_drift=false` for data that is already detrended.

`correct_bias=true` (default) applies the Monte-Carlo-measured unbias
correction `σ_unbiased = σ_raw / √B`, where `B = b0 + b1·(τ/T)` is the
τ/T-linear bias from [`bias_correction`](@ref)`(…, :mhtot, …)`; the
coefficients are measured with the default drift removal in the loop, so the
calibration matches what this function computes. Pass `correct_bias=false`
for the raw kernel value. EDF uses the measured `_coeff_mhtot` fit. See the
"MHTOTDEV bias and EDF" theory page for the measurement methodology.
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


# ──────────────────────────────────────────────────────────────────────
# ── stab/api/mtie.jl ────────────────────────────────────────────────────

# api/mtie.jl — User wrapper for Maximum Time Interval Error

"""
$(SIGNATURES)

Maximum Time Interval Error (MTIE) for a phase record, per ITU-T G.810.
For each averaging factor `m`, returns the largest peak-to-peak phase
excursion observed in any sliding window of `m+1` samples (spanning
`τ = m·τ₀`).

MTIE is a σ_x quantity (units of seconds), reported as a single
deterministic envelope rather than a statistical σ — there is no
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


# ──────────────────────────────────────────────────────────────────────
# ── stab/api/pdev.jl ────────────────────────────────────────────────────

# api/pdev.jl — User wrapper for parabolic deviation

"""
$(SIGNATURES)

Parabolic deviation σ_PDEV(τ) for a phase record, per Vernotte–Lenczner–
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
# m-max (`(N − 1) ÷ 2`, see `_default_m_values`).
pdev(data::PhaseData;     kwargs...) = pdev(data, _default_m_values(length(data.x), :pdev); kwargs...)
pdev(data::FrequencyData; kwargs...) = pdev(data, _default_m_values(length(data.y), :pdev); kwargs...)

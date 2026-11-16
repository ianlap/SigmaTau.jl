module SigmaTauRecipesBaseExt

using SigmaTau
using SigmaTau: StabilityResult, StabilitySuite, SpectralResult,
                DynamicStabilityResult, HatResult, ci_lower, ci_upper, estimators,
                PowerLawFit, HoldoverResult, HoldoverPrediction, SteeringResult,
                SteeringPrediction, TimescaleResult, KalmanResult, powerlaw_sigma,
                _estimator_label
using RecipesBase

# Tick positions at integer powers of 10 spanning the positive finite values
# in `vals` (any iterator of numbers). Passed as default `xticks`/`yticks` so
# a narrow-range log axis keeps integer-decade labels instead of falling back
# to fractional powers (10^2.5-style). Returns `:auto` when no positive finite
# value exists, deferring to the backend.
function _decade_ticks(vals)
    lo, hi = Inf, -Inf
    for v in vals
        (isfinite(v) && v > 0) || continue
        lo = min(lo, v)
        hi = max(hi, v)
    end
    lo <= hi || return :auto
    return 10.0 .^ (floor(log10(lo)):ceil(log10(hi)))
end

# Ordinate values of a result: the curve plus its CI bounds (empty when
# computed with `ci=false`), so default y ticks span the error bars too.
_ordinate_values(r::StabilityResult) = Iterators.flatten((r.dev, ci_lower(r), ci_upper(r)))

# Single deviation curve on log-log axes. When CI bounds are present they render
# as error bars by default, or as a filled band when the `ci_band` attribute is
# set (e.g. `plot(result; ci_band=true)`). `ci_band` is consumed here so it does
# not leak to the backend as an unknown attribute.
@recipe function f(res::StabilityResult)
    xscale --> :log10
    yscale --> :log10
    xlabel --> "Averaging Time τ (s)"
    ylabel --> _estimator_label(res.deviation_type)
    label  --> _estimator_label(res.deviation_type)
    seriestype := :path
    markershape --> :circle
    markersize  --> 3
    xticks --> _decade_ticks(res.tau)
    yticks --> _decade_ticks(_ordinate_values(res))
    minorgrid --> true
    gridalpha --> 0.2
    minorgridalpha --> 0.05

    ci_band = pop!(plotattributes, :ci_band, false)
    if !isempty(res.ci)
        lo = res.dev .- ci_lower(res)
        hi = ci_upper(res) .- res.dev
        if ci_band
            ribbon := (lo, hi)
        else
            yerror := (lo, hi)
        end
    end

    return res.tau, res.dev
end

# Estimators in seconds: TDEV and the time errors.
const _TIME = Set(e.name for e in estimators() if e.axis === :sigma_x)

# Label of curve `i` of a collection: the clock's name for a hat, the
# estimator otherwise.
_curve_label(h::HatResult, i) = h.names[i]
_curve_label(results, i) = _estimator_label(results[i].deviation_type)

# A suite draws as the vector of its results.
@recipe f(suite::StabilitySuite) = suite.results

# Overlay several results (e.g. to compare clocks): the σ_y estimators on one
# pair of axes and those in seconds on a second beside it, each with decade
# ticks spanning its own curves.
@recipe function f(results::AbstractVector{StabilityResult})
    xscale --> :log10
    yscale --> :log10
    xlabel --> "Averaging Time τ (s)"
    xticks --> _decade_ticks(Iterators.flatten(r.tau for r in results))
    minorgrid --> true
    gridalpha --> 0.2
    minorgridalpha --> 0.05
    time = [r.deviation_type in _TIME for r in results]
    split = any(time) && !all(time)
    if split
        layout --> (1, 2)
        size --> (800, 360)
        left_margin --> (5, :mm)        # room for the axis titles beside two panels
        bottom_margin --> (5, :mm)
    end
    for (i, r) in enumerate(results)
        same_axis = results[time .== time[i]]
        @series begin
            if split
                subplot := 1 + time[i]
            end
            label --> _curve_label(results, i)
            ylabel --> (time[i] ? "Deviation (s)" : "Deviation")
            yticks --> _decade_ticks(Iterators.flatten(_ordinate_values(q) for q in same_axis))
            r
        end
    end
end

# Dynamic deviation (DADEV / DHDEV) in the SP1065 presentation: one σ(τ)
# curve per window time, stacked along the time axis as 3-D lines. Plots'
# 3-D axes do not support log scales, so log10(τ) and log10(σ) are the
# plotted coordinates and the decade ticks are labeled in plain-text
# `1e<k>` form. Non-positive and NaN cells (windows that cannot support the
# averaging factor) map to NaN and terminate the curve.
@recipe function f(res::DynamicStabilityResult)
    title --> _estimator_label(res.deviation_type)
    z = map(v -> (isfinite(v) && v > 0) ? log10(v) : NaN, res.dev)
    logτ = log10.(res.tau)
    kx = floor(Int, minimum(logτ)):ceil(Int, maximum(logτ))
    zfin = filter(isfinite, z)
    kz = isempty(zfin) ? (0:0) : (floor(Int, minimum(zfin)):ceil(Int, maximum(zfin)))
    xlabel --> "Averaging Time τ (s)"
    ylabel --> "Time t (s)"
    xticks --> (Float64.(kx), ["1e$k" for k in kx])
    zticks --> (Float64.(kz), ["1e$k" for k in kz])
    zguide --> "σ"
    # Legend gets noisy past a handful of windows; label only small maps.
    labelled = length(res.t) <= 8
    for (i, tc) in enumerate(res.t)
        @series begin
            seriestype := :path3d
            label := labelled ? "t = $(round(tc; sigdigits=4)) s" : ""
            logτ, fill(tc, length(logτ)), z[i, :]
        end
    end
    return nothing
end

# Spectral density on frequency axes. S_y / S_x are power densities → log-log;
# ℒ(f) is already in dBc/Hz → log frequency axis with a linear (dB) ordinate,
# so only its frequency axis gets the decade-tick treatment.
@recipe function f(res::SpectralResult)
    xscale --> :log10
    xlabel --> "Frequency f (Hz)"
    seriestype := :path
    # Drop the DC bin (f = 0): undefined on a log frequency axis. `Sy`/`Sx`
    # carry it at index 1; `L` already excludes it in the API.
    keep = res.freq .> 0
    xticks --> _decade_ticks(res.freq[keep])
    gridalpha --> 0.2
    minorgridalpha --> 0.05
    if res.spectral_type === :L
        label  --> "ℒ(f)"
        ylabel --> "ℒ(f) (dBc/Hz)"
        xminorgrid --> true
    else
        yscale --> :log10
        label  --> string(res.spectral_type)
        ylabel --> (res.spectral_type === :Sy ? "S_y(f) (1/Hz)" : "S_x(f) (s²/Hz)")
        yticks --> _decade_ticks(res.psd[keep])
        minorgrid --> true
    end
    return res.freq[keep], res.psd[keep]
end

# ── Models, holdover, steering, time scales and filters ──────────────────

# Noise name of exponent α of S_y(f) ∝ f^α.
_noise_name(α::Integer) = String(SigmaTau._NOISE_SYMBOLS[α + 5])

# Level of a quantile, for a legend.
_level(c) = (v = round(100c; digits=1); string(isinteger(v) ? Int(v) : v, " %"))

# A power-law fit: the measured points, the model through them, and each
# noise type of the model alone, dashed. A spectrum of a sampled record is
# fitted as Σ h_α g^α with g = sin(πfτ₀)/(πτ₀), so its components are too.
@recipe function f(fit::PowerLawFit)
    spectrum = fit.source in (:Sy, :Sx, :Sphi, :L)
    g = isnan(fit.tau0) ? fit.x : sin.(π .* fit.x .* fit.tau0) ./ (π * fit.tau0)
    xscale --> :log10
    yscale --> :log10
    xlabel --> (spectrum ? "Frequency f (Hz)" : "Averaging Time τ (s)")
    ylabel --> (spectrum ? "S_y(f) (1/Hz)" : _estimator_label(fit.source))
    xticks --> _decade_ticks(fit.x)
    yticks --> _decade_ticks(Iterators.flatten((fit.y, fit.fit)))
    minorgrid --> true
    gridalpha --> 0.2
    minorgridalpha --> 0.05
    @series begin
        seriestype := :scatter
        label --> "measured"
        fit.x, fit.y
    end
    @series begin
        label --> "fit"
        fit.x, fit.fit
    end
    for (α, h) in zip(fit.alphas, fit.h)
        h > 0 || continue
        @series begin
            linestyle --> :dash
            label --> _noise_name(α)
            fit.x, spectrum ? h .* g .^ α :
                powerlaw_sigma(Dict(α => h), fit.x; deviation=fit.source, tau0=fit.tau0)
        end
    end
end

# Time error against holdover time on log-log axes, with the limit as a
# horizontal line when one was asked for.
@recipe function f(h::Union{HoldoverResult,HoldoverPrediction})
    xscale --> :log10
    yscale --> :log10
    xlabel --> "Holdover Time τ (s)"
    ylabel --> "Time Error (s)"
    xticks --> _decade_ticks(h.tau)
    minorgrid --> true
    gridalpha --> 0.2
    minorgridalpha --> 0.05
    curves = h isa HoldoverResult ?
        [("RMS", h.rms, :solid), (_level(h.confidence) * " of trials within", h.quantile, :solid),
         ("worst", h.worst, :dot), ("Kalman σ", h.claimed, :dash)] :
        [[("σ", h.sigma, :solid), (_level(h.confidence) * " within", h.quantile, :solid)];
         [(_noise_name(α), sqrt.(h.budget[:, j]), :dash) for (j, α) in enumerate(h.alphas)]]
    filter!(c -> !isempty(c[2]), curves)        # `claimed` is empty for a least-squares fit
    yticks --> _decade_ticks(Iterators.flatten(c[2] for c in curves))
    for (name, y, style) in curves
        @series begin
            seriestype := :path
            linestyle --> style
            label --> name
            h.tau, y
        end
    end
    if isfinite(h.limit)
        @series begin
            seriestype := :hline
            linestyle --> :dashdot
            color --> :gray
            label --> "limit"
            [h.limit]
        end
    end
end

# The steady state of a loop: the Allan deviation of the steered clock.
@recipe function f(p::SteeringPrediction)
    xscale --> :log10
    yscale --> :log10
    xlabel --> "Averaging Time τ (s)"
    ylabel --> "ADEV"
    label --> "steered, predicted"
    xticks --> _decade_ticks(p.tau)
    yticks --> _decade_ticks(p.adev)
    minorgrid --> true
    p.tau, p.adev
end

# One panel of a time series: a line per column of `y`, and a band of
# ± `band` about the first when `band` is not empty.
struct _Panel
    ylabel::String
    y::Matrix{Float64}
    labels::Vector{String}
    band::Vector{Float64}
end
_Panel(ylabel, y::AbstractVector, label; band=Float64[]) = _Panel(ylabel, reshape(y, :, 1), [label], band)
_Panel(ylabel, y::AbstractMatrix, labels) = _Panel(ylabel, y, labels, Float64[])

# Panels stacked over one time axis.
struct _Panels
    t::Vector{Float64}
    panels::Vector{_Panel}
end

@recipe function f(p::_Panels)
    n = length(p.panels)
    layout --> (n, 1)
    link --> :x
    size --> (700, 220n + 60)
    left_margin --> (5, :mm)        # room for the axis titles
    bottom_margin --> (5, :mm)
    for (i, panel) in enumerate(p.panels), j in axes(panel.y, 2)
        @series begin
            subplot := i
            seriestype := :path
            ylabel --> panel.ylabel
            xlabel --> (i == n ? "Time (s)" : "")
            label --> panel.labels[j]
            if j == 1 && !isempty(panel.band)
                ribbon := panel.band
                fillalpha --> 0.25
            end
            p.t, panel.y[:, j]
        end
    end
end

# A steered record: its phase, and the frequency correction in force.
@recipe function f(r::SteeringResult)
    legend --> false
    t = (0:length(r.correction)-1) .* r.steered.tau0
    return _Panels(t, [_Panel("Phase (s)", r.steered.x, "steered"),
                       _Panel("Frequency Correction", r.correction, "correction")])
end

# A time scale: each clock against the scale and, for the algorithms that
# have them, the clocks' weights.
@recipe function f(r::TimescaleResult)
    names = r.comparisons.names
    panels = [_Panel("Clock − Scale (s)", r.x, names)]
    all(isnan, r.weights) || push!(panels, _Panel("Weight", r.weights, fill("", length(names))))
    return _Panels(r.t, panels)
end

# A Kalman filter: the phase and the frequency it estimates, with a band of
# one standard deviation. The frequency includes the Markov frequency states.
@recipe function f(r::KalmanResult)
    np = SigmaTau._npoly(r.model)
    freq = [2; np+1:np+length(SigmaTau._fm_taus(r.model))]
    σ(c) = [sqrt(max(sum(view(r.P, c, c, k)), 0.0)) for k in axes(r.P, 3)]
    legend --> false
    return _Panels(r.t, [_Panel("Phase (s)", r.x[:, 1], "phase"; band=σ([1])),
                         _Panel("Frequency", vec(sum(view(r.x, :, freq); dims=2)), "frequency"; band=σ(freq))])
end

end # module

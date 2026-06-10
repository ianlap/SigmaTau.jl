module SigmaTauRecipesBaseExt

using SigmaTau: StabilityResult, StabilitySuite, SpectralResult,
                DynamicStabilityResult, ci_lower, ci_upper
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
    ylabel --> uppercase(string(res.deviation_type))
    label  --> uppercase(string(res.deviation_type))
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

# Overlay every result in a suite on one set of axes (compare deviations).
# Decade ticks span the union of the per-result ranges, set here so all
# series share one axis policy.
@recipe function f(suite::StabilitySuite)
    xscale --> :log10
    yscale --> :log10
    xlabel --> "Averaging Time τ (s)"
    ylabel --> "Deviation"
    xticks --> _decade_ticks(Iterators.flatten(r.tau for r in suite.results))
    yticks --> _decade_ticks(Iterators.flatten(_ordinate_values(r) for r in suite.results))
    minorgrid --> true
    gridalpha --> 0.2
    minorgridalpha --> 0.05
    for r in suite.results
        @series begin
            label --> uppercase(string(r.deviation_type))
            r
        end
    end
end

# Overlay an arbitrary collection of results (e.g. compare clocks).
@recipe function f(results::AbstractVector{StabilityResult})
    xscale --> :log10
    yscale --> :log10
    xlabel --> "Averaging Time τ (s)"
    ylabel --> "Deviation"
    xticks --> _decade_ticks(Iterators.flatten(r.tau for r in results))
    yticks --> _decade_ticks(Iterators.flatten(_ordinate_values(r) for r in results))
    minorgrid --> true
    gridalpha --> 0.2
    minorgridalpha --> 0.05
    for r in results
        @series begin
            label --> uppercase(string(r.deviation_type))
            r
        end
    end
end

# Dynamic deviation map (DADEV / DHDEV): log10(σ) over (t, τ). Default
# rendering is a heatmap with a log-scale τ axis; pass
# `seriestype = :path3d` for the SP1065-style dynamic-deviation view — one
# σ(τ) curve per window time, stacked along the time axis as 3-D lines.
# Plots' 3-D axes do not support log scales, so the 3-D branch plots
# log10(τ) and log10(σ) as the coordinates and labels the decade ticks in
# plain-text `1e<k>` form — safe under both GR and the PGFPlotsX/lualatex
# docs build. The map is stored windows × taus; the heatmap z matrix is
# rows-along-y, so it is transposed there. Non-positive and NaN cells
# (windows that cannot support the averaging factor) map to NaN and render
# blank (heatmap) or terminate the curve (3-D lines).
@recipe function f(res::DynamicStabilityResult)
    title --> uppercase(string(res.deviation_type))
    z = map(v -> (isfinite(v) && v > 0) ? log10(v) : NaN, res.dev)

    st = get(plotattributes, :seriestype, :heatmap)
    if st === :path3d || st === :path3D
        logτ = log10.(res.tau)
        kx = floor(Int, minimum(logτ)):ceil(Int, maximum(logτ))
        zfin = filter(isfinite, z)
        kz = isempty(zfin) ? (0:0) :
             (floor(Int, minimum(zfin)):ceil(Int, maximum(zfin)))
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
    else
        xlabel --> "Time t (s)"
        ylabel --> "Averaging Time τ (s)"
        colorbar_title --> "log10 σ"
        yscale --> :log10
        yticks --> _decade_ticks(res.tau)
        @series begin
            seriestype := :heatmap
            label := ""
            res.t, res.tau, permutedims(z)
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

end # module

module SigmaTauRecipesBaseExt

using SigmaTau: StabilityResult, StabilitySuite
using RecipesBase

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

    ci_band = pop!(plotattributes, :ci_band, false)
    if !isempty(res.ci_lower) && !isempty(res.ci_upper)
        lo = res.dev .- res.ci_lower
        hi = res.ci_upper .- res.dev
        if ci_band
            ribbon := (lo, hi)
        else
            yerror := (lo, hi)
        end
    end

    return res.tau, res.dev
end

# Overlay every result in a suite on one set of axes (compare deviations).
@recipe function f(suite::StabilitySuite)
    xscale --> :log10
    yscale --> :log10
    xlabel --> "Averaging Time τ (s)"
    ylabel --> "Deviation"
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
    for r in results
        @series begin
            label --> uppercase(string(r.deviation_type))
            r
        end
    end
end

end # module

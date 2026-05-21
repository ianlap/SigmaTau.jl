module SigmaTauRecipesBaseExt

using SigmaTau: StabilityResult
using RecipesBase

@recipe function f(res::StabilityResult)
    xscale --> :log10
    yscale --> :log10
    xlabel --> "Averaging Time τ (s)"
    ylabel --> uppercase(string(res.deviation_type))
    label  --> uppercase(string(res.deviation_type))
    seriestype := :path

    if !isempty(res.ci_lower) && !isempty(res.ci_upper)
        yerror := (res.dev .- res.ci_lower, res.ci_upper .- res.dev)
    end

    return res.tau, res.dev
end

end # module

# PlotRecipes.jl
# We will define @recipe functions for StabilityResult here.

# using RecipesBase

# @recipe function f(res::StabilityResult)
#     xscale --> :log10
#     yscale --> :log10
#     xlabel --> "Averaging Time τ (s)"
#     ylabel --> string(res.deviation_type)
#     yerror --> (res.dev .- res.ci_lower, res.ci_upper .- res.dev)
#     
#     return res.tau, res.dev
# end

# # Computing the Allan deviation
#
# `adev` returns a fully-populated [`StabilityResult`](@ref): the σ_y(τ)
# point estimates, a power-law-noise classification per τ, the
# equivalent number of degrees of freedom, and χ² confidence bounds —
# all in one call. There's no separate "compute the noise type" step or
# "compute the CI" step; they're part of the default output.
#
# This tutorial:
#
# 1. Runs `adev` on the WPM + RWFM fixture from
#    [Phase data and frequency data](01_phase_data.md).
# 2. Walks through every field of the returned `StabilityResult`.
# 3. Renders σ_y(τ) on log-log axes with χ² error bars via the bundled
#    `RecipesBase` extension.

using SigmaTau
using Random

# ## Fixture
#
# Same deterministic seed as the previous tutorial so the figures
# line up.

Random.seed!(20260509)

τ₀ = 1.0
N  = 4096
wpm  = randn(N) .* 1e-9
rwfm = cumsum(cumsum(randn(N) .* 1e-12))
x    = wpm .+ rwfm
pd   = PhaseData(x, τ₀)

# ## Run `adev` on an octave-spaced τ-grid
#
# The averaging factor `m` selects which τ values get sampled — the
# package multiplies by `τ₀` internally to produce `result.tau`. An
# octave grid (`m = 1, 2, 4, …`) covers a wide range with few points
# and matches the convention in NIST SP1065.

m_grid = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
result = adev(pd, m_grid)

# ## What's in the result?
#
# Every relevant quantity is on the returned struct — there's no
# follow-up call needed for noise classification, EDF, or confidence
# bounds.

result.tau          # τ values (s)

#-

result.dev          # σ_y(τ) point estimates

#-

result.noise_type   # power-law classification per τ
                    # (:WHPM, :FLPM, :WHFM, :FLFM, :RWFM)

#-

result.edf          # equivalent degrees of freedom
                    # (Greenhall–Riley formula, used by the χ² CI)

#-

result.ci_lower     # lower χ² confidence bound (default 95%)

#-

result.ci_upper     # upper χ² confidence bound

#-

result.deviation_type  # which deviation kernel produced this result

# Skipping CI computation (when you just want the point estimates) is
# a single keyword:
#
# ```julia
# adev(pd, m_grid; calc_ci = false)
# ```
#
# `result.ci_lower`, `result.ci_upper`, and `result.edf` come back
# empty in that mode.

# ## The plot recipe
#
# `SigmaTau` ships a `RecipesBase` extension
# (`ext/SigmaTauRecipesBaseExt.jl`); loading any `Plots`-compatible
# backend brings in a log-log τ–σ recipe with χ² error bars sourced
# from `result.ci_lower / ci_upper`.

using Plots

plot(result;
     title  = raw"Overlapping Allan deviation (WPM + RWFM mix)",
     xlabel = raw"Averaging time \(\tau\) (s)",
     ylabel = raw"\(\sigma_y(\tau)\)",
     legend = :topright,
     lw     = 1.5)

# ## Reading the curve
#
# The slope of σ_y(τ) on log-log axes encodes the dominant noise type.
# At short τ the WPM contribution dominates and the curve falls as
# τ⁻¹; past a knee, the RWFM contribution takes over and rises as
# τ^(1/2). The `result.noise_type` vector tracks this regime change
# automatically; in the test suite,
# `test/stab/runtests.jl` (`Stable32 cross-validation`) verifies these
# σ_y values match Stable32's reference output to four-five
# significant figures across the same fixture.

# # Three-cornered hat: separating individual clock noises from pairwise comparisons
#
# A comparison of two clocks measures the sum of their noises. With three
# clocks and the three pairwise difference records, each clock's own noise
# can be recovered, provided the three noises are independent: the pairwise
# Allan variances obey
#
# ```math
# \sigma^2_{ij}(\tau) = \sigma^2_i(\tau) + \sigma^2_j(\tau),
# ```
#
# three equations in three unknowns, so for example
# ``\sigma^2_1 = \tfrac{1}{2}(\sigma^2_{12} + \sigma^2_{13} - \sigma^2_{23})``.
# `nch` solves this system, and its N-clock generalization, from the
# pairwise `StabilityResult`s. Because the inputs are synthesized here, the
# recovery can be checked against each clock's true σ_y(τ).

using SigmaTau
using Random

# ## Three independent clocks
#
# White-FM clocks of similar grade (the classical hat is well conditioned
# only when no single clock dominates), from independent random streams:

N  = 8192
τ₀ = 1.0
m_values = unique(round.(Int, exp10.(range(0, log10(N ÷ 128); length = 12))))

clock(σ₁, seed) = noise_gen(PhaseData, N, τ₀; sigma1 = Dict(0 => σ₁), rng = Xoshiro(seed)).x

x1 = clock(1.00e-12, 51)
x2 = clock(1.05e-12, 152)
x3 = clock(1.10e-12, 253)

# ## Pairwise comparisons
#
# What a lab actually measures is the three differences. `ci = false`
# because the hat has no confidence model of its own; it works on the
# point estimates.

r12 = adev(PhaseData(x1 .- x2, τ₀), m_values; ci = false)
r23 = adev(PhaseData(x2 .- x3, τ₀), m_values; ci = false)
r31 = adev(PhaseData(x3 .- x1, τ₀), m_values; ci = false)

# ## Separation
#
# `nch` takes the pairs in cyclic order A−B, B−C, C−A (the sign of each
# difference is irrelevant) and returns one result per clock. A τ where the
# inversion goes negative is reported as `NaN`: an invalid estimate, and the
# usual sign that the independence assumption is failing.

clock1, clock2, clock3 = nch(r12, r23, r31)

# ## Check against the truth

truth1 = adev(PhaseData(x1, τ₀), m_values; ci = false)
truth2 = adev(PhaseData(x2, τ₀), m_values; ci = false)
truth3 = adev(PhaseData(x3, τ₀), m_values; ci = false)

println("τ (s)   recovered/true   clock 1   clock 2   clock 3")
for i in (1, length(m_values) ÷ 2, length(m_values))
    println(lpad(round(Int, clock1.tau[i]), 5), "                  ",
            join((rpad(round(c.dev[i] / t.dev[i]; digits = 3), 10)
                  for (c, t) in ((clock1, truth1), (clock2, truth2), (clock3, truth3)))))
end

# Ratios scatter around 1 and the scatter grows with τ, as the number of
# independent differences `N − 2m` shrinks. A longer record, or a fourth
# clock to overdetermine the system, tightens it.

using Plots

p = plot(truth1; label = "clock 1 true", ylabel = "σ_y(τ)", legend = :topright, lw = 1.5)
plot!(p, clock1; label = "clock 1 recovered", ls = :dash)
plot!(p, truth2; label = "clock 2 true", lw = 1.5)
plot!(p, clock2; label = "clock 2 recovered", ls = :dash)
plot!(p, truth3; label = "clock 3 true", lw = 1.5)
plot!(p, clock3; label = "clock 3 recovered", ls = :dash)

# ## When it breaks
#
# Correlated noise (a shared reference, common temperature or supply
# coupling) biases the separation and tends to show up as negative
# variances at specific τ. One clock much noisier than the others makes the
# two quiet recoveries ill-conditioned. Both are diagnosed the same way:
# add a fourth, independent clock and check that each clock's recovered σ
# is invariant to the choice of triplet. For N > 3 clocks, fill an
# upper-triangular matrix of pairwise results and call `nch(pairwise)`.

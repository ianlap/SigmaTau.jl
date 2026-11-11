# # Three clocks and the three-cornered hat
#
# A comparison of two clocks shows the noise of both, and no measurement
# shows one clock alone. Three clocks compared in pairs can be separated:
# by the three-cornered hat, from the deviations of the three comparisons
# ([Gray and Allan (1974)](@cite gray-1974-three-cornered-hat);
# [SP1065 §10.14](@cite riley-2008-sp1065)), or by the Groslambert
# covariance, from the three records themselves
# ([Vernotte et al. (2016)](@cite vernotte-2016-gcov)). This tutorial
# simulates three clocks, so that the answer is known, separates them both
# ways, and shows where each method stops working.

using SigmaTau
using Random

# ## Three clocks and their comparisons
#
# The clocks have white FM noise of ``σ_y = 1.0``, ``1.2`` and
# ``1.5 × 10^{-12}`` at 1 s, each from its own random stream. A measurement
# gives the three differences, sampled at the same instants:

τ₀ = 1.0
N  = 16_384
σ₁ = [1.0e-12, 1.2e-12, 1.5e-12]          # white FM: σ_y at 1 s

clock(σ, seed) =
    noise_gen(PhaseData, N, τ₀; sigma1 = Dict(0 => σ), rng = Xoshiro(seed)).x
a, b, c = clock(σ₁[1], 1), clock(σ₁[2], 2), clock(σ₁[3], 3)

ab = PhaseData(a .- b, τ₀)
bc = PhaseData(b .- c, τ₀)
ca = PhaseData(c .- a, τ₀)
nothing #hide

# ## The hat
#
# For independent clocks the Allan variance of a difference is the sum of
# the variances of its two clocks, ``σ_{ab}^2 = σ_a^2 + σ_b^2``. The three
# comparisons give three such sums, and each clock follows from them:
#
# ```math
# σ_a^2 = \tfrac{1}{2} \left( σ_{ab}^2 + σ_{ca}^2 - σ_{bc}^2 \right)
# ```
#
# ([SP1065 §10.14, Eqs. (68) and (69)](@cite riley-2008-sp1065)).
# [`ncornered_hat`](@ref) takes the results for A−B, B−C and C−A and returns
# those of the clocks A, B and C:

m = 2 .^ (0:10)
τ = m .* τ₀

hat = ncornered_hat(adev(ab, m), adev(bc, m), adev(ca, m))

using Plots

plot(hat; label = ["clock A" "clock B" "clock C"], legend = :bottomleft)
plot!(τ, σ₁' ./ sqrt.(τ); ls = :dash, color = :gray, label = ["true" "" ""])

# The three curves follow the noise put into the three clocks. The results
# of any one estimator serve, as long as its variance adds for independent
# clocks: [`mdev`](@ref) or [`hdev`](@ref) in place of `adev`.
#
# The intervals are wider than those of a direct measurement. A separated
# variance keeps a fraction of the degrees of freedom that the clock would
# have against a perfect reference: 0.4 for three clocks of equal stability,
# and less for the most stable of three unequal ones
# ([SP1065 §10.14, Eq. (70)](@cite riley-2008-sp1065)). For clock A that
# equation gives 0.22. A simulation also has the clock alone, and the EDF of
# the two results can be compared directly. The columns are ``τ`` and the
# fraction that the hat keeps:

alone = adev(PhaseData(a, τ₀), m)
[τ hat[1].edf ./ alone.edf][1:4, :]

# ## Where the hat fails
#
# The hat subtracts one measured variance from the sum of two others, and
# the difference can come out negative: for a clock much more stable than
# the other two, for clocks that are correlated, and where the data are too
# few, most often at long ``τ``. SP1065 advises the method for clocks of
# similar stability, and not in place of a low-noise reference
# ([§10.14](@cite riley-2008-sp1065)). Where the variance is negative the
# deviation is `NaN`. Here clock A gives way to a clock four times more
# stable than B. The columns are ``τ``, the hat's estimate of that clock
# and its true deviation:

q  = clock(0.3e-12, 4)
qb = PhaseData(q .- b, τ₀)
cq = PhaseData(c .- q, τ₀)

hat_q = ncornered_hat(adev(qb, m), adev(bc, m), adev(cq, m))[1]
[τ hat_q.dev 0.3e-12 ./ sqrt.(τ)]

# Up to 64 s the estimate lies within about a third of the true value. From
# 128 s on there is none.
#
# ## The noise of the instruments
#
# Each record also holds the noise of the instrument that measured it, and
# the hat counts it as noise of the two clocks of that record. Here every
# record gets white PM of 2 ps rms, the noise of its own counter:

rng = Xoshiro(10)
counter() = 2e-12 .* randn(rng, N)

ab_n = PhaseData(ab.x .+ counter(), τ₀)
bc_n = PhaseData(bc.x .+ counter(), τ₀)
ca_n = PhaseData(ca.x .+ counter(), τ₀)
nothing #hide

# [`groslambert`](@ref) works on the records in place of their deviations.
# It multiplies the second differences of the two records that share a
# clock and averages the product
# ([Vernotte et al. (2016), Eq. (7)](@cite vernotte-2016-gcov)). The noise
# of a counter is in one record only, so it averages out of the product,
# and the clock remains
# ([Groslambert et al. (1981), Eqs. (17) to (19)](@cite groslambert-1981-crosscorrelations);
# [Vernotte et al. (2016), Eqs. (13) and (15)](@cite vernotte-2016-gcov)).
# The columns are ``τ``, the hat and the covariance for clock A, and its
# true deviation:

hat_n = ncornered_hat(adev(ab_n, m), adev(bc_n, m), adev(ca_n, m))
gcov  = groslambert(ab_n, bc_n, ca_n, m)

[τ hat_n[1].dev gcov[1].dev σ₁[1] ./ sqrt.(τ)][1:6, :]

# The counters more than double the hat's estimate at 1 s, and the
# covariance stays within about a tenth of the clock:

plot([hat_n[1], gcov[1]]; label = ["hat" "covariance"], legend = :bottomleft)
plot!(τ, σ₁[1] ./ sqrt.(τ); ls = :dash, color = :gray, label = "clock A, true")

# Without instrument noise the three records sum to zero, and the two
# estimates are then the same number:

groslambert(ab, bc, ca, m)[1].dev ≈ hat[1].dev

# Both methods need records taken at the same time. SP1065 asks it of the
# hat, and comparisons made one after another can give a wrong or a
# negative variance
# ([Groslambert et al. (1981), Eqs. (4) to (6)](@cite groslambert-1981-crosscorrelations)).
#
# ## A clock that is not resolved
#
# Both functions report a ``χ^2`` interval on the equivalent degrees of
# freedom of the separated variance. It holds for a clock that dominates
# its records. For a clock more stable than the others it becomes too
# narrow once the clock is no longer resolved, and below 1 EDF none is
# given. For three clocks, `groslambert` with `interval = :klts` reports
# the interval of [Lantz et al. (2019)](@cite lantz-2019-klts) in its
# place. Where the records do not exclude a clock without noise, its lower
# limit is 0 and the upper limit is the result. For the stable clock of
# [Where the hat fails](@ref), at 95 %, the columns are ``τ``, the two
# limits and the true deviation:

klts = groslambert(qb, bc, cq, m; interval = :klts, confidence = 0.95)[1]
[τ ci_lower(klts) ci_upper(klts) 0.3e-12 ./ sqrt.(τ)]

# From 8 s on the records set an upper limit only. At 2 s the interval
# misses the true value: the docstring of [`groslambert`](@ref) gives the
# coverage of both intervals in simulation.
#
# ## More than three clocks
#
# With ``N`` clocks the pairs fill a matrix: of results for
# [`ncornered_hat`](@ref), of records for [`groslambert`](@ref). One entry
# for each pair is enough
# ([SP1065 §10.14, Eqs. (71) and (72)](@cite riley-2008-sp1065)). A fourth
# clock beside the first three:

x = [a, b, c, clock(1.3e-12, 5)]
pairs = Matrix{StabilityResult}(undef, 4, 4)
for i in 1:4, j in i+1:4
    pairs[i, j] = adev(PhaseData(x[i] .- x[j], τ₀), m)
end
hat4 = ncornered_hat(pairs)

[τ hat[1].edf hat4[1].edf][end-3:end, :]

# The columns are ``τ`` and the EDF of clock A from three clocks and from
# four. The fourth clock adds three comparisons, and with them degrees of
# freedom where the curve has fewest.
#
# ## Clocks read against one reference
#
# A counter with a channel for each clock gives the clocks against one
# reference, and not the pairs. [`ClockComparisons`](@ref) holds such
# readings, and both functions take it. They form the difference of every
# pair, in which the reference cancels, and count the reference as one
# clock more, which comes last. Clocks B and C read against A are the three
# clocks of the first section:

cc = ClockComparisons([b .- a c .- a], τ₀; names = ["B", "C"], reference = "A")

ncornered_hat(cc, m)[3].dev ≈ hat[1].dev

# The command [`sigmatau hat`](../reference/cli.md#hat) does the same from
# the columns of a file.
#
# ## Which to use
#
# | At hand | Use |
# |:---|:---|
# | clocks read against one reference | either function on a [`ClockComparisons`](@ref) |
# | the deviations of each pair | [`ncornered_hat`](@ref) |
# | the difference records | [`groslambert`](@ref), which keeps the noise of the instruments out |
# | three records, and a clock that may not be resolved | [`groslambert`](@ref) with `interval = :klts` |
# | two measurements of one clock | [`cross_adev`](@ref) ([Two-channel measurements](../guide/cross.md)) |

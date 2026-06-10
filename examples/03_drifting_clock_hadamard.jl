# # Characterizing a drifting clock with the Hadamard family
#
# A rubidium standard or GPSDO in its first weeks after power-on has
# linear frequency drift on top of its white-FM noise floor — the
# physics package is still settling toward thermal equilibrium. Drift
# at the 1e-10-per-day level is common in early life. On a record like
# this, ADEV reports the drift, not the clock: the deterministic ramp
# swamps the random noise at long τ. The Hadamard family exists for
# exactly this record.
#
# This tutorial:
#
# 1. Synthesises one day of 1 Hz data — calibrated WHFM noise plus an
#    explicit linear frequency drift.
# 2. Runs `adev` and shows the drift ramp masking the noise floor.
# 3. Runs `hdev` and shows the third difference removing the drift.
# 4. Introduces `mhdev` and `htdev` — the time-deviation analog,
#    original to this package.
# 5. Truncates to a three-hour record and compares `mhdev` against
#    `mhtotdev` degrees of freedom and confidence intervals at long τ.
# 6. Closes with a summary of which estimator to reach for.

using SigmaTau
using Random

# ## One day of a drifting clock
#
# [`noise_gen`](@ref) synthesises a calibrated power-law mixture; here a
# single component, white frequency modulation (α = 0) at
# σ_y(1 s) = 1e-11 — a typical rubidium / GPSDO class. `noise_gen` has
# no drift parameter, because deterministic drift is not noise: we add
# the term ``D \, t`` to the fractional-frequency series by hand.

Random.seed!(20260609)

τ₀ = 1.0                    # 1 Hz sampling
N  = 86_400                 # one day of data
σ1 = 1.0e-11                # WHFM level: σ_y(τ = 1 s)
D  = 2.0e-10 / 86_400       # drift rate, 2e-10 per day, in 1/s

fd_noise = noise_gen(FrequencyData, N, τ₀; sigma1 = Dict(0 => σ1))

t  = (0:N-1) .* τ₀
y  = fd_noise.y .+ D .* t   # drift injected explicitly
fd = FrequencyData(y, τ₀)

# A common octave τ-grid for every estimator in this tutorial, from
# 1 s out to 8192 s (≈ 2.3 h):

m_values = [2^k for k in 0:13]
taus     = m_values .* τ₀

# ## ADEV: the drift ramp
#
# For a pure linear frequency drift, the Allan deviation is
#
# ```math
# \sigma_y(\tau) \;=\; \frac{D\,\tau}{\sqrt{2}},
# ```
#
# where ``D`` is the drift rate in fractional frequency per second and
# ``\tau`` the averaging time. The second difference of phase does not
# annihilate a term linear in ``t`` in ``y(t)``; it passes it through
# as a residual proportional to ``\tau``, so the drift appears as a
# rising τ⁺¹ ramp. The ramp crosses the WHFM floor
# ``\sigma_1 \tau^{-1/2}`` where the two terms are equal:

τ_cross = (sqrt(2) * σ1 / D)^(2 / 3)
println("drift overtakes the WHFM floor near τ ≈ ", round(τ_cross; digits = 1), " s")

#-

r_adev = adev(fd, m_values)

using Plots

plot(r_adev;
     label  = "ADEV",
     xlabel = raw"Averaging time \(\tau\) (s)",
     ylabel = raw"\(\sigma_y(\tau)\)",
     title  = "ADEV of a drifting clock: the ramp masks the floor",
     legend = :bottomleft,
     lw     = 1.5)
plot!(taus, σ1 ./ sqrt.(taus); ls = :dash, label = raw"WHFM floor \(\sigma_1\tau^{-1/2}\)")
plot!(taus, D .* taus ./ sqrt(2); ls = :dot, label = raw"drift ramp \(D\tau/\sqrt{2}\)")

# Past the crossover the curve tracks the drift line, not the clock's
# noise. Stable32 users handle this by detrending first (the Drift
# function) and re-running; the Hadamard family makes that step
# unnecessary.

# ## HDEV: the third difference removes the drift
#
# `hdev` replaces ADEV's second difference of phase with a third
# difference. A linear term in ``y(t)`` is quadratic in phase, and a
# third difference annihilates any quadratic — so constant frequency
# drift contributes exactly zero, with no detrending pass over the
# data. See [the Allan-family theory page](../theory/allan_family.md)
# for the kernel definition.

r_hdev = hdev(fd, m_values)

plot(r_adev;
     label  = "ADEV",
     xlabel = raw"Averaging time \(\tau\) (s)",
     ylabel = raw"\(\sigma_y(\tau)\)",
     title  = "ADEV vs HDEV on the same drifting record",
     legend = :bottomleft,
     lw     = 1.5)
plot!(r_hdev; label = "HDEV", lw = 1.5)
plot!(taus, σ1 ./ sqrt.(taus); ls = :dash, label = raw"WHFM floor \(\sigma_1\tau^{-1/2}\)")

# HDEV follows the true τ⁻¹/² noise floor across the whole grid. The
# clock under the drift is recovered without touching the data.

# ## MHDEV and HTDEV
#
# `mhdev` adds an m-point phase average inside the third-difference
# window, exactly as MDEV does inside ADEV's. The averaging changes
# the WPM slope from −1 to −3/2 while leaving FPM at −1, so MHDEV
# distinguishes white-PM from flicker-PM where HDEV (like ADEV)
# cannot — and it keeps the third difference's drift immunity. On this
# WHFM-dominated record MHDEV simply parallels HDEV; its value shows
# on records where a drifting clock also has short-τ PM noise of
# ambiguous color.

r_mhdev = mhdev(fd, m_values)

# `htdev` is the time-deviation analog: it rescales MHDEV by
#
# ```math
# \sigma_{x,\mathrm{HT}}(\tau) \;=\; \frac{\tau}{\sqrt{10/3}}\,\sigma_{y,\mathrm{MH}}(\tau),
# ```
#
# where ``\sigma_{y,\mathrm{MH}}`` is MHDEV; HTDEV has units of seconds
# (it is a σ_x quantity). HTDEV is to MHDEV what TDEV is to MDEV — a
# time-error summary built on the drift-immune third-difference kernel —
# and it is original to this package: SP1065, IEEE 1139-2022, and
# NBS-TN-1337 do not define it (see
# [the HTDEV section of the theory page](../theory/allan_family.md)).

r_htdev = htdev(fd, m_values)

println("HTDEV is the τ/√(10/3) rescale of MHDEV — max deviation from the identity: ",
        maximum(abs.(r_htdev.dev .- r_mhdev.dev .* r_mhdev.tau ./ sqrt(10 / 3))))

# A few rows of the time-error budget for this clock, drift excluded
# by construction:

for i in (1, 5, 9, 13)
    println("τ = ", lpad(round(Int, r_htdev.tau[i]), 5), " s   σ_x = ",
            round(r_htdev.dev[i]; sigdigits = 3), " s")
end

# ## Short records: MHTOTDEV
#
# One day was generous. Suppose the record stops after three hours —
# a warm-up check before a deployment, say. At the longest reachable
# τ, `mhdev` has very few independent third differences left, so its
# χ² confidence interval is wide. `mhtotdev` extends each analysis
# subsegment by detrended reflection (the Greenhall total-variance
# methodology), keeping every subsegment in play at long τ and
# raising the equivalent degrees of freedom. Both its bias correction
# (`correct_bias = true` by default) and its EDF model are
# Monte-Carlo-calibrated against synthesized known-noise — MHTOTDEV is
# original to this package and has no published coefficient tables —
# and that calibration is valid only for `τ/τ₀ ≥ 16`, the same floor
# the Hadamard-total literature quotes. See
# [the total-family theory page](../theory/total_family.md) and
# [MHTOTDEV bias and EDF](../theory/mhtotdev_bias_edf.md) for the
# methodology.
#
# Unlike its parent `mhdev`, the total kernel is not drift-immune on
# its own — the boundary extension re-admits residual drift — so
# `mhtotdev` removes the record's least-squares frequency drift itself
# before analyzing (`remove_drift = true`, the default; per-window
# drift removal was evaluated in both the phase and frequency domains
# and measurably damages the statistic, so the global removal is the
# design). Since drift is deterministic, removing it discards no
# statistical information; it is the same practice SP1065 recommends
# before any total estimator. One piece of hygiene remains for the
# *comparison* below: `mhdev`'s noise identification runs on the
# record as given and would classify the longest-τ points as redder
# than the true noise on a record drifting this strongly. Detrending
# the record once lets both estimators' χ² intervals classify the
# underlying noise:

N_short  = 10_800   # three hours at 1 Hz
fd_short = detrend(FrequencyData(y[1:N_short], τ₀))   # noise-ID hygiene for the mhdev comparison
m_short  = [16, 32, 64, 128, 256, 512, 1024, 2048]   # all ≥ the τ/τ₀ ≥ 16 floor

r_mh  = mhdev(fd_short, m_short)
r_mht = mhtotdev(fd_short, m_short)

# The octave-by-octave comparison, with the default 68.3% (1σ)
# confidence level:

ci_width(r, i) = (r.ci[i].hi - r.ci[i].lo) / r.dev[i]

println("τ (s)    T/τ    edf mhdev  edf mhtotdev  rel CI mhdev  rel CI mhtotdev")
for i in eachindex(m_short)
    Ttau = (N_short - 1) * τ₀ / r_mh.tau[i]
    println(lpad(round(Int, r_mh.tau[i]), 5), "  ",
            lpad(round(Ttau; digits = 1), 6), "  ",
            lpad(round(r_mh.edf[i];  digits = 1), 9), "  ",
            lpad(round(r_mht.edf[i]; digits = 1), 12), "  ",
            lpad(round(ci_width(r_mh,  i); digits = 2), 12), "  ",
            lpad(round(ci_width(r_mht, i); digits = 2), 15))
end

# Read the table from the bottom up. From τ = 16 s through τ = 1024 s
# the extension raises the equivalent degrees of freedom by a steady
# ≈ 30 % and tightens the relative confidence interval accordingly.
# The gain collapses only at the kernel's reach limit (`mhtotdev`
# needs at least one subsegment, `N − 4m + 1 ≥ 1`, so τ cannot exceed
# T/4, where `T` is the record duration): at the last octave both
# estimators are down to an EDF of 2–3, and no extension trick
# manufactures much independent information that close to the record
# length. The relative gain is uniform, but it buys the most where
# the unextended interval is widest — the last usable decade, which
# is exactly where a warm-up check needs its answer. Plotted with CI
# bands:

plot(r_mh;
     label   = "MHDEV (3 h record)",
     ci_band = true,
     xlabel  = raw"Averaging time \(\tau\) (s)",
     ylabel  = raw"\(\sigma_y(\tau)\)",
     title   = "MHDEV vs MHTOTDEV confidence on a short record",
     legend  = :bottomleft,
     lw      = 1.5)
plot!(r_mht; label = "MHTOTDEV (3 h record)", ci_band = true, lw = 1.5)

# ## Which estimator, when
#
# All four Hadamard-family estimators reject linear frequency drift.
# Within the family:

println("""
Symptom                                              Reach for
---------------------------------------------------  ---------
drift contaminates ADEV at long τ                    hdev
drift + ambiguous PM noise (WPM vs FPM) at short τ   mhdev
time-error budget (σ_x) for a drifting clock         htdev
drifting clock, record too short for long-τ CI       mhtotdev
""")

# | Symptom | Reach for |
# |---|---|
# | Drift contaminates ADEV at long τ | `hdev` |
# | Drift plus ambiguous PM noise (WPM vs FPM) at short τ | `mhdev` |
# | Time-error budget (σ_x) for a drifting clock | `htdev` |
# | Drifting clock on a record too short for long-τ confidence | `mhtotdev` |
#
# For records without drift the Allan-family counterparts (`adev`,
# `mdev`, `tdev`, `mtotdev`) are slightly more efficient — tighter CI
# per τ on the same record — so the Hadamard family is the tool for
# the drifting clock, not the default for every clock.

# # Reading your data
#
# Measurement records arrive as plain text: Stable32 `.DAT` exports, counter
# logs, CSV dumps from a data-acquisition system. This tutorial is a cookbook
# for getting each of those into a [`PhaseData`](@ref) or
# [`FrequencyData`](@ref) record. Once the data is in one of those two
# containers, everything downstream — `adev`, `mdev`, the total family — is
# identical.
#
# Covered here:
#
# 1. Phase versus frequency records, and which one your instrument produces.
# 2. Two-column files (time, phase) with `read_phase` defaults.
# 3. Single-column files: `time_col = 0` and an explicit `tau0`.
# 4. A counter CSV with MJD time tags, parsed by hand.
# 5. Gaps and outliers: `NaN`, `fillgaps`, and the MAD test.
# 6. Removing offsets and drift with `detrend`.
# 7. Writing results out for plotting elsewhere.
#
# Everything below writes only to a temporary directory and runs as-is with
# `julia --project=examples examples/04_reading_your_data.jl`.

using SigmaTau
using Random
using Statistics

Random.seed!(20260609)

# ## Phase records and frequency records
#
# A phase record ``x[k]`` is the time difference between the device under test
# and the reference, in seconds, sampled every ``\tau_0`` seconds. A
# fractional-frequency record ``y[k]`` is the dimensionless frequency offset
# averaged over each interval, ``y[k] = (x[k+1] - x[k]) / \tau_0``. Which one
# you have depends on the instrument: time-interval counters and dual-mixer
# time-difference systems measure phase; most frequency counters report
# averaged fractional frequency.
#
# Every deviation function in `SigmaTau` accepts either container and converts
# internally, so you rarely convert by hand. When you do need the conversion
# explicitly, it is one line each way:
#
# ```julia
# y = diff(x) ./ τ₀       # phase → fractional frequency
# x = cumsum(y) .* τ₀     # fractional frequency → phase (absolute offset lost)
# ```
#
# The absolute offset is lost in the round trip because differencing destroys
# the constant term; none of the deviations depend on it.

# ## A two-column file: `read_phase` with defaults
#
# The most common layout is two whitespace- or comma-separated columns: a time
# tag and a phase value. We synthesise one — a white-FM phase record sampled
# at 1 s — and write it to a temporary directory.

dir = mktempdir()

τ₀ = 1.0
N  = 2048
x  = cumsum(randn(N)) .* 1e-10        # white FM: phase is a random walk

twocol = joinpath(dir, "twocol.DAT")
open(twocol, "w") do io
    for k in 1:N
        println(io, (k - 1) * τ₀, "  ", x[k])
    end
end

# `read_phase` with no keyword arguments expects exactly this layout: column 1
# is time (`time_col = 1`), column 2 is phase in seconds (`value_col = 2`).
# The sample interval is inferred from the time column as the *median* of
# `diff(t)` — the median rather than the mean, so a few missing rows do not
# bias the inferred cadence.

pd = read_phase(twocol)

length(pd.x), pd.tau0

# The delimiter is auto-detected from the file extension: `,` for `.csv`,
# tab for `.tsv`, whitespace for everything else. Pass `delim` explicitly to
# override.

# ## A single-column file: the classic stumbling block
#
# Stable32 itself usually writes a *single* column of phase values, preceded
# by a short header. The repository carries one such file as a validation
# fixture — an 8192-point simulated phase record produced by Stable32's noise
# generator. Its header looks like this:

datpath = joinpath(pkgdir(SigmaTau), "test", "fixtures", "validation",
                   "stable32gen.DAT")
foreach(println, readlines(datpath)[1:10])

# With no time column there is nothing to infer ``\tau_0`` from, and
# `read_phase` throws an `ArgumentError` rather than guessing. You must say
# three things: skip the header (`header = 10`), there is no time column
# (`time_col = 0`), and the phase values are therefore in column 1
# (`value_col = 1`) — and supply `tau0` yourself. Here the header records it
# on the `Tau:` line, so we parse it from there:

tau_line = first(filter(startswith("Tau:"), readlines(datpath)))
τ₀_dat   = parse(Float64, split(tau_line)[2])

pd_dat = read_phase(datpath; header = 10, time_col = 0, value_col = 1,
                    tau0 = τ₀_dat)

length(pd_dat.x), pd_dat.tau0

# Forgetting `value_col = 1` is the usual failure mode: the default
# `value_col = 2` points past the only column and the reader raises an error
# instead of silently reading the wrong thing.

# ## A counter CSV with MJD time tags: no special reader needed
#
# Counter logs often carry Modified Julian Date time tags and phase in
# engineering units (nanoseconds, picoseconds). We synthesise a 600-row log
# sampled every 10 s, with phase in nanoseconds:

csvpath = joinpath(dir, "counter_log.csv")
mjd0    = 60800.0
τ_csv   = 10.0
M       = 600
xns     = cumsum(randn(M)) .* 0.05    # phase in ns

open(csvpath, "w") do io
    println(io, "mjd,phase_ns")
    for k in 1:M
        println(io, mjd0 + (k - 1) * τ_csv / 86400, ',', xns[k])
    end
end

# You do not need a special reader for this — plain Julia is enough. Parse
# the two columns, convert MJD days to seconds, nanoseconds to seconds, and
# construct the `PhaseData` directly:

rows = split.(readlines(csvpath)[2:end], ',')
mjd  = parse.(Float64, first.(rows))
xs   = parse.(Float64, last.(rows)) .* 1e-9     # ns → s
t    = (mjd .- mjd[1]) .* 86400.0               # days → s

pd_csv = PhaseData(xs, round(t[2] - t[1]; digits = 3))

length(pd_csv.x), pd_csv.tau0

# The `round` absorbs floating-point jitter in the time tags: a `Float64` MJD
# near 60800 resolves only ~0.6 µs (`eps(60800.0) * 86400` seconds), so the
# sample spacings recovered from raw MJD tags are not exactly 10.0 s.
#
# `read_phase` can read the same file — but mind the units. If you let it
# infer `tau0` from the MJD column, you get the cadence in *days*:

read_phase(csvpath; header = 1, scaling = 1e-9).tau0   # ≈ 1.157e-4 — days!

# That call succeeds silently and every τ in your results would be wrong by a
# factor of 86400. Pass `tau0` in seconds explicitly (and `scaling = 1e-9`
# for the ns → s conversion):

pd_csv2 = read_phase(csvpath; header = 1, scaling = 1e-9, tau0 = 10.0)
pd_csv2.x ≈ pd_csv.x

# ## Gaps and outliers
#
# The deviation estimators expect an equispaced, gap-free record. The kernels
# do not skip `NaN` samples — a single `NaN` propagates into every σ_y(τ)
# estimate whose sums touch it:

xg = copy(x)                  # the 2048-point record from above
xg[1001:1012] .= NaN          # a 12-sample dropout

adev(PhaseData(xg, τ₀), [8]; ci = false).dev

# [`fillgaps`](@ref) imputes the missing samples using Howe and
# Schlossberger's reflect-and-filter algorithm (PTTI 2009): points adjacent
# to each gap are mirrored, inverted, low-pass filtered, and joined with an
# endpoint-matched linear ramp. The point of that machinery is physical: the
# fill preserves the local noise character, so the σ_y(τ) curve is not
# distorted by the imputed segment the way a straight-line or zero fill would
# distort it.

pd_filled = fillgaps(PhaseData(xg, τ₀))
adev(pd_filled, [8]; ci = false).dev

# For raw files where dropouts appear as *missing rows* rather than `NaN`
# placeholders, pass `fillgaps = true` to `read_phase` / `read_frequency`:
# the record is snapped to a uniform grid (spacing = the minimum time step),
# absent samples become `NaN`, and the same fill runs. This requires a time
# column.
#
# `SigmaTau` has no automatic outlier detector. Standard practice
# [riley-2008-sp1065](@cite) is the median-absolute-deviation test on the
# frequency record: flag samples far from the median, set them to `NaN`, and
# treat them as gaps. Here we inject a counter glitch and remove it:

y      = diff(pd_filled.x) ./ τ₀
y[700] = 3e-9                          # a 30σ glitch

med = median(y)
mad = median(abs.(y .- med))
bad = abs.(y .- med) .> 5 * mad / 0.6745
count(bad)

#-

y[bad] .= NaN
fd_clean = fillgaps(FrequencyData(y, τ₀))
adev(fd_clean, [8]; ci = false).dev

# where `mad / 0.6745` rescales the median absolute deviation to be an
# estimate of the standard deviation for normally distributed data, and 5 is
# the deglitching threshold in those units (Stable32's default for its
# outlier check).

# ## Removing offsets and drift with `detrend`
#
# A constant frequency offset appears as a linear trend in phase; a linear
# frequency drift appears as a linear trend in ``y[k]`` (quadratic in
# ``x[k]``). Deterministic drift inflates the long-τ end of the σ_y(τ) curve,
# masking the random-walk noise floor underneath, so the convention is to
# remove the drift first and characterise the residual noise.
#
# [`detrend`](@ref) does this on either container and returns a new record,
# leaving the original untouched. `method` is one of `:linear` (least-squares
# straight line, the default), `:endpoint` (line through the first and last
# samples), `:mean` (subtract the mean), or `:none`.

ydrift   = randn(4096) .* 1e-12 .+ 5e-16 .* (0:4095)   # white FM + drift
fd_drift = FrequencyData(ydrift, 1.0)
fd_flat  = detrend(fd_drift)                            # method = :linear

adev(fd_drift, [512]; ci = false).dev[1],
adev(fd_flat,  [512]; ci = false).dev[1]

# The drift contribution dominated σ_y(512 s) before detrending and is gone
# after. The same operation is available at load time through the `detrend`
# keyword of `read_phase` / `read_frequency`, e.g.
# `read_frequency(path; detrend = :linear)`.
#
# One caution: the detrend modes fit at most a straight line. Applied to
# *phase* data, `:linear` removes a constant frequency offset — not a
# frequency drift. To remove a linear frequency drift, detrend the frequency
# representation (`detrend(FrequencyData(diff(pd.x) ./ pd.tau0, pd.tau0))`).

# ## Writing results out
#
# A [`StabilityResult`](@ref) is a plain struct — the vectors `tau`, `dev`,
# `noise_type`, `ci_lower`, `ci_upper`, `edf`, plus a `deviation_type` tag —
# so exporting it for a colleague's plotting tool is a short loop. We run `adev` on the Stable32 fixture from earlier
# and write a tab-separated table:

result  = adev(pd_dat)        # default octave τ-grid, CI included
outpath = joinpath(dir, "adev_results.tsv")

open(outpath, "w") do io
    println(io, "tau_s\tadev\tci_lower\tci_upper\tnoise")
    for k in eachindex(result.tau)
        println(io, result.tau[k], '\t', result.dev[k], '\t',
                result.ci_lower[k], '\t', result.ci_upper[k], '\t',
                result.noise_type[k])
    end
end

foreach(println, readlines(outpath)[1:5])

# Any spreadsheet, gnuplot script, or Stable32 user can take it from there.
#
# For round-tripping *within* `SigmaTau`, prefer the native format:
# [`save_result`](@ref) writes a self-describing tab-delimited file that
# [`load_result`](@ref) reads back into an identical `StabilityResult`
# (and [`save_suite`](@ref) / [`load_suite`](@ref) do the same for a whole
# [`stability`](@ref) session):

nativepath = joinpath(dir, "adev_native.tsv")
save_result(nativepath, result)
load_result(nativepath).dev ≈ result.dev

# ## Where to go next
#
# - [Phase data and frequency data](01_phase_data.md) — the two containers in
#   detail, and the phase ↔ frequency identity.
# - [Computing the Allan deviation](02_compute_adev.md) — every field of the
#   `StabilityResult` you just wrote to disk.

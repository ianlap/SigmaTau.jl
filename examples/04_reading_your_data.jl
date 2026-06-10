# # Reading and saving your data
#
# Measurement records arrive as plain text: Stable32 `.DAT` exports, counter
# logs, CSV dumps from a data-acquisition system. SigmaTau's IO surface for
# all of them is small: [`read_phase`](@ref) and [`read_frequency`](@ref)
# bring raw records into a [`PhaseData`](@ref) or [`FrequencyData`](@ref)
# container, [`fillgaps`](@ref), [`remove_outliers`](@ref), and
# [`detrend`](@ref) clean them up, and [`save`](@ref) writes data records and
# finished analyses to plain text files that `read_phase` /
# [`load_result`](@ref) / [`load_suite`](@ref) read back. This tutorial walks
# that whole path.
#
# Covered here:
#
# 1. Phase versus frequency records, and which reader to call.
# 2. Two-column files with the `read_phase` defaults.
# 3. Single-column Stable32 `.DAT` files: `header`, `time_col = 0`,
#    `value_col = 1`, and an explicit `tau0`.
# 4. A counter CSV in engineering units: `scaling`, and the `tau0` units trap.
# 5. The escape hatch: constructing a `PhaseData` by hand when no reader
#    keyword fits.
# 6. Data preparation: `fillgaps` for dropouts, `find_outliers` /
#    `remove_outliers` for glitches, `detrend` for offsets and drift.
# 7. Saving and sharing: `save` for data records, `save_result` /
#    `save_suite` for one deviation or a whole session, and the loaders.
#
# Everything below writes only to a temporary directory and runs as-is with
# `julia --project=examples examples/04_reading_your_data.jl`.

using SigmaTau
using Random

Random.seed!(20260609)

# ## Phase records and frequency records
#
# A phase record ``x[k]`` is the time difference between the device under test
# and the reference, in seconds, sampled every ``\tau_0`` seconds. A
# fractional-frequency record ``y[k]`` is the dimensionless frequency offset
# averaged over each interval. Which one you have depends on the instrument:
# time-interval counters and dual-mixer time-difference systems measure phase;
# most frequency counters report averaged fractional frequency.
#
# Call the reader that matches the file: `read_phase` returns a `PhaseData`,
# `read_frequency` a `FrequencyData`, and the two take identical keyword
# arguments. Every deviation function accepts either container and converts
# internally — see [the phase-data tutorial](01_phase_data.md) for the
# conversion identity.

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

# ## A counter CSV in engineering units
#
# Counter logs often carry Modified Julian Date (MJD) time tags and phase in
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

# `read_phase` handles this directly: `header = 1` skips the column-name row,
# and `scaling = 1e-9` multiplies the values into seconds. The one thing the
# reader cannot do for you is convert *time* units. The MJD column is in
# days, so the inferred cadence comes out in days too:

read_phase(csvpath; header = 1, scaling = 1e-9).tau0   # ≈ 1.157e-4 — days!

# That call succeeds silently, and every τ in your results would be wrong by
# a factor of 86400. Whenever the time column is not in seconds, pass `tau0`
# in seconds explicitly — it overrides the inference:

pd_csv = read_phase(csvpath; header = 1, scaling = 1e-9, tau0 = 10.0)

length(pd_csv.x), pd_csv.tau0

# ## The escape hatch: constructing the container by hand
#
# `read_phase` requires the whole table to parse as numbers. A log with a
# text column — a lock-status flag, say — defeats it, and no keyword
# combination will help:

logpath = joinpath(dir, "counter_status.log")
open(logpath, "w") do io
    for k in 1:M
        println(io, mjd0 + (k - 1) * τ_csv / 86400, "  OK  ", xns[k])
    end
end

# For a genuinely foreign format like this, parse it with plain Julia and
# construct the `PhaseData` directly — the containers are ordinary structs,
# and everything downstream is identical no matter how the record got in:

rows   = split.(readlines(logpath))
xs     = [parse(Float64, r[3]) * 1e-9 for r in rows if r[2] == "OK"]
pd_log = PhaseData(xs, τ_csv)

pd_log.x ≈ pd_csv.x

# Three lines is typical. Reach for this only when the file truly fits no
# `read_phase` keyword combination; for anything column-shaped the reader is
# shorter and validates as it goes.

# ## Gaps and outliers: `fillgaps`
#
# The deviation estimators expect an equispaced, gap-free record, and the
# kernels do not skip `NaN` samples — a single `NaN` propagates into every
# σ_y(τ) estimate whose sums touch it:

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
# The same machinery handles outliers. [`find_outliers`](@ref) implements
# standard practice [riley-2008-sp1065](@cite) — the median-absolute-deviation
# test on the frequency record, flagging samples more than `nsigma` scaled
# MADs from the median:

y      = diff(pd_filled.x) ./ τ₀
y[700] = 3e-9                          # inject a counter glitch
fd_glitchy = FrequencyData(y, τ₀)

find_outliers(fd_glitchy)

# The scale factor 0.6745 inside the test rescales the median absolute
# deviation into an estimate of the standard deviation for normally
# distributed data, and the default `nsigma = 5` is the deglitching threshold
# in those units (Stable32's default for its outlier check). On a `PhaseData`
# the same call tests the *first differences* — a phase glitch appears as a
# step in frequency — and flags the phase sample(s) where each step lands.
#
# [`remove_outliers`](@ref) completes the idiom: flag, set `NaN`, and impute
# with the same Howe fill as `fillgaps`, so the replacement preserves the
# noise character:

fd_clean = remove_outliers(fd_glitchy)

adev(fd_clean, [8]; ci = false).dev[1]

# ## Removing offsets and drift: `detrend`
#
# A constant frequency offset appears as a linear trend in phase; a linear
# frequency drift appears as a linear trend in ``y[k]``. Deterministic drift
# inflates the long-τ end of the σ_y(τ) curve, masking the random-walk noise
# floor underneath, so the convention is to remove the drift first and
# characterise the residual noise.
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
# after. The same operation is available at read time through the `detrend`
# keyword — here via `read_frequency`, which takes the same keywords as
# `read_phase` and returns a `FrequencyData`:

driftpath = joinpath(dir, "drift.txt")
write(driftpath, join(ydrift, '\n'))

fd_flat2 = read_frequency(driftpath; time_col = 0, value_col = 1, tau0 = 1.0,
                          detrend = :linear)
fd_flat2.y ≈ fd_flat.y

# One caution: the detrend modes fit at most a straight line. Applied to
# *phase* data, `:linear` removes a constant frequency offset — not a
# frequency drift. To remove a linear frequency drift, detrend the frequency
# representation (`detrend(FrequencyData(diff(pd.x) ./ pd.tau0, pd.tau0))`).

# ## Saving a data record: `save`
#
# The trip out of SigmaTau matters too: generate a reference record with
# [`noise_gen`](@ref) — or clean up a measured one — and hand the samples to
# software outside SigmaTau. [`save`](@ref) writes any `PhaseData` /
# `FrequencyData` as a plain two-column text file (sample time in seconds,
# value) behind a short comment header recording the provenance and the
# sample interval:

wfm      = noise_gen(PhaseData, 1024, τ₀; sigma1 = Dict(0 => 1e-12))
datafile = joinpath(dir, "white_fm.txt")
save(datafile, wfm)

foreach(println, readlines(datafile)[1:5])

# The `source` line reads `user` for records constructed in the session and
# the originating path for records that came from `read_phase` /
# `read_frequency` — every record carries that provenance in its `source`
# field, and `detrend`, `fillgaps`, and `remove_outliers` preserve it. On the
# way back in, `read_phase` skips the comment header automatically and picks
# `tau0` up from it, so the file round-trips with no keyword arguments at
# all:

wfm2 = read_phase(datafile)
wfm2.x ≈ wfm.x, wfm2.tau0 == wfm.tau0, wfm2.source == datafile

# ## Saving one result: `save_result` and `load_result`
#
# A finished analysis should leave the session as a file, not a screenshot.
# [`save_result`](@ref) writes a [`StabilityResult`](@ref) to a
# self-describing tab-delimited text file; [`load_result`](@ref) reads it
# back into an identical record. We run `adev` on the Stable32 fixture from
# earlier and save it:

result     = adev(pd_dat)        # default octave τ-grid, CI included
resultpath = joinpath(dir, "adev_results.tsv")
save_result(resultpath, result)

foreach(println, readlines(resultpath)[1:8])

# The file is a plain table: two comment lines record the deviation type and
# whether confidence intervals were computed, then a header row and one row
# per τ with `tau`, `dev`, `noise_type`, `ci_lower`, `ci_upper`, `edf`, and
# `neff`. Any spreadsheet, gnuplot script, or Stable32 user can take it
# as-is; there is no binary format to escape from.
#
# `load_result` reconstructs the `StabilityResult`, confidence intervals and
# all:

result2 = load_result(resultpath)
result2.dev ≈ result.dev, result2.edf ≈ result.edf

# ## Saving a whole session: `save_suite` and `load_suite`
#
# [`stability`](@ref) computes several deviations from one record in a single
# call and returns a [`StabilitySuite`](@ref). [`save_suite`](@ref) writes
# the whole session — every result plus the metadata needed to reproduce it —
# to one file, and [`load_suite`](@ref) reads it back. `source_file` is
# optional provenance, recorded in the header:

suite     = stability(pd_dat)    # adev, mdev, hdev, mhdev by default
suitepath = joinpath(dir, "session.tsv")
save_suite(suitepath, suite; source_file = "stable32gen.DAT")

foreach(println, readlines(suitepath)[1:14])

# The metadata header records the package version, a timestamp, the input
# kind, τ₀, the number of samples, the confidence level, and the deviation
# set; each deviation then follows as its own block in the same six-column
# table that `save_result` writes. Loading it back returns a suite you can
# index by deviation symbol:

suite2 = load_suite(suitepath)
suite2[:hdev].dev ≈ suite[:hdev].dev

# ## Where to go next
#
# - [Phase data and frequency data](01_phase_data.md) — the two containers in
#   detail, and the phase ↔ frequency identity.
# - [Computing the Allan deviation](02_compute_adev.md) — every field of the
#   `StabilityResult` you just wrote to disk.

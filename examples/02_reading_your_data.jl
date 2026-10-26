# # Reading and saving data
#
# Records arrive as plain text: Stable32 `.DAT` exports, counter logs, CSV
# dumps. `read_phase` and `read_frequency` bring them in, `fillgaps`,
# `remove_outliers` and `detrend` clean them, and `save` writes data and
# results back out. This page is the whole IO surface.

using SigmaTau
using Random
Random.seed!(20260609)

dir = mktempdir()   # everything below writes here

# ## Two-column files
#
# The default layout is a time column and a value column, whitespace, comma
# (`.csv`) or tab (`.tsv`) separated. The sample interval is inferred from the
# time column, so no keywords are needed:

τ₀ = 1.0
x  = cumsum(randn(2048)) .* 1e-10      # white FM phase record
twocol = joinpath(dir, "twocol.DAT")
open(twocol, "w") do io
    for k in eachindex(x)
        println(io, (k - 1) * τ₀, "  ", x[k])
    end
end

pd = read_phase(twocol)

# If the values are in engineering units, `scaling` converts them
# (`scaling = 1e-9` for nanoseconds). If the time column is not in seconds,
# pass `tau0` explicitly; it overrides the inference.

# ## Stable32 `.DAT` files
#
# Stable32 writes a 10-line text header and a single column of values. Say
# so once: skip the header, no time column, values in column 1, and the
# sample interval (the header's `Tau:` line). The repository's validation
# fixture is such a file:

datpath = joinpath(pkgdir(SigmaTau), "test", "fixtures", "validation",
                   "stable32gen.DAT")
foreach(println, readlines(datpath)[1:10])

#-

pd_dat = read_phase(datpath; header = 10, time_col = 0, value_col = 1, tau0 = 1.0)

# `read_frequency` takes the same keywords and returns a `FrequencyData`.
# Every deviation function accepts either container.

# ## Cleaning
#
# The kernels need an equispaced, gap-free record; one `NaN` poisons every
# σ(τ) whose sums touch it. `fillgaps` imputes missing samples with the
# Howe–Schlossberger reflect-and-filter method, which preserves the local
# noise character instead of flattening it:

xg = copy(x)
xg[1001:1012] .= NaN
pd_gap = PhaseData(xg, τ₀)

adev(pd_gap, [8]; ci = false).dev[1], adev(fillgaps(pd_gap), [8]; ci = false).dev[1]

# `remove_outliers` applies Stable32's median-absolute-deviation test
# (default 5 scaled MADs on the frequency record, or on the first differences
# of a phase record) and imputes the flagged samples the same way:

y = diff(x) ./ τ₀
y[700] = 3e-9                         # a counter glitch
fd = FrequencyData(y, τ₀)

find_outliers(fd), adev(remove_outliers(fd), [8]; ci = false).dev[1]

# `detrend` removes a mean, an endpoint line or a least-squares line
# (`method = :mean`, `:endpoint`, `:linear`), or a quadratic from phase data.
# On a frequency record `:linear` removes a frequency drift:

ydrift = randn(4096) .* 1e-12 .+ 5e-16 .* (0:4095)
fd_drift = FrequencyData(ydrift, 1.0)

adev(fd_drift, [512]; ci = false).dev[1], adev(detrend(fd_drift), [512]; ci = false).dev[1]

# All three are also keywords on the readers: `read_phase(path; fillgaps = true,
# detrend = :linear)` cleans on the way in. `fillgaps = true` needs a time
# column, since it is how missing rows are detected.

# ## Saving
#
# `save` writes a `PhaseData` / `FrequencyData` as a two-column text file with
# a short comment header; `read_phase` reads it back with no keywords, since
# the header carries τ₀:

wfm = noise_gen(PhaseData, 1024, τ₀; sigma1 = Dict(0 => 1e-12))
datafile = save(joinpath(dir, "white_fm.txt"), wfm)
foreach(println, readlines(datafile)[1:5])

#-

read_phase(datafile).x ≈ wfm.x

# The same `save` writes a `StabilityResult` or a `StabilitySuite` as a
# self-describing tab-separated table that any spreadsheet can open;
# `load_result` and `load_suite` read them back:

result = adev(pd_dat)
resultpath = save(joinpath(dir, "adev.tsv"), result)
foreach(println, readlines(resultpath)[1:6])

#-

load_result(resultpath).dev ≈ result.dev

#-

suite = stability(pd_dat)              # adev, mdev, hdev, mhdev
suitepath = save(joinpath(dir, "session.tsv"), suite; source_file = "stable32gen.DAT")
load_suite(suitepath)[:hdev].dev ≈ suite[:hdev].dev

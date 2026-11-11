# # Reading and saving data
#
# A record usually arrives as a text file: columns of numbers, a table that
# a logger wrote, or the export of an instrument or of another program.
# This tutorial writes three small files, reads each one back as a record,
# looks at what the reader made of it, and saves a record and its results.
# [Data and files](../reference/data.md) is the reference for every reader
# and format.

using SigmaTau
using Dates
using Random

dir = mktempdir()            # the files of this page go here
clock = noise_gen(PhaseData, 2048, 1.0; sigma1 = Dict(0 => 1e-11), rng = Xoshiro(1))
nothing #hide

# ## Columns of numbers
#
# The plainest file has a time column and a column of samples, separated by
# blanks, by commas, or by tabs in a `.tsv` file. This one holds the time
# and the phase of the simulated clock, both in seconds:

file = joinpath(dir, "clock.txt")
open(file, "w") do io
    for (k, x) in enumerate(clock.x)
        println(io, k - 1, "  ", x)
    end
end
foreach(println, readlines(file)[1:3])

# [`read_phase`](@ref) needs no keywords for it. The sample interval is the
# step of the time column:

p = read_phase(file)
p.tau0, length(p)

# The file is read line by line, so a long log costs the memory of its
# samples and no more. Keywords cover the other layouts: `value_col` and
# `time_col` choose the columns, `header` skips lines of text, `scaling`
# converts the unit of the samples, and a file with a single column takes
# its sample interval from `tau0`. [`read_frequency`](@ref) reads
# fractional frequency in the same way.
#
# ```julia
# read_phase("log.txt"; value_col = 4, header = 1, scaling = 1e-9)  # nanoseconds
# read_frequency("samples.txt"; tau0 = 0.1)                         # one column
# ```
#
# ## Missing readings
#
# A logger that drops readings leaves a file with rows missing. The readers
# place every sample by its time, so a missing row becomes a missing sample,
# `NaN`, and the samples on either side of it keep their distance. This file
# lacks twelve rows:

gapfile = joinpath(dir, "gaps.txt")
open(gapfile, "w") do io
    for (k, x) in enumerate(clock.x)
        k in 1001:1012 || println(io, k - 1, "  ", x)
    end
end

g = read_phase(gapfile)
length(g), count(isnan, g.x)

# The Allan and Hadamard deviations, their modified and time forms,
# [`mtie`](@ref) and [`tierms`](@ref) leave the missing samples out, and a
# result counts the terms that were used:

adev(p, [8]).neff, adev(g, [8]).neff

# The total, Thêo, parabolic and dynamic estimators and the spectra need a
# complete record. [Timestamps and gaps](@ref) lists what every function
# does with missing samples, and
# [Gaps, outliers and drift](../guide/preparing.md) weighs filling them
# against analyzing around them.
#
# ## A table from a logger
#
# [`read_data`](@ref) reads files that carry more than numbers: a line of
# column names, units, timestamps, or the header of an instrument. This
# table holds the offset of a 1 PPS output in nanoseconds, beside a
# temperature:

pps = noise_gen(PhaseData, 600, 1.0; sigma1 = Dict(2 => 1e-9), rng = Xoshiro(2))
table = joinpath(dir, "gpsdo.csv")
open(table, "w") do io
    println(io, "timestamp,PPS offset (ns),temperature (C)")
    for (k, x) in enumerate(pps.x)
        t = DateTime(2026, 3, 1, 12) + Second(k - 1)
        println(io, t, "Z,", round(1e9 * x; digits = 3), ",31.5")
    end
end
foreach(println, readlines(table)[1:3])

# The reader finds the format, takes the time column and the clock column
# by their names, and converts the samples to seconds by the unit in the
# name. The timestamps give the sample interval and the UTC time of the
# first sample:

q = read_data(table)
detect_format(table), q.tau0, start_time(q), q.x[1:2]

# The same call reads the `.phd` and `.frd` files of Stable32, the `.tim`
# files of TimeLab, the readings of a counter, the BIPM files of
# UTC − UTC(k) and RINEX clock files. The docstring of [`read_data`](@ref)
# states how each format is read, how a column is chosen when a table has
# several, and where a file can be read wrongly.
#
# ## A first look
#
# [`summarize`](@ref) reports what a record holds before any analysis: its
# length and span, the number of missing samples, the indices of its
# outliers and the slope of its fractional frequency, per day. `pairs`
# lists the fields one to a line:

pairs(summarize(g))

# [Gaps, outliers and drift](../guide/preparing.md) goes through the
# preparation of a record. `read_phase` and `read_frequency` can also
# remove a trend and fill the gaps as they read, with the keywords
# `detrend` and `fillgaps`.
#
# ## Saving a record
#
# [`save`](@ref) writes a record as text, under comment lines that hold
# the file it came from, its sample interval and its start time:

saved = save(joinpath(dir, "pps.txt"), q)
foreach(println, readlines(saved)[1:6])

# The readers take the last two from those lines, so the file reads back as
# it was written, without keywords:

back = read_phase(saved)
back.x == q.x, start_time(back) == start_time(q)

# ## Saving results
#
# A [`StabilityResult`](@ref) is saved as a table with tab-separated
# columns, which a spreadsheet opens, under comment lines that describe it.
# The file begins with those lines and the names of the columns:

r = adev(q)
resultfile = save(joinpath(dir, "adev.tsv"), r)
foreach(println, readlines(resultfile)[1:7])

# [`load`](@ref) returns what `save` wrote, whether a record, a result or a
# [`StabilitySuite`](@ref):

load(resultfile).dev == r.dev

#-

suite = stability(q)
suitefile = save(joinpath(dir, "suite.tsv"), suite; source_file = "gpsdo.csv")
load(suitefile)[:mdev].dev == suite[:mdev].dev

# The command line and the explorer read the same files
# ([The record](../reference/cli.md#The-record);
# [Opening a file](../explorer/data.md#Opening-a-file)).

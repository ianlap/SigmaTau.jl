# Data and files

Read records from text files and write records, results and suites back out.
The readers accept whitespace-, comma- (`.csv`) and tab-separated (`.tsv`)
columns, skip `#` comment lines and Stable32 headers, and take τ₀ from the
time column, a `# tau0:` header or a Stable32 `Tau:` line.
[Reading and saving data](../tutorials/02_reading_your_data.md) walks through
them.

## Reading

```@docs
read_phase
read_frequency
```

## Timestamps and gaps

A time column places the samples. The grid step is the mean of the
regular steps, rounded to the shortest value that keeps every timestamp
within a quarter step of its slot, so timestamps rounded to a coarse
resolution (MJD to six decimals is 0.0864 s) still give the true τ₀. Each
step between timestamps is rounded to a whole number of grid steps, so
jittered timestamps keep every reading in order and a reading missing from
the log leaves a `NaN` sample in its slot instead of joining the samples
around it end to end. [`adev`](@ref), [`mdev`](@ref), [`hdev`](@ref),
[`mhdev`](@ref), [`tdev`](@ref), [`htdev`](@ref), [`mtie`](@ref) and
[`tierms`](@ref) skip the terms that read a missing sample; the total
deviations, Theo1 and the spectra return `NaN` until the gaps are filled
with [`fillgaps`](@ref) (or `fillgaps=true` on the readers). A record with
no missing timestamps is read exactly as written, and a column whose steps
fit no grid (not increasing, or steps a quarter or more off a whole multiple)
leaves the samples in file order. In a table, a row with no reading in the
clock column is a missing reading too, so a column logged every 5 minutes
in a 30-second table keeps its 5-minute step.

ISO 8601 timestamps (`2026-11-08T21:41:03.5Z`, with a `±hh:mm` offset or
none for UTC) and MJD time tags also give the record its UTC start time,
which [`save`](@ref) writes as `# start:` and the readers read back.

The readers stream the file line by line and keep only the columns they
return, so a multi-day log at a high rate costs the memory of its samples.

```@docs
start_time
```

## Instrument files

[`read_data`](@ref) reads the files instruments and their software write,
detecting the layout: TimeLab `.tim` logs (which is how the Microchip
53100A, 5120A/5125A and TimePod are usually recorded), Stable32 `.phd` /
`.frd`, Keysight 53230A-style counter readings, and logs with named
columns such as GPSDO exports. `sigmatau --input-format` and the explorer's
import dialog use the same readers.

```@docs
read_data
detect_format
INPUT_FORMATS
```

## Checking a record

```@docs
summarize
```

## Writing and loading

```@docs
save
load
```

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

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

## Writing and loading

```@docs
save
save_result
load_result
save_suite
load_suite
```

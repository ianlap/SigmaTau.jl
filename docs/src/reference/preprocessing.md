# Preprocessing

Clean a record before the analysis: remove a trend, fill gaps and replace
outliers. Each function returns a new record and leaves its input untouched.
[`read_phase`](@ref) and [`read_frequency`](@ref) can apply `detrend` and
`fillgaps` while reading.

```@docs
detrend
fillgaps
find_outliers
remove_outliers
```

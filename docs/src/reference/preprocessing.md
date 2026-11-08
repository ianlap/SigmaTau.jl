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

## Frequency drift

[`drift_fit`](@ref) fits Stable32's drift models to the frequency and
reports the aging per day with its standard error. Compare the fits and the
residuals before removing anything; `detrend` and the CLI's
`--remove-drift` remove a polynomial fitted to the phase.

```@docs
drift_fit
DriftFit
```

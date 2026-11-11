# Preprocessing

These functions prepare a record for the analysis: they remove a trend,
fill missing samples, and find or replace outliers. Each returns a new
record and leaves its argument as it was.
[Gaps, outliers and drift](../guide/preparing.md) says when each step is
called for and what it costs.

[`read_phase`](@ref) and [`read_frequency`](@ref) apply the first two
while they read, with the keywords `detrend` and `fillgaps`, and
[`read_data`](@ref) the second.

```@docs
detrend
fillgaps
find_outliers
remove_outliers
```

## Frequency drift

[`drift_fit`](@ref) fits a model to the fractional frequency of a record
and reports the drift rate with its standard error. It measures the drift
and removes nothing: [`detrend`](@ref) removes a line or a parabola fitted
to the samples.

```@docs
drift_fit
DriftFit
```

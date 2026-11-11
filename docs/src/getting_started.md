# Getting started

This tutorial takes a set of samples to a stability plot and a saved
result, in the Julia REPL. It uses the package and, for the plot, `Plots`:

```julia-repl
pkg> add SigmaTau Plots
```

## A record

An analysis starts from a record: the samples and the interval ``τ_0``
between them, in seconds. Phase samples, in seconds, make a
[`PhaseData`](@ref), and fractional-frequency samples a
[`FrequencyData`](@ref). Here the samples are simulated, 4096 values of
white frequency noise:

```@example start
using SigmaTau, Random

y = 1e-12 .* randn(Xoshiro(1), 4096)   # fractional frequency
f = FrequencyData(y, 1.0)              # sampled every second
```

The estimators take a record of either kind and never a bare vector, so the
sample interval always travels with the data.
[Conventions](explanation/conventions.md) defines phase, frequency and the
other quantities used here.

## A deviation

[`adev`](@ref) computes the overlapping Allan deviation:

```@example start
r = adev(f)
```

Each row is one averaging time ``τ = mτ_0``. With no second argument the
averaging factors ``m`` are spaced by octaves, up to the longest the record
supports. A named grid or the factors themselves can be given instead:

```@example start
adev(f, Decade)                 # m = 1, 10, 100, 1000
adev(f, [1, 3, 10, 30, 100])    # these factors
nothing # hide
```

The result is a [`StabilityResult`](@ref) whose fields are vectors with one
entry for each averaging time:

```@example start
r.tau[1], r.dev[1], r.noise_type[1], r.ci[1]
```

`noise_type` is the noise the package identified at that averaging time,
here white FM, and `ci` is the confidence interval that follows from it, at
68.3 % by default. Both can be set in the call:

```@example start
adev(f, [1, 10, 100]; confidence = 0.95, noise = :WHFM)
```

[Noise identification](explanation/noise_identification.md) and
[Confidence intervals](explanation/confidence.md) explain how the two
columns are obtained and how far to trust them. `ci=false` leaves both out.

Every other estimator is called the same way: [`mdev`](@ref),
[`hdev`](@ref), [`totdev`](@ref) and the rest of the
[estimator reference](reference/estimators.md).
[Choosing an estimator](guide/choosing.md) says which to use when.

## Several deviations

[`stability`](@ref) runs a set of estimators on one record and returns a
[`StabilitySuite`](@ref), indexed by name:

```@example start
suite = stability(f; devs = (:adev, :mdev, :hdev, :totdev))
suite[:mdev].dev[1:3]
```

## A plot

Loading `Plots` turns on the plot recipes. A result or a suite plots on
logarithmic axes with its confidence intervals:

```@example start
using Plots
plot(suite)
```

All four fall as ``τ^{-1/2}``, the slope of white FM
([Power-law noise](explanation/noise.md)).
[Plotting and tables](guide/plotting.md) has the other plots: bands,
several clocks on one axis, spectra and dynamic stability.

## Saving

[`save`](@ref) writes a record, a result or a suite as text that any
program can read, and [`load`](@ref) reads it back:

```@example start
file = joinpath(mktempdir(), "adev.tsv")
save(file, r)
load(file).dev == r.dev
```

## A file of your own

[`read_phase`](@ref) and [`read_frequency`](@ref) read columns of text, and
[`read_data`](@ref) reads the files of instruments and of other programs:

```julia
p = read_phase("clock.txt")                       # time and phase columns
p = read_phase("clock.txt"; value_col = 3, scaling = 1e-9)   # nanoseconds
p = read_data("clock.tim")                        # a TimeLab file
```

[Reading and saving data](tutorials/02_reading_your_data.md) goes through
the readers, and [Gaps, outliers and drift](guide/preparing.md) through the
checks a measured record needs before its deviations mean anything.

## Next

- [A drifting clock and the Hadamard family](tutorials/01_drifting_clock_hadamard.md)
  shows what drift does to the Allan deviation and what to use instead.
- [Choosing an estimator](guide/choosing.md) maps the question asked of a
  record onto an estimator.
- [The explorer](explorer/index.md) does the same analysis without code.

# Getting started

This page takes you from a vector of numbers to a plotted, saved stability
analysis in the REPL. It assumes SigmaTau is installed (see [Home](index.md)).

## 1. Wrap your samples in a record

Every analysis starts from a record: the samples plus the sample interval
τ₀ in seconds. Phase (time-error) samples in seconds go in a
[`PhaseData`](@ref); fractional-frequency samples go in a
[`FrequencyData`](@ref).

```@example start
using SigmaTau, Random
Random.seed!(1)

τ₀ = 1.0
y  = 1e-12 .* randn(4096)          # white FM, fractional frequency
x  = τ₀ .* cumsum(y)               # the same record as phase, in seconds

p = PhaseData(x, τ₀)
f = FrequencyData(y, τ₀)
nothing # hide
```

Wrapping is the one way in: the estimators take a record, never a bare
vector, so the sample interval always travels with the data. Both kinds of
record work everywhere; frequency records are integrated to phase internally.

To read a file instead, use [`read_phase`](@ref) or
[`read_frequency`](@ref). They handle two-column text files, CSV and TSV,
Stable32 `.DAT` exports and files that [`save`](@ref) wrote, and they infer τ₀
from the time column or the header:

```julia
p = read_phase("clock.dat")                          # time and phase columns
p = read_phase("counter.csv"; value_col=3, scaling=1e-9)   # ns → s
p = read_phase("phase.txt"; time_col=0, value_col=1, tau0=1.0)   # one column, no time
```

[Reading and saving data](tutorials/02_reading_your_data.md) covers the
readers, gaps and outliers in detail.

## 2. Compute a deviation

Each estimator has three call forms:

```@example start
adev(p)                    # default: octave-spaced τ up to the estimator's limit
adev(p, Decade)            # a named τ grid (a TauMode)
adev(p, 1:10)              # any integer range or vector of averaging factors m
r = adev(p, [1, 2, 4, 8, 16, 32, 64])   # τ = m·τ₀
```

The result is a [`StabilityResult`](@ref), shown above as a table. Its
fields are plain vectors with one entry per τ: `r.tau` (seconds), `r.dev`,
`r.noise_type`, `r.ci` (the `(lo, hi)` bounds), `r.edf` and `r.neff`.

```@example start
r.dev[1], r.noise_type[1], r.ci[1]
```

`r.edf` holds the equivalent degrees of freedom behind each interval and
`r.neff` the number of terms averaged. Intervals default to 68.3 %
(`DEFAULT_CONFIDENCE`); pass `confidence=0.95` for wider ones, or
`ci=false` to skip the intervals when you only want the deviation (the
total estimators still identify the noise for their bias correction unless
you also pass `correct_bias=false`).

```@example start
adev(p; confidence=0.95).ci[1]
```

## 3. Compute several at once

[`stability`](@ref) runs a set of estimators on one record and returns a
[`StabilitySuite`](@ref), indexed by estimator name:

```@example start
suite = stability(p; devs=(:adev, :mdev, :hdev, :totdev))
suite[:hdev].dev[1:3]
```

Not sure which estimators you need? See
[Choosing an estimator](guide/choosing.md).

## 4. Plot

Loading `Plots` turns on the plot recipes: `plot(result)` draws σ(τ) on
log-log axes with confidence-interval bars, and `plot(suite)` overlays the
suite.

```@example start
using Plots
plot(suite)
```

[Plotting and tables](guide/plotting.md) shows the other styles: bands,
overlays of several clocks, spectra and dynamic maps.

## 5. Save

[`save`](@ref) writes records, results and suites as plain text that any
tool can read; [`read_phase`](@ref) reads records back and [`load`](@ref)
reads results and suites.

```@example start
dir = mktempdir()
save(joinpath(dir, "adev.tsv"), r)
load(joinpath(dir, "adev.tsv")).dev ≈ r.dev
```

## Where next

- [A drifting clock and the Hadamard family](tutorials/01_drifting_clock_hadamard.md):
  why drift spoils ADEV and what to use instead.
- [Noise, EDF and confidence](explanation/noise_and_confidence.md): what the
  noise types and intervals mean.
- [Reference](reference/estimators.md): every estimator and its keywords.

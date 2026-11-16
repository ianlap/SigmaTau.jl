# Plotting and tables

The package draws nothing itself. Its plot recipes load when `Plots` does,
or any other package built on `RecipesBase`, and a result prints as a
table and converts to one.

```@example plotting
using SigmaTau, Plots, Random

# white FM of 1e-12 at 1 s and random-walk FM, 8192 samples at 1 s
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 3e-14),
              rng = Xoshiro(3))
r = adev(p)
nothing # hide
```

## A deviation

`plot(r)` draws a [`StabilityResult`](@ref) on logarithmic axes, with a bar
for the confidence interval at each averaging time:

```@example plotting
plot(r)
```

The attributes of `Plots` apply as usual, and `ci_band = true` draws the
interval as a band:

```@example plotting
plot(r; ci_band = true, fillalpha = 0.25, lw = 1.5, title = "Clock A",
     legend = :bottomleft)
```

The interval is the one the call computed: `adev(p; confidence = 0.95)`
widens it and `adev(p; ci = false)` leaves it out. The limits are also
vectors, `ci_lower(r)` and `ci_upper(r)`. [`mtie`](@ref) and
[`tierms`](@ref) have no interval and plot as plain curves.

## Several curves

A [`StabilitySuite`](@ref) plots its deviations on one pair of axes:

```@example plotting
plot(stability(p); legend = :bottomleft)
```

So does a vector of results, which is the way to compare clocks. `label`
takes one name for each:

```@example plotting
q = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 3e-12), rng = Xoshiro(6))
plot([adev(p), adev(q)]; label = ["clock A" "clock B"], legend = :bottomleft)
```

The clocks that [`ncornered_hat`](@ref) and [`groslambert`](@ref) separate
come as a [`HatResult`](@ref), whose curves take the names of the clocks.
`plot!(adev(q))` adds a curve to a figure that exists.

The time deviations ([`tdev`](@ref), [`htdev`](@ref), [`ttotdev`](@ref))
and the time errors ([`mtie`](@ref), [`tierms`](@ref)) are in seconds.
Beside the σ_y estimators they get axes of their own:

```@example plotting
plot(stability(p; devs = (:adev, :mdev, :tdev)); legend = :bottomleft)
```

Such a figure has two panels, so `plot!` cannot add a mix of the two kinds
to a figure of one: draw them together, or give `plot!` one kind.

## A reference line

The slope of a curve names the noise
([Power-law noise](../explanation/noise.md) has the table). A line of a
given slope is an ordinary series, anchored here at the first point:

```@example plotting
line = r.dev[1] .* (r.tau ./ r.tau[1]) .^ (-1/2)

plot(r; legend = :bottomleft)
plot!(r.tau, line; ls = :dash, color = :gray, label = "τ^(-1/2), white FM")
```

## Spectra

[`psd_y`](@ref), [`psd_x`](@ref) and [`ssb_phase_noise`](@ref) return a
[`SpectralResult`](@ref), which plots against a logarithmic frequency axis:

```@example plotting
plot(psd_y(p))
```

``ℒ(f)`` from `ssb_phase_noise(p; f_carrier = 10e6)` is in dBc/Hz and gets
a linear vertical axis.

## Dynamic stability

[`dadev`](@ref) and [`dhdev`](@ref) return the deviation in a window that
slides along the record. The plot is the one of
[SP1065 §5.2.20](@cite riley-2008-sp1065): a curve for each window, set
out along a time axis. The record here is white FM whose level rises tenfold
halfway through, which a deviation of the whole record would average over:

```@example plotting
rng = Xoshiro(7)
y = [1e-11 .* randn(rng, 4096); 1e-10 .* randn(rng, 4096)]
d = dadev(FrequencyData(y, 1.0); window = 1024, step = 512)

plot(d; camera = (40, 30), xlabel = "τ (s)", ylabel = "t (s)", size = (700, 480))
```

`camera` turns the view, and the short axis labels keep clear of the tick
labels. The result holds the map itself, `d.dev`, with a row for each
window center in `d.t` and a column for each averaging time in `d.tau`,
for any other picture of it:

```@example plotting
heatmap(d.t, d.tau, log10.(d.dev'); yscale = :log10, xlabel = "time (s)",
        ylabel = "τ (s)", title = "log10 of the Allan deviation")
```

## Models, holdover, steering and time scales

A [`PowerLawFit`](@ref) plots the measured points, the model through them
and, dashed, each noise type of the model alone:

```@example plotting
fit = powerlaw_fit(r)
plot(fit; legend = :bottomleft)
```

A [`HoldoverResult`](@ref) plots the RMS time error, the error that a
share of the trials stay within and the largest, against holdover time,
with the limit asked for as a line. A [`HoldoverPrediction`](@ref) plots
the error the noise model predicts and what each noise type adds to it:

```@example plotting
plot(holdover(p; train = 512, tau = [16, 64, 256, 1024], limit = 1e-9);
     legend = :topleft)
```

A [`KalmanResult`](@ref) plots the phase and the frequency that the filter
estimates, over time, each with a band of one standard deviation. The
frequency includes the flicker states of the model:

```@example plotting
plot(kalman_filter(p, clock_model(fit)))
```

A [`TimescaleResult`](@ref) plots each clock against the scale and, for
the algorithms that weigh the clocks, their weights. A
[`SteeringResult`](@ref) plots the steered phase and the frequency
correction over time, and a [`SteeringPrediction`](@ref) the Allan
deviation of the steered clock.

```@example plotting
ms = [ClockModel(; tau0 = 1.0, q1 = q, q2 = 1e-6q) for q in (1e-24, 4e-24, 1.6e-23)]
cc, _ = simulate_ensemble(ms, 4096; reference = 0, rng = Xoshiro(1))
plot(timescale(cc; algorithm = :at1))
```

## Tables

A result prints as a table, one row for each averaging time:

```@example plotting
adev(p, [1, 10, 100, 1000])
```

A result or a suite is also a table in the sense of `Tables.jl`, with the
columns `deviation_type`, `tau`, `dev`, `noise_type`, `ci_lower`,
`ci_upper`, `edf` and `neff`, and so is a [`HatResult`](@ref), with a
`clock` column in front. Any package that reads such tables takes one
directly:

```julia
using DataFrames, CSV

df = DataFrame(r)
CSV.write("adev.csv", stability(p))        # one row for each deviation and τ
```

[`save`](@ref) needs none of these packages. It writes a result or a suite
as tab-separated text, with the settings of the call in comment lines, and
[`load`](@ref) reads the file back.

## Saving a figure

`savefig` takes the format from the extension. `size` is in pixels and
`dpi` sets the resolution of a PNG, while PDF and SVG are vector formats:

```julia
fig = plot(r; ci_band = true, size = (800, 500), dpi = 300)
savefig(fig, "adev.png")
savefig(fig, "adev.pdf")
```

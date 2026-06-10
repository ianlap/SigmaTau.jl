# Plotting

SigmaTau ships its plot recipes in a package extension
(`ext/SigmaTauRecipesBaseExt.jl`). The extension loads automatically when
`RecipesBase` is in the session — in practice, the moment you load `Plots`
(or any other package that loads `RecipesBase`; a backend package alone, such
as `GR`, does not). SigmaTau itself takes no plotting dependency,
so headless batch analysis never pays for a graphics stack.

This page is a cookbook: each section is a recipe you can paste into the REPL.

## Setup

```julia
using Plots, SigmaTau
```

`Plots` defaults to the GR backend, which needs no extra installation. The
recipes are backend-agnostic; the figures on this page are rendered through
the docs build's backend (PGFPlotsX when a LaTeX engine is available, GR
otherwise).

## A single result

`plot(result)` on a [`StabilityResult`](reference/types.md) draws the σ–τ
curve as a line on log-log axes:

```@example plotting
using Plots, SigmaTau, Random
Random.seed!(42)

# Synthetic record: WFM at σ_y(1 s) = 1e-12 plus RWFM at 3e-14,
# 8192 phase samples at τ₀ = 1 s — a fall-then-rise curve.
p = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 1e-12, -2 => 3e-14))
r = adev(p)

plot(r)
```

The recipe sets, as overridable defaults:

- `xscale = :log10`, `yscale = :log10`,
- `xlabel = "Averaging Time τ (s)"`,
- `ylabel` and legend `label` equal to the deviation name (`"ADEV"` here,
  taken from `result.deviation_type`).

When the result carries confidence bounds (`ci=true`, the default), the
recipe attaches asymmetric error bars: the lower offset is
`dev .- ci_lower` and the upper offset is `ci_upper .- dev`, so each bar
spans exactly `[ci_lower, ci_upper]`. The bars are asymmetric because the
χ² confidence interval is asymmetric about the point estimate.

Every default yields to a standard `Plots` attribute passed at the call
site. The one exception is the series type, which the recipe fixes to a
line path — add markers with `marker = :circle` rather than
`seriestype = :scatter`:

```@example plotting
plot(r; lw = 1.5, marker = :circle, ms = 3,
     title = "Overlapping Allan deviation", legend = :bottomleft)
```

## Confidence-interval styles

The recipe supports one custom attribute, `ci_band`:

| Call                          | CI rendering                       |
|-------------------------------|------------------------------------|
| `plot(r)`                     | error bars (default)               |
| `plot(r; ci_band = true)`     | filled band (ribbon)               |
| `plot(adev(p; ci = false))`   | none — the CI vectors are empty    |

```@example plotting
plot(r; ci_band = true, fillalpha = 0.3)
```

`ci_band` is consumed by the recipe, so it never reaches the backend as an
unknown attribute. The band edges are `ci_lower` and `ci_upper`, the same
bounds the error bars span; only the rendering changes.

The bounds are plain vectors on the result, so you can also draw them
yourself when you want full control over the style:

```@example plotting
plot(r.tau, r.dev;
     xscale = :log10, yscale = :log10,
     xlabel = "Averaging Time τ (s)", label = "ADEV")
plot!(r.tau, r.ci_lower; ls = :dot, color = :gray, label = "68.3 % CI")
plot!(r.tau, r.ci_upper; ls = :dot, color = :gray, label = "")
```

The default confidence level is `DEFAULT_CONFIDENCE = 0.683` (1σ); pass
`confidence = 0.95` to the deviation call to widen the bounds. `mtie` has no
published EDF model, so its results carry no CI and always plot as a bare
curve.

## Overlaying several deviations

`plot(suite)` on a [`StabilitySuite`](reference/types.md) overlays every
result on one set of axes, one labelled series per deviation, with the
generic y-label `"Deviation"`:

```@example plotting
suite = stability(p)          # :adev, :mdev, :hdev, :tdev by default
plot(suite; legend = :bottomleft)
```

!!! note "Mixed units in the default suite"
    The default suite includes `tdev`, which is a *time* deviation in
    seconds, while ADEV/MDEV/HDEV are dimensionless fractional-frequency
    deviations. The overlay is still useful for shape comparison, but the
    ordinate mixes units. Index the suite (`suite[:tdev]`) to plot it on its
    own axes.

A `Vector{StabilityResult}` overlays the same way — useful for comparing
clocks rather than deviations: `plot([adev(p1), adev(p2)])`. Both series
then inherit the same default label (`"ADEV"`), so for clock comparisons the
manual form with explicit labels is usually clearer:

```@example plotting
p2 = noise_gen(PhaseData, 8192, 1.0; sigma1 = Dict(0 => 3e-12))

plot(adev(p);   label = "clock A", legend = :bottomleft)
plot!(adev(p2); label = "clock B")
```

Every `plot!` call adds to the existing figure; the user-supplied `label`
overrides the recipe's default.

## Log-log conventions and slope guides

Because the log scales are recipe *defaults*, the usual `Plots` attributes
restyle the axes. Decade ticks come from powers of ten:

```@example plotting
plot(r; xticks = 10.0 .^ (0:3),
     xlabel = "τ (s)", ylabel = "σ_y(τ)", legend = :bottomleft)
```

On log-log axes the slope of σ_y(τ) encodes the dominant noise type
[riley-2008-sp1065](@cite); α is the spectral exponent of the
fractional-frequency PSD, `S_y(f) ∝ f^α`:

| Noise type            | α  | ADEV slope |
|-----------------------|---:|:-----------|
| White PM (`:WHPM`)    |  2 | τ⁻¹        |
| Flicker PM (`:FLPM`)  |  1 | ≈ τ⁻¹      |
| White FM (`:WHFM`)    |  0 | τ^(−1/2)   |
| Flicker FM (`:FLFM`)  | −1 | τ⁰         |
| Random-walk FM (`:RWFM`) | −2 | τ^(+1/2) |

ADEV cannot separate white PM from flicker PM — both fall as roughly τ⁻¹ —
which is one reason `mdev` (slope τ^(−3/2) under white PM) is computed
alongside it.

A reference slope line is an ordinary series anchored at a point on the
curve. For the τ^(−1/2) guide of white FM, anchored at the first point:

```@example plotting
guide = r.dev[1] .* (r.tau ./ r.tau[1]) .^ (-1/2)

plot(r; legend = :bottomleft)
plot!(r.tau, guide; ls = :dash, color = :gray, label = "τ^(-1/2) (WHFM)")
```

The synthetic record above follows the guide at short τ (the white-FM
regime) and peels upward at long τ where the random-walk FM component takes
over.

## Spectral densities

The same extension plots [`SpectralResult`](reference/types.md)s from `Sy`,
`Sx`, and `L`. `S_y(f)` and `S_x(f)` are power densities and render on
log-log axes with the units in the y-label; `ℒ(f)` is already in dBc/Hz, so
it gets a log frequency axis with a linear (dB) ordinate. The recipe drops
the DC bin, which is undefined on a log frequency axis:

```@example plotting
plot(Sy(p))
```

`L` needs the carrier frequency to convert from S_x: pass `f_carrier`
(in Hz) to the `L` call itself, e.g. `plot(L(p; f_carrier = 10e6))`.

## Publication export

`size` (pixels) and `dpi` are standard `Plots` attributes; `savefig` picks
the format from the file extension. `dpi` affects only raster formats such
as PNG — PDF and SVG are vector formats and ignore it:

```julia
fig = plot(r; ci_band = true, size = (800, 500), dpi = 300,
           lw = 1.5, legend = :bottomleft)

savefig(fig, "adev.png")   # raster, 300 dpi
savefig(fig, "adev.pdf")   # vector
```

For LaTeX documents, switch the backend to PGFPlotsX (`Plots.pgfplotsx()`,
requires a TeX engine) and the figure fonts will match the surrounding
text — this is what the docs build itself does.

## Where to next

- [Tutorial 2: Computing the Allan deviation](tutorials/02_compute_adev.md) —
  the recipe applied to a worked WPM + RWFM record.
- [Tutorial 6: Three-cornered hat](tutorials/06_three_cornered_hat.md) — a
  fully manual log-log plot for derived (non-`StabilityResult`) quantities.
- [Theory: Overview](theory/overview.md) — what the slopes mean physically.

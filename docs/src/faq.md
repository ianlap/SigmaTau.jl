# FAQ and Troubleshooting

Short answers to the questions that come up most often, with pointers into the
rest of the documentation. If your question is about Stable32 equivalence,
start with [Migrating from Stable32](migration_from_stable32.md).

## I get "Package Plots not found" when running an example

The example scripts live in their own Julia environment. The package's runtime
`Project.toml` deliberately does not depend on `Plots` or `PGFPlotsX` — they
are visualisation-only — so a clean checkout run under `--project=.` fails
with `Package Plots not found`. Run the examples under the `examples` project
instead:

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'   # first run only
julia --project=examples examples/00_julia_for_metrologists.jl
```

`examples/Project.toml` carries `Plots` and `PGFPlotsX` and points at the
package source via a path entry, so the scripts always run against your
checkout.

## How do I read a single-column .DAT file?

`read_phase` defaults to the common two-column layout (`time_col=1`,
`value_col=2`), inferring the sample interval from the time column. A bare
single-column `.DAT` has neither, so you must say so explicitly:

```julia
p = read_phase("clock.DAT"; time_col=0, value_col=1, tau0=1.0)
```

`time_col=0` means "no time column" and `value_col=1` reads the lone column.
`tau0` is required here because with no timestamps there is nothing to infer
the sample interval from; `read_phase` throws an `ArgumentError` if you omit
it. `read_frequency` takes the same keywords.

## Why do my confidence intervals differ slightly from Stable32?

The deviation values themselves agree — the core kernels are validated against
Stable32 to the 5-significant-figure precision of its published outputs (see
[Validation](validation/methodology.md) and the
[Stable32 comparison tables](validation/stable32.md)). The intervals around
them can differ for two reasons.

First, SigmaTau identifies the dominant noise type independently at every τ
(lag-1 ACF with a B1/R(n) fallback), and by default it quadratically detrends
each decimated subseries before classifying; on a borderline record a τ point
can land in a different power-law family than Stable32 assigns, which changes
the EDF and hence the interval.
(`identify_noise(x, m_values; detrend=false)` reproduces Stable32's
no-detrend convention point-for-point.)

Second, for the total family above α = 0 the published EDF coefficient
tables run out, and the two programs substitute different EDF models there.
For the core estimators the EDF model is the same: SigmaTau uses the
Greenhall–Riley spectral-sum formulas [greenhall-2003-edf-stability](@cite),
which Stable32 has also used since v1.41 — wherever the two programs agree
on the noise type, the interval bounds agree to within 0.05 % (see
[Validation: Stable32](validation/stable32.md)).

## When should I use the Hadamard family instead of Allan?

Use `hdev` / `mhdev` / `htdev` when the record contains linear frequency
drift, or noise redder than random-walk FM. The Hadamard kernel is a third
difference of phase rather than a second difference, and the third difference
annihilates a linear-in-`t` term in `y(t)` — so constant frequency drift drops
out of the statistic instead of adding a spurious `+τ` slope, which is exactly
what happens to ADEV on a drifting rubidium or cesium record
[riley-2008-sp1065](@cite). The same extra difference order keeps the variance
integral convergent down to `α = −4`, where the Allan family diverges for
`α ≤ −3`. For a drift-free record in the usual `−2 ≤ α ≤ +2` range, ADEV
suffices and has slightly tighter intervals. See
[Theory: Allan family](theory/allan_family.md).

## Why are MHTOTDEV's confidence intervals "Monte-Carlo-calibrated"?

MHTOTDEV is novel to SigmaTau: no other library computes it and no published
reference gives its bias factor or its EDF. Every other total estimator
inherits a published coefficient table; MHTOTDEV has none, so SigmaTau
measured both quantities by Monte Carlo against synthesized known-α noise,
following the methodology of the NIST total-variance papers
[howe-2000-tothvar-ptti](@cite). The fitted coefficients ship in the package
(the EDF table and the `bias_correction(…, :mhtot, …)` entries), and the full
procedure — estimators, fit forms, validity window, and the resulting
coefficient table — is documented in
[Theory: MHTOTDEV bias and EDF](theory/mhtotdev_bias_edf.md) so the
calibration is reproducible.

## What does the τ/τ₀ ≥ 16 note on the total estimators mean?

The EDF coefficient models for the modified and Hadamard totals are stated to
be valid only for `τ/τ₀ ≥ 16`, where `τ₀` is the sample interval; Howe et al.
note their Hadamard-total EDF model "should be used only if … τ/τ₀ ≥ 16"
[howe-2000-tothvar-ptti](@cite), and the SP1065 MTOTDEV coefficient formula
carries the same floor. Below that floor the boundary extension that
defines a total estimator confers no real advantage over the plain estimator,
and the coefficient fits degrade. The σ values themselves are valid at every
τ — the window applies to the EDF and bias models, not the deviation.
SigmaTau respects it: for MTOTDEV at `m < 16` the EDF falls back to the
modified-Allan Greenhall–Riley value (which matches Stable32's reported
degrees of freedom there), and the MHTOTDEV coefficients were fit only over
the `τ/τ₀ ≥ 16` window.

## Which Julia version do I need?

Julia 1.11 or later. The `[compat]` section of `Project.toml` declares
`julia = "1.11"`, and the examples environment requires the same.

## How do I choose a τ grid?

Three ways. Calling a deviation with data only uses the default octave-spaced
grid (`1, 2, 4, 8, …`), the same spacing Stable32 defaults to. Passing a
`Vector{Int}` of averaging factors gives full control (`τ = m·τ₀`, where `m`
is the averaging factor and `τ₀` the sample interval). Passing a `TauMode`
selector picks a named spacing:

```julia
adev(p)                            # octave (default)
adev(p, [1, 2, 4, 8, 16, 32, 64])  # explicit averaging factors
adev(p, Decade)                    # decade-spaced
```

The selectors are `AllTaus`, `Octave`, `HalfOctave`, `QuarterOctave`,
`Decade`, and `HalfDecade`; `tau_values(mode, N, kernel)` returns the
underlying integer grid. All grids are clamped to the kernel's algorithmic
m-max for the record length. Octave spacing is the usual choice: denser grids
(`AllTaus`) look smoother but adjacent points are strongly correlated, so they
add little independent information.

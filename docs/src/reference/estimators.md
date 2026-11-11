# Estimators

Every estimator takes a [`PhaseData`](@ref) or a [`FrequencyData`](@ref)
record and returns a [`StabilityResult`](@ref). The averaging factors
``m``, with ``τ = mτ_0``, come second, in one of three forms:

```julia
adev(data)                  # the octave grid, up to the limit of the estimator
adev(data, Decade)          # another grid, a TauMode
adev(data, [1, 10, 100])    # the averaging factors, as a vector or a range
```

The estimators share their keywords, shown here with their defaults. The
docstring of [`adev`](@ref) describes the common ones in full:

| Keyword | Taken by | Effect |
|:--------|:---------|:-------|
| `ci=true` | all | Identify the noise type at each ``τ`` and report the equivalent degrees of freedom and the confidence interval. `false` returns the deviation and the number of terms, and from `totdev`, `mtotdev`, `htotdev` and `mhtotdev` the noise type that their bias factor rests on |
| `confidence=0.683` | all | The level of the interval |
| `noise=:auto` | all but `mtie` and `tierms` | The noise type behind the interval and the bias factor: identified at each ``τ``, or stated |
| `correct_bias=true` | the total estimators and `theo1` | The bias factors divided out: those of this package, `:published`, `:stable32`, or none with `false` |
| `remove_drift=true` | `hdev`, `mhdev`, `htdev`, `htotdev`, `mhtotdev` | Identify the noise type on the record less its frequency drift. `false` identifies it on the record as given, as Stable32 does |
| `overlapping=true` | `adev`, `hdev` | `false` gives the non-overlapped estimator |

`mtie` and `tierms` report no interval and accept `ci` and `confidence`
without effect.

[Choosing an estimator](../guide/choosing.md) compares the estimators.
[Noise identification](../explanation/noise_identification.md),
[Confidence intervals](../explanation/confidence.md) and
[Bias correction](../explanation/bias.md) explain what the keywords
select. Some estimators leave out the terms that a missing sample enters
and the others need a complete record:
[Timestamps and gaps](@ref) lists them.

## Allan family

```@docs
adev
mdev
tdev
```

## Hadamard family

```@docs
hdev
mhdev
htdev
```

## Total estimators

```@docs
totdev
mtotdev
ttotdev
htotdev
mhtotdev
```

## Thêo1 and ThêoH

```@docs
theo1
theoh
```

## Parabolic deviation

```@docs
pdev
```

## Time interval error

```@docs
mtie
tierms
```

## Telecom masks

The ITU-T recommendations limit the wander of a clock, or of a
synchronization interface, by masks on MTIE and TDEV over the observation
interval ``τ``. [`mask_check`](@ref) compares an [`mtie`](@ref) or a
[`tdev`](@ref) result with one of them, and [`TELECOM_MASKS`](@ref) lists
the masks with the measurement that the recommendations prescribe. The
command line checks a result with [`analyze --mask`](cli.md), and the
explorer draws the limits on its
[Stability](../explorer/stability.md#Masks) page.

```@docs
TELECOM_MASKS
telecom_mask
mask_check
```

## Dynamic stability

[`dadev`](@ref) and [`dhdev`](@ref) compute a deviation inside a window
that slides along the record and return the map ``σ_y(t, τ)`` as a
[`DynamicStabilityResult`](@ref).
[Dynamic stability](../guide/plotting.md#Dynamic-stability) draws it.

```@docs
dadev
dhdev
```

## Three-cornered hat

[`ncornered_hat`](@ref) and [`groslambert`](@ref) separate the deviations
of three or more clocks that are measured against each other. The tutorial
[Three clocks and the three-cornered hat](../tutorials/03_three_cornered_hat.md)
works through both.

```@docs
ncornered_hat
groslambert
```

## Two channels

[`cross_adev`](@ref) estimates the Allan deviation of a device from two
simultaneous records that share the device and nothing else: the noise of
the two channels, or of the two references, averages toward zero.
[Two-channel measurements](../guide/cross.md) says when the method
applies and what it cannot remove. The cross-spectra are on the
[Spectra](spectra.md#Cross-spectra) page.

```@docs
cross_adev
CrossStabilityResult
```

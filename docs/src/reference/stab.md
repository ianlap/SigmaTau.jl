# SigmaTau — Stability analysis

## Allan family deviations

```@docs
adev
mdev
hdev
mhdev
tdev
htdev
```

## Total family deviations

```@docs
totdev
mtotdev
ttotdev
htotdev
mhtotdev
```

## MTIE and PDEV

```@docs
mtie
pdev
```

## Tau grids and suite API

```@docs
TauMode
tau_values
DEFAULT_DEVIATIONS
stability
```

The exported [`TauMode`](@ref) instances are `AllTaus`, `Octave`,
`HalfOctave`, `QuarterOctave`, `Decade`, and `HalfDecade`.

## Noise identification

```@docs
identify_noise
noise_gen
```

## Spectral densities

```@docs
Sy
Sx
L
```

## EDF, bias, and confidence intervals

```@docs
calculate_edf
bias_correction
confidence_intervals
DEFAULT_CONFIDENCE
```

## I/O and preprocessing

```@docs
read_phase
read_frequency
detrend
fillgaps
save_result
load_result
save_suite
load_suite
```

## Advanced / research kernels

These are not exported, but remain reachable as `SigmaTau._adev_core` and
friends for benchmarking, parity validation against reference implementations,
and research use. The underscore prefix signals that the calling convention is
not stability-guaranteed across minor versions; prefer the wrappers above for
application code.

```@docs
SigmaTau._adev_core
SigmaTau._mdev_core
SigmaTau._tdev_core
SigmaTau._hdev_core
SigmaTau._mhdev_core
SigmaTau._totdev_core
SigmaTau._mtotdev_core
SigmaTau._htotdev_core
SigmaTau._mhtotdev_core
SigmaTau._mtie_core
SigmaTau._pdev_core
SigmaTau._kernel_m_max
SigmaTau._default_m_values
```

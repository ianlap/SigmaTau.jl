# SigmaTau.Stab — Stability analysis

## Allan family deviations

```@docs
adev
mdev
hdev
mhdev
tdev
htdev
ldev
```

## Total family deviations

```@docs
totdev
mtotdev
ttotdev
htotdev
mhtotdev
```

## Noise identification

```@docs
identify_noise
```

## EDF, bias, and confidence intervals

```@docs
calculate_edf
bias_correction
confidence_intervals
```

## I/O

```@docs
save_result
load_result
```

## Advanced / research kernels

These are exported for benchmarking, parity validation against reference
implementations, and research use. The underscore prefix signals that the
calling convention is not stability-guaranteed across minor versions; prefer
the wrappers above for application code.

```@docs
_adev_core
_mdev_core
_tdev_core
_hdev_core
_mhdev_core
_totdev_core
_mtotdev_core
_htotdev_core
_mhtotdev_core
_mtie_core
_pdev_core
```

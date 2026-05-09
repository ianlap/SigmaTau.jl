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

## Internal kernels

These are exported for advanced use and benchmarking. Most users should
prefer the wrappers above.

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
```

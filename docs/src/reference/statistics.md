# Noise and statistics

Synthesize power-law noise, identify it, and compute the EDF, bias factors
and confidence intervals that every estimator applies when `ci=true`.
[Noise, EDF and confidence](../explanation/noise_and_confidence.md) explains
the method; these functions are the building blocks, for when you need them
outside an estimator call.

## Noise

```@docs
noise_gen
identify_noise
```

## EDF, bias and confidence intervals

```@docs
DEFAULT_CONFIDENCE
calculate_edf
bias_correction
confidence_intervals
```

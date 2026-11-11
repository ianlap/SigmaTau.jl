# Noise and statistics

An estimator called with `ci=true` identifies the noise type at each
averaging time, computes the equivalent degrees of freedom and the
confidence limits, and corrects the bias of a total estimator. The
functions that do each step are exported for use on their own, beside
[`noise_gen`](@ref), which simulates the records that this documentation
tries them on.

| Function | Step |
|:---------|:-----|
| [`noise_gen`](@ref) | Simulates power-law noise ([Simulated noise](../explanation/noise.md#Simulated-noise)) |
| [`identify_noise`](@ref) | Names the noise type at each ``τ`` ([Noise identification](../explanation/noise_identification.md)) |
| [`noise_id_nnls`](@ref) | Fits an amount of every noise type ([Fitting a model to a record](../explanation/noise.md#Fitting-a-model-to-a-record)) |
| [`calculate_edf`](@ref) | Gives the equivalent degrees of freedom ([Confidence intervals](../explanation/confidence.md#Equivalent-degrees-of-freedom)) |
| [`bias_correction`](@ref) | Gives the bias factor of a total estimator ([Bias correction](../explanation/bias.md)) |
| [`confidence_intervals`](@ref) | Gives the confidence limits ([Confidence intervals](../explanation/confidence.md)) |

## Noise

```@docs
noise_gen
identify_noise
noise_id_nnls
```

## EDF, bias and confidence intervals

```@docs
DEFAULT_CONFIDENCE
calculate_edf
bias_correction
confidence_intervals
```

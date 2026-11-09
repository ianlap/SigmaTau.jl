# Spectra

One-sided Welch estimates of the frequency and phase spectra and the
single-sideband phase noise, returned as a [`SpectralResult`](@ref).
`S_x(f) = S_y(f) / (2πf)²` and `ℒ(f) = ½ (2π f_carrier)² S_x(f)`.
[Plotting](../guide/plotting.md#Spectral-densities) shows them with
power-law slope guides.

```@docs
psd_y
psd_x
ssb_phase_noise
```

## Power-law model

The power-law model `S_y(f) = Σ h_α f^α` links the two domains.
[`powerlaw_fit`](@ref) fits the h_α to a record's σ(τ) or to a spectrum,
including a datasheet's ℒ(f) points. [`powerlaw_sigma`](@ref) and
[`powerlaw_psd`](@ref) turn a model into σ(τ) or into S_y(f), S_x(f) or
ℒ(f).

```julia
f = powerlaw_fit([1, 10, 100, 1e3, 1e4], [-95, -125, -140, -150, -152];
                 quantity=:L, f_carrier=10e6)
powerlaw_sigma(f, [1.0, 10.0]; tau0=1e-3)    # ADEV sampled at 1 ms
```

```@docs
powerlaw_fit
PowerLawFit
powerlaw_sigma
powerlaw_psd
```

## Clock model

[`clock_model`](@ref) turns a power-law fit into the state-space
[`ClockModel`](@ref) the Kalman filters run on: white PM becomes the
measurement noise, white, random-walk and random-run FM drive the phase,
frequency and drift states. [The clock model](../explanation/clock_model.md)
derives the mapping and the matrices.

```julia
m = clock_model(powerlaw_fit(p; alphas=[2, 0, -2]))
state_transition(m, 1.0), process_noise(m, 1.0)
powerlaw_sigma(m, [1.0, 100.0])    # σ(τ) the model stands for
```

```@docs
ClockModel
clock_model
state_transition
process_noise
```

## Two channels

[`cross_psd_y`](@ref) and [`cross_psd_x`](@ref) average the cross-spectrum
of two simultaneous records of one device; [`cross_adev`](@ref) is the
matching cross Allan variance. [Two-channel measurements](../guide/cross.md)
shows when and how to use them.

```@docs
cross_psd_y
cross_psd_x
CrossSpectralResult
cross_adev
CrossStabilityResult
```

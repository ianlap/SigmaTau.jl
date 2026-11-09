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

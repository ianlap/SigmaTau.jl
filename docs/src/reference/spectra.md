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

# Spectra

One-sided Welch estimates of the frequency and phase spectra and the
single-sideband phase noise, returned as a [`SpectralResult`](@ref).
`S_x(f) = S_y(f) / (2πf)²` and `ℒ(f) = ½ (2π f_carrier)² S_x(f)`.
[Plotting](../guide/plotting.md#Spectral-densities) shows them with
power-law slope guides.

```@docs
Sy
Sx
L
```

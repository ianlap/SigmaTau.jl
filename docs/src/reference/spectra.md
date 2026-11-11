# Spectra and noise models

## Spectra

``S_y(f)``, the spectrum of the fractional frequency, is a one-sided Welch
estimate. ``S_x(f)`` of the time fluctuations and the phase noise ``ℒ(f)``
of a carrier are derived from it, through the response of the first
difference that takes phase to frequency. Each is returned as a
[`SpectralResult`](@ref). [Conventions](../explanation/conventions.md#Spectra)
gives their units and the relations between them, and
[Plotting and tables](../guide/plotting.md#Spectra) draws them.

```@docs
psd_y
psd_x
ssb_phase_noise
integrated_jitter
```

## Power-law model

The model ``S_y(f) = \sum_α h_α f^α`` links the spectra and the
deviations. [`powerlaw_fit`](@ref) fits the ``h_α`` to the deviation of a
record or to a spectrum, and [`powerlaw_sigma`](@ref) and
[`powerlaw_psd`](@ref) give the deviation and the spectra of a model.
[Power-law noise](../explanation/noise.md) has the conversions.

The spectrum can be a table of points, such as the phase noise on a data
sheet:

```@example spectra
using SigmaTau

# ℒ(f) in dBc/Hz of a 10 MHz oscillator, at offsets from 1 Hz to 10 kHz
f = powerlaw_fit([1, 10, 100, 1e3, 1e4], [-100, -130, -150, -158, -160];
                 quantity = :L, f_carrier = 10e6)
```

The Allan deviation that this model gives at 1 ms, 10 ms and 0.1 s, in a
record sampled every millisecond:

```@example spectra
powerlaw_sigma(f, [0.001, 0.01, 0.1]; tau0 = 1e-3)
```

The table says nothing of the noise below 1 Hz, which sets the deviation
at averaging times of a second and more, so the model is evaluated only
where the table supports it.

```@docs
powerlaw_fit
PowerLawFit
powerlaw_sigma
powerlaw_psd
```

## A model built by hand

[`noise_fit`](@ref) builds the model one noise type at a time. Each
component has the slope of its noise type and is fitted over a chosen
range of the ``τ`` grid to what the other components leave. It serves
where the noise types and their ranges are known, and where
[`powerlaw_fit`](@ref) spends a component on a bump of the curve.

```@example spectra
using Random

# white PM, white FM and random-walk FM
p = noise_gen(PhaseData, 8192, 1.0; rng = Xoshiro(5),
              sigma1 = Dict(2 => 2e-12, 0 => 1e-12, -2 => 3e-14))
m = noise_fit(mhdev(p))
add_component!(m, :WHPM, 1:2)      # white PM, fitted at the first two τ
add_component!(m, :WHFM, 5:8)
```

[`noise_residual`](@ref) returns the variance that the model leaves at
each ``τ``, and [`undo!`](@ref) takes the last component back.

```@docs
noise_fit
NoiseFit
add_component!
undo!
noise_residual
```

## Clock model

[`clock_model`](@ref) turns a power-law model into the
[`ClockModel`](@ref) on which the [Kalman filter](kalman.md) and the
[time scales](ensemble.md) run. White, random-walk and random-run FM drive
the phase, the frequency and the drift, and PM noise becomes the variance
of a reading. By default flicker FM has no state and is left out. With
`flicker = true` flicker FM and flicker PM each become a set of Markov
states. [The clock model](../explanation/clock_model.md) derives the
mapping.

```@example spectra
clock_model(m)
```

```@docs
ClockModel
clock_model
state_transition
process_noise
```

## Cross-spectra

[`cross_psd_y`](@ref) and [`cross_psd_x`](@ref) average the cross-spectrum
of two simultaneous records that share a device and nothing else, so that
the noise of the two channels averages toward zero.
[Two-channel measurements](../guide/cross.md) says when the method applies
and what it cannot remove. [`cross_adev`](@ref) is its counterpart in the
time domain.

```@docs
cross_psd_y
cross_psd_x
CrossSpectralResult
```

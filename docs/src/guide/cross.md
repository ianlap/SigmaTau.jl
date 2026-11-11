# Two-channel measurements

A measurement of a device against a reference holds the noise of both, and
that of the measuring channel. When the device is measured twice at the
same time, through two channels or against two references, the two records
share the device and nothing else. The average of their product keeps the
shared part, and the rest averages toward zero. The method is the
cross-spectrum in the frequency domain
([Baudiquez et al., 2020](@cite baudiquez-2020-cross-spectrum-statistics))
and the Allan covariance in the time domain
([Vernotte et al., 2016](@cite vernotte-2016-gcov)).

Two arrangements use it:

- **Two references.** The device is measured against two references at
  once, which gives the records `a = device − ref1` and
  `b = device − ref2`. The cross estimates show the device even where each
  reference is noisier than it.
- **Two channels.** One signal is split into the two inputs of an
  instrument. The cross estimates show the source below the noise of each
  channel. The same reference in both inputs shows the floor of the
  instrument.

## The cross Allan deviation

The example simulates a device with white FM noise against two references
whose noise is ten times larger in power:

```@example cross
using SigmaTau, Random

N = 65536
wfm(h, seed) = noise_gen(FrequencyData, N, 1.0; h = Dict(:WHFM => h),
                         rng = Xoshiro(seed)).y
device, ref1, ref2 = wfm(1e-24, 2), wfm(1e-23, 3), wfm(1e-23, 4)
a = FrequencyData(device .- ref1, 1.0)
b = FrequencyData(device .- ref2, 1.0)

x = cross_adev(a, b, [1, 4, 16, 64, 256])
[x.tau x.dev x.dev_a]
```

The second column is the cross Allan deviation and the third the Allan
deviation of channel `a` alone, which its reference sets. The device was
simulated at ``σ_y(1\,\mathrm{s}) = 7.1 × 10^{-13}``, falling as
``τ^{-1/2}``:

```@example cross
using Plots
plot(x.tau, [x.dev_a x.dev sqrt.(1e-24 ./ 2x.tau)]; xscale = :log10,
     yscale = :log10, marker = :circle, xlabel = "τ (s)", ylabel = "σ_y(τ)",
     label = ["channel a" "cross" "device, as simulated"])
```

[`cross_adev`](@ref) returns a [`CrossStabilityResult`](@ref). Its `covar`
is the cross Allan variance and `dev` the square root. `uncertainty` is
the standard uncertainty of `covar`, which is also the level that two
records with nothing in common would scatter about:

```@example cross
round.(x.uncertainty ./ x.covar; digits = 2)
```

At 64 s and 256 s the uncertainty is half the estimate, which is a quarter
of the deviation, and the point at 64 s lies that far below the simulated
line. Where the uncertainty exceeds the estimate, the estimate can be
negative: `covar` keeps the sign and `dev` is `NaN`. A longer record lowers
the uncertainty as the inverse square root of its length.

## The cross-spectrum

[`cross_psd_y`](@ref) and [`cross_psd_x`](@ref) do the same for the
spectrum. They average the product of the two channels' transforms over
the segments of Welch's method:

```@example cross
using Statistics
s = cross_psd_y(a, b; nperseg = 512)
band = 20:200                    # 0.04 Hz to 0.39 Hz
s.nseg, mean(s.psd[band]), mean(s.psd_a[band]), mean(s.floor[band])
```

Over the band the cross-spectrum reads the ``10^{-24}`` per hertz of the
device, where channel `a` alone reads ``1.1 × 10^{-23}``. `floor` is the
level that the unshared noise leaves after `nseg` averages. A bin whose
`psd` is not well above it is not resolved.

The floor falls as the inverse square root of the number of segments, so
the segment length trades resolution in frequency for depth. The default,
a quarter of the record, gives 7 segments and a floor above this device.
Give `nperseg` when the aim is to see below the channels.

## What the method cannot remove

- **Noise the channels share.** A common reference, a common supply or the
  temperature of the room appears in both records and is counted as the
  device.
- **Unequal timing.** The two records must be simultaneous, with the same
  sample interval and length and no gaps. [`fillgaps`](@ref) fills a
  record with gaps first.

With three devices compared in pairs, [`groslambert`](@ref) applies the
same covariance to give each one its own deviation, and
[`ncornered_hat`](@ref) does it from the three deviations. The tutorial
[Three clocks and the three-cornered hat](../tutorials/03_three_cornered_hat.md)
works through both.

## From the command line

```bash
sigmatau cross a.dat --with b.dat                  # cross Allan deviation, CSV
sigmatau cross a.dat --with b.dat --format json    # adds the cross-spectrum
```

`--with` names the second channel, which is read with the options of the
first ([Command line](cli.md)).

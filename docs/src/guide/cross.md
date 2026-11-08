# Two-channel measurements

A single measurement of a device against a reference sees the noise of both
plus the noise of the measuring channel. When the device is measured twice at
the same time, through two channels or against two references, the noise the
two records share is the device's, and everything else is independent between
them. Averaging the product of the two records keeps the shared part and
drives the rest toward zero. Two setups use this:

- **Two references.** Measure the device against references R₁ and R₂
  simultaneously: `a = DUT − R₁`, `b = DUT − R₂`. The cross estimates give the
  device's noise even where each reference is noisier than it.
- **One signal, two channels.** Split one source into two inputs of a
  dual-channel instrument. The cross estimates give the source's noise below
  the channels' own floor; feeding the same reference into both shows the
  instrument floor itself.

```@example cross
using SigmaTau, Random
N = 1 << 16
dut = noise_gen(FrequencyData, N, 1.0; h=Dict(:WHFM => 1e-24), rng=Xoshiro(2)).y
r1  = noise_gen(FrequencyData, N, 1.0; h=Dict(:WHFM => 1e-23), rng=Xoshiro(3)).y
r2  = noise_gen(FrequencyData, N, 1.0; h=Dict(:WHFM => 1e-23), rng=Xoshiro(4)).y
a, b = FrequencyData(dut .- r1, 1.0), FrequencyData(dut .- r2, 1.0)

x = cross_adev(a, b, [1, 4, 16, 64])
[x.tau x.dev x.uncertainty x.dev_a]
```

The device alone has σ_y(1 s) ≈ 7.1e-13 (white FM, h₀ = 1e-24); either
channel reads 2.4e-12 (`dev_a`), set by its reference. `cross_adev` recovers
7.1e-13 at τ = 1 s and follows the device down the τ⁻¹ᐟ² line. `uncertainty` is the 1σ
statistical uncertainty of `covar` (the cross variance) from the record
itself, and `covar` can come out negative where the shared noise is below it;
`dev` is `NaN` there.

```@example cross
s = cross_psd_y(a, b; nperseg=512)
(s.nseg, sum(s.psd[20:200]) / 181, sum(s.floor[20:200]) / 181)
```

`cross_psd_y` returns the real part of the averaged cross-spectrum (`psd`),
each channel's own spectrum (`psd_a`, `psd_b`) and the residual left by
unshared noise after `nseg` averages, `floor = √(S_aa S_bb / nseg)`. The
floor falls as the square root of the number of averages: four times the
data halves it.

From the command line, `sigmatau cross a.dat --with b.dat` writes the cross
Allan deviation as CSV, and `--format json` adds the cross-spectrum.

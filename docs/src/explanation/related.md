# Related software

Other programs compute the statistics of this package. This page describes
the ones it has been compared with, or that cover the same ground, in the
versions current in November 2026. The descriptions follow each program's
own documentation, and the comparisons state what was compared and how
closely the results agree.

## Stable32

Stable32 is the Windows program of Hamilton Technical Services, whose
W. J. Riley wrote the NIST handbook that this package follows
([Riley, 2008](@cite riley-2008-sp1065)). It was sold from 1982 until 2017
and then given to the IEEE UFFC Society for free distribution. The current
version is 1.62, of 2016 ([www.wriley.com](https://www.wriley.com/)).

Its manual ([Hamilton Technical Services, 2008](@cite hamilton-2008-stable32-manual))
lists the Allan and Hadamard deviations in their normal and overlapping
forms, the modified Allan and time deviations, the total, modified total,
time total and Hadamard total deviations, Thêo1, ThêoBR and ThêoH, MTIE
and TIE rms. Each run identifies the noise type and reports χ² confidence
limits. The program also plots and edits data, removes drift and outliers,
computes spectra and dynamic stability, and separates three clocks with
the three-cornered hat. It has functions this package lacks: frequency
jump detection, filtering and a histogram.

The cross-check suite compares this package with Stable32's result files
on its own test record, the 1000-point test suite of
[SP1065 §12.4](@cite riley-2008-sp1065), and on a flicker PM record from
its noise generator, with the options that select Stable32's conventions.
The deviations agree within ``5 × 10^{-5}``, one unit of the last of the
five digits Stable32 prints. On the test record the noise types are the
same at every averaging time compared, and the confidence limits agree
within 1 part in 1000, apart from the cases listed in the
[cross-check report](crosscheck.md).

This package departs from Stable32 where the literature or a simulation
supports a different choice, in the noise identification, the bias
factors and the handling of drift.
[Coming from Stable32](../guide/from_stable32.md) lists each difference
and the option that reproduces Stable32.

## allantools

allantools is a Python package by Anders E. E. Wallin and contributors
under the LGPL, version 3 or later
([Wallin et al., 2024](@cite wallin-2024-allantools)). Version 2024.6 has
the functions `adev`, `oadev`, `mdev`, `tdev`, `hdev`, `ohdev`, `pdev`,
`totdev`, `mtotdev`, `ttotdev`, `htotdev`, `theo1`, `mtie` and `tierms`,
the gap-tolerant `gradev` and the Groslambert codeviation `gcodev`. Each
returns the deviation with the number of terms and an error estimate of
the deviation over the square root of that number. The equivalent degrees
of freedom, the χ² interval and the lag-1 noise identification are
separate functions (`edf_greenhall`, `confidence_interval`,
`autocorr_noise_id`). The total deviations are returned without bias
correction. The package also has real-time estimators, noise generators
and classes for datasets and plots.

The cross-check suite runs allantools 2024.6 on seven power-law records
of 4096 points and on two reference records, and compares every function
above except `gcodev` over a grid of averaging factors. The deviations
agree within ``10^{-9}`` (``10^{-8}`` for the modified Allan and time
deviations), and so do the EDF of `edf_greenhall` and the χ² factors. The
lag-1 estimates of the noise exponent agree within ``10^{-6}``. The two
packages differ in the EDF they give the total deviations; the
[cross-check report](crosscheck.md) lists the cases. A
continuous-integration job installs allantools and AllanDeviations.jl and
repeats the comparison on every pull request.

[Coming from allantools](../guide/from_allantools.md) maps its calls onto
this package.

## AllanDeviations.jl

AllanDeviations.jl is a Julia package by Julien Kluge under the MIT
license ([Kluge, 2022](@cite kluge-2022-allandeviations)). Version 0.3.0
computes the Allan, modified Allan, Hadamard, time and total deviations
and MTIE (`allandev`, `mallandev`, `hadamarddev`, `timedev`, `totaldev`,
`mtie`), the Allan and Hadamard deviations with or without overlap. Each
returns the deviation, the number of terms and an error estimate of the
deviation over the square root of that number.

The cross-check suite compares all six on the seven power-law records.
The deviations agree within ``10^{-9}`` (``10^{-8}`` for the modified
Allan and time deviations), and MTIE is identical.

## SigmaTheta

SigmaTheta is a set of command-line programs in C by François Vernotte,
Université de Franche-Comté, under the CeCILL license
([Vernotte, 2025](@cite vernotte-2025-sigmatheta)). Its programs compute
the Allan, modified Allan, Hadamard and parabolic deviations, give the
Allan and modified Allan deviations 95 % and 68 % confidence intervals,
and find the dominant noise type at each ``τ`` from asymptotes fitted to
the deviation curve. Others remove linear drift and simulate noise
records. For three clocks they compute the Groslambert codeviation and
its confidence interval by the KLTG and KLTS methods
([Vernotte and Lantz, 2019](@cite vernotte-2019-tch-gcov-uncertainty);
[Lantz et al., 2019](@cite lantz-2019-klts)).

[`groslambert`](@ref) computes the KLTS interval with `interval=:klts`,
and SigmaTheta is its authors' implementation. On three simulated clocks
of 4096 points the Groslambert codeviations of the two programs agreed to
every printed digit, and given the same inputs the KLTS bounds agreed
within 2 %. The comparison was made once, against the version of April
2025, and is not part of the test suite. The suite does hold [`pdev`](@ref)
to nine values on the 1000-point test suite that the allantools tests give,
computed there with the `PDev` program of SigmaTheta 3.0 (October 2023).
They agree within ``10^{-11}``.

## ntpstats

ntpstats is a Python toolkit for validating and studying network time
synchronization by Thiago de Freitas, under the MIT license
([de Freitas, 2026](@cite de-freitas-2026-ntpstats)). Beside its tools for
NTP and PTP, version 3.9.0 has the Allan, modified Allan, time and
Hadamard deviations, the total, modified total, time total and Hadamard
total deviations, Thêo1, ThêoBR and ThêoH, MTIE and TIE rms, with χ²
intervals, lag-1 noise identification, a power-law fit, spectra, an
N-cornered hat and masks. Its documentation reports agreement with
Stable32 1.62.

This package has not been compared with ntpstats.

## Estimators by program

The estimators that each program documents, in the versions named above:

| Deviation | SigmaTau | Stable32 | allantools | AllanDeviations.jl | SigmaTheta | ntpstats |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| Allan | ● | ● | ● | ● | ● | ● |
| Modified Allan | ● | ● | ● | ● | ● | ● |
| Time | ● | ● | ● | ● | | ● |
| Hadamard | ● | ● | ● | ● | ● | ● |
| Modified Hadamard, Hadamard time | ● | | | | | |
| Parabolic | ● | | ● | | ● | |
| Total | ● | ● | ● | ● | | ● |
| Modified total, time total, Hadamard total | ● | ● | ● | | | ● |
| Modified Hadamard total | ● | | | | | |
| Thêo1 | ● | ● | ● | | | ● |
| ThêoBR, ThêoH | ● | ● | | | | ● |
| MTIE | ● | ● | ● | ● | | ● |
| TIE rms | ● | ● | ● | | | ● |
| Three-cornered hat or Groslambert covariance | ● | ● | ● | | ● | ● |

In the last row Stable32 has the classical three-cornered hat, and the
other programs have the Groslambert covariance as well.

The modified Hadamard deviation is defined by
[Bregni and Jmoda (2006)](@cite bregni-2006-mhvar) and described in
[SP1065 §5.2.10](@cite riley-2008-sp1065). The modified Hadamard total
deviation is constructed in this package, after an outline by
[Vernotte and Howe (2000)](@cite vernotte-2000-total-structure)
([MHTOTDEV bias and EDF](@ref)). None of the programs above computes
either, so their tests rest on identities, on direct evaluation of the
defining sums and, for the total, on the simulations described on that
page.

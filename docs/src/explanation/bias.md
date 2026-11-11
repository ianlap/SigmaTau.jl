# Bias correction

The total estimators and Thêo1 gain their degrees of freedom by reusing the
record, and they pay for it with a bias: the expected value of each differs
from that of the variance it stands in for, its *parent*, by a factor that
depends on the noise type. The package divides the factor out by default.
This page lists the factors and their sources, shows how well the corrected
estimators agree with their parents, and gives the options that reproduce
other conventions.

## The factors

The bias factor is the ratio of expected values,

```math
B = \frac{E[\text{total variance}]}{E[\text{parent variance}]} ,
```

and an estimator divides its deviation by ``\sqrt{B}``. The parents are the
Allan variance for [`totdev`](@ref), the modified Allan variance for
[`mtotdev`](@ref) and [`ttotdev`](@ref), the Hadamard variance for
[`htotdev`](@ref) and the modified Hadamard variance for
[`mhtotdev`](@ref). [`bias_correction`](@ref) returns ``B``:

```@example bias
using SigmaTau, Printf

five = [:WHPM, :FLPM, :WHFM, :FLFM, :RWFM]
seven = [five; :FWFM; :RRFM]
function row(name, kind, m, noises)
    B = bias_correction(noises, kind, fill(m, length(noises)), 1025)
    println(rpad(name, 20), join(map(b -> @sprintf("%7.3f", b), B)))
end
println(rpad("", 20), join(lpad.(seven, 7)))
row("totdev, τ = T/8", :totvar, 128, five)
row("totdev, τ = T/2", :totvar, 512, five)
row("mtotdev, ttotdev", :mtot, 64, five)
row("htotdev", :htot, 64, seven)
row("mhtotdev", :mhtot, 64, seven)
```

A value below 1 means that the uncorrected estimator reads low. The total
deviation reads low on flicker FM and random-walk FM and high on PM noise,
both in proportion to ``τ/T``. The modified and Hadamard totals read low by
a fixed factor from ``τ = 16τ_0``, and the Hadamard total reads high on PM
noise.

| Estimator | Noise | Factor | Source |
|:---|:---|:---|:---|
| `totdev` | flicker FM, random-walk FM | ``1 - a\,τ/T`` with ``a = 1/(3\ln 2)`` and ``3/4`` | [Greenhall et al. (1999), Eq. (25)](@cite greenhall-1999-totvar); [SP1065 §5.11](@cite riley-2008-sp1065) |
| `totdev` | white FM | 1 | the same |
| `totdev` | white PM, flicker PM | ``1 + (8/3)⌊m/2⌋/(N-2)``, exact; ``1 + (0.593 + 0.068\ln m)\,τ/T`` | this package: [Total estimators on PM noise](@ref) |
| `mtotdev`, `ttotdev` | white PM to random-walk FM | 0.995, 0.851, 0.771, 0.717, 0.679 | [Howe and Peppler (2001), Table 2](@cite howe-2001-total-estimators) |
| `htotdev` | white FM to random-run FM | the same five values | [Howe et al. (2005), Table I](@cite howe-2005-tothvar) |
| `htotdev` | white PM, flicker PM | 1.10 to 1.38 and 1.02 to 1.30, by ``m`` | this package: [Total estimators on PM noise](@ref) |
| `mhtotdev` | white PM to random-run FM | 0.78 to 0.90 | this package: [MHTOTDEV bias and EDF](@ref) |

!!! note "SP1065 Table 11"
    [SP1065 §5.12](@cite riley-2008-sp1065) lists modified total bias
    factors of 1.06 to 1.31 and describes them as the ratio of the expected
    modified total variance to the modified Allan variance. The estimator
    reads low, so that ratio is below 1, and those values are not used.

## Averaging factors below 16

The published factors of the modified and Hadamard totals are stated for
``τ ≥ 16τ_0``. Below that the expectation depends on the averaging factor,
and the package divides by its exact value at each ``m`` from 2 to 15,
computed as on the page [Total estimators on PM noise](@ref). For the
modified total:

```@example bias
ms = [1, 2, 3, 4, 5, 8, 16]
println(rpad("m", 6), join(lpad.(ms, 7)))
for n in five
    B = bias_correction(fill(n, length(ms)), :mtot, ms, 1025)
    println(rpad(n, 6), join(map(b -> @sprintf("%7.3f", b), B)))
end
```

The Hadamard total has the same rows two noise types down, from white FM to
random-run FM. White PM dips at odd ``m``, by 6 % at ``m = 3``, and white
FM is 6 % above its limit at ``m = 2``. The two flicker rows are the mean
of two flicker models, each within 2.4 % of it at ``m = 2`` and within 1 %
from ``m = 5``.

At ``m = 1`` nothing is approximate. The modified total variance is half
the modified Allan variance on any record, and the modified Hadamard total
is 124/225 of the modified Hadamard variance. The total and Hadamard total
equal the Allan and Hadamard variances term for term.

## Accuracy of the corrected estimators

In simulation (1024 points, 2000 records per noise type) with the noise
type stated, the corrected variance of each of the four total estimators
averages within 2 % of its parent up to ``τ = T/8`` and within 4 % at the
longest averaging times, on every noise type.

A default call identifies the noise type, and the factor follows the
label. With 64 averaged points or more the result is the same within 2 %,
and down to 8 averaged points within 5 %. With fewer the label is
uncertain ([Noise identification](noise_identification.md)), and a record
given the wrong label is divided by the wrong factor. Mean of the corrected
variance over the parent variance at ``τ = T/4``, with 4 averaged points
left:

| Estimator | white PM | flicker PM | white FM | flicker FM | random-walk FM | flicker-walk FM | random-run FM |
|:---|---:|---:|---:|---:|---:|---:|---:|
| `totdev` | 1.00 | 1.19 | 1.06 | 0.96 | 0.96 | | |
| `mtotdev` | 1.00 | 1.10 | 1.04 | 0.99 | 0.98 | | |
| `htotdev` | 1.00 | 1.25 | 1.00 | 0.99 | 1.01 | 1.01 | 1.00 |
| `mhtotdev` | 1.04 | 1.04 | 1.02 | 1.01 | 1.01 | 0.99 | 1.00 |

These are variances; the deviation is off by half as much. Flicker PM
stands out. With fewer than 6 averaged points it is labeled white FM in
most records, so the total and the Hadamard total keep the bias that they
have on PM noise, and the modified total is divided by the factor of white
FM. White FM reads high in the total and the modified total, for the
records labeled flicker or random-walk FM. The modified Hadamard total is
the least exposed, because its factors lie within 15 % of each other, and
its 4 % on white PM is there with the noise type stated as well. State the
noise type at the longest averaging times when it is known:

```julia
htotdev(p; noise = :RWFM)
```

## Thêo1

Thêo1 is unbiased against the Allan variance on white FM only
([Howe (2006), §4](@cite howe-2006-theoh)). [`theo1`](@ref) removes the
bias in the form called ThêoBR: it multiplies the Thêo1 variance by the
mean ratio of the Allan variance to Thêo1 at the same averaging time,
measured on the record itself at averaging times from ``9τ_0`` up to half
the record ([SP1065 §5.2.16, Eq. (33)](@cite riley-2008-sp1065)). No noise
type is needed. [`theoh`](@ref) uses the corrected value above a fifth of
the record and the Allan deviation below.

The method has three limits. The first two show in simulation on records
of one noise type (1024 points, 2000 records per noise type).

**It is one number.** On PM noise the bias varies with ``m`` throughout
([Howe (2006), Table 3](@cite howe-2006-theoh)): on white PM the
definition gives an expected value of ``H_{m/2}`` times the Allan variance
at the same ``τ``, where ``H_n`` is the ``n``th harmonic number, 2.45 at
``m = 12`` and 5.9 at ``m = 400``. One ratio cannot follow it. On these
records the corrected variance is 0.48 of the Allan variance at ``9τ_0``
on white PM and 0.55 on flicker PM, right near a fifth of the record, and
1.2 times the Allan variance at ``T/2`` on both. Do not use Thêo1 where PM
noise dominates.

**It is an estimate.** The last terms of the mean are ratios at averaging
times near half the record, where the Allan variance rests on a few
differences, and each counts as much as the first. On FM noise the ratio
scatters by 26 to 32 % in variance, at every record length tried from 256
to 16 384 points, and the corrected deviation carries that scatter at
every averaging time. Its rms error is 14 to 15 % at ``9τ_0``, where that
of the Allan deviation is 5 to 7 %. It is the more precise of the two only
from about a tenth of the record on: 22 to 34 % against 27 to 35 % at
``T/5``, and 40 to 60 % against 62 to 64 % at ``T/2``. The corrected
variance averages 6 to 9 % below the Allan variance between ``48τ_0`` and
``T/5``, and 8 to 14 % above it at ``T/2``. The interval does not include
the scatter of the ratio ([Confidence intervals](confidence.md)).

**It assumes one noise type.** One ratio is applied at every averaging
time. Its size depends on the noise type: in the fit of
[Howe (2006), Table 3](@cite howe-2006-theoh), reprinted as
[SP1065 Table 12](@cite riley-2008-sp1065), it tends to 2.70 on
random-walk FM and to 1.87 on flicker FM, is 1 on white FM, and falls with
``m`` toward 0.14 on flicker PM and 0.09 on white PM. Where one noise type
gives way to another within the record, the mean mixes the ratios of both,
and the corrected values at the long averaging times come out low. Mean of
the corrected variance over the true Allan variance at ``τ = 0.4T``, in
simulation (2048 points, 1000 records per case), for two noise types whose
Allan variances cross at ``τ_c``:

| Noise types | ``τ_c = 0.01T`` | ``0.03T`` | ``0.1T`` | ``0.3T`` |
|:---|:---:|:---:|:---:|:---:|
| White PM, then white FM | 0.90 (0.69) | 0.85 (0.59) | 0.84 (0.64) | 0.91 (0.87) |
| White PM, then flicker FM | 0.94 (0.69) | 0.82 (0.37) | 0.66 (0.25) | 0.78 (0.64) |
| White FM, then flicker FM | 0.97 (0.85) | 0.94 (0.78) | 0.92 (0.74) | 0.93 (0.79) |
| White FM, then random-walk FM | 1.04 (0.91) | 0.99 (0.71) | 0.88 (0.50) | 0.86 (0.62) |

The error grows toward the longest averaging times: where white PM gives
way to flicker FM at ``0.1T`` to ``0.3T``, the corrected variance is 0.46
to 0.55 of the true one at ``0.7T``. On such a record, take the long
averaging times from [`totdev`](@ref) or [`htotdev`](@ref), which reach
half the record with a factor for the noise type at each averaging time.

The values in parentheses are those of `correct_bias=:published`, which
averages the ratio up to a tenth of the record only. That is the range of
the sum in the papers that define the method
([Howe (2006), Eq. (6)](@cite howe-2006-theoh);
[Taylor and Howe (2010), Eq. (1)](@cite taylor-2010-fast-theobr)). SP1065
and [IEEE Std 1139-2022, Eq. (A.37)](@cite ieee1139-2022-definitions)
print half the record, the range that Stable32 uses and the default here.
The two fail in different places. The shorter range measures the ratio
where the first noise type holds and applies it where the second does, and
the corrected variance falls to a quarter of the true one in the table. On
a record of one FM noise type it is the better of the two: its ratio
scatters by 7 to 14 %, the corrected variance averages within about 2 % of
the Allan variance from ``48τ_0`` up, and from there on the rms error of
the deviation is below that of the Allan deviation, with 28 to 53 % at
``T/2``.

`correct_bias=:stable32` takes the default range on at most 3000 samples
of the record, as Stable32 does ([`theo1`](@ref) has the details). On
records of up to 3000 samples it agreed with the default within 0.3 % in
simulation. `correct_bias=false` returns Thêo1 as defined.

## Other conventions

The factors derived in this package are in no other program, so
`correct_bias` also takes the convention of each reference:

| `correct_bias` | `totdev` | `mtotdev`, `ttotdev` | `htotdev` | `mhtotdev` | `theo1` |
|:---|:---|:---|:---|:---|:---|
| `true`, the default | every noise type | exact below ``16τ_0`` | every noise type, exact below ``16τ_0`` | corrected | ThêoBR, ratio up to half the record |
| `:published` | FM only | the published factor at every ``τ`` | FM only, the published factor at every ``τ`` | corrected | ThêoBR, ratio up to a tenth of the record |
| `:stable32` | FM only | uncorrected | FM only, the published factor at every ``τ`` | corrected | ThêoBR as by default, on at most 3000 samples |
| `false` | uncorrected | uncorrected | uncorrected | uncorrected | uncorrected |

`:published` takes each correction as the paper that derived it prints
it, `:stable32` is what Stable32 applies to each estimator, and `false` is
what allantools returns. The modified Hadamard total exists in neither program and has no
published factor.

The setting changes the deviation and the interval that scales with it. It
does not change the noise identification or the EDF, so where this package
and Stable32 find a different noise type the factors still differ
([Coming from Stable32](../guide/from_stable32.md)). The command line names
the four `--bias exact`, `published`, `stable32` and `none`, the explorer
lists them under **Bias correction**, and Python takes `True`, `False`,
`"published"` and `"stable32"`.

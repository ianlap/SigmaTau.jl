# Dead time

A counter that measures for a time ``τ_0`` and then waits before it measures
again leaves *dead time* between its readings. The readings are frequency
averages over ``τ_0`` that start every ``T_0 > τ_0``, they cannot be summed
to phase, and their two-sample variance is not the Allan variance. For
power-law noise the two differ by a known factor, given by the bias
functions of [Barnes (1969)](@cite barnes-1969-bias-functions) and
[Barnes and Allan (1990)](@cite barnes-1990-dead-time). This page defines the
functions, corrects an Allan deviation with them, and says where the
correction holds.

## The bias functions

Write ``σ^2(N, T, τ)`` for the sample variance, with divisor ``N - 1``, of
``N`` frequency averages over ``τ`` that start every ``T``, and
``r = T/τ``. The Allan variance is the expected ``σ^2(2, τ, τ)``. On noise
whose Allan variance goes as ``τ^μ``:

| Function | Ratio of expected variances | Use |
|:---|:---|:---|
| [`bias_b1`](@ref)`(N, r, μ)` | ``σ^2(N, T, τ)`` over ``σ^2(2, T, τ)`` | ``N`` samples against two |
| [`bias_b2`](@ref)`(r, μ)` | ``σ^2(2, T, τ)`` over ``σ^2(2, τ, τ)`` | dead time against none |
| [`bias_b3`](@ref)`(N, M, r, μ)` | ``σ^2(N, M, T, τ)`` over ``σ^2(N, T, τ)`` | dead time spread among ``M`` averaged readings against the same dead time in one block |

In ``B_3`` the samples are averages of ``M`` readings, so that ``τ = Mτ_0``
and ``T = MT_0``. The exponent is ``μ = -2`` for white and flicker PM, −1
for white FM, 0 for flicker FM and 1 for random-walk FM, and the functions
take any ``-2 ≤ μ ≤ 2`` or the noise type.

The two-sample deviation of ``M``-reading averages, which is what
[`adev`](@ref) computes from the readings at averaging factor ``M``,
estimates ``\sqrt{B_2(r, μ)\,B_3(2, M, r, μ)}`` times ``σ_y(Mτ_0)``. For
100 s measured once an hour, ``r = 36``:

```@example deadtime
using SigmaTau, Printf

noises = [:WHPM, :WHFM, :FLFM, :RWFM]
println(rpad("M", 6), join(lpad.(noises, 8)))
for M in (1, 4, 16, 64)
    println(rpad(M, 6), join(map(n -> @sprintf("%8.3f", sqrt(bias_b2(36, n) * bias_b3(2, M, 36, n))), noises)))
end
```

White FM needs no correction. On white PM the uncorrected deviation reads
low at ``M = 1`` and high beyond it, and falls as ``τ^{-1/2}``, the slope of
white FM. On flicker FM it reads high by an amount that fades once
``M > r``, and on random-walk FM by a factor that settles there, where
``B_3`` tends to 2/3
([Barnes and Allan (1990), Eq. (8)](@cite barnes-1990-dead-time)).

## Correcting an Allan deviation

Give the readings to [`FrequencyData`](@ref) with the averaging time of one
reading as `tau0`, and the time not measured between readings as
`dead_time`. The noise type has to be stated. Here a flicker FM record is
read for 16 s out of every 64 s:

```@example deadtime
using Random

y = noise_gen(FrequencyData, 2^17, 1.0; sigma1 = Dict(:FLFM => 1e-12),
              rng = MersenneTwister(7)).y
readings = FrequencyData([sum(@view y[i:i+15]) / 16 for i in 1:64:length(y)-15], 16.0)

ms = [1, 2, 4, 8, 16]
uncorrected = adev(readings, ms; ci = false).dev
corrected = adev(readings, ms; dead_time = 48, noise = :FLFM, ci = false).dev
whole = adev(FrequencyData(y, 1.0), 16 .* ms; ci = false).dev
println("   τ (s)   uncorrected     corrected  whole record")
for i in eachindex(ms)
    @printf("%8d  %12.3e  %12.3e  %12.3e\n", 16ms[i], uncorrected[i], corrected[i], whole[i])
end
```

The result is ``σ_y(τ)`` at ``τ = m\,τ_0``, the time measured, not at the
``m\,T_0`` that the readings span. On the command line the same correction
is `sigmatau analyze --freq --devs adev --tau0 16 --dead-time 48 --noise
FLFM`, where `--tau0` has to be given: the time column of such a file holds
the period ``T_0``. In Python it is the `dead_time` keyword of `st.adev`.

## Where the correction holds

The functions are those of the two-sample variance on one power-law noise
type, and the correction is offered where that describes the data
([SP1065 §5.15](@cite riley-2008-sp1065)). Stable32 draws the same line:
it corrects the normal and overlapping Allan deviations when the noise type
is selected by hand
([Stable32 manual, "Dead Time"](@cite hamilton-2008-stable32-manual)).

- [`adev`](@ref) alone takes `dead_time`, in its overlapping and
  non-overlapped forms. No such functions exist for the modified, Hadamard,
  total or Thêo estimators, and [`deviation`](@ref) refuses the keyword for
  them. [`stability`](@ref), `sigmatau report` and the explorer do not apply
  it.
- The noise type is not identified. The slope of the uncorrected deviation
  is not that of the noise, as white PM shows above, and the lag-1 test
  assumes adjacent samples. State it with `noise`.
- The factor is that of the stated type at each ``τ``. A record with two
  noise types of comparable size is corrected for only one of them.
- The interval is the uncorrected one scaled by the same factor. Its EDF is
  that of a record of the same length without dead time, for which no
  published correction exists, so take it as approximate.
- A record with `NaN` gaps is refused, and so is phase data, which has no
  dead time.
- Flicker PM shares ``μ = -2`` with white PM, and the functions are those of
  white PM. In simulation (``2^{20}`` samples, readings of 8 or 16 samples,
  ``r`` = 1.5, 4 and 36) the corrected deviation on flicker PM read 3 to 6 %
  high at ``m = 1`` and 20 to 25 % low at ``m = 16``. On white PM, white
  FM, flicker FM and random-walk FM it was within 1 % of the deviation of
  the whole record at every ``m`` up to 16.

The functions describe ideal averages of a continuous process. A reading
that spans a single step of a simulated record is not one: with one sample
to a reading the same simulation gave a random-walk FM deviation 15 % low
at ``m = 1``.

## Agreement with the published tables

The functions reproduce 2935 of 2958 entries read from the tables of the
two notes to the four figures printed: ``B_2`` complete, ``B_1`` at
``r`` = 0.001, 0.01, 1.01 and 2, and ``B_3`` at ``r`` = 0.01, 0.03, 1.01, 2
and 32. Fifteen entries differ by one unit in the last figure. Eleven of
those are in the ``B_2`` table of the 1969 note, at places where the 1990
note prints the value computed here.

The remaining eight are in the white PM row at ``r = 0.01``. There an average
starts exactly where the one 100 places before it ends, the two share a
phase sample, and the published expressions give, for instance,
``B_3(2, 128, 0.01, -2) = 150``. The tables print the value that holds at
any other ratio, ``M = 128``, and
[Barnes and Allan (1990), Table 1](@cite barnes-1990-dead-time) lists
``B_1(N, r, -2) = 1`` and ``B_3(2, M, r, -2) = M`` without that exception.
The functions return the value of the expressions. The case needs
overlapping averages, ``r = 1/k``, and does not arise with dead time.

!!! note "Two printed values"
    [SP1065, Fig. 24](@cite riley-2008-sp1065) gives ``\sqrt{B_2} = 1.92``
    for flicker FM at ``r = 36``; the expression gives 1.915.
    [Barnes and Allan (1990), Table 1](@cite barnes-1990-dead-time) lists
    ``B_1(N, 1, μ) = 1`` for ``μ < 0`` and a sum over ``n^μ`` for
    ``μ > 0``. Those are the limits of ``B_1`` at large ``r``, not its
    values at ``r = 1``, which the line above it in the table gives.

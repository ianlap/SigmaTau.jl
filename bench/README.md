# Benchmarks

Three scripts time the same estimators, on the same samples and the same
averaging factors, in SigmaTau.jl, in
[allantools](https://github.com/aewallin/allantools) and in
[AllanDeviations.jl](https://github.com/JulienKluge/AllanDeviations.jl).

## Running them

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t auto bench/bench_sigmatau.jl 25000
pip install allantools
python bench/bench_allantools.py 25000
```

`bench_sigmatau.jl` draws a record of 25 000 phase samples, times
SigmaTau.jl on it, and writes the record to `bench/results/input_25000.txt`
and the timings to `bench/results/sigmatau_25000.json`.
`bench_allantools.py` times allantools on that record, writes
`bench/results/allantools_25000.json` and prints the tables of this page as
Markdown. The two JSON files in the repository hold the timings of the
tables below, and a run replaces them.

```bash
# AllanDeviations.jl, run between the two, adds its columns to the first table
julia --project=tools/validation/allandeviations -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/validation/allandeviations bench/bench_allandeviations.jl 25000

# a phase record of your own: the last column of each line, sampled at --tau0 seconds
julia --project=. -t auto bench/bench_sigmatau.jl --input phase.txt --tau0 1 --tag mine
python bench/bench_allantools.py --input phase.txt --tau0 1 --tag mine

# the tables again, from the saved timings
python bench/bench_allantools.py 25000 --report-only
```

A first argument other than `25000` sets the length of the record, and all
three scripts take it. `bench_sigmatau.jl` takes the largest averaging
factor as a second argument, 512 by default, and the other two read the
factors from its file.

## Method

Every library gets the same phase samples at τ₀ = 1 s and the averaging
factors 1, 2, 4, …, 512, and computes all ten in one call. `adev` and
`hdev` are the overlapping estimators, `oadev` and `ohdev` in allantools.

- **Record.** `randn(Xoshiro(20260601), N) .* 1e-9`: white PM, the same
  samples on every run.
- **Kernels.** SigmaTau.jl is called with `ci=false`, and the total
  deviations also with `correct_bias=false`. Each library then computes the
  deviations and nothing else: no noise identification, no bias correction
  and no confidence interval.
- **Default call.** `bench_sigmatau.jl` also times each estimator as it is
  called without keywords. That call identifies the noise at every
  averaging factor, corrects the bias of the total deviations, and computes
  the EDF and the confidence intervals. allantools has no call that does
  this work: the error it returns is the deviation divided by the square
  root of the number of terms. Its time is the same in both tables.
- **Threads.** SigmaTau.jl is started with every thread of the machine
  (`-t auto`), which `totdev`, `mtotdev` and `htotdev` use. Its other four
  kernels run on one thread, as allantools and AllanDeviations.jl do
  throughout.
- **Repeats.** A SigmaTau.jl kernel is called once, so that compilation is
  not counted, and then timed as the fastest of 20 calls, or of 3 for
  `mtotdev` and `htotdev`. The default call is timed once, after one call
  that compiles it. An allantools kernel is timed as the fastest of 5
  calls, and `mtotdev` and `htotdev`, which take minutes, once. An
  AllanDeviations.jl kernel is timed as the fastest of 20 calls, after one
  that compiles it. On a record given with `--input` every kernel is timed
  once.
- **Printing.** Times have three significant figures. A speed-up is the
  ratio of the two times as printed: a whole number from 10 up, and one
  decimal below.

## Results

Tables 1 and 2 are one run on a workstation with an AMD Ryzen 5 9600X
(6 cores, 12 threads) under Fedora, otherwise idle: SigmaTau.jl 0.5.0 at
commit `9c6acb7`, Julia 1.12.6 with 12 threads, Python 3.14.4, numpy 2.4.4
and allantools 2024.06. The record has 25 000 samples. The two scripts were
run three times, and the tables are those of the run with the smallest sum
of SigmaTau.jl times. They have no columns for AllanDeviations.jl, whose
script was added after the run.

**Table 1.** Kernels.

| Kernel | SigmaTau (12 threads) | allantools | Speed-up |
|---|---:|---:|---:|
| adev | 13.5 µs | 180 µs | 13× |
| mdev | 34.1 µs | 860 µs | 25× |
| hdev | 19.0 µs | 233 µs | 12× |
| tdev | 32.2 µs | 860 µs | 27× |
| totdev | 7.33 µs | 520 µs | 71× |
| mtotdev | 26.3 ms | 108 s | 4,106× |
| htotdev | 27.2 ms | 115 s | 4,228× |

**Table 2.** The default call of SigmaTau.jl, against the allantools times
of Table 1.

| Kernel | SigmaTau kernel-only | SigmaTau default call | allantools | Speed-up (allantools ÷ default call) |
|---|---:|---:|---:|---:|
| adev | 13.5 µs | 66.8 µs | 180 µs | 2.7× |
| mdev | 34.1 µs | 125 µs | 860 µs | 6.9× |
| hdev | 19.0 µs | 110 µs | 233 µs | 2.1× |
| tdev | 32.2 µs | 87.2 µs | 860 µs | 9.9× |
| totdev | 7.33 µs | 64.5 µs | 520 µs | 8.1× |
| mtotdev | 26.3 ms | 26.8 ms | 108 s | 4,030× |
| htotdev | 27.2 ms | 29.3 ms | 115 s | 3,925× |

For the first five estimators the noise identification and the intervals
take 50 to 90 µs, which is more than the kernel takes, and the ratio falls
to between 2 and 10. For `mtotdev` and `htotdev` the kernel is nearly the
whole call, and the ratio stays near 4000.

## Reading the tables

**They are for version 0.5.0.** The tables are not timed again at each
release. Noise identification does more work since that version, and the
default call of the first five estimators takes about twice as long: about
0.4 ms on `main` at `7674356`, against about 0.2 ms for version 0.5.0, each
the fastest of 50 calls on one 4-core machine with Julia 1.11. With the
allantools times unchanged, the first five ratios of Table 2 would be about
half of those printed.

**The ratios are those of one machine.** `totdev`, `mtotdev` and `htotdev`
use every thread in SigmaTau.jl and one in allantools, so their ratios
change with the number of cores. All the ratios change with the processor
and with the version of Julia.

**The modified totals.** Both libraries do work in proportion to N·m at
each averaging factor of `mtotdev` and `htotdev`. allantools runs the loop
over the subsequences, and the loop over the terms of each, in Python.
SigmaTau.jl compiles both loops, computes half of the terms of a
subsequence and takes the other half from the mirror symmetry of its
extension, and divides the subsequences among the threads.

Stable32 is compared for its results, in the
[validation](https://ianlap.github.io/SigmaTau.jl/dev/explanation/validation/)
of the package, and is not timed.

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
(6 cores, 12 threads) under Fedora 44 (Linux 6.19.14): SigmaTau.jl 0.6.0 at
commit `6565d93`, Julia 1.12.6 with 12 threads, Python 3.14.4, numpy 2.5.4,
allantools 2024.06 and AllanDeviations.jl 0.3.0. It was idle apart from one
desktop process using about a third of one core. The record has 25 000
samples. The three scripts were run three times, and the tables
are those of the run with the smallest sum of SigmaTau.jl times.

**Table 1.** Kernels.

| Kernel | SigmaTau (12 threads) | allantools | Speed-up | AllanDeviations.jl | Speed-up |
|---|---:|---:|---:|---:|---:|
| adev | 14.9 µs | 177 µs | 12× | 91.4 µs | 6.1× |
| mdev | 34.9 µs | 863 µs | 25× | 114 µs | 3.3× |
| hdev | 20.3 µs | 234 µs | 12× | 95.0 µs | 4.7× |
| tdev | 36.2 µs | 857 µs | 24× | 114 µs | 3.1× |
| totdev | 7.36 µs | 530 µs | 72× | 129 µs | 18× |
| mtotdev | 15.9 ms | 108 s | 6,792× | — | — |
| htotdev | 15.8 ms | 115 s | 7,278× | — | — |

**Table 2.** The default call of SigmaTau.jl, against the allantools times
of Table 1.

| Kernel | SigmaTau kernel-only | SigmaTau default call | allantools | Speed-up (allantools ÷ default call) |
|---|---:|---:|---:|---:|
| adev | 14.9 µs | 159 µs | 177 µs | 1.1× |
| mdev | 34.9 µs | 230 µs | 863 µs | 3.8× |
| hdev | 20.3 µs | 168 µs | 234 µs | 1.4× |
| tdev | 36.2 µs | 248 µs | 857 µs | 3.5× |
| totdev | 7.36 µs | 234 µs | 530 µs | 2.3× |
| mtotdev | 15.9 ms | 17.6 ms | 108 s | 6,136× |
| htotdev | 15.8 ms | 16.0 ms | 115 s | 7,188× |

For the first five estimators the noise identification and the intervals
take 140 to 230 µs, which is more than the kernel takes, and the ratio falls
to between 1.1 and 3.8. For `mtotdev` and `htotdev` the kernel is nearly the
whole call, and the ratio stays between 6000 and 7200.

**Thêo1.** On allantools' public phase record of a Cs 5071A against a
hydrogen maser (`tests/Cs5071A/5071A_phase.txt.gz` in the allantools
repository, downloaded and not stored here), N = 556 990 at τ₀ = 1 s,
`theo1` at the octave factors up to m = 524 288 takes 1.44 s as a kernel
(`ci=false, correct_bias=false`) and 14.0 s as the default call, each the
fastest of three calls on the same machine and commit with 12 threads. The
default call includes the ThêoBR bias ratio, the mean over the ⌊N/6 − 3⌋
ladder of AVAR/Thêo1 terms out to half the record, with Thêo1 computed at 48
factors along it and interpolated between them, as well as the noise
identification and the intervals. allantools is not timed here.

## Reading the tables

**They are for one commit.** The tables are not timed again at each
release, and they hold for `main` at `6565d93` only.

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

# Benchmarks

Two scripts time the same deviation kernels, on the same samples and the same
τ grid, in SigmaTau and in [allantools](https://github.com/aewallin/allantools).

```bash
# synthetic record (default)
julia --project=. -t auto bench/bench_sigmatau.jl 25000      # writes results/input_25000.txt, sigmatau_25000.json
pip install allantools
python bench/bench_allantools.py 25000                       # writes allantools_25000.json, prints the tables

# a phase record from a file (last column = phase, sampled at --tau0 seconds)
julia --project=. -t auto bench/bench_sigmatau.jl --input path/to/phase.txt --tau0 1 --tag measured
python bench/bench_allantools.py --input path/to/phase.txt --tau0 1 --tag measured

# reprint the tables from saved JSON without timing
python bench/bench_allantools.py 25000 --report-only
```

## Method

Both libraries get the same phase samples, τ₀ = 1 s and the octave grid
`m = 1, 2, 4, …, 512`, with one call per kernel over the whole grid.
SigmaTau is called with `ci=false` and raw (uncorrected) totals, so both sides
compute the same statistic. Every kernel is compiled (Julia) or imported
(Python) and called once before timing. SigmaTau runs with `-t auto`.
allantools is single-threaded numpy, since it has no parallel path.

- **Synthetic record**: N = 25 000 white-PM phase samples from one fixed seed
  (`Xoshiro(20260601)`); this is the same record on every run. Each SigmaTau
  kernel is timed as the minimum of 20 calls (3 for mtotdev and htotdev).
  Each allantools kernel is timed as the minimum of 5 calls (a single call for
  mtotdev and htotdev, which take about two minutes each). The whole pair of
  scripts was run three times. The table reports the run with the smallest
  summed SigmaTau time, as one whole run with no kernels mixed across runs.
- **Measured record**: a 406 763-sample phase record at τ₀ = 1 s (not in the
  repository), timed with `--input`, one timed call per kernel in each library.
- **End-to-end default call**: the same synthetic run also times SigmaTau's
  default call once per kernel (`ci=true`, bias-corrected totals: noise
  identification, EDF and χ² intervals included), after one compiling call.
  allantools has no counterpart, so its column is the same kernel time as in
  Table 1.

Times are printed to three significant figures. Each speed-up is the ratio of
the two times printed on its row, given as a whole number from 10× up and to
one decimal below 10×. `bench_allantools.py` produces the tables below
verbatim. The JSON files in `results/` are the raw timings behind Tables 1 to 3.

## Recorded results

Tables 1 to 3 come from a 6-core / 12-thread AMD Ryzen 5 9600X workstation
(Fedora, Linux 6.19), SigmaTau 0.5.0, Julia 1.12.6 with `-t auto` = 12
threads, Python 3.14.4, numpy 2.4.4 and allantools 2024.06. The machine
was otherwise idle. They were timed at commit 36379ab (`src/` identical
to `main` at 23634b0), before PR #57 added a detrending pass to the
MDEV/MHDEV prefix sums.

**Table 1.** Synthetic record, N = 25 000, one fixed-seed record. Each time is
the minimum of repeated calls, from the best of three runs (see Method).

| Kernel | SigmaTau (12 threads) | allantools | Speed-up |
|---|---:|---:|---:|
| adev | 13.3 µs | 179 µs | 13× |
| mdev | 31.1 µs | 858 µs | 28× |
| hdev | 18.6 µs | 234 µs | 13× |
| tdev | 28.4 µs | 862 µs | 30× |
| totdev | 7.04 µs | 511 µs | 73× |
| mtotdev | 27.1 ms | 107 s | 3,948× |
| htotdev | 27.6 ms | 114 s | 4,130× |

**Table 2.** Measured phase record, N = 406 763, one timed call per kernel.

| Kernel | SigmaTau (12 threads) | allantools | Speed-up |
|---|---:|---:|---:|
| adev | 266 µs | 8.27 ms | 31× |
| mdev | 825 µs | 20.5 ms | 25× |
| hdev | 341 µs | 14.2 ms | 42× |
| tdev | 492 µs | 20.3 ms | 41× |
| totdev | 155 µs | 12.3 ms | 79× |
| mtotdev | 444 ms | 1,814 s | 4,086× |
| htotdev | 457 ms | 1,932 s | 4,228× |

**Table 3.** End-to-end default call, synthetic record, N = 25 000, same run
as Table 1. One timed call per kernel. The speed-up divides the allantools
time by the SigmaTau default-call time.

| Kernel | SigmaTau kernel-only | SigmaTau default call | allantools | Speed-up (allantools ÷ default call) |
|---|---:|---:|---:|---:|
| adev | 13.3 µs | 62.6 µs | 179 µs | 2.9× |
| mdev | 31.1 µs | 115 µs | 858 µs | 7.5× |
| hdev | 18.6 µs | 111 µs | 234 µs | 2.1× |
| tdev | 28.4 µs | 83.2 µs | 862 µs | 10× |
| totdev | 7.04 µs | 64.8 µs | 511 µs | 7.9× |
| mtotdev | 27.1 ms | 31.1 ms | 107 s | 3,441× |
| htotdev | 27.6 ms | 28.4 ms | 114 s | 4,014× |

The default call adds per-τ noise identification, EDF and χ² intervals (and
the bias correction for the totals), which the kernel-only ratios in Tables 1
and 2 leave out. For the five fast kernels that work costs about 50 to 90 µs
per call at this N. It dominates the default-call time of the five fast kernels and cuts their ratio
to 2 to 10×. For mtotdev and htotdev the kernel dominates and the
end-to-end ratio stays in the thousands. allantools returns only a simple
σ/√n error estimate, with no noise identification or χ² intervals, so Table 3
compares SigmaTau doing more work against allantools doing less.

### Earlier re-run: 4-core cloud container

This is an independent run of the earlier version of these scripts on a
4-core cloud container. It used Julia 1.12.7 with 4 threads, Python 3.13 and
allantools 2024.06, with N = 25 000 and the minimum of repeats as in Table 1.
The CPU model and numpy version were not recorded and no JSON was kept. The
speed-ups were recomputed from the printed times under the rule above.

| Kernel  | SigmaTau (4 threads) | allantools | Speed-up |
|---------|---------------------:|-----------:|---------:|
| adev    | 44 µs    | 480 µs     | 11×      |
| mdev    | 102 µs   | 1.58 ms    | 15×      |
| hdev    | 63 µs    | 709 µs     | 11×      |
| tdev    | 183 µs   | 1.60 ms    | 8.7×     |
| totdev  | 28 µs    | 2.05 ms    | 73×      |
| mtotdev | 102 ms   | 167 s      | 1,637×   |
| htotdev | 96 ms    | 179 s      | 1,865×   |

The modified-total gap is algorithmic as well as compiled-versus-interpreted.
allantools evaluates each subsegment's reflection in Python-level loops.
SigmaTau streams the extension through a per-thread buffer and splits the
subsegment sum across threads. Absolute times depend on the machine, while
the ratios are the portable result. Stable32 is validated for numerical
agreement but is not timed here.

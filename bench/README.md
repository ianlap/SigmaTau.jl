# Benchmarks

Two scripts time the same deviation kernels, on the same samples and the same
τ grid, in SigmaTau and in [allantools](https://github.com/aewallin/allantools).

```bash
julia --project=. -t auto bench/bench_sigmatau.jl 25000      # writes results/input_25000.txt
pip install allantools
python bench/bench_allantools.py 25000                       # prints the comparison table
```

Method: white-PM phase, τ₀ = 1 s, octave grid `m = 1 … 512`, one call per
kernel over the whole grid, `ci=false` and raw (uncorrected) totals so both
libraries compute the same statistic. Each kernel is compiled or imported
before timing and the minimum of several runs is reported. SigmaTau runs
with `-t auto`; allantools is single-threaded numpy, as it has no parallel
path.

## Recorded results

12-core Linux workstation, Julia 1.12 / Python 3.14, allantools 2024.06.

Synthetic record, N = 25 000, mean of 30 realizations:

| Kernel  | SigmaTau | allantools | Speedup |
|---------|---------:|-----------:|--------:|
| adev    | 15 µs    | 393 µs     | 20×     |
| mdev    | 56 µs    | 1.13 ms    | 20×     |
| hdev    | 25 µs    | 341 µs     | 14×     |
| tdev    | 51 µs    | 991 µs     | 19×     |
| totdev  | 254 µs   | 684 µs     | 2.7×    |
| mtotdev | 30.3 ms  | 107.3 s    | 3,543×  |
| htotdev | 31.5 ms  | 114.3 s    | 3,632×  |

Measured phase record, N = 406 763, single call per kernel:

| Kernel  | SigmaTau | allantools | Speedup |
|---------|---------:|-----------:|--------:|
| adev    | <1 ms    | 9 ms       | 23×     |
| mdev    | 1 ms     | 21 ms      | 28×     |
| hdev    | <1 ms    | 15 ms      | 42×     |
| tdev    | 1 ms     | 21 ms      | 30×     |
| totdev  | 5 ms     | 16 ms      | 3.2×    |
| mtotdev | 0.469 s  | 1 803.6 s  | 3,842×  |
| htotdev | 0.455 s  | 1 928.0 s  | 4,233×  |

Independent re-run of the same scripts on a 4-core cloud container (Julia
1.12.7, Python 3.13, allantools 2024.06, N = 25 000, minimum of repeats):

| Kernel  | SigmaTau (4 threads) | allantools | Speedup |
|---------|---------------------:|-----------:|--------:|
| adev    | 44 µs    | 480 µs     | 11×     |
| mdev    | 102 µs   | 1.58 ms    | 15×     |
| hdev    | 63 µs    | 709 µs     | 11×     |
| tdev    | 183 µs   | 1.60 ms    | 9×      |
| totdev  | 28 µs    | 2.05 ms    | 74×     |
| mtotdev | 102 ms   | 167 s      | 1,641×  |
| htotdev | 96 ms    | 179 s      | 1,860×  |

The modified-total gap is algorithmic as well as compiled-versus-interpreted:
allantools evaluates each subsegment's reflection in Python-level loops, while
SigmaTau streams the extension through a per-thread buffer and splits the
subsegment sum across threads. Absolute times depend on the machine; the
ratios are the portable result. Stable32 is validated for numerical agreement
but is not timed here.

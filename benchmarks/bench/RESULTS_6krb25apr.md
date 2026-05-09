# Bench results — `6krb25apr.txt` (m capped at 512)

Run on remote workstation (12-core Linux, Julia 1.12.6 / Python 3.14.4, allantools 2024.06).
Repo at `a4a93e0` (post total-family perf work).

Both kernels called **one-shot** (single fn call passing all 10 taus). No allantools
internals were touched — we use the public `at.<kernel>(x, rate, data_type, taus)` API.

## Comparison table

```
=== 6krb25apr.txt (N=406763, τ₀=1 s, 10 m values 1..512) ===
SigmaTau threads: 12   |   allantools threads: 1
Memory: SigmaTau = Julia @timed .bytes (cumulative alloc);
        allantools = peak ΔRSS via psutil polling. Not directly comparable.

  kernel     SigmaTau (s)  allantools (s)  speedup    SigmaTau alloc   allantools ΔRSS
  ------------------------------------------------------------------------------------
  adev            0.000           0.009     23.2×       0.0 MiB           2.1 MiB
  mdev            0.001           0.021     27.6×      11.2 MiB           0.0 MiB
  hdev            0.000           0.015     41.7×       0.0 MiB           3.1 MiB
  tdev            0.001           0.021     29.6×      11.2 MiB           0.0 MiB
  totdev          0.005           0.016      3.2×      55.9 MiB           6.2 MiB
  mtotdev         0.469        1803.561  3,841.9×      47.0 MiB           0.1 MiB
  htotdev         0.455        1927.985  4,233.4×      50.1 MiB           6.1 MiB
  ------------------------------------------------------------------------------------
  total wall      0.932        3731.627  4,003.5×
```

## Headline

- **mtotdev: 0.469 s vs 1804 s → 3,842× faster.**
- **htotdev: 0.455 s vs 1928 s → 4,233× faster.**
- For totdev specifically the speedup is much smaller (3.2×) — that kernel is the cheap one of the total family on both sides; the perf work paid off mostly in the *modified*-total kernels.

## Notes on the numbers

- **SigmaTau small-kernel timings (adev/mdev/hdev/tdev) are at-or-below `@elapsed` resolution** (printf-rounded to 0.000–0.001 s; raw values from JSON are 13–810 µs). For quotable per-kernel numbers on these, run with `BenchmarkTools.@btime` from the Mac side. The speedup column rounds to 23–42× but at this granularity the *ratio* is the only meaningful figure.
- **allantools peak ΔRSS is often near zero** for the heavy kernels (mtotdev 0.1 MiB, mdev/tdev 0 MiB) because allantools writes into a small fixed accumulator and the input array dominates baseline RSS. SigmaTau's "alloc" is cumulative `@timed .bytes`, which counts every transient allocation including reused buffers; not directly comparable to peak ΔRSS.
- **Why m_max=512?** Per the prior per-m profile, allantools `mtotdev` follows `t(m) ≈ 3.385 + 1.615·m s`. m=512 alone is ~890 s; m=1024 doubles to ~1780 s, and m=131072 alone would be ~58 hours single-threaded. m=512 is the cap that keeps total allantools wall under ~1 hour while still exercising the full algorithm shape.

## Bench machinery (left on remote, not pushed)

- `benchmarks/bench/results_sigmatau.json` — full SigmaTau JSON (per-kernel time + alloc + gctime, run meta)
- `benchmarks/bench/results_allantools.json` — allantools JSON (per-kernel time + peak ΔRSS, run meta)
- `benchmarks/bench/render_table.py` — merges the two JSONs into the table above; supports `--kernels k1,k2,...` filter
- `benchmarks/bench/Project.toml` — bench env (path-source dev'd `SigmaTauStability` + `SigmaTauBase` + `JSON`); created because root `Project.toml`'s `[workspace] members = [...]` block isn't accepted by Julia 1.12.6 Pkg
- `benchmarks/bench/bench_sigmatau.jl` — added `@timed` per-kernel memory capture, JSON output, CLI entry, optional 3rd arg `m_max`
- `benchmarks/bench/bench_allantools.py` — added `RSSPeakSampler` (25 ms-poll psutil), one-shot mode (default) with optional `--per-m` diagnostic mode, `--m-max` flag, JSON output
- `benchmarks/bench/bench_*.log` — raw run logs
- `.bench-venv/` — Python venv (allantools 2024.06, numpy 2.4.4, psutil 7.2.2)

## Reproduce

```bash
# SigmaTau (~1 s wall + warmup, 12 threads, m capped at 512)
julia --project=benchmarks/bench -t auto \
  benchmarks/bench/bench_sigmatau.jl \
  reference/clock_data/6krb25apr.txt \
  benchmarks/bench/results_sigmatau.json \
  512

# allantools (~62 min wall, single-thread, m capped at 512)
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .bench-venv/bin/python benchmarks/bench/bench_allantools.py \
  reference/clock_data/6krb25apr.txt \
  --out benchmarks/bench/results_allantools.json \
  --m-max 512

# Combined table (filter to the 7 kernels both implement)
.bench-venv/bin/python benchmarks/bench/render_table.py \
  --kernels adev,mdev,hdev,tdev,totdev,mtotdev,htotdev \
  benchmarks/bench/results_sigmatau.json \
  benchmarks/bench/results_allantools.json
```

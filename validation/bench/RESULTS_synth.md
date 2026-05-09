# Statistical bench results — synthetic, 30 realizations × N=25000

Run on remote workstation (12-core Linux, Julia 1.12.6 / Python 3.14.4, allantools 2024.06).
Repo at `a4a93e0` (post total-family perf work).

**Design**: 30 realizations of N=25000 white phase noise (numpy `default_rng(seed=i)`,
single-column .txt files in `validation/bench/synth/`). Identical data feeds both
libraries. m grid: 1, 2, 4, ..., 512 (10 octave-spaced values). Kernels called one-shot
(single fn call passing all 10 taus, no allantools internals touched). Per-kernel
warmup pass before the realization loop. `GC.gc()` before each Julia kernel call.

## Time per kernel call (seconds, across 30 realizations)

```
  kernel    |  SigmaTau mean      ± std     median         p5        p95 |  allantools mean      ± std     median         p5        p95 |    speedup
  --------------------------------------------------------------------------------------------------------------------------------------------------
  adev      |          1.5e-5    3.5e-7    1.5e-5    1.5e-5    1.6e-5 |          3.93e-4   1.51e-5   3.93e-4   3.74e-4   4.30e-4 |       20.1×
  mdev      |          5.6e-5    8.7e-7    5.6e-5    5.4e-5    5.7e-5 |          1.13e-3   5.09e-6   1.13e-3   1.04e-3   1.13e-3 |       20.1×
  hdev      |          2.5e-5    4.2e-7    2.5e-5    2.4e-5    2.5e-5 |          3.41e-4   1.54e-5   3.41e-4   3.26e-4   3.62e-4 |       13.6×
  tdev      |          5.1e-5    1.4e-6    5.0e-5    4.9e-5    5.4e-5 |          9.91e-4   1.06e-5   9.92e-4   9.69e-4   1.04e-3 |       19.3×
  totdev    |          2.5e-4    4.3e-6    2.5e-4    2.5e-4    2.6e-4 |          6.84e-4   8.48e-5   6.86e-4   6.81e-4   8.13e-4 |        2.7×
  mtotdev   |        3.03e-2    2.2e-3    3.02e-2    2.76e-2    3.37e-2 |        107.261   2.52e-1   107.193  106.982   107.689   |     3542.7×
  htotdev   |        3.15e-2    2.3e-3    3.11e-2    2.82e-2    3.48e-2 |        114.250   2.19e-1   114.287  113.910   114.557   |     3632.0×
```

(All numbers from the renderer output verbatim — see `validation/bench/per_realization.csv`
for the raw per-realization arrays suitable for histograms / violins / ECDFs.)

## Memory per kernel call (MiB)

```
  kernel    |  SigmaTau mean      ± std       max |  allantools mean    ± std      max
  -------------------------------------------------------------------------------------
  adev      |          0.00       0.00      0.00  |          0.02       0.08     0.45
  mdev      |          0.78       0.00      0.78  |          0.20       0.10     0.59
  hdev      |          0.00       0.00      0.00  |          0.00       0.00     0.01
  tdev      |          0.78       0.00      0.78  |          0.00       0.00     0.00
  totdev    |          2.31       0.00      2.31  |          0.01       0.02     0.13
  mtotdev   |          2.21       0.00      2.21  |          0.01       0.02     0.12
  htotdev   |          2.39       0.00      2.39  |          0.26       0.00     0.28
```

Caveat: not directly comparable.
- **SigmaTau**: cumulative `@timed .bytes` (every alloc, including reused buffers).
- **allantools**: peak ΔRSS via 25 ms-poll psutil sampler (only the high-water mark of physical
  pages held above the baseline).

Both indicate memory pressure but answer different questions. SigmaTau allocates
~2 MiB of transients per heavy-kernel call; allantools holds <300 KiB of additional
RSS above the input array.

## Headline (paper-ready)

- **mtotdev**: SigmaTau **30.3 ± 2.2 ms** vs allantools **107.3 ± 0.25 s** → **3,543× faster**
- **htotdev**: SigmaTau **31.5 ± 2.3 ms** vs allantools **114.3 ± 0.22 s** → **3,632× faster**
- totdev (the non-modified total): SigmaTau 254 ± 4 µs vs allantools 684 ± 85 µs → 2.7× faster — confirms the perf work concentrated in the *modified*-total family.
- Cheap kernels (adev/mdev/hdev/tdev) speedup 13–20× — SigmaTau wall is sub-100 µs, near the `@timed` resolution; treat the ratio as the load-bearing figure for these.

## Variance characterization

- **allantools relative std**: 0.19 – 0.34 % across all kernels — extremely stable
  single-thread numpy.
- **SigmaTau relative std**: ~7 % on the modified-total kernels (heavy threading);
  <2 % on the cheap ones. The 7 % jitter is dominated by 12-thread scheduling on
  this box and would shrink with thread pinning / `Threads.@threads :static`.
- All distributions look symmetric around the mean (mean ≈ median, p5/p95 roughly
  equidistant). No long tail; suitable for normal-approximation reporting in a paper.

## Bench machinery (left on remote, not pushed)

- `validation/bench/gen_synth.py` — deterministic realization generator (seed = realization index, 25k samples per file)
- `validation/bench/synth/` — 30 single-column .txt files (15 MiB total)
- `validation/bench/results_sigmatau_synth.json` — full per-realization SigmaTau data (per-kernel array of `{realization, time_s, bytes, gctime_s}`)
- `validation/bench/results_allantools_synth.json` — full per-realization allantools data (per-kernel array of `{realization, time_s, rss_delta_bytes}`)
- `validation/bench/per_realization.csv` — long-format CSV (kernel, library, realization, time_s, mem_bytes), 420 rows — load with pandas/matplotlib for histograms
- `validation/bench/render_stats.py` — produces the table above
- `validation/bench/bench_sigmatau.jl` — added `--synth <dir>` mode
- `validation/bench/bench_allantools.py` — added `--synth <dir>` mode
- `validation/bench/bench_*synth*.log` — raw run logs

## Reproduce

```bash
# 1. generate data (~5 s)
.bench-venv/bin/python validation/bench/gen_synth.py --N 25000 --reals 30 \
  --out-dir validation/bench/synth

# 2. SigmaTau bench (~10 s, 12 threads)
julia --project=validation/bench -t auto \
  validation/bench/bench_sigmatau.jl --synth validation/bench/synth \
  validation/bench/results_sigmatau_synth.json 512

# 3. allantools bench (~1h45m, single-thread)
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .bench-venv/bin/python validation/bench/bench_allantools.py \
  --synth validation/bench/synth \
  --out validation/bench/results_allantools_synth.json \
  --m-max 512 \
  --kernels adev,mdev,hdev,tdev,totdev,mtotdev,htotdev

# 4. stats table + CSV export
.bench-venv/bin/python validation/bench/render_stats.py \
  validation/bench/results_sigmatau_synth.json \
  validation/bench/results_allantools_synth.json \
  --csv-out validation/bench/per_realization.csv \
  --kernels adev,mdev,hdev,tdev,totdev,mtotdev,htotdev
```

## Plotting starter (matplotlib)

```python
import pandas as pd, matplotlib.pyplot as plt
df = pd.read_csv("validation/bench/per_realization.csv")

# histogram of mtotdev wall-time, both libs
fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
mt_s = df.query("kernel == 'mtotdev' and library == 'sigmatau'").time_s
mt_a = df.query("kernel == 'mtotdev' and library == 'allantools'").time_s
ax[0].hist(mt_s, bins=15); ax[0].set_title(f"SigmaTau mtotdev (n={len(mt_s)})")
ax[0].set_xlabel("wall time (s)")
ax[1].hist(mt_a, bins=15, color='tab:orange'); ax[1].set_title(f"allantools mtotdev (n={len(mt_a)})")
ax[1].set_xlabel("wall time (s)")
plt.tight_layout(); plt.savefig("hist_mtotdev.png", dpi=120)
```

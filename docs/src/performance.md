# Performance

This page records timings against
[allantools](https://github.com/aewallin/allantools) 2024.06 on the same input
records and tau grids.

!!! note "What is and isn't benchmarked here"
    The numbers below compare SigmaTau against **allantools** on identical
    input. **Stable32** is a closed-source Windows application; SigmaTau is
    cross-validated against it for *numerical agreement* (see
    [Validation](validation/methodology.md)), but it is not part of these
    speed comparisons. Treat "faster than Stable32" as unverified — the
    claim here is specifically *faster than allantools*.

All runs were on a 12-core Linux workstation, Julia 1.12.6 / Python 3.14.4,
allantools 2024.06. Each kernel was called one-shot (a single call passing the
full τ grid) through each library's public API. SigmaTau ran multithreaded
(`-t auto`, 12 threads); allantools ran single-threaded numpy, as it has no
parallel path. The modified-total kernels are where SigmaTau's threading and
sliding-window reductions pay off most.

## Synthetic, 30 realizations × N = 25 000

White phase noise, identical single-column inputs fed to both libraries, τ
grid `1, 2, 4, …, 512` (10 octave-spaced averaging factors). Mean wall time per
kernel call across 30 realizations:

| Kernel  | SigmaTau (mean) | allantools (mean) | Speedup |
|---------|----------------:|------------------:|--------:|
| adev    | 15 µs           | 393 µs            | 20.1×   |
| mdev    | 56 µs           | 1.13 ms           | 20.1×   |
| hdev    | 25 µs           | 341 µs            | 13.6×   |
| tdev    | 51 µs           | 991 µs            | 19.3×   |
| totdev  | 254 µs          | 684 µs            | 2.7×    |
| mtotdev | 30.3 ms         | 107.3 s           | 3,543×  |
| htotdev | 31.5 ms         | 114.3 s           | 3,632×  |

The cheap kernels (adev / mdev / hdev / tdev) finish in well under 100 µs,
near the timer resolution. The more reliable comparison is the ratio. The
largest difference is in the *modified*-total family: `mtotdev` and `htotdev`
are roughly **3,500–3,600×** faster. `totdev` (the non-modified total) gains a
more modest 2.7×.

## Real record, N = 406 763

A real phase record (`6krb25apr.txt`, τ₀ = 1 s), same τ grid capped at m = 512,
single one-shot call per kernel:

| Kernel  | SigmaTau | allantools | Speedup |
|---------|---------:|-----------:|--------:|
| adev    | <1 ms    | 9 ms       | 23.2×   |
| mdev    | 1 ms     | 21 ms      | 27.6×   |
| hdev    | <1 ms    | 15 ms      | 41.7×   |
| tdev    | 1 ms     | 21 ms      | 29.6×   |
| totdev  | 5 ms     | 16 ms      | 3.2×    |
| mtotdev | 0.469 s  | 1 803.6 s  | 3,842×  |
| htotdev | 0.455 s  | 1 928.0 s  | 4,233×  |

End to end, the full seven-kernel sweep runs in **0.93 s** versus allantools'
**3 731.6 s** (just over an hour) — about **4,000×** on this record. The
modified-total kernels dominate that total; on allantools `mtotdev` alone is
~30 minutes at this length, which is why the comparison caps the averaging
factor at 512.

## Caveats

- **Memory is reported on different bases** and isn't directly comparable:
  SigmaTau's figure is cumulative `@timed` bytes (every transient allocation,
  including reused buffers, ~2 MiB per heavy-kernel call), while allantools' is
  peak ΔRSS from a sampling poller (<300 KiB above the input array). Both
  indicate pressure; they answer different questions.
- **SigmaTau's run-to-run variance is higher on the threaded kernels** (~7 %
  relative std on the modified-total family from 12-thread scheduling jitter,
  versus <0.35 % for single-thread numpy). Distributions are symmetric with no
  long tail; the jitter shrinks with thread pinning.
- Numbers are from one workstation. Absolute times will vary with hardware and
  thread count; the *ratios* are the portable result.

## Reproducing

The figures above are a recorded snapshot; the dedicated benchmark harness is
no longer kept in-tree. To reproduce them, time each public deviation through a
warm-started one-shot call (pass the full τ grid in a single invocation) on
identical single-column inputs, running SigmaTau with `-t auto` and allantools
single-threaded — the method described at the top of this page.

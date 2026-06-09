---
name: Bug report
about: Report incorrect behavior, a crash, or a numerical discrepancy
title: ""
labels: bug
assignees: ""
---

## What happened

Describe the bug.

## Minimal reproduction

Use the smallest snippet that triggers it. Prefer a self-contained generator
(e.g. `noise_gen` or `randn` with a fixed seed) over an attached data file:

```julia
using SigmaTau
# ...
```

## Expected vs. actual

- **Expected:** expected behavior.
- **Actual:** observed behavior. Paste the full error and stacktrace if there is
  one.

## Numerical discrepancies

If a deviation value disagrees with a reference, please include:

- The exact call you made.
- The value SigmaTau returned.
- The reference value and its source (Stable32 / allantools / hand calc),
  including the τ (or averaging factor `m`) at which they differ.

## Environment

- SigmaTau.jl version (or commit):
- Julia version (`julia --version`):
- OS:

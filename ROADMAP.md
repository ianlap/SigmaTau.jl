# SigmaTau Roadmap

Last updated: 2026-06-07.

SigmaTau has two active planning tracks:

- [RELEASE.md](RELEASE.md) — the immediate Julia `0.4.0` release checklist.
- [PYTHON_PORT.md](PYTHON_PORT.md) — the future Python package design, gated on
  a frozen Julia reference release.

## Product Principles

- Math correctness outranks API novelty, speed claims, and backwards
  compatibility.
- Public APIs should return structured objects, not loose tuples.
- Defaults should match time-and-frequency convention, but every meaningful
  policy choice should be explicit and serializable.
- Avoid hidden global mutable state. Per-call options are easier to reproduce
  than process-wide configuration.
- Keep plotting optional. Core analysis must not depend on plotting packages.
- Keep IO boring and inspectable. TSV round trips are appropriate for scientific
  workflows and cross-language validation.
- Performance optimizations must preserve parity tests. No `fastmath` or
  reordered reductions unless the numerical tolerance change is deliberate and
  documented.

## Current Direction

Keep the Julia user surface flat for the `0.4.0` release:

```julia
using SigmaTau

data = read_phase("clock.dat"; tau0=1.0, time_col=0)
suite = stability(data; devs=(:adev, :mdev, :hdev, :pdev), taus=Octave)
save_suite("clock_stability.tsv", suite)
```

The current source layout is serviceable. A post-release mechanical cleanup can
revisit directory names and ownership boundaries, but it should not block the
tag.

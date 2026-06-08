# SigmaTau.jl 0.4.0 Release Checklist

Last updated: 2026-06-07.

This release is `0.4.0`, not `0.3.x`: it adds exported symbols and removes the
deprecated `ldev` alias.

## Status

| Gate | Status | Notes |
|------|--------|-------|
| Full package tests | Done | `Pkg.test()` passed locally: 2360 tests. |
| Stable32 fixtures | Done | Exercised by `test/stab/runtests.jl`. |
| allantools fixtures | Done | Exercised by `test/stab/allantools_cross_validation.jl`. |
| README API shape | Done | Current exported symbols and keywords match the documented quickstart. |
| License | Done | MIT license is present. |
| Registry compat | Done | Every dependency has `[compat]`; `Distributions` is widened to `0.25`. |
| Threaded CI coverage | Done | CI includes a Julia 1.11 / Ubuntu / `JULIA_NUM_THREADS=auto` leg. |
| Public docs wiring | Done | Missing `@docs` entries were added for suites, tau grids, defaults, and suite IO. |
| Docs build | Done | The docs env instantiates locally; `docs/make.jl` builds with GR fallback when no TeX engine is installed. |
| Tables.jl extension | Done | `StabilityResult` and `StabilitySuite` implement row tables via weakdep extension. |
| Project version bump | Done | `Project.toml` is set to `0.4.0`. |

## Immediate Release Work

The docs build and full test suite were re-verified on 2026-06-07: `Pkg.test()`
passed 2360/2360, and `julia --project=docs docs/make.jl` built cleanly (inventory
version 0.4.0, GR fallback). Only the release mechanics remain:

1. Commit the release-prep diff.
2. Tag and publish.

Re-run both checks once more immediately before tagging if any code lands after
this checkpoint.

## API Decisions

Keep these decisions scoped to `0.4.0`:

- Keep the flat `using SigmaTau; adev(...)` public surface.
- Keep named deviation functions (`adev`, `mdev`, `hdev`, `pdev`, `mtie`, etc.).
- Keep `stability` as the batch API.
- Remove the deprecated `ldev` alias; `htdev` is the canonical function.
- Do not add mutable global defaults for confidence level or tau grids.
- Support Tables.jl through a weakdep extension, not a hard dependency.

Potential v1 polish, not `0.4.0` blockers:

- Add `analyze(...)` as a friendlier alias for `stability(...)`.

## Validation Contract

Every estimator should have at least one anchor:

- External reference: Stable32, allantools, NIST examples, or published tables.
- Legacy reference: MATLAB-era or prior SigmaTau kernels under `test/`.
- Internal identity: e.g. `tdev = mdev * tau / sqrt(3)`, `pdev(tau0) = adev(tau0)`.
- Monte Carlo contract: only when no external implementation exists, with a
  reproducible harness, seed, artifact, and documented validity window.

Tolerance policy:

- Direct legacy/kernel parity tests target `rtol = 1e-12` where platform SIMD
  order allows it.
- The broader all-noise-type kernel sweep uses `rtol = 1e-11` for known
  cross-platform SIMD reduction drift.
- Stable32 CSV comparisons tolerate Stable32's printed precision.
- Total-family comparisons separate raw kernel parity from bias-corrected API
  policy.
- Multithreaded total-family kernels use fixed-order chunk reductions, but
  exact bits may still vary across thread counts because the number of chunks
  follows the configured thread count. Treat observed cross-thread drift at the
  few-ULP level as acceptable unless it exceeds the documented test tolerance.

## Post-Release Cleanup

The internal layout can be revisited after `0.4.0` if the payoff is clear. A
possible mechanical layout is:

```text
src/
  data/
  deviations/core/
  deviations/api/
  stats/
  noise/
  spectral/
  io/
```

Do not do that restructure in the release-prep commit.

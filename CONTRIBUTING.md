# Contributing to SigmaTau.jl

Contributions are welcome: bug reports, examples, documentation fixes,
benchmarks, validation data, and new stability estimators all help.

SigmaTau.jl is a Julia package for time-and-frequency stability analysis. The
main goals are accurate numerical results, clear APIs, and reproducible
validation against established references.

## Getting set up

You need Julia 1.11 or newer.

```bash
git clone https://github.com/ianlap/SigmaTau.jl
cd SigmaTau.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` is not tracked; Julia resolves dependencies from
`Project.toml`.

## Running tests

Run the full test suite with:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests live under `test/`. The stability tests include Stable32 and
allantools fixture checks from `test/fixtures/validation/`, plus legacy-kernel
parity tests for estimators without an external implementation.

For smaller changes, it is fine to run a focused test file first, then run the
full suite before opening a PR.

## Numerical changes

For changes to deviation kernels, EDF/confidence intervals, noise
identification, or bias correction, please include a test or reference that
shows why the new result is correct. Published formulas should cite the source
(for example NIST SP1065, Greenhall/Riley, or IEEE 1139).

If you are reporting a numerical discrepancy, please include:

- the input data or a small generator;
- the exact SigmaTau call;
- the value SigmaTau returned;
- the reference value and where it came from.

## Documentation

The manual is built with Documenter.jl from `docs/`. Runnable examples live in
`examples/` and are rendered into tutorials with Literate.jl during the docs
build.

## Pull requests

Small, focused PRs are easiest to review. Please include tests or a short note
about what was checked. If a change affects the public API, update the manual
and add a brief entry under `## [Unreleased]` in `CHANGELOG.md`.

By contributing, you agree that your contributions are licensed under the
project's [MIT License](LICENSE).

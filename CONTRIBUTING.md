# Contributing

Bug reports, validation data, benchmarks and new estimators are welcome.

```bash
git clone https://github.com/ianlap/SigmaTau.jl
cd SigmaTau.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

For a change to a kernel, EDF model, bias correction or noise identifier,
include a test or a reference that shows the new result is right, and cite the
source (NIST SP1065, Greenhall & Riley 2003, IEEE 1139-2022, or the paper).
When reporting a numerical discrepancy, include the input (or a seeded
generator), the exact call, the value returned, and the reference value with
its source.

Public API changes need a docstring update and a line under `Unreleased` in
`CHANGELOG.md`. Contributions are licensed under the project's MIT license.

# Contributing to SigmaTau.jl

Thanks for your interest in contributing. SigmaTau.jl is a single, flat Julia
package for time-and-frequency stability analysis (Allan / Modified Allan /
Hadamard / Total deviation families, MTIE, parabolic deviation, lag-1 noise
identification, Greenhall–Riley EDF/χ² confidence intervals, and calibrated
power-law noise generation). Bug reports, new deviations, validation fixtures,
documentation, and benchmarks are all welcome.

This guide covers how to get set up, the conventions the codebase follows, and
how changes are verified. For a per-component map of the package, see
[`project_overview.md`](project_overview.md).

## Getting set up

You need Julia 1.11 or newer.

```bash
git clone https://github.com/ianlap/SigmaTau.jl
cd SigmaTau.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` is gitignored; `instantiate` resolves it from `Project.toml`.

For interactive development, [Revise.jl](https://github.com/timholy/Revise.jl)
hot-patches source edits without restarting Julia (struct-field changes,
`Project.toml` edits, and new `@eval`'d definitions still need a fresh session).

## Running the tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite lives under `test/` and runs five sub-suites (`types/`, `stab/`,
`io/`, `umbrella_smoke.jl`, `recipes.jl`). Please add or update tests for any
behavior you change, and make sure the full suite passes before opening a PR.

Two periodic quality checks (not part of `Pkg.test()`):

```julia
using Aqua; Aqua.test_all(SigmaTau)      # method ambiguities, stale deps, piracy
using JET;  report_package(SigmaTau)     # type instabilities across the package
```

## Numerical correctness

Deviation kernels are cross-validated against three independent references:

- **Stable32** (W. Riley) — the de facto industry desktop reference.
- **allantools** (A. Wallin) — Python second reference.
- **AllanLab** (MATLAB) — third reference, locked-in fixtures.

Fixtures live under `reference/validation/` (read-only). After any change to a
core kernel (`src/stab/core/*.jl`) run the parity testsets in
`test/stab/runtests.jl` (Stable32 cross-validation, `allantools_cross_validation.jl`,
and the rtol = 1e-12 `legacy_kernels.jl` contract) and confirm they still pass.

Do not invent new χ²/EDF formulas. Confidence-interval and EDF math cites its
source — NIST SP1065 (Riley & Howe), Greenhall & Riley 2003, IEEE Std
1139-2022. New statistical expressions need a citation.

## Code conventions

- **Core / API split is firm.** Core kernels (`_adev_core`, …, in
  `src/stab/core/`) take `Vector{Float64}` and return raw arrays. Public API
  functions (`adev`, …, in `src/stab/api/`) take `PhaseData` / `FrequencyData`
  and return a `StabilityResult`. Never collapse the two layers.
- **`StabilityResult` fields stay non-parametric `Vector{Float64}`.** Do not
  parameterize them.
- **The `edf` / CI fields are empty when `calc_ci=false` and populated when
  `calc_ci=true`.** Preserve this contract.
- **Plot recipes live only in `ext/SigmaTauRecipesBaseExt.jl`** (a `RecipesBase`
  package extension). Do not add plotting code to `src/`.
- **Every exported function has a docstring.** 4-space indent, no trailing
  whitespace, match the surrounding style.

### Adding a new deviation

1. Add the pure kernel to the relevant `src/stab/core/*.jl` as
   `_yourdev_core(x::Vector{Float64}, m_values::Vector{Int}, tau0::Float64)`.
2. Add **both** a `PhaseData` and a `FrequencyData` public entry point in the
   matching `src/stab/api/*.jl`; the `FrequencyData` method delegates via
   `_freq_to_phase`. Add the zero-arg and `TauMode` convenience methods too.
3. Wire it into the `stability` router (`src/stab/api/suite.jl`) if it should be
   reachable from the compute-all entry point.
4. Add validation tests; if no external reference exists, validate against
   legacy kernels (rtol = 1e-12) and internal consistency.
5. Export it from `src/SigmaTau.jl` and add it to the umbrella smoke test.

## Documentation, CHANGELOG, and TODO

Every shipped code change should, in the same commit:

1. Remove the matching item from [`TODO.md`](TODO.md) if one exists.
2. Add a terse, past-tense entry under `## [Unreleased]` in
   [`CHANGELOG.md`](CHANGELOG.md) (Keep-a-Changelog format; tag breaking
   changes). Mark pure-docs/typo changes as not warranting an entry in the
   commit body.
3. Refresh [`project_overview.md`](project_overview.md) if the change alters the
   public surface (new exported function, new type, new role).

User-facing docs are built with Documenter.jl from `docs/`; runnable tutorials
under `examples/` render via Literate.jl.

## Commit and PR style

- Commit messages: imperative mood, ≤ 72-character subject, body explains *why*.
- CHANGELOG entries and commit messages use a terse, factual voice — past tense,
  no marketing language, no emoji.
- Open a pull request against `main`. The PR template includes a short checklist
  (tests pass, CHANGELOG/TODO updated, docstrings present).

## Reporting bugs and proposing features

Use the [issue templates](https://github.com/ianlap/SigmaTau.jl/issues/new/choose).
For a numerical discrepancy, please include the input data (or a generator),
the call you made, the value you got, and the reference value you expected.

By contributing, you agree that your contributions are licensed under the
project's [MIT License](LICENSE).

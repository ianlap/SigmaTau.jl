# Contributing

Bug reports, validation data, benchmarks and new estimators are welcome, as
issues or pull requests at <https://github.com/ianlap/SigmaTau.jl>.

## Reporting a numerical difference

Give the input, or a seeded generator of it, the exact call, the value
returned, and the reference value with its source.

## Changing the code

```bash
git clone https://github.com/ianlap/SigmaTau.jl
cd SigmaTau.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

A change to a kernel, an EDF model, a bias factor or the noise
identification comes with a test or a reference that shows the new result
is right, and it cites its source: NIST SP1065, IEEE Std 1139-2022, or the
paper. A source that the documentation cites has a record in
`docs/src/refs.bib`.

A change to the public interface comes with its docstring and a line under
`Unreleased` in `CHANGELOG.md`.

## Tests

| Suite | Command |
|:---|:---|
| The package, with the cross-check against other programs | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| The cross-check alone | `julia --project=. -e 'using Test; include("test/validation/runtests.jl")'` |
| The command line | `julia --project=cli cli/test/runtests.jl` ([`cli/README.md`](cli/README.md)) |
| The explorer | [`app/web/README.md`](app/web/README.md) |
| The Python package | `pip install -e "python/[dev]"`, then `pytest python/tests -q` |

## Documentation

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl          # writes docs/build
```

The build runs the examples of the pages and the doctests of the
docstrings, turns the scripts of `examples/` into the tutorials, and fails
on an exported name without a docstring.

A docstring follows the layout of those in `src/`: the signature with the
return type, one sentence that says what the function does, the
definition, then `# Arguments`, `# Keywords`, `# Returns` and `# Examples`.
A statement that rests on a source cites it with the section or the
equation, as in `[SP1065 §5.2.4, Eq. (11)](@cite riley-2008-sp1065)`.

## License

Contributions are licensed under the MIT license of the project.

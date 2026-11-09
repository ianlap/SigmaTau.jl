# sigmatau command line

`SigmaTauCLI` is the `sigmatau` program: one binary with the subcommands
`analyze`, `check`, `spectrum`, `dynamic`, `cross`, `drift`, `report`, `live`,
`completion` and `version`. The explorer (`app/web`, `app/desktop`) runs it with `--format json`. Usage, the
CSV columns, the JSON schema and the exit codes are in
[docs/src/guide/cli.md](../docs/src/guide/cli.md).

```
cli/
  Project.toml        SigmaTau from this repository, FFTW from the trim-safe fork
  src/SigmaTauCLI.jl  entry point, errors, exit codes
  src/options.jl      the option table (parser, help and completions read it)
  src/record.jl       reading or synthesising the record
  src/commands.jl     analyze, check, spectrum, dynamic, cross, drift
  src/report.jl       report: the self-contained HTML page
  src/live.jl         live: instrument or simulator readings to running σ(τ)
  src/output.jl       JSON and CSV writers
  src/completion.jl   bash, zsh, fish and PowerShell scripts
  test/runtests.jl    tests (in process, and against a build with SIGMATAU_EXE)
```

## Build

```bash
julia --project=cli -e 'using Pkg; Pkg.instantiate()'
juliac --output-exe sigmatau --bundle cli/build --trim=safe --experimental cli
```

The build takes a few minutes and must pass `--trim=safe` with no errors.
It depends on the FFTW.jl fork `ianlap/FFTW.jl` (branch `juliac-trim`) with
the `eager_library = true` preference, both set in `Project.toml`: stock
FFTW.jl dispatches dynamically when it frees plans and loads its library
lazily, which trimming rejects or which crashes at the first FFT.

The code follows a few rules so it stays trim-safe: options live in one
concrete struct, estimators are picked through the library's registry (an
if-chain), keyword arguments are concrete at every call (no
`Union{Nothing,…}` keywords, so optional ones get their own call sites), no
`do`-block closures over reassigned locals, no `showerror`, no JSON package,
and all output is rendered into an `IOBuffer` before it is written.

## Test

```bash
julia --project=cli cli/test/runtests.jl
SIGMATAU_EXE=cli/build/bin/sigmatau julia --project=cli cli/test/runtests.jl
```

The second form also runs the compiled binary and checks that it answers
byte for byte like the in-process call.

Measured on a 4-core Linux container (Julia 1.13.3): the build takes about
a minute once the packages are precompiled, the executable is 15 MB in a
132 MB bundle, `version` returns in 24 ms, six deviations on 8192 points in
29 ms, and `analyze --spectrum` with the default four estimators on 65 536
points in 64 ms.

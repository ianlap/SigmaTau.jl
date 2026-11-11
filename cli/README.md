# The `sigmatau` command line

`SigmaTauCLI` is the source of `sigmatau`, the package as one compiled
program. `sigmatau help` lists its commands. The
[guide](https://ianlap.github.io/SigmaTau.jl/dev/guide/cli/) covers the
common tasks, and the
[reference](https://ianlap.github.io/SigmaTau.jl/dev/reference/cli/) lists
every command and option, the columns of the CSV output, the JSON schema
and the exit codes. The explorer (`app/web`, `app/desktop`) runs this
program for every analysis it shows.

```text
cli/
  Project.toml        SigmaTau from this repository, FFTW from a fork
  src/SigmaTauCLI.jl  entry point, errors, exit codes
  src/options.jl      the table of options, read by the parser, the help and completion
  src/record.jl       reads or simulates the record
  src/commands.jl     the analysis commands
  src/report.jl       the HTML report
  src/live.jl         readings from an instrument or the simulator
  src/output.jl       the JSON and CSV writers
  src/completion.jl   the bash, zsh, fish and PowerShell scripts
  test/runtests.jl    the tests
```

## Build

With Julia 1.13 or later. On Julia 1.12 the build stops: a program that
1.12 trims cannot start a task, so the threaded loops of the total
deviations and Thêo1 would fail in it. With juliaup, `juliaup add 1.13`
and `julia +1.13` in the first two lines leave the default version as it is.

```bash
julia -e 'using Pkg; Pkg.Apps.add("JuliaC")'      # installs ~/.julia/bin/juliac
julia --project=cli -e 'using Pkg; Pkg.instantiate()'
juliac --output-exe sigmatau --bundle cli/build --trim=safe --experimental cli
cli/build/bin/sigmatau version
```

The build takes a few minutes. The folder `cli/build` holds the program
and the Julia runtime libraries it loads, about 150 MB in all.

`Project.toml` takes FFTW.jl from the fork `ianlap/FFTW.jl`, branch
`juliac-trim`, with the preference `eager_library = true`. The registered
FFTW.jl dispatches dynamically when it frees a plan and loads its library
on first use. Under trimming these fail the build or crash the program at
its first FFT.

## Code that trims

The build must pass `--trim=safe` without errors. The code in `cli/src`
keeps to these rules for it:

- the options live in one concrete struct;
- an estimator is chosen through the registry of the package, which is an
  if-chain;
- keyword arguments have a concrete type at every call, so an optional
  keyword gets a call of its own in place of a `Union{Nothing,…}`;
- no `do` block closes over a local variable that is assigned again;
- no `showerror`, and no JSON package: `output.jl` writes the JSON;
- output is written to an `IOBuffer` first and to the stream at the end.

## Test

```bash
julia --project=cli cli/test/runtests.jl
SIGMATAU_EXE=cli/build/bin/sigmatau julia --project=cli cli/test/runtests.jl
```

The first form runs the commands inside the Julia process. The second
also runs the compiled program and requires its output to match that of
the same call in the process: the same text, with every number equal
within a relative difference of 10⁻¹².

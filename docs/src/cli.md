# Command-line interface

SigmaTau ships a non-interactive command-line interface that wraps the same
numerical kernels used from Julia, so phase/frequency stability analysis can be
driven from bash, PowerShell, Python, MATLAB, CI pipelines, or a desktop
frontend. There is no interactive/TUI behavior: results go to stdout (or a
file), diagnostics go to stderr, and the process returns a meaningful exit code.

## Running it

The entry point is `bin/sigmatau.jl`, run through Julia with the package project
active:

```bash
julia --project=. bin/sigmatau.jl --help
julia --project=. bin/sigmatau.jl version
julia --project=. bin/sigmatau.jl analyze clock.DAT --phase --tau0 1 \
      --devs adev,mdev,hdev --tau octave --format csv --out results.csv
```

The same logic is callable from Julia as `SigmaTau.cli_main(args)::Cint` (it
returns the exit code and never calls `exit`, so it is safe to embed and is
ready for a future PackageCompiler build). In the examples below, substitute
`julia --project=. bin/sigmatau.jl` wherever `sigmatau` appears.

## Commands

| Command | Purpose |
|---------|---------|
| `sigmatau analyze FILE [options]` | compute stability deviations from a data file |
| `sigmatau completion SHELL` | print a shell completion script (`bash`/`zsh`/`fish`/`powershell`) |
| `sigmatau version` | print the SigmaTau version |
| `sigmatau help` | print usage (also `--help`, `-h`, or no arguments) |

## `analyze`

`FILE` is a whitespace/CSV/TSV table; by default column 1 is the timestamp and
column 2 is the value. Leading `#` comment lines (and a `# tau0:` header written
by `save`) are skipped automatically.

| Option | Meaning |
|--------|---------|
| `--phase` | treat FILE as phase data (default) |
| `--freq` | treat FILE as fractional-frequency data |
| `--tau0 FLOAT` | base sample interval τ₀ in seconds |
| `--time-col N` | 1-based time column (`0` = no time column) |
| `--value-col N` | 1-based value column |
| `--devs LIST` | comma-separated deviations (default `adev,mdev,hdev,mhdev`) |
| `--tau GRID` | named grid: `octave` (default), `decade`, `alltaus`, `halfoctave`, `quarteroctave`, `halfdecade` |
| `--taus LIST` | explicit comma list of integer averaging factors m (τ = m·τ₀); mutually exclusive with `--tau` |
| `--ci BOOL` | compute confidence intervals (default `true`) |
| `--confidence FLOAT` | CI confidence level in (0,1) (default `0.683`, the 1-σ convention) |
| `--format FMT` | output format: `csv` (default) or `json` |
| `--out PATH` | write output to PATH (default: stdout) |
| `--header N` | non-comment header lines to skip |
| `--delim CHAR` | field delimiter (default: auto from extension; `\t` for a tab) |
| `--scaling FLOAT` | multiply samples by this factor (e.g. `1e-12` for ps→s) |
| `--detrend MODE` | detrend before analysis: `none`/`mean`/`endpoint`/`linear` |
| `--quiet` | suppress informational messages |
| `--debug` | print stack traces on error |

Valid deviations for `--devs`:

```text
adev mdev tdev hdev mhdev htdev
totdev mtotdev ttotdev htotdev mhtotdev
mtie tierms pdev theo1 theoh
```

The dynamic deviations `dadev`/`dhdev` are **not** available in `analyze`: they
produce a time-resolved σ_y(t, τ) map (a waterfall / 3-D surface of σ versus τ
over time) that does not fit the per-τ table schema. Use the `dadev`/`dhdev`
Julia API for those. They are still offered in shell completion as known names.

### Output schema

**CSV** — one row per (deviation, τ), with the fixed header:

```text
source,quantity,deviation,tau,dev,ci_lower,ci_upper,edf,neff,noise_type
```

`source` is the input path, `quantity` is `phase` or `frequency`, `tau` is in
seconds. When confidence intervals are off (`--ci false`) or undefined for a
deviation, `ci_lower`/`ci_upper`/`edf`/`noise_type` are written **blank** (never
`NaN` or `nothing`); `neff` is always present.

**JSON** — a stable, versioned schema (NaN/missing serialize as `null`):

```json
{
  "schema_version": 1,
  "source": "clock.DAT",
  "analysis": {
    "quantity": "phase",
    "tau0": 1.0,
    "n": 1024,
    "ci": true,
    "confidence": 0.683,
    "tau_mode": "Octave",
    "deviations": ["adev", "mdev"]
  },
  "results": [
    {
      "deviation": "adev",
      "points": [
        {"tau": 1.0, "dev": 1.2e-11, "ci_lower": 1.0e-11, "ci_upper": 1.4e-11,
         "edf": 12.0, "neff": 1023, "noise_type": "WHFM"}
      ]
    }
  ]
}
```

### Exit codes

| Code | Meaning |
|------|---------|
| `0` | success |
| `1` | usage error (bad option, missing FILE, unknown deviation, invalid τ grid/list) |
| `2` | input/read error (bad path, non-numeric content) |
| `3` | analysis error (e.g. record too short for the requested grid) |
| `4` | output/write error |

Ordinary user mistakes print a one-line message to stderr without a Julia stack
trace; pass `--debug` to see the trace. A downstream consumer closing the pipe
early (`sigmatau analyze … | head`) exits cleanly.

## Shell completion

`completion SHELL` prints a static completion script that completes subcommands,
options, known option values (`octave`/`decade`, `csv`/`json`, `true`/`false`),
and deviation names — **while preserving your shell's native file-path
completion** for `FILE` and path options like `--out` (all files and
directories, no extension filtering). It performs no analysis or file reading.

```bash
# bash — write to the completion dir, or source it for the current shell
sigmatau completion bash > /usr/local/etc/bash_completion.d/sigmatau
source <(sigmatau completion bash)

# zsh — install onto a directory in $fpath, then restart the shell
sigmatau completion zsh > "${fpath[1]}/_sigmatau"

# fish
sigmatau completion fish > ~/.config/fish/completions/sigmatau.fish

# PowerShell — add to your profile
sigmatau completion powershell | Out-String | Invoke-Expression
```

(Substitute `julia --project=. bin/sigmatau.jl` for `sigmatau` until an
installed shim exists. Defining a shell alias or wrapper named `sigmatau` makes
the completion scripts apply directly.)

## Examples

```bash
# Phase .DAT with a single value column, decade grid, JSON to a file
sigmatau analyze clock.DAT --time-col 0 --value-col 1 --tau0 1 \
  --devs adev,hdev --tau decade --format json --out clock.json

# Fractional-frequency CSV, explicit averaging factors, no CIs, to stdout
sigmatau analyze freq.csv --freq --tau0 10 --taus 1,2,4,8,16 --ci false

# Phase data in picoseconds, scaled to seconds, with a linear detrend
sigmatau analyze raw.txt --time-col 0 --value-col 1 --tau0 1 \
  --scaling 1e-12 --detrend linear --devs adev,mdev,tdev
```

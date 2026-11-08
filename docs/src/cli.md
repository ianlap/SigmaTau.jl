# Command line

`sigmatau` is one compiled program with a subcommand per task. It runs the
same library code as the Julia API, starts in about 30 ms, and needs no Julia
install. The stability explorer runs the same binary with `--format json`.

```bash
sigmatau analyze clock.dat --tau0 1 --devs adev,mdev,tdev
sigmatau analyze freq.txt --freq --m 1,10,100 --format json --out result.json
sigmatau check clock.dat --time-col 1 --value-col 3 --fillgaps
sigmatau spectrum clock.dat --nperseg 4096
sigmatau dynamic clock.dat --dev dhdev --window 4096
```

## Getting it

Build it from the repository with Julia 1.12 or newer and
[JuliaC](https://github.com/JuliaLang/JuliaC.jl) (`pkg> app add JuliaC`):

```bash
julia --project=cli -e 'using Pkg; Pkg.instantiate()'
juliac --output-exe sigmatau --bundle cli/build --trim=safe --experimental cli
cli/build/bin/sigmatau version
```

The bundle (`cli/build/`, about 130 MB) holds the executable and the Julia
runtime libraries it needs; keep them together when you copy it. Without a
build, `julia --project=cli -m SigmaTauCLI ARGS…` runs the same commands, but
each call pays Julia's start-up.

## Commands

| Command | Does |
| --- | --- |
| `sigmatau analyze FILE` | stability deviations σ(τ), with confidence intervals, EDF and noise type per τ |
| `sigmatau check FILE` | read the record and summarise it: samples, span, gaps, 5σ outliers, linear drift per day |
| `sigmatau spectrum FILE` | Welch power spectral densities S_y(f) and S_x(f) (Hann taper, half overlap) |
| `sigmatau dynamic FILE` | dynamic ADEV or HDEV: σ(τ) in sliding windows |
| `sigmatau completion SHELL` | a completion script for bash, zsh, fish or PowerShell |
| `sigmatau version` | the SigmaTau version |
| `sigmatau help [COMMAND]` | usage, or one command's options (also `COMMAND --help`) |

`FILE` is a whitespace, comma or tab separated table. By default column 1 is
time and column 2 the value. Lines starting with `#` and a Stable32 header
block are skipped, and a `# tau0:` line (as `save` writes) sets τ₀. Every
analysis command also takes `--synth N` in place of `FILE`, a synthetic
white FM plus random-walk FM record for trying things out.

## Options

The record (all four analysis commands):

| Option | Meaning |
| --- | --- |
| `--phase` / `--freq` | the value column is phase in seconds (default) or fractional frequency |
| `--tau0 SECONDS` | sample interval; default from the file header or the time column |
| `--time-col N`, `--value-col N` | 1-based columns; `--time-col 0` means no time column |
| `--reference-col N` | subtract this column from the value column (several clocks against one reference) |
| `--header N` | extra header lines to skip |
| `--delim CHAR` | field delimiter; default from the extension, `'\t'` for a tab |
| `--scaling X` | multiply the values by X (`1e-9` for nanoseconds) |
| `--detrend MODE` | `none` (default), `mean`, `endpoint` or `linear`, applied by the reader to the column as read |
| `--remove-drift MODE` | `none` (default), `linear` (frequency offset) or `quadratic` (offset and linear frequency drift): a least-squares fit to the phase, over the non-gap samples, removed after reading |
| `--fillgaps` | put the record on an even time grid and fill the gaps (needs a time column) |
| `--remove-outliers` | replace 5σ outliers in frequency and fill them |
| `--synth N`, `--seed S` | synthetic record instead of `FILE` (seed 42) |

`analyze`:

| Option | Meaning |
| --- | --- |
| `--devs LIST` | comma list (default `adev,mdev,hdev,mhdev`) of `adev mdev tdev hdev mhdev htdev totdev mtotdev ttotdev htotdev mhtotdev theo1 theoh pdev mtie tierms` |
| `--grid GRID` | `octave` (default), `half-octave`, `quarter-octave`, `decade`, `half-decade` or `all`; `halfoctave` and the like work too |
| `--m LIST` | explicit averaging factors, τ = m·τ₀, in place of `--grid` |
| `--confidence P` | interval level (default 0.683, one sigma) |
| `--no-ci` | skip intervals, EDF and noise identification |
| `--no-bias-correction` | leave the total estimators and Thêo1 uncorrected |
| `--spectrum` | add the spectrum to the JSON output |
| `--nperseg N` | Welch segment length for `--spectrum` |

`spectrum`: `--nperseg N` (default: the largest power of two ≤ N/4, at least 16).

`dynamic`: `--dev dadev|dhdev` (default `dadev`), `--grid GRID`,
`--window N` (default: the largest power of two ≤ N/8, at least 16) and
`--step N` (default half the window).

Output (all commands): `--format csv|json` (default `csv`), `--out PATH`
in place of stdout, `--series N` to include up to N record points in the
JSON, and `--quiet` to drop the "wrote PATH" note.

## CSV

`analyze` writes one row per deviation and τ:

```text
source,quantity,deviation,tau,dev,ci_lower,ci_upper,edf,neff,noise_type
```

`tau` is in seconds. Without intervals (`--no-ci`, or a τ where they are
undefined) `ci_lower`, `ci_upper`, `edf` and `noise_type` are blank; a value
is never written as `NaN`. `check` writes
`source,quantity,n,tau0,span_s,gaps,outliers,drift_per_day`, `spectrum`
writes `f,sy,sx` (without the DC bin), and `dynamic` writes `t,tau,dev`, one
row per window and τ, with `t` the window centre in seconds.

## JSON

Every command answers with one object in a versioned schema. A field is
only added within a schema version, never renamed or removed; non-finite
numbers are `null`.

```json
{
  "schema_version": 1,
  "sigmatau": "0.5.0",
  "command": "analyze",
  "source": "clock.dat",
  "record": {
    "quantity": "phase", "n": 16384, "tau0": 1.0, "span_s": 16384.0,
    "gaps": 0, "outliers": 2, "outlier_t": [812.0, 9031.0], "outlier_x": [3.1e-9, 3.4e-9],
    "drift_per_day": 3.1e-14, "removed_offset": null, "removed_drift_per_day": null,
    "series": { "t": [0.0, 4.0], "x": [0.0, 1.2e-12] }
  },
  "analysis": {
    "quantity": "phase", "tau0": 1.0, "n": 16384, "tau_mode": "octave",
    "ci": true, "confidence": 0.683, "bias_correction": true,
    "deviations": ["adev", "tdev"],
    "detrend": "none", "remove_drift": "none", "fillgaps": false, "remove_outliers": false, "reference_col": 0
  },
  "results": [
    { "deviation": "adev", "label": "ADEV", "axis": "sigma_y",
      "points": [ { "tau": 1.0, "dev": 1.2e-11, "ci_lower": 1.1e-11, "ci_upper": 1.3e-11,
                    "edf": 16382.0, "neff": 16382, "noise_type": "WHFM" } ] }
  ],
  "spectrum": { "nperseg": 4096, "noverlap": 2048, "f": [], "sy": [], "sx": [] },
  "elapsed_ms": 41.7
}
```

- `record` is in every answer; `series` only with `--series N`, as an evenly
  strided subset of the phase. At most 500 outliers are listed in
  `outlier_t`/`outlier_x`; `outliers` is the full count.
- `axis` is `sigma_x` (seconds) for TDEV, HTDEV, TTOTDEV, MTIE and TIE rms,
  and `sigma_y` for the rest. With `--m`, `tau_mode` is `"m"` and `m` lists
  the factors.
- `spectrum` comes from `sigmatau spectrum`, or `analyze --spectrum`. Its
  arrays leave out the DC bin.
- `sigmatau dynamic` adds `"dynamic": {"deviation", "window", "step", "t",
  "tau", "dev"}`, where `dev` has one row per window (centre time `t`, in
  seconds) and one column per τ.
- On failure, with `--format json`, stdout holds
  `{"schema_version": 1, "sigmatau", "command", "source", "error": "…", "exit_code": n}`.

## Exit codes

| Code | Meaning |
| --- | --- |
| 0 | success |
| 1 | usage: unknown command or option, bad value, missing FILE |
| 2 | input: missing file, unreadable or non-numeric content, bad column |
| 3 | analysis: record too short, gaps where the spectrum needs none |
| 4 | output: `--out` cannot be written |

Errors print one line, `sigmatau: …`, on stderr. A consumer that closes the
pipe early (`sigmatau analyze … | head`) is not an error.

## Shell completion

`sigmatau completion SHELL` prints a static script that completes commands,
options, deviation names and option values, and leaves file names to the
shell.

```bash
source <(sigmatau completion bash)                                  # bash, this shell
sigmatau completion bash > /usr/local/etc/bash_completion.d/sigmatau  # bash, always
sigmatau completion zsh > "${fpath[1]}/_sigmatau"                     # zsh
sigmatau completion fish > ~/.config/fish/completions/sigmatau.fish   # fish
sigmatau completion powershell | Out-String | Invoke-Expression      # PowerShell ($PROFILE)
```

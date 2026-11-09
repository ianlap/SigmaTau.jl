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
|:---|:---|
| `sigmatau analyze FILE` | stability deviations σ(τ), with confidence intervals, EDF and noise type per τ |
| `sigmatau check FILE` | read the record and summarise it: samples, span, gaps, 5σ outliers, linear drift per day |
| `sigmatau spectrum FILE` | Welch power spectral densities S_y(f) and S_x(f) (Hann taper, half overlap) |
| `sigmatau dynamic FILE` | dynamic ADEV or HDEV: σ(τ) in sliding windows |
| `sigmatau drift FILE` | fit a frequency drift model: aging per day with its uncertainty |
| `sigmatau cross FILE --with FILE2` | two channels measuring one device: cross Allan deviation and cross-spectrum ([two-channel measurements](cross.md)) |
| `sigmatau live SOURCE` | running σ(τ) from an instrument or the simulator, one JSON line per update |
| `sigmatau report FILE` | a self-contained HTML report: record, settings, σ(τ) plots, per-τ tables and the command line (prints to PDF) |
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
|:---|:---|
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
| `--range A:B` | analyse only seconds A to B from the first sample (`3600:` and `:86400` leave one end open); applied before the preprocessing |
| `--synth N`, `--seed S` | synthetic record instead of `FILE` (seed 42) |

`analyze`:

| Option | Meaning |
|:---|:---|
| `--devs LIST` | comma list (default `adev,mdev,hdev,mhdev`) of `adev mdev tdev hdev mhdev htdev totdev mtotdev ttotdev htotdev mhtotdev theo1 theoh pdev mtie tierms` |
| `--grid GRID` | `octave` (default), `half-octave`, `quarter-octave`, `decade`, `half-decade` or `all`; `halfoctave` and the like work too |
| `--m LIST` | explicit averaging factors, τ = m·τ₀, in place of `--grid` |
| `--tau-min S`, `--tau-max S` | leave out averaging times below or above S seconds |
| `--confidence P` | interval level (default 0.683, one sigma) |
| `--no-ci` | skip intervals, EDF and noise identification |
| `--no-bias-correction` | leave the total estimators and Thêo1 uncorrected |
| `--spectrum` | add the spectrum to the JSON output |
| `--nperseg N` | Welch segment length for `--spectrum` |
| `--fit` | add the power-law noise model fitted to MHDEV to the JSON output |

`spectrum`: `--nperseg N` (default: the largest power of two ≤ N/4, at least 16).

`dynamic`: `--dev dadev|dhdev` (default `dadev`), `--grid GRID`,
`--window N` (default: the largest power of two ≤ N/8, at least 16) and
`--step N` (default half the window).

`drift`: `--model linear|quadratic|log|diffusion` (default `linear`); see
[`drift_fit`](@ref).

`cross`: `--with FILE2`, the second channel, read with the same options as
`FILE`; `--grid` or `--m` for the τ values and `--nperseg N` for the
cross-spectrum. Both records must have the same length and τ₀ and no gaps.

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
row per window and τ, with `t` the window centre in seconds. `drift` writes
one row, `source,quantity,model,a,a_err,b,b_err,c,c_err,drift_per_day,drift_per_day_err,rms,n`,
with `c` blank for the linear model. `cross` writes
`tau,covar,dev,uncertainty,dev_a,dev_b,neff`, one row per τ (see
[`cross_adev`](@ref)).

## JSON

Every command answers with one object in a versioned schema. A field is
only added within a schema version, never renamed or removed; non-finite
numbers are `null`.

```json
{
  "schema_version": 1,
  "sigmatau": "0.6.0",
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
- `powerlaw` comes from `analyze --fit`: the fitted model
  (`deviation`, `alphas`, `h`, `noise_types`; see [`powerlaw_fit`](@ref)),
  its σ(τ) for each ADEV, MDEV, TDEV, HDEV and MHDEV result in `curves`, and
  `sigma_at`, the model's ADEV at 1 s, 10 s, 100 s, 1000 s, 10⁴ s and one day
  (those within half the record). A record the fit cannot use gives
  `{"error": "…"}` there instead.
- `sigmatau dynamic` adds `"dynamic": {"deviation", "window", "step", "t",
  "tau", "dev"}`, where `dev` has one row per window (centre time `t`, in
  seconds) and one column per τ.
- `sigmatau drift` adds `"drift": {"model", "coef", "stderr",
  "drift_per_day", "drift_per_day_err", "rms", "n"}`, and with `--series N`
  a `series` of `t`, `y` and `fit`.
- `sigmatau cross` adds `"with"` (the second file), `"cross_adev": {"tau",
  "covar", "dev", "uncertainty", "dev_a", "dev_b", "neff"}` and
  `"cross_spectrum": {"nperseg", "nseg", "f", "sy", "floor", "sy_a", "sy_b"}`
  without the DC bin.
- On failure, with `--format json`, stdout holds
  `{"schema_version": 1, "sigmatau", "command", "source", "error": "…", "exit_code": n}`.

## Live measurement

`sigmatau live SOURCE` reads one reading at a time and keeps ADEV, MDEV,
HDEV and MHDEV up to date with [`LiveStability`](@ref), at constant work per
reading. `SOURCE` is

- `sim`, a simulated clock (white PM, white FM and random-walk FM) that gives
  one phase reading per τ₀ (`--fast` runs it flat out);
- `tcp://HOST:PORT`, an LXI instrument's SCPI socket (usually port 5025),
  sent `--query` (default `READ?`) before every reading;
- a path: a serial port set up with `stty` first, or a log file, read to its
  end. Without `--query` the instrument is expected to stream.

`tcp://` needs macOS or Linux. A reading of `9.91E+37` (SCPI's missed gate)
stops the run with exit code 2: a running estimate cannot bridge a gap.

The first number on each line is the reading: phase in seconds by default,
fractional frequency with `--freq`, or a counter frequency in Hz with `--hz`
(against `--nominal`, by default the first reading to 4 significant
figures). `--tau0` is the time between readings (default 1 s), `--devs`
the deviations (default all four) and `--m` or `--grid` the τ values.

```bash
sigmatau live tcp://192.168.1.20:5025 --hz --nominal 10e6 --tau0 1 --devs adev,mdev
stty -F /dev/ttyUSB0 9600 raw && sigmatau live /dev/ttyUSB0 --freq --tau0 0.1
```

Every `--every N` readings (default 1) it writes one JSON object on a line:

```
{"schema_version": 1, "sigmatau", "command": "live", "source",
 "quantity": "phase" | "frequency" | "hz", "tau0", "nominal", "n",
 "readings": [the readings since the last line],
 "tau": […], "results": [{"deviation": "adev", "dev": […], "neff": […]}, …],
 "done": false, "elapsed_ms"}
```

`dev` is `null` at τ values that do not have two terms yet. It stops after
`--count N` readings or when the source ends, with a last line that has
`"done": true`; until then it runs until stopped. The values carry no
confidence intervals: run `analyze` on the saved readings for those.
## Report

`sigmatau report` takes the same record and analysis options as `analyze`
and writes HTML instead of CSV or JSON (so `--format` does not apply):

```sh
sigmatau report clock.dat --devs adev,mdev,tdev --out clock-report.html
```

The file has no outside references: styles and SVG plots are inline. It
lists the record (samples, τ₀, span, gaps, outliers, drift, preprocessing),
the settings, one plot per axis (σ_y and σ_x) with the confidence band, a
table per estimator (τ, σ, interval, EDF, noise type, terms) and the command
that produced it.

## Exit codes

| Code | Meaning |
|:---|:---|
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

# Command line

`sigmatau-cli` runs a SigmaTau analysis from the shell and prints one JSON
document on standard output. It is the compiled analysis core behind the
[explorer](explorer.md), so it starts in about 30 ms with no Julia session,
and any language that can run a program and parse JSON can use it.

## Build

The program is compiled with `juliac --trim`; you need Julia 1.12 or newer and
JuliaC (`pkg> app add JuliaC`):

```bash
julia --project=app/web/cli -e 'using Pkg; Pkg.instantiate()'
juliac --output-exe sigmatau-cli --bundle app/web/cli/build --trim=safe --experimental app/web/cli
```

The executable is `app/web/cli/build/bin/sigmatau-cli`. Without a build,
`julia --project=app/web/cli -m SigmaTauCLI [flags]` runs the same code,
paying Julia's start-up on each call.

## Examples

```bash
# ADEV, MDEV, HDEV and MHDEV of a phase file, octave τ grid, 68.3 % intervals
sigmatau-cli --input clock.dat

# Frequency data with no time column, sampled every 10 s, at 95 %
sigmatau-cli --input freq.txt --kind frequency --time-col 0 --tau0 10 --confidence 0.95

# Clock in column 3 against the reference in column 4, values in ns
sigmatau-cli --input pair.csv --value-col 3 --reference-col 4 --scaling 1e-9

# Total estimators on a decade grid, with the spectrum and a dynamic ADEV map
sigmatau-cli --input clock.dat --devs totdev,mtotdev,htotdev --grid decade --spectrum --dynamic dadev

# A synthetic white FM + random-walk FM record, to try it out
sigmatau-cli --synth 8192 --seed 1
```

## Options

| Flag | Meaning | Default |
|:---|:---|:---|
| `--input PATH` | read a record with SigmaTau's readers | |
| `--synth N`, `--seed S` | synthetic white FM + random-walk FM record of N samples | seed 42 |
| `--kind phase\|frequency` | what the value column holds | `phase` |
| `--tau0 S` | sample interval in seconds; otherwise from the time column or the header | |
| `--time-col I`, `--value-col I` | 1-based columns; time column 0 means none | 1, 2 |
| `--reference-col J` | subtract this column from the value column (several clocks in one file) | none |
| `--header N` | header lines to skip; `#` lines and a Stable32 header are skipped anyway | 0 |
| `--scaling X` | multiply the values (1e-9 for ns) | 1 |
| `--detrend none\|mean\|endpoint\|linear` | remove a trend first | `none` |
| `--fillgaps` | put the record on an even time grid and fill gaps (needs a time column) | off |
| `--remove-outliers` | remove 5σ outliers before the analysis | off |
| `--devs a,b,…` | any of adev mdev tdev hdev mhdev htdev totdev mtotdev ttotdev htotdev mhtotdev theo1 theoh pdev mtie tierms; `none` reads and checks the record only | adev,mdev,hdev,mhdev |
| `--grid octave\|half-octave\|quarter-octave\|decade\|half-decade\|all` | τ grid | `octave` |
| `--confidence P` | confidence level of the intervals | 0.683 |
| `--no-ci` | skip noise identification and intervals | |
| `--no-bias-correction` | raw total estimators and Thêo1 | |
| `--series N` | at most N record points in the output | 4000 |
| `--spectrum` | add Welch S\_y(f) and S\_x(f) (Hann window, half overlap) | off |
| `--nperseg N` | Welch segment length | largest power of 2 ≤ N/4 |
| `--dynamic dadev\|dhdev` | add a dynamic ADEV or HDEV map | off |
| `--window W`, `--step S` | dynamic window and step, in samples | N/8, W/2 |

## Output

On success the exit code is 0 and the output looks like this (arrays
shortened):

```json
{
  "sigmatau": "0.5.0",
  "record": { "n": 16384, "tau0": 1.0, "kind": "phase", "source": "clock.dat",
              "series": { "t": [0.0, 4.0, ...], "x": [0.0, 1.2e-12, ...] },
              "check": { "samples": 16384, "span_s": 16384.0, "gaps": 0, "outliers": 2,
                         "outlier_t": [...], "outlier_x": [...], "drift_per_day": 3.1e-14 } },
  "options": { "grid": "octave", "ci": true, "confidence": 0.683, "bias_correction": true },
  "results": [
    { "dev": "adev", "axis": "sigma_y", "tau": [1.0, 2.0, ...], "sigma": [...],
      "ci_lo": [...], "ci_hi": [...], "edf": [...], "noise": ["WHFM", ...], "neff": [...] }
  ],
  "spectrum": { "nperseg": 4096, "f": [...], "sy": [...], "sx": [...] },
  "dynamic": { "dev": "adev", "window": 2048, "step": 1024, "t": [...], "tau": [...],
               "map": [[...], ...] },
  "elapsed_ms": 41.7
}
```

- `results` has one entry per estimator, in the order of `--devs`. `axis` is
  `sigma_x` (seconds) for tdev, htdev, ttotdev, mtie and tierms and
  `sigma_y` otherwise. `ci_lo`, `ci_hi`, `edf` and `noise` are empty with
  `--no-ci`.
- `check` is the record check: NaN gaps, points past 5σ (at most 500 listed)
  and the linear frequency drift per day.
- `spectrum` and `dynamic` appear only when asked for. The spectrum leaves
  out the DC bin and needs a record without gaps. In `dynamic`, `t` holds the
  window centres in seconds and `map` has one row per window, one column per τ.
- Non-finite numbers are written as `null`.

On failure the exit code is 2 and the output is `{"error": "message"}`.
Running with no arguments prints a usage line and exits with 1.

## From Python

```python
import json, subprocess

out = subprocess.run(
    ["sigmatau-cli", "--input", "clock.dat", "--devs", "adev,hdev"],
    capture_output=True, text=True,
)
doc = json.loads(out.stdout)
if out.returncode != 0:
    raise RuntimeError(doc["error"])
adev = doc["results"][0]
print(adev["tau"], adev["sigma"])
```

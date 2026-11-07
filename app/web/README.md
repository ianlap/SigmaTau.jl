# SigmaTau web explorer

A browser front end for SigmaTau over a small compiled analysis core. It is
the successor to the Makie explorer in `app/`: plain HTML, CSS and Plotly.js,
so it is easy to restyle and can later be wrapped in a desktop window (Tauri)
without changes.

```
app/web/
  cli/        SigmaTauCLI: the analysis core, built with juliac --trim
  ui/         static page: index.html, theme.css, app.css, js/
  serve.py    local server: serves ui/ and runs the core (stdlib Python)
```

## Run it

```bash
python3 app/web/serve.py --open          # http://127.0.0.1:9385
```

`serve.py` looks for the core in this order: `$SIGMATAU_CLI`, then
`app/web/cli/build/bin/sigmatau-cli`, then `julia --project=app/web/cli -m
SigmaTauCLI` (no build needed, but each call pays Julia's start-up).
Plotly loads from `ui/vendor/plotly.min.js` if present, else from jsDelivr;
`python3 app/web/serve.py --fetch-plotly` saves a local copy for offline use.

## Build the core

Julia 1.12 or newer, with JuliaC installed (`pkg> app add JuliaC`):

```bash
julia --project=app/web/cli -e 'using Pkg; Pkg.instantiate()'
juliac --output-exe sigmatau-cli --bundle app/web/cli/build --trim=safe --experimental app/web/cli
```

Measured on a 4-core Linux container (Julia 1.13.3): an 11.6 MB executable in a
128 MB bundle, 32 ms start-up, six deviations on 8192 points in 42 ms. A clean
`--trim=safe` build needs a patched FFTW.jl (its plan cleanup dispatches
dynamically); until that is upstreamed, build with `--trim=unsafe-warn`, whose
output was identical in testing. The synthetic record in the core avoids FFTW.

## Core interface

`sigmatau-cli` takes flags and prints one JSON document on stdout.

| Flag | Meaning | Default |
| --- | --- | --- |
| `--input PATH` | read a record with SigmaTau's readers | |
| `--synth N`, `--seed S` | synthetic white FM + random-walk FM record | seed 42 |
| `--kind phase\|frequency` | what the value column holds | `phase` |
| `--tau0 S` | sample interval; else from the file | |
| `--time-col I`, `--value-col I` | 1-based columns; time 0 = none | 1, 2 |
| `--header N` | non-comment header lines to skip | 0 |
| `--scaling X` | multiply values (1e-9 for ns) | 1 |
| `--detrend none\|mean\|endpoint\|linear` | | `none` |
| `--devs a,b,…` | any of adev mdev tdev hdev mhdev htdev totdev mtotdev ttotdev htotdev mhtotdev theo1 theoh pdev mtie tierms | adev,mdev,hdev,mhdev |
| `--grid octave\|half-octave\|quarter-octave\|decade\|half-decade\|all` | τ grid | `octave` |
| `--confidence P` / `--no-ci` | interval level, or none | 0.683 |
| `--no-bias-correction` | for the total estimators and Thêo1 | |
| `--series N` | at most N record points in the output | 4000 |

Success (exit 0):

```json
{
  "sigmatau": "0.5.0",
  "record": { "n": 16384, "tau0": 1.0, "kind": "phase", "source": "synthetic",
              "series": { "t": [0.0, 4.0, ...], "x": [0.0, 1.2e-12, ...] } },
  "options": { "grid": "octave", "ci": true, "confidence": 0.683, "bias_correction": true },
  "results": [
    { "dev": "adev", "axis": "sigma_y", "tau": [1.0, 2.0, ...], "sigma": [...],
      "ci_lo": [...], "ci_hi": [...], "edf": [...], "noise": ["WHFM", ...], "neff": [...] }
  ],
  "elapsed_ms": 41.7
}
```

`axis` is `sigma_x` (seconds) for tdev, htdev, ttotdev, mtie and tierms.
`ci_lo`, `ci_hi`, `edf` and `noise` are empty with `--no-ci`. Non-finite
numbers are written as `null`. Failure (exit 2): `{"error": "message"}`.

`serve.py` exposes it as `POST /api/analyze` with a JSON body
`{source: {type: "file", name, content} | {type: "synth", n, seed}, kind,
tau0, time_col, value_col, header, scaling, detrend, devs, grid, confidence
(null for none), bias_correction, series}` and returns the core's JSON.
`GET /api/health` reports whether a core was found.

## Front end

- `ui/theme.css` holds every color, font and size; restyle there.
- `ui/js/api.js` is the only code that talks to the core, and
  `ui/js/platform.js` the only code that opens or saves files. A Tauri
  shell replaces those two and keeps the rest.
- Pages: Data (clocks and the import dialog), Stability (estimators,
  options, σ(τ) plot, per-τ table), Spectrum and Dynamic (placeholders until
  the core returns those).
- Several clocks can be loaded; Stability shows them overlaid (one line
  style per clock) or side by side with shared axes.

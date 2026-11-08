# SigmaTau web explorer

A browser front end for SigmaTau over a small compiled analysis core: plain
HTML, CSS and Plotly.js, so it is easy to restyle, and the same pages run in
the Tauri desktop window in `app/desktop/`.

```
app/web/
  ui/         static page: index.html, theme.css, app.css, js/
  serve.py    local server: serves ui/ and runs the core (stdlib Python)

app/desktop/  the Tauri window: the same pages, no Python (see its README)
cli/          the analysis core: the sigmatau command line (see its README)
```

## Run it

```bash
python3 app/web/serve.py --open          # http://127.0.0.1:9385
```

`serve.py` looks for the core in this order: `$SIGMATAU_CLI`, then
`cli/build/bin/sigmatau`, then `julia --project=cli -m SigmaTauCLI` (no
build needed, but each call pays Julia's start-up).
Plotly loads from `ui/vendor/plotly.min.js` if present, else from jsDelivr;
`python3 app/web/serve.py --fetch-plotly` saves a local copy for offline use.

## The core

The pages run the `sigmatau` command line (`cli/`); build it as its README
says. Measured on a 4-core Linux container (Julia 1.13.3): a 15 MB
executable in a 132 MB bundle, 24 ms to start, four estimators on 16 384
points in under 70 ms.

`ui/js/core.js` builds the command and its flags, always with
`--format json`, and reshapes the answer (the versioned schema in
`docs/src/cli.md`) into the columns the pages plot. `serve.py` exposes it as
`POST /api/analyze` with the body

```json
{ "args": ["analyze", "--devs", "adev,mdev", "--format", "json"],
  "file": { "name": "clock.dat", "content": "0 0\n1 1.2e-9\n…" } }
```

It checks that the first argument is `check`, `analyze`, `spectrum` or
`dynamic` and that there is no `--out`, saves `file` (or nothing, for
`--synth`) in a temporary folder, puts its path after the command and
returns the binary's JSON: status 200 on success and 422 with
`{"error": …}` on failure. `GET /api/health` reports whether a core was
found. The desktop shell takes the same request.

## Front end

- `ui/theme.css` holds every color, font and size; restyle there.
- `ui/js/core.js` is the only code that talks to the core, and
  `ui/js/platform.js` the only code that opens or saves files. Inside the
  desktop window (`app/desktop`) those two call its commands instead.
- Layout: a side panel (clocks, estimators, options, Run analysis) and five
  pages. Data has the file card, the import dialog (pick the time and value
  columns, record type, units, τ₀; gaps flagged from the whole file), the
  record check and the phase plot with outliers marked. Stability has the
  σ(τ) plot with the per-τ table, or "Compare at τ" when several clocks are
  loaded. Spectrum shows S_y(f) and ℒ(f) (or S_x(f)) with log-binned
  averaging and power-law slope guides; carrier, view and averaging apply
  without a rerun. Dynamic shows a DADEV or DHDEV map (or a 3D waterfall),
  computed on its own when the page opens or its settings change; click a
  cell to pin it. With several clocks Spectrum overlays them and Dynamic
  draws one map per clock on a shared colour scale. Help explains the pages,
  estimators, noise types and file formats.
- With one clock the panel holds the current page's settings; with several
  it keeps the clock list and the page settings move to a toolbar.
- Analysis runs only on Run analysis (F5); changed inputs mark the results
  out of date in the panel and the status bar.
- Several clocks can be loaded, with no limit. With more than one, colour
  means clock (eight colours, then they repeat) and line style means
  estimator; the eye button hides a clock.
- One file can hold several clocks: in the import dialog each column is
  Time, Clock or Ignore, and every Clock column becomes a clock named from its
  header. "Measured against" records the shared reference. The Reference
  select in the Clocks list re-references the clocks of that file to the
  chosen one (x_i - x_R); clocks without a shared reference are greyed out,
  and the reference itself leaves the plots. With "Remember for files like
  this" (on by default) a file with the same extension and header loads
  without the dialog; Import settings… reopens it.
- Keys: F5 or Ctrl+Enter runs, 1 to 4 switch pages, ? opens Help, a file dropped anywhere opens the
  import dialog.

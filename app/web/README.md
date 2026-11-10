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
`docs/src/guide/cli.md`) into the columns the pages plot. `serve.py` exposes it as
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

The Live page runs `sigmatau live` for as long as a measurement lasts:
`POST /api/live/start` with `{"args": ["live", SOURCE, flags…]}` starts it
(stopping any other) and answers `{"id": n}`; `GET /api/live/poll?id=n&since=k`
answers `{"lines": [JSON lines from k on], "next", "running", "error"}`; and
`POST /api/live/stop` ends it. The desktop shell has the same three as the
commands `live_start`, `live_poll` and `live_stop`.

## Tests

```bash
SIGMATAU_EXE=cli/build/bin/sigmatau julia --project=cli cli/test/runtests.jl
SIGMATAU_CLI=cli/build/bin/sigmatau python3 -m unittest discover -s app/web/test -v
```

The first runs the CLI tests against the library and checks that the compiled
binary prints the same output; the second sends `serve.py` the page's
requests. The `Web core` workflow builds with `--trim=safe` and runs both on
every change to `src/`, `cli/` or `app/web/`.

## Front end

- `ui/theme.css` holds every color, font and size; restyle there.
- `ui/js/core.js` is the only code that talks to the core, and
  `ui/js/platform.js` the only code that opens or saves files. Inside the
  desktop window (`app/desktop`) those two call its commands instead.
- Layout: a side panel (clocks, estimators, options, Run analysis) and six
  pages. Data has the file card, the import dialog (pick the time and value
  columns, record type, units, τ₀; gaps flagged from the whole file), the
  record check and the phase plot with outliers marked. Stability has the
  σ(τ) plot with the per-τ table, or "Compare at τ" when several clocks are
  loaded. Spectrum shows S_y(f) and ℒ(f) (or S_x(f)) with log-binned
  averaging and power-law slope guides; carrier, view and averaging apply
  without a rerun. Dynamic shows a DADEV or DHDEV map (or a 3D waterfall),
  computed on its own when the page opens or its settings change; click a
  cell to pin it. With several clocks Spectrum overlays them and Dynamic
  draws one map per clock on a shared colour scale. Holdover shows the time
  error after a training stretch against holdover time (RMS, a chosen share
  of the trials, the worst trial, the noise model's prediction and √6·TDEV
  for reference) and the time within a limit, also computed on its own. Live reads an
  instrument (LXI socket, serial port or file) or the simulator and redraws
  σ(τ) and the frequency as readings arrive; "Analyse readings" loads them
  as a clock for the full analysis. Help explains the pages,
  estimators, noise types and file formats.
- With one clock the panel holds the current page's settings; with several
  it keeps the clock list and the page settings move to a toolbar (Live
  keeps its settings in the panel).
- The round button left of Help picks the theme: the system's, light, or
  dark (Espresso); the choice is remembered.
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
- Run analysis turns into Cancel while a run is going; F5, Ctrl+Enter or Esc
  stops it (F5 restarts it when the inputs changed meanwhile), and `serve.py`
  stops the core process when the page hangs up. F5 never reloads the page.
- Every clock can be renamed: in the import dialog (one name per Clock column)
  or in the name field on the Data page. Synthesize can make several clocks
  at once (`--clocks`, up to 8), each noisier than the last by a factor of two,
  all against one perfect reference; the Ensemble page uses them. The Synthesis
  page sets each clock's noise (clock-model coefficients or h_α, converted both
  ways) and previews the MDEV it predicts (`sigmatau model`).
- Keys: F5 or Ctrl+Enter runs, 1 to 7 switch pages, ? opens Help, a file dropped anywhere opens the
  import dialog.

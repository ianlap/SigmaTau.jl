# The explorer in a browser

The explorer is a set of static pages, in HTML, CSS and JavaScript with
Plotly, over the [`sigmatau`](../../cli/README.md) program, which does
every computation. `serve.py` serves the pages to a browser and runs the
program for them. The desktop app in [`app/desktop`](../desktop/README.md)
shows the same pages in a window of its own, without Python. The
[manual](https://ianlap.github.io/SigmaTau.jl/dev/explorer/) describes the
pages and how to work with them.

```text
app/web/
  ui/         the pages: index.html, theme.css, app.css, js/
  serve.py    the local server, on the Python standard library alone
  test/       tests of the server, of the import code and of the pages in Chromium
```

## Run

```bash
python3 app/web/serve.py --open          # http://127.0.0.1:9385
```

`serve.py` looks for the program in this order: the path in the
environment variable `SIGMATAU_CLI`, then `cli/build/bin/sigmatau`, then
`julia --project=cli -m SigmaTauCLI`, which needs no build and pays the
start-up of Julia on every call. Plotly is loaded from
`ui/vendor/plotly.min.js` when the file exists, and from jsDelivr
otherwise. `python3 app/web/serve.py --fetch-plotly` saves that file, for
use without a network.

## How the pages reach the program

`ui/js/core.js` builds each command with its options, always with
`--format json`, and turns the answer, in the
[JSON schema](https://ianlap.github.io/SigmaTau.jl/dev/reference/cli/) of
the program, into the columns that the pages plot. `serve.py` takes the
command as `POST /api/analyze` with the body

```json
{ "args": ["analyze", "--devs", "adev,mdev", "--format", "json"],
  "file": { "name": "clock.dat", "content": "0 0\n1 1.2e-9\n…" } }
```

The first argument must be `check`, `analyze`, `spectrum`, `dynamic`,
`drift`, `holdover`, `steer`, `timescale` or `model`, and `--out` is refused. The
server saves `file` in a temporary folder and puts its path after the
command. A request for a simulated record has no `file`, and `--synth N`
must then follow the command directly, so that no request names a path on
the machine of the server. The answer is the JSON of the program, with
status 200 when the command succeeded and 422 when it failed. A request of
another form is answered with status 400 and an `error` member.
`GET /api/health` reports whether a program was found. When the page
closes the connection, as its **Cancel** button does, the server stops the
program.

The server answers 403 to a `POST` whose content type is not
`application/json`, and on `127.0.0.1` to a request addressed to another
host name. A page of another site can therefore not call it through the
browser.

The Live page keeps `sigmatau live` running for as long as a measurement
lasts:

| Request | Effect |
|:---|:---|
| `POST /api/live/start` with `{"args": ["live", SOURCE, …]}` | Starts a measurement, stops any other, and answers `{"id": n}` |
| `GET /api/live/poll?id=n&since=k` | Answers `{"lines", "next", "running", "error"}`, where `lines` holds the JSON lines of the program from line `k` on |
| `POST /api/live/stop` | Ends the measurement |

A server that listens on another address than `127.0.0.1` starts a
measurement only from `sim` or a `tcp://` source: a device or a file would
be read on the machine of the server.

The desktop app answers the same requests through commands of its own.

## Tests

```bash
SIGMATAU_CLI=cli/build/bin/sigmatau python3 -m unittest discover -s app/web/test -v
node --test app/web/test/import.test.mjs app/web/test/help.test.mjs
SIGMATAU_CLI=cli/build/bin/sigmatau node --test app/web/test/explorer.e2e.mjs
```

The first sends `serve.py` the requests of the pages. The second tests
what the import dialog makes of the top of a file, and the addresses and
the search of the help. The third drives the pages in Chromium against the
compiled program, and needs Playwright (`npm install playwright`, then
`npx playwright install chromium`). The *Web core* workflow builds the
program with `--trim=safe` and runs all three, with the tests of `cli/`,
on every change to `src/`, `cli/` or `app/web/`.

## The code of the pages

- `ui/theme.css` holds every color, font and size, for the light and the
  dark theme. Restyle there.
- `ui/js/core.js` is the only module that runs the program, and
  `ui/js/platform.js` the only one that opens and saves files. Inside the
  desktop app the two call its commands in place of the server.
- `ui/js/state.js` holds the state of the application: one object, changed
  through `update()`, after which the listeners draw again.
- `ui/js/main.js` is the entry point: routing between the pages, loading
  records, the run, the export menu and the keyboard shortcuts.
- `ui/js/pages/` has one module for each page, `ui/js/panel.js` the side
  panel, `ui/js/importer.js` the import dialog and `ui/js/plot.js` the
  look that every plot shares.
- `ui/js/session.js` writes and reads session files, and `ui/js/report.js`
  writes the HTML report.
- `ui/js/help-content.js` is the help: the
  [manual](../../docs/src/explorer/), as `julia docs/explorer_help.jl`
  writes it. Change the pages of the manual and run that script. The build
  of the documentation stops when the file is behind them. `ui/js/help.js`
  finds and searches its sections, `ui/js/helpmenu.js` is the menu under
  **?** and `ui/js/pages/help.js` the page.

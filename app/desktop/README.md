# SigmaTau Explorer (desktop)

A Tauri window around the web explorer in `app/web/`. The window is the OS
webview (WebKit on macOS, WebView2 on Windows, WebKitGTK on Linux) showing
`app/web/ui/` unchanged, and a small Rust shell runs the compiled analysis
core in place of `serve.py`. The Python server stays as the development path
for working on the pages in a browser.

```
app/desktop/
  src/main.rs         commands: health, analyze, cancel, save_file
  tauri.conf.json     window, bundle; ships app/web/cli/build/ as core/
  capabilities/       webview permissions (core window events only)
  icons/              app-icon.png and the sizes generated from it
```

## Build

1. Build the core (see `app/web/README.md`), which leaves the juliac bundle in
   `app/web/cli/build/`.
2. Optional, for offline use: `python3 app/web/serve.py --fetch-plotly`
   puts Plotly into `ui/vendor/`, which the app then ships.
3. With Rust and the Tauri CLI (`cargo install tauri-cli --version "^2"`):

```bash
cd app/desktop
cargo tauri build      # target/release/bundle/macos/SigmaTau Explorer.app
cargo tauri dev        # run from the source tree
```

`cargo build` alone also works and gives a bare executable that finds the
core in the source tree (debug builds) or through `$SIGMATAU_CLI`.

## How it finds the core

`$SIGMATAU_CLI`, then `core/bin/sigmatau-cli` in the app's resources, then
(debug builds only) `app/web/cli/build/bin/sigmatau-cli`.

## Commands

The pages detect the shell through `window.__TAURI__` in `ui/js/core.js` and
`ui/js/platform.js`; nothing else in the front end knows about it.

| Command | Does |
| --- | --- |
| `health()` | `{ ok, core }`, as `GET /api/health` |
| `analyze(id, req)` | the request of `POST /api/analyze`; returns the core's JSON |
| `cancel(id)` | kills the core process of a running `analyze` (the pages abort superseded runs) |
| `save_file(name, data)` | save dialog, then writes the bytes; false if cancelled |

Each `analyze` writes the record to a fresh temporary folder, runs the core
with a 10-minute limit and removes the folder. Set `SIGMATAU_DESKTOP_LOG=1`
to print page-load times and each core call (flags, outcome, milliseconds)
on stderr.

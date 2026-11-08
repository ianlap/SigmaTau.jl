# SigmaTau Explorer (desktop)

A Tauri window around the web explorer in `app/web/`. The window is the OS
webview (WebKit on macOS, WebView2 on Windows, WebKitGTK on Linux) showing
`app/web/ui/` unchanged, and a small Rust shell runs the compiled analysis
core in place of `serve.py`. The Python server stays as the development path
for working on the pages in a browser.

```
app/desktop/
  src/main.rs         commands: health, analyze, cancel, save_file
  tauri.conf.json     window and bundle settings
  bundle.sh           cargo tauri build, then copies the core into the app
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
app/desktop/bundle.sh          # target/release/bundle/macos/SigmaTau Explorer.app
cd app/desktop && cargo tauri dev   # run from the source tree
```

`bundle.sh` copies the core in after Tauri bundles, keeping its symlinks
(Tauri's own resource copy turns them into files and the app grows from
114 MB to 289 MB), then signs the app ad hoc. It does the macOS step only;
on other systems put `app/web/cli/build/` next to the executable as `core/`.

Measured on an Intel MacBook Pro (macOS 15, Julia 1.13.3): a 5.7 MB shell,
a 114 MB app (30 MB zipped), the page loaded about 1.1 s after launch, and
the first analysis back 0.1 s after that.

`cargo build` alone also works and gives a bare executable that finds the
core in the source tree (debug builds) or through `$SIGMATAU_CLI`.

## How it finds the core

`$SIGMATAU_CLI`, then `core/bin/sigmatau-cli` in the app's resources
(`Contents/Resources/core` on macOS), then
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

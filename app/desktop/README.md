# The explorer as a desktop app

A [Tauri](https://tauri.app) window around the pages of
[`app/web`](../web/README.md). The window is the webview of the system,
WebKit on macOS and WebKitGTK on Linux, and it shows `app/web/ui/`
unchanged. A small Rust program runs the
[`sigmatau`](../../cli/README.md) program in place of `serve.py`.
[Install and run](https://ianlap.github.io/SigmaTau.jl/dev/explorer/install/)
in the manual says how to install a released app.

```text
app/desktop/
  src/main.rs          the commands that the pages call
  tauri.conf.json      the window and the bundle (tauri.linux.conf.json: Linux settings)
  bundle.sh            builds the app with the program inside, signs it and writes dist/
  entitlements.plist   hardened-runtime exceptions for a Developer ID build
  capabilities/        the permissions of the webview
  icons/               app-icon.png and the sizes made from it
```

## Build

1. Build the program as [`cli/README.md`](../../cli/README.md) says, which
   leaves it in `cli/build/`.
2. For use without a network, `python3 app/web/serve.py --fetch-plotly`
   saves Plotly to `app/web/ui/vendor/`, and the app then carries it.
3. With Rust and the Tauri CLI (`cargo install tauri-cli --version "^2"`):

```bash
app/desktop/bundle.sh               # writes app/desktop/dist/
cd app/desktop && cargo tauri dev   # runs from the source tree
```

On macOS `bundle.sh` makes `SigmaTau Explorer.app` and puts a `.dmg` and a
`.zip` of it in `dist/`. It copies the program into the app after Tauri
has bundled it, which keeps the symbolic links among the Julia libraries.
Tauri's own copy of resources turns the links into files and more than
doubles the size of the app. On Linux it copies Tauri's `.deb` and
`.AppImage` to `dist/`. There the program goes in as a Tauri resource
(`tauri.linux.conf.json`), so those packages carry the libraries twice.
The Linux build needs the development packages of WebKitGTK:
`libwebkit2gtk-4.1-dev libgtk-3-dev librsvg2-dev
libayatana-appindicator3-dev patchelf`.

`cargo build` alone gives a bare executable. It finds the program in the
source tree, in a debug build, or through `SIGMATAU_CLI`.

### Signing on macOS

Without a Developer ID the app is signed ad hoc. It runs on the Mac that
built it. On another Mac, macOS stops the first launch of a downloaded
copy, and the manual says how to allow it.

With a Developer ID, set `APPLE_SIGNING_IDENTITY` to the
`Developer ID Application: …` identity in the keychain. `bundle.sh` then
signs every library of the program and the app with the hardened runtime
(`entitlements.plist`). It also notarizes and staples the app when
`APPLE_ID`, `APPLE_PASSWORD` (an app-specific password) and `APPLE_TEAM_ID`
are set.

### Releases

`.github/workflows/DesktopRelease.yml` builds the app for macOS, on Apple
silicon and on Intel, and for Linux, with a `.tar.gz` of the `sigmatau`
program for each. On a release tag `v…`, the one TagBot pushes, it adds
the files to that GitHub release. A tag such as `app-v0.1.0` makes a draft
release of its own, and a run started by hand uploads the files as
artifacts of the run. The header of the workflow lists the secrets that
turn on Developer ID signing.

The launcher of the AppImage points `LD_LIBRARY_PATH` at the libraries of
the AppImage. The app therefore starts the program without that variable,
so that the program loads its own Julia runtime.

## How the app finds the program

In this order: the path in the environment variable `SIGMATAU_CLI`, then
`core/bin/sigmatau` among the resources of the app
(`Contents/Resources/core` on macOS), then, in a debug build only,
`cli/build/bin/sigmatau`.

## Commands

The pages detect the app through `window.__TAURI__` in `ui/js/core.js` and
`ui/js/platform.js`. No other code of the pages knows of it.

| Command | Effect |
|:---|:---|
| `health()` | Answers `{ ok, core }`, as `GET /api/health` does |
| `analyze(id, req)` | Takes the request of `POST /api/analyze` and returns the JSON of the program |
| `cancel(id)` | Stops the program of a running `analyze` |
| `save_file(name, data)` | Opens a save dialog and writes the bytes; `false` when the dialog is canceled |
| `open_url(url)` | Shows an `https` address in the browser of the system: the links from the help to the documentation. An address with a space, a quote or a character outside printable ASCII is refused |
| `live_start(args)`, `live_poll(id, since)`, `live_stop()` | The three requests of the Live page |

Each `analyze` writes the record to a new temporary folder, runs the
program with a limit of 10 minutes and removes the folder. With
`SIGMATAU_DESKTOP_LOG=1` the app prints the load time of the page and
every call of the program, with its arguments, its outcome and its
duration, on standard error.

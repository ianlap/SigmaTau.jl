# SigmaTau Explorer (desktop)

A Tauri window around the web explorer in `app/web/`. The window is the OS
webview (WebKit on macOS, WebView2 on Windows, WebKitGTK on Linux) showing
`app/web/ui/` unchanged, and a small Rust shell runs the compiled analysis
core in place of `serve.py`. The Python server stays as the development path
for working on the pages in a browser.

```
app/desktop/
  src/main.rs         commands: health, analyze, cancel, save_file
  tauri.conf.json     window and bundle settings (tauri.linux.conf.json adds Linux's)
  bundle.sh           builds the app with the core inside, signs it, writes dist/
  entitlements.plist  hardened-runtime exceptions for a Developer ID build
  capabilities/       webview permissions (core window events only)
  icons/              app-icon.png and the sizes generated from it
```

## Build

1. Build the core (see `cli/README.md`), which leaves the juliac bundle in
   `cli/build/`.
2. Optional, for offline use: `python3 app/web/serve.py --fetch-plotly`
   puts Plotly into `ui/vendor/`, which the app then ships.
3. With Rust and the Tauri CLI (`cargo install tauri-cli --version "^2"`):

```bash
app/desktop/bundle.sh               # writes app/desktop/dist/
cd app/desktop && cargo tauri dev   # run from the source tree
```

On macOS `bundle.sh` makes `SigmaTau Explorer.app` and puts a `.dmg` and a
`.zip` of it in `dist/`. It copies the core in after Tauri bundles, keeping
its symlinks (Tauri's own resource copy turns them into files and the app
grows from 114 MB to 289 MB). On Linux it copies Tauri's `.deb` and
`.AppImage` to `dist/`; there the core goes in as a Tauri resource
(`tauri.linux.conf.json`), so those packages carry the duplicated libraries.
The Linux build needs the WebKitGTK development packages
(`libwebkit2gtk-4.1-dev libgtk-3-dev librsvg2-dev libayatana-appindicator3-dev patchelf`).

### Signing (macOS)

Without a Developer ID the app is signed ad hoc. It runs on the Mac that
built it; on another Mac, macOS blocks the first launch of a downloaded copy.
Right-click the app and choose Open, then Open again (on macOS 15, open it
once, then allow it under System Settings → Privacy & Security), or clear the
download flag:

```bash
xattr -dr com.apple.quarantine "/Applications/SigmaTau Explorer.app"
```

With a Developer ID, set `APPLE_SIGNING_IDENTITY` to the
"Developer ID Application: …" identity in your keychain. `bundle.sh` then
signs every library of the core and the app with the hardened runtime
(`entitlements.plist`), and notarizes and staples it when `APPLE_ID`,
`APPLE_PASSWORD` (an app-specific password) and `APPLE_TEAM_ID` are set too.

### Releases

`.github/workflows/DesktopRelease.yml` builds macOS (Apple silicon and Intel)
and Linux packages, and a `.tar.gz` of the `sigmatau` command line for each.
On a release tag `v…` (the one TagBot pushes) it adds the files to that
GitHub release; a tag like `app-v0.1.0` makes a draft release of its own; a
manual run uploads them as workflow artifacts, listed in the run's summary.
The workflow's header lists the secrets that turn on Developer ID signing.

The AppImage runs on most x86_64 distributions; it needs FUSE 2
(`sudo dnf install fuse` on Fedora, `sudo apt install libfuse2` on Ubuntu),
or run it with `--appimage-extract-and-run`. Its launcher points
`LD_LIBRARY_PATH` at the AppImage's own libraries, so the app starts the
core without it and the core loads its own Julia runtime.

Measured on an Intel MacBook Pro (macOS 15, Julia 1.13.3): a 5.7 MB shell,
a 114 MB app (30 MB zipped), the page loaded about 1.1 s after launch, and
the first analysis back 0.1 s after that.

`cargo build` alone also works and gives a bare executable that finds the
core in the source tree (debug builds) or through `$SIGMATAU_CLI`.

## How it finds the core

`$SIGMATAU_CLI`, then `core/bin/sigmatau` in the app's resources
(`Contents/Resources/core` on macOS), then
(debug builds only) `cli/build/bin/sigmatau`.

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
to print page-load times and each core call (arguments, outcome, milliseconds)
on stderr.

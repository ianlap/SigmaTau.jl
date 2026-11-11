# Install and run

The explorer is a set of pages and the [`sigmatau`](../guide/cli.md)
program, which does every computation. The desktop app holds both. From a
checkout of the repository, a short script serves the pages to a browser
and runs the program for them.

## Desktop app

A release carries the app for macOS, on Apple silicon and on Intel, and for
Linux on x86-64:

| System | File | Install |
|:-------|:-----|:--------|
| macOS | `SigmaTau-Explorer-<version>-macos-<arch>.dmg` | Open it and drag the app to Applications |
| Linux | the `.deb` package | `sudo apt install ./<file>.deb` |
| Linux | the `.AppImage` | `chmod +x <file>.AppImage`, then run it |

The app needs no Julia, no Python and no network. It takes over 100 MB on
disk, most of it the Julia runtime of the program inside.

macOS stops the first launch of a downloaded copy that Apple has not
notarized. Right-click the app and choose **Open**, then **Open** again. On
macOS 15, open it once and then allow it under System Settings, Privacy &
Security. Clearing the download flag does the same:

```bash
xattr -dr com.apple.quarantine "/Applications/SigmaTau Explorer.app"
```

The AppImage needs FUSE 2 (`sudo dnf install fuse` on Fedora,
`sudo apt install libfuse2` on Ubuntu), or the option
`--appimage-extract-and-run`.

## In a browser, from the repository

With Python 3 and a checkout of the repository:

```bash
python3 app/web/serve.py --open
```

The script serves the pages at `http://127.0.0.1:9385` and opens them in the
default browser. It uses the Python standard library alone.

| Option | Effect |
|:-------|:-------|
| `--open` | Open the pages in the default browser |
| `--port N` | Serve on port `N` in place of 9385 |
| `--host ADDRESS` | Listen on another address in place of `127.0.0.1` |
| `--fetch-plotly` | Save the Plotly library to `app/web/ui/vendor` and exit |

The script looks for the program in this order:

1. the path in the environment variable `SIGMATAU_CLI`;
2. `cli/build/bin/sigmatau`, where the
   [build](../guide/cli.md#Install) leaves it;
3. `julia --project=cli -m SigmaTauCLI`, which needs Julia 1.13 or later and
   the packages of `cli` installed
   (`julia --project=cli -e 'using Pkg; Pkg.instantiate()'`).

The third way needs no build, but every request then waits for Julia to
start and to compile the program, about half a minute each time. With the
compiled program the first analysis of the explorer, four deviations with
the spectra and the noise model of 16384 samples, returns in a fraction of
a second.

!!! warning "The script has no access control"
    Whoever reaches its address can run the program on this computer, on a
    record sent with the request or on a simulated one, and can start a
    measurement. Keep the default address, `127.0.0.1`, which only this
    computer reaches, unless everyone on the network is trusted. On another
    address the [Live](live.md) page reads from the simulator and from the
    network only, and not from a serial port or a file, so that no request
    reads a file of this computer.

The pages draw with [Plotly](https://plotly.com/javascript/), which they
load from a content delivery network unless `--fetch-plotly` has saved a
copy. Without a network and without that copy there are no plots.

## Building the desktop app

With the program built in `cli/build`, Rust and the Tauri command line
(`cargo install tauri-cli --version "^2"`):

```bash
python3 app/web/serve.py --fetch-plotly    # optional: plots without a network
app/desktop/bundle.sh                      # writes app/desktop/dist/
```

On macOS the script writes the app with a `.dmg` and a `.zip` of it, and on
Linux a `.deb` and an `.AppImage`. The Linux build needs the WebKitGTK
development packages. `app/desktop/README.md` in the repository lists them
and covers signing and notarization on macOS.

`cd app/desktop && cargo tauri dev` runs the app from the source tree
without bundling it.

# Using the explorer

The SigmaTau explorer is a point-and-click front end for the whole workflow:
load a record, check it, pick estimators, and read σ(τ), the spectra and the
dynamic deviation, for one clock or several. It runs the same library code as
the REPL, compiled into a small analysis core, so its numbers are SigmaTau's
numbers.

It comes in two forms that show the same pages:

- **Desktop app** (`app/desktop/`): a native window built with Tauri. No
  Julia, Python or browser needed once it is built.
- **Browser** (`app/web/`): the pages served by a small Python script on
  `http://127.0.0.1:9385`. This is the quickest way to run it from a
  checkout and the development path for the pages.

## Run it in a browser

From a checkout of the repository, with Python 3:

```bash
python3 app/web/serve.py --open
```

`serve.py` needs the analysis core. If you have not built it (see below), it
falls back to running the core through Julia, which works without a build
but pays Julia's start-up on every analysis. `--port N` changes the port and
`--fetch-plotly` saves a local copy of Plotly for offline use.

## Build the analysis core

The pages run the [`sigmatau` command line](cli.md) with `--format json`, so
the core is that program, compiled with `juliac --trim`. You need Julia 1.12
or newer and JuliaC (`pkg> app add JuliaC`):

```bash
julia --project=cli -e 'using Pkg; Pkg.instantiate()'
juliac --output-exe sigmatau --bundle cli/build --trim=safe --experimental cli
```

The build takes a few minutes and leaves a bundle of about 130 MB in
`cli/build/`. `serve.py` looks for the core in `$SIGMATAU_CLI`, then
`cli/build/bin/sigmatau`, then falls back to `julia --project=cli -m
SigmaTauCLI`. The core starts in about 25 ms and runs four estimators on
16 384 points in under 70 ms.

## Build the desktop app

With the core built, Rust, and the Tauri CLI
(`cargo install tauri-cli --version "^2"`):

```bash
app/desktop/bundle.sh               # macOS: target/release/bundle/macos/SigmaTau Explorer.app
cd app/desktop && cargo tauri dev   # or run it from the source tree
```

On macOS the app is about 114 MB (30 MB zipped) and opens in about a second.
`app/desktop/README.md` covers the other platforms and how the window finds
the core.

## A first analysis

1. **Data.** Open a phase or frequency file, or drop one anywhere in the
   window. The import dialog shows the first lines of the file: mark each
   column as Time, Clock or Ignore, and set the record type (phase or
   frequency), the units and τ₀ if there is no time column. Or synthesize a
   test record.
2. **Check the record.** The record check lists the samples, span, gaps,
   points past 5σ and the frequency drift, next to the phase plot with
   outliers marked. Detrending, gap filling and outlier removal are on the
   same page.
3. **Pick estimators.** In the side panel, **+ Add estimator** opens the
   list; [Choosing an estimator](choosing.md) explains each one. Choose the
   τ grid and the confidence level under Options.
4. **Run analysis** (F5 or Ctrl+Enter). Nothing recomputes on its own: after
   you change the record, the estimators or the options, the Run button and
   the status bar turn amber until you run again.

## The pages

| Page | Shows |
|:---|:---|
| **Data** | The file, the import settings, the record check and the phase plot. |
| **Stability** | σ(τ) on log-log axes with confidence intervals as bars or a band. With one clock, the per-τ table (τ, deviation, interval, noise type, EDF); with several, "Compare at τ" with the most stable clock in bold. Click the plot to move τ. |
| **Spectrum** | S\_y(f) and either ℒ(f) at a carrier you choose or S\_x(f), Welch estimates with log-binned averaging and power-law slope guides. Changes here apply without a rerun. |
| **Dynamic** | A DADEV or DHDEV map σ(t, τ), or a 3D waterfall: does the stability change during the record? Recomputes when its settings change; click a cell to pin it. |
| **Live** | Readings from an instrument or the simulator as they arrive, with the running σ(τ); see below. |
| **Help** | The pages, estimators, noise types, file formats and keyboard shortcuts. |

The **File** menu saves the current plot as PNG or SVG and every number as
TSV or JSON. **Report (HTML)** writes one self-contained page with the
analysis settings, each record (source file and its SHA-256, samples, τ₀,
span, gaps, outliers, drift, preprocessing), the σ_y and σ_x plots, a table
per clock and estimator, and the `sigmatau` command lines that reproduce
it; print it from a browser for a PDF.

**Save session** writes the clocks, their file contents, import settings and
analysis settings to one JSON file. **Open session…**, or dropping a session
file on the window, replaces what is loaded and reruns the analysis.

## Several clocks

Add clocks with **+ Add clock**, or mark several columns of one file as Clock
in the import dialog; there is no limit. With more than one clock, colour
means clock and line style means estimator, and the eye button hides a
clock. Stability and Spectrum overlay the clocks; Dynamic draws one map per
clock on a shared colour scale.

Clocks read from one file usually share the measurement's reference
("Measured against" in the import dialog). **Reference** in the Clocks list
re-references them to one of the clocks, x\_i − x\_R, and the reference itself
leaves the plots. This is how you look at clock B against clock A when both
were measured against a common house reference.

## Files

Text files with one sample per line: `.dat`, `.txt`, `.csv` (commas), `.tsv`
(tabs), anything else split on spaces. Lines starting with `#` are comments;
Stable32 headers ending in `# Header End` are skipped and their `Tau:` line
gives τ₀. `NaN` marks a missing sample. Import settings are remembered for
files with the same extension and header (turn this off with "Remember for
files like this"); **Import settings…** on the Data page reopens the dialog.

## Keyboard shortcuts

| Key | Action |
|:---|:---|
| F5 or Ctrl+Enter | Run the analysis |
| 1, 2, 3, 4, 5 | Data, Stability, Spectrum, Dynamic, Live |
| ? | Help |
| Esc | Close a dialog |

On a plot, hover for values, drag to zoom and double-click to reset.

## Part of a record, and comparing runs

On the Data page, **Select range** turns a drag across the phase plot into
the analysis range: Stability, Spectrum and Dynamic then use only that
stretch (the record check and the plot still show the whole record, with the
range shaded). **whole record** clears it. Each clock keeps its own range.
**τ range** in the panel's options leaves out averaging times outside the
limits.

**Pin run** on the Stability page keeps the current curves on the plot,
thin and dashed, so the next run (another range, drift removal or τ grid)
can be compared with it. Hover a pin for the settings it was made with; ×
removes it.

## Live measurement

The **Live** page runs [`sigmatau live`](cli.md) and plots what it streams.
In the panel, pick the source (the simulator, an instrument's SCPI socket
by address, or a serial port or log file by path), what each reading is
(phase in seconds, fractional frequency, or a counter frequency in Hz
against a nominal), τ₀ and the estimators (ADEV, MDEV, HDEV, MHDEV), then
**Start**. The upper plot shows the fractional frequency so far, the lower
one σ(τ) updated at every reading. **Stop** ends the run; **Load the
readings as a clock** then adds them to the Data page and runs the full
analysis with confidence intervals. The running estimates carry no
intervals.

**Fit** draws the power-law noise model, S_y(f) = Σ h_α f^α fitted to MHDEV
(see [`powerlaw_fit`](@ref)), as a dotted line beside each ADEV, MDEV, TDEV,
HDEV and MHDEV curve. The panel under the per-τ table lists the h_α and the
model's ADEV at round averaging times (1 s to one day), the numbers a
datasheet quotes.

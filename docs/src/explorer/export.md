# Exports, sessions and reports

The **File** menu saves the plots on screen, the numbers of the page in
view, a report of the analysis and the whole session. In a browser each file is a
download. The desktop app asks where to save it.

| Entry | File | Holds |
|:------|:-----|:------|
| **Plot as PNG** | `sigmatau-<page>.png` | The plots of the page, at twice their size on screen |
| **Plot as SVG** | `sigmatau-<page>.svg` | The same plots as vector graphics |
| **… as TSV** | `sigmatau-<table>.tsv` | The numbers of the page in view, as a table ([Results](#Results)) |
| **… as JSON** | `sigmatau-<table>.json` | Everything the analysis core returned for them |
| **Stability report (HTML, prints to PDF)** | `sigmatau-report.html` | The analysis of the Stability page as one page that can be printed |
| **Save session…** | `sigmatau-session.json` | The clocks, their files and the settings |
| **Open session…** | | Opens a saved session |

## Plots

The image is the plot as it is on screen, with its zoom and with any pinned
runs. In the light theme the margin around the axes is transparent. In the
dark theme the image keeps the dark card behind the plot.

On a page with several plots, such as the two spectra of the Spectrum page
or the maps of several clocks on the Dynamic page, one image holds every
plot that is shown, placed as on screen. In a browser, the camera button
that appears on hovering over a plot saves that plot alone as a PNG.

## Results

The two results entries name what they write, which follows the page in
view. Data, Synthesis, Live and Help have no results of their own, and the
entries are greyed out there.

| Page | Entry | File | One line for each |
|:-----|:------|:-----|:------------------|
| Stability | **Stability results** | `sigmatau-results` | clock, estimator and averaging time |
| Spectrum | **Spectra** | `sigmatau-spectrum` | clock and Fourier frequency |
| Dynamic | **Dynamic maps** | `sigmatau-dynamic` | clock, window and averaging time |
| Holdover | **Holdover table** | `sigmatau-holdover` | clock and holdover time |
| Steering | **Steering table** | `sigmatau-steering` | clock and averaging time |
| Ensemble, Time scale | **Time scale table** | `sigmatau-timescale` | averaging time of the scale, then clock and averaging time |
| Ensemble, N-cornered hat | **Hat results** | `sigmatau-hat` | clock, estimator and averaging time |

The TSV file is a table for a spreadsheet or a script, with tabs between
the columns. Those of Stability:

| Column | Holds |
|:-------|:------|
| `clock` | The name of the clock |
| `estimator` | The estimator, such as `adev` |
| `tau_s` | The averaging time in seconds |
| `sigma` | The deviation |
| `ci_lo`, `ci_hi` | The limits of the confidence interval. Empty without intervals |
| `edf` | The equivalent degrees of freedom |
| `noise` | The noise type, such as `WFM` |
| `neff` | The number of terms in the estimate |

The JSON file is the output of [`sigmatau analyze`](../reference/cli.md#JSON)
for each clock, with the name of the clock added: the record check, the
settings of the run, the deviations, the spectra, the noise model and the
mask checks. With several clocks it is a list. Both files hold every clock
with results, hidden clocks included.

The tables of the other pages hold the clocks the page shows, with the
columns of the command line's CSV output for the same command
([Command-line reference](../reference/cli.md)), units in the names of the
columns:

| Table | Columns |
|:------|:--------|
| Spectra | `clock`, `f_hz`, `sy`, `sx`, and `L_dbc_hz` when the lower plot shows ℒ(f) |
| Dynamic maps | `clock`, `deviation`, `t_s` (the center of the window), `tau_s`, `dev` |
| Holdover table | `clock`, `train_s`, `tau_s`, `rms_s`, `quantile_s`, `worst_s`, `mean_s`, `trials`, `limit_s`, `time_to_limit_s`, `claimed_s`, and the `reference_s` and `predicted_s` curves of the plot |
| Steering table | `clock`, `tau_s`, `free`, `steered`, `predicted` |
| Time scale table | `clock`, `against`, `tau_s`, `adev`, `weight`, `weight_mean`, `error_rms_s`: the scale against the reference, then each clock against the scale or the reference as plotted |
| Hat results | `clock`, `estimator`, `tau_s`, `sigma`, `ci_lo`, `ci_hi`, `edf`, `neff`, for each clock separated |

Their JSON files hold the output of `holdover`, `steer`, `dynamic`,
`timescale` or `hat` that the page drew, with the clocks named as on the
page. A table not computed yet, or computed for other settings, is not
written: the page computes it once it is in view.

## The report

The report is a single HTML file with its figures inside, which a browser
prints to PDF. It records what was analyzed and how:

| Section | Holds |
|:--------|:------|
| Analysis | The estimators, the τ grid and limits, the confidence level, the bias correction, the noise type and the reference clock |
| Records | For each clock: the file and its SHA-256 hash, the record type, the samples, the sample interval, the span, the gaps, the outliers, the drift and the preparation |
| Stability | The plots of the Stability page for the frequency and the time estimators, on a white background |
| Results | A table for each clock and estimator: averaging time, deviation, interval, EDF, noise type and number of terms |
| Reproduce | The `sigmatau analyze` command of each clock |

The report states the settings on screen, so it is refused while the
results are out of date: run the analysis first. Hidden clocks are left
out.

The commands under Reproduce name each file without its folder and print
CSV. Run in the folder of the file with the same version of the program,
they return the deviations of the report. The hash identifies the file that
was analyzed. A gzip file is named and hashed as unpacked, the form in
which the program reads it.

## Sessions

**Save session…** writes one JSON file that holds:

- every clock, with the content of its file, its import settings, its
  preparation, its analysis range and its noise model if built by hand;
- the reference clock;
- the estimators and options of the side panel;
- the display settings of the Stability, Spectrum and Dynamic pages and
  the settings of the Holdover and Steering pages.

It holds no results. **Open session…**, or the file dropped on the window,
puts these clocks and settings in the place of those on screen and runs the
analysis. Pinned runs and the settings of the Synthesis, Ensemble and Live
pages are not saved.

A session file contains the records themselves. It is as large as they are,
and it hands the whole analysis to whoever receives it.

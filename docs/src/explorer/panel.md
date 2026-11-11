# The side panel

The side panel holds what an analysis is made of: the clocks, the
estimators and the options, with the **Run analysis** button below them.
The Data, Stability, Spectrum, Dynamic, Holdover and Steering pages share
it. The Synthesis, Ensemble and Live pages show their own settings in its
place.

With one clock, the Spectrum, Dynamic, Holdover and Steering pages put
their settings in the panel too, in place of the estimators and options.
With several clocks the panel keeps its content, and those four pages show
their settings in a bar above the plots.

## Clocks

The list of clocks appears when more than one is loaded. Clocks are added
by **+ Add clock**, which opens a file, by marking several columns of one
file as clocks in the [import dialog](data.md#The-import-dialog), by
simulating several on the [Data](data.md#Simulated-records) or
[Synthesis](synthesis.md) page, and by loading the readings of a
[live measurement](live.md).

Each row shows the color of the clock on every plot, its name, where it
comes from and the number of its samples. There are eight colors, used
again from the ninth clock on.

| Control | Action |
|:--------|:-------|
| A click on the row | Shows that clock on the Data page |
| ◉ | Hides the clock: it leaves the plots, the tables and the report. A run still analyzes it, and the TSV and JSON exports hold its results. ○ shows it again |
| × | Removes the clock. Removing one of several clocks that were simulated together simulates the others again |

### Reference

Clocks read from one file were measured against one reference at the same
instants. **Reference** exchanges that reference for one of the clocks: the
record of every other clock of the file becomes its readings less those of
the reference clock, ``x_i - x_R``, and the reference clock leaves the
plots. This shows clock B against clock A when both were measured against a
third. **As measured** returns to the records as read.

Clocks from other files share neither the reference nor the instants. They
cannot be chosen as the reference, and they stay as measured when one is
chosen.

## Estimators

The list holds the estimators that **Run analysis** computes. With one
clock each estimator has a color. With several clocks the color is that of
the clock, and each estimator has a line style. × removes an estimator.

**+ Add estimator** opens the list of all sixteen, with a field that
searches their names and descriptions:

| Group | Estimators |
|:------|:-----------|
| Allan family | ADEV, MDEV, HDEV, MHDEV, PDEV |
| Total | TOTDEV, MTOTDEV, HTOTDEV, MHTOTDEV |
| Thêo | THEO1, THEOH |
| Time | TDEV, HTDEV, TTOTDEV, MTIE, TIE rms |

The first three groups are fractional frequency deviations ``σ_y(τ)`` and
the last is time deviations ``σ_x(τ)``, in seconds. The Stability page
shows one kind at a time.
[Estimators and noise types](estimators.md) says what each one gives, which
of them need a record without gaps and which take seconds on a long record.

## Options

The options apply to every clock and every estimator. With several clocks
they are folded under the heading **Options**, which shows the grid and the
confidence level.

| Option | Choices | Effect |
|:-------|:--------|:-------|
| **τ grid** | Octave, Half octave, Quarter octave, Decade, Half decade, All τ | The averaging times: 1, 2, 4, … sample intervals for the octave grid, and every one that the record allows for All τ |
| **Confidence intervals** | On or off, at 68.3 %, 90 %, 95 % or 99 % | The level of the [confidence intervals](../explanation/confidence.md). Off gives the deviations alone, without intervals, EDF or noise types |
| **τ range (s)** | A smallest and a largest averaging time | Leaves out the averaging times outside the limits. An empty field is no limit |
| **Bias correction** | Exact, Published tables, As Stable32, None | The [bias correction](../explanation/bias.md) of TOTDEV, MTOTDEV, HTOTDEV, MHTOTDEV, TTOTDEV and Thêo1 |
| **Noise type** | Identified at each τ, From the fitted model, or one of the seven noise types | The noise type behind the intervals and the bias factors |

The explorer starts with the octave grid, intervals at 68.3 %, the exact
bias correction and the noise type identified at each τ.

The four bias corrections are those of the keyword `correct_bias` of the
estimators: **Exact** is `true`, **Published tables** is `:published`,
**As Stable32** is `:stable32` and **None** is `false`.

**Identified at each τ** takes the noise type from
[noise identification](../explanation/noise_identification.md). **From the
fitted model** takes the strongest term, at each averaging time, of the
power-law model fitted to the record, the one that **Fit** shows on the
[Stability](stability.md#The-noise-model) page. A stated type holds for
every averaging time, for a record whose noise is known.

## Running

**Run analysis** computes every estimator in the list for every clock but
the reference, together with the spectra, the noise model and the mask
checks. F5 and Ctrl+Enter do the same. The button turns amber when the
record, the estimators or the options have changed since the last run, and
it reads **Cancel** during a run. The line below it gives the number of
estimators, the number of samples or clocks and the time the last run took.

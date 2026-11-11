# Command-line reference

```text
sigmatau COMMAND [FILE] [options]
```

A call runs one command on one file. An option that takes a value is
written `--tau0 1` or `--tau0=1`. `sigmatau help COMMAND`, or `--help`
after a command, prints the options of that command.
[Command line](../guide/cli.md) shows the common tasks.

The samples of output on this page are written by the program when the
documentation is built. They show its numbers to four significant digits,
where the program writes them in full.

```@setup cliref
using SigmaTau
include(joinpath(pkgdir(SigmaTau), "docs", "transcripts.jl"))
```

## Commands

| Command | Result |
|:---|:---|
| `analyze FILE` | deviations against averaging time, each with its confidence interval, EDF and noise type |
| `check FILE` | a summary of the record: samples, span, gaps, outliers and drift |
| `spectrum FILE` | the Welch spectra ``S_y(f)`` and ``S_x(f)`` |
| `dynamic FILE` | the Allan or Hadamard deviation in a window that slides along the record |
| `drift FILE` | the fit of a frequency drift model |
| `holdover FILE` | the time error of the clock left to run on a frequency learned from the record, against holdover time |
| `steer FILE` | the Allan deviation of the record steered to its reference in replay, beside that of the record and of the steady state that the clock model predicts |
| `cross FILE --with FILE2` | the cross Allan deviation and cross-spectrum of two channels |
| `timescale FILE --value-cols LIST` | a time scale from several clocks read against one reference |
| `model --synth N --noise-model R:Q1:Q2:Q3` | the deviations that a clock model predicts |
| `report FILE` | the analysis of `analyze` as one HTML file |
| `live SOURCE` | deviations kept up to date while readings arrive |
| `completion SHELL` | a completion script for `bash`, `zsh`, `fish` or `powershell` |
| `version` | the version of the package |
| `help [COMMAND]` | the commands, or the options of one |

## The record

Every command but `live` and `model` reads a record, and takes the options
of this section.

### Reading

| Option | Meaning |
|:---|:---|
| `--input-format FMT` | layout of the file: `columns`, `stable32`, `timelab`, `counter`, `table`, `bipm_utc` or `rinex_clk`. By default `auto`, which detects it as [`read_data`](@ref) does |
| `--phase`, `--freq` | the samples are phase in seconds, or fractional frequency. A file of columns is phase unless `--freq` is given. The other layouts tell which they hold, and either flag overrides them |
| `--tau0 SECONDS` | sample interval. By default the one in the header of the file, else the step of the time column |
| `--time-col N`, `--value-col N` | columns of the time and of the samples, counted from 1. By default 1 and 2. `--time-col 0` reads a file that has no time column |
| `--reference-col N` | column to subtract from the value column, for clocks measured against one reference in one file. For files of columns only |
| `--header N` | number of lines to skip after the leading `#` lines, such as a line of column names |
| `--delim CHAR` | character between the fields, `'\t'` for a tab. By default a comma where the first line of data has one, a tab where a table or a `.tsv` file has one, and blanks otherwise |
| `--scaling X` | factor applied to the samples, `1e-9` for nanoseconds |
| `--column NAME` | the column of a table with named columns, or the clock of a RINEX clock file, by name |
| `--nominal HZ` | nominal frequency, for readings in hertz. By default their mean, rounded to four significant figures |
| `--keep-zeros` | read the zeros of a Stable32 `.phd` or `.frd` file as samples. By default a zero there is a missing sample, as in Stable32, except the first and the last sample of a phase file |

A file of columns is read as [`read_phase`](@ref) or
[`read_frequency`](@ref) reads it: lines that begin with `#` are skipped,
and a `# tau0:` line among them gives the sample interval. The other
layouts are read as [`read_data`](@ref) reads them. A table whose first
line names its columns is read by name, and a unit in the name of the
column is applied, so `phase_ns` needs no `--scaling`. With `--time-col`,
`--value-col` or `--header` given, such a table is read by position
instead, and `--header 1` skips its line of names. The name of the time
column still sets its unit, so a column named `MJD` is in days either way.

The program holds every record as phase. It integrates a frequency file of
``N`` samples to ``N + 1`` phase samples, and reports that number. A
frequency file with missing samples must be filled with `--fillgaps`,
whatever the analysis.

### Preparing

| Option | Meaning |
|:---|:---|
| `--fillgaps` | fill the missing samples, as [`fillgaps`](@ref) does: the readings that the time column shows to be absent, and those written as `NaN` |
| `--range A:B` | keep seconds `A` to `B`, counted from the first sample. Either end can be left open: `3600:` or `:86400` |
| `--detrend MODE` | remove `mean`, `endpoint` or `linear` from the samples as the file holds them, by the methods of [`detrend`](@ref). By default `none` |
| `--remove-outliers` | replace the outliers and fill their places, as [`remove_outliers`](@ref) does |
| `--remove-drift MODE` | fit a line (`linear`) or a parabola (`quadratic`) to the phase by least squares and remove it. By default `none` |

The steps apply in the order of the table. `--detrend linear` removes a
frequency offset from a phase file and a linear frequency drift from a
frequency file. `--remove-drift` works on the phase whatever the file
holds: `linear` removes a frequency offset, and `quadratic` an offset and a
linear frequency drift. Its fit leaves out the missing samples and those
that `--fillgaps` put in, and the JSON output reports what it removed.

### Simulated records

| Option | Meaning |
|:---|:---|
| `--synth N` | in place of `FILE`, a simulated record of `N` phase samples: white FM of ``10^{-12}`` at one sample interval, and a random walk of the frequency in steps of ``10^{-14}`` |
| `--seed S` | seed of the simulation. By default 42 |
| `--noise-model R:Q1:Q2:Q3` | simulate a clock of this [`ClockModel`](@ref) instead: the variance `R` of its white PM in s², and the diffusion coefficients of its white FM, random-walk FM and random-run FM. Given several times, it simulates several clocks |
| `--clocks K` | simulate `K` clocks, from 2 to 64, with white FM and random-walk FM, each a factor of 2 above the last in deviation |
| `--clock I` | which of several simulated clocks the command reads. By default 1 |

The sample interval of a simulated record is `--tau0`, or 1 s.

## Output

| Option | Meaning |
|:---|:---|
| `--format FMT` | `csv`, the default, or `json` |
| `--out PATH` | write to a file. By default the output goes to standard output |
| `--series N` | include samples of the record in the JSON output, about `N`, for a plot |
| `--quiet` | no note on standard error that the file was written |

A value that does not exist is an empty field in CSV and `null` in JSON,
never the text `NaN`. `model` and `live` write JSON only, and `report`
writes HTML.

### JSON

Each answer is one object on one line. Its `schema_version` is 1, and
within a version members are added but not renamed or removed. The samples
on this page are indented, and `…` stands for the rest of an array and for
the members left out. This is the answer of `analyze`:

```@example cliref
sigmatau_json("analyze --synth 1024 --devs adev,tdev --m 1,8 --series 2 --format json") # hide
```

`record` describes the record after its preparation. Every command that
reads a record gives it, except `timescale`:

| Member | Meaning |
|:---|:---|
| `quantity` | `"phase"` or `"frequency"`, what the file holds |
| `n`, `tau0`, `span_s` | number of phase samples, sample interval, and their product, in seconds |
| `start` | UTC time of the first sample, where the file gives one, else `""` |
| `gaps` | number of missing samples |
| `outliers`, `outlier_t`, `outlier_x` | number of outliers that [`find_outliers`](@ref) flags by default, and the time and phase of the first 500 |
| `drift_per_day` | least-squares slope of the fractional frequency, per day |
| `removed_offset`, `removed_drift_per_day` | the frequency offset and the drift that `--remove-drift` removed |
| `series` | with `--series N`, the times and the phase of the samples kept: every sample of a record of `N` or fewer. Of a longer record, the first sample, the last, and from each of `(N − 2) ÷ 3` equal stretches the lowest sample, the highest and the first missing one, which is `null`, so that a plot spans the record and shows every outlier and every gap |

`analysis` repeats the settings of the call, for `analyze`, `check`,
`spectrum`, `dynamic` and `drift`. In it, `tau_mode` is the name of the
grid, or `"m"` when `--m` gave the factors, which `m` then lists, and
`range` holds the two ends of `--range`.

Each entry of `results` is one deviation. Its `axis` is `"sigma_x"` for
the deviations in seconds (TDEV, HTDEV, TTOTDEV, MTIE and TIE rms) and
`"sigma_y"` for the others.

When a command fails with `--format json`, the object holds `error` and
`exit_code` in place of the results. `live` and `model`, which write JSON
only, answer this way without the option:

```@example cliref
sigmatau_json("analyze missing.dat --format json") # hide
```

## `analyze`

| Option | Meaning |
|:---|:---|
| `--devs LIST` | deviations, separated by commas. By default `adev,mdev,hdev,mhdev` |
| `--grid GRID` | grid of averaging factors: `octave` (the default), `half-octave`, `quarter-octave`, `decade`, `half-decade` or `all` |
| `--m LIST` | the averaging factors themselves, `1,2,5,10`, in place of a grid |
| `--tau-min S`, `--tau-max S` | leave out the averaging times below or above `S` seconds |
| `--confidence P` | confidence level of the intervals. By default 0.683 |
| `--no-ci` | no intervals, EDF or noise types |
| `--noise TYPE` | noise type behind the intervals and the bias factors: `auto` (the default, identified at each averaging time), `fit` (the dominant term at each averaging time of the model that [`powerlaw_fit`](@ref) fits to the record), or one of `WHPM`, `FLPM`, `WHFM`, `FLFM`, `RWFM`, `FWFM`, `RRFM` for every averaging time |
| `--bias MODE` | bias correction of the total deviations and Thêo1: `exact` (the default), `published`, `stable32` or `none`, for `correct_bias = true`, `:published`, `:stable32` and `false` ([Bias correction](../explanation/bias.md)). `--no-bias-correction` is `--bias none` |

The deviations are `adev`, `mdev`, `hdev`, `mhdev`, `pdev`, `totdev`,
`mtotdev`, `htotdev`, `mhtotdev`, `theo1`, `theoh`, `tdev`, `htdev`,
`ttotdev`, `mtie` and `tierms`, each the function of that name
([Estimators](estimators.md)). A grid stops at the longest averaging time
of each deviation. The factors of `--m` go to every deviation, and one that
a deviation cannot reach gives a row without a value. The record needs
eight samples or more.

The CSV output has a row for each deviation and averaging time:

```text
source,quantity,deviation,tau,dev,ci_lower,ci_upper,edf,neff,noise_type
```

`tau` is in seconds, and `neff` is the number of terms averaged.

Four options add a member to the JSON output and need `--format json`:

| Option | Member |
|:---|:---|
| `--spectrum` | `spectrum`, as the `spectrum` command gives it. `--nperseg N` sets the segment length. For a record with gaps the member is `{"error": "…"}`, and the deviations are computed all the same |
| `--fit` | `powerlaw`, the power-law model that [`powerlaw_fit`](@ref) fits to the modified Hadamard deviation of the record |
| `--component NOISE:A:B` | `powerlaw`, built one noise type at a time with [`noise_fit`](@ref): `NOISE`, a name such as `WHFM` or a value of ``α``, is fitted between `A` and `B` seconds to what the components before it leave. Either end can be left open. Repeat the option for each component, in order |
| `--mask NAME` | `masks`, the limits of an ITU-T wander mask with the check of each MTIE and TDEV result against them. `NAME` is one of the [`TELECOM_MASKS`](@ref) (`G811`, `G813_opt1`, `G813_opt2`, `G8262_opt1`, `G8262_opt2`, `G823_PRC`, `G823_SSU`, `G823_SEC`, `G823_PDH`), or `all` |

```@example cliref
sigmatau_json("analyze --synth 1024 --component WHFM::8 --component RWFM:32: --format json"; only = "powerlaw") # hide
```

`alphas`, `h` and `noise_types` are the model
``S_y(f) = \sum_α h_α f^α``. `mode` is `"fit"` with `--fit` and
`"components"` with `--component`, which adds `components` and `residual`,
the deviation that the components leave. `curves` holds the deviations of
the model at the averaging times of each ADEV, MDEV, TDEV, HDEV and MHDEV
result. `sigma_at` is its Allan deviation at 1 s, 10 s, 100 s, 1000 s,
``10^4`` s and one day, as far as half the record. `clock_model` is the
[`clock_model`](@ref) of the fit. Where the record or a component cannot
be fitted, `powerlaw` holds an `error` instead.

```@example cliref
sigmatau_json("analyze --synth 1024 --devs mtie,tdev --mask G811 --format json"; only = "masks") # hide
```

`mtie` and `tdev` are the limits of the mask, in seconds, at the averaging
times `tau`, from 0.1 s to ``10^6`` s, and `null` outside the range of the
recommendation. Each entry of `checks` is one result against the mask, as
[`mask_check`](@ref) returns it.

## `check`

`check` takes the options of the record and of the output. Its CSV output
is one row:

```text
source,quantity,n,tau0,span_s,gaps,outliers,drift_per_day
```

The columns are the members of `record` in the JSON output, which holds
nothing else. The record needs three samples or more.

## `spectrum`

| Option | Meaning |
|:---|:---|
| `--nperseg N` | length of the segments. By default the largest power of two up to a quarter of the record, and at least 16 |

The spectra are those of [`psd_y`](@ref) and [`psd_x`](@ref), with a Hann
taper and segments that overlap by half. The record needs 16 samples or
more and no gaps. The CSV columns are `f,sy,sx`, without the row of zero
frequency.

```@example cliref
sigmatau_json("spectrum --synth 1024 --format json"; only = "spectrum") # hide
```

## `dynamic`

| Option | Meaning |
|:---|:---|
| `--dev NAME` | `dadev`, the default, or `dhdev` |
| `--window N` | length of the window, in samples. By default the largest power of two up to an eighth of the record, and at least 16 |
| `--step N` | step of the window, in samples. By default half the window |
| `--grid GRID` | grid of averaging factors within the window. By default `octave` |

The result is that of [`dadev`](@ref) or [`dhdev`](@ref). The CSV columns
are `t,tau,dev`, with a row for each window and averaging time. `t` is the
center of the window, in seconds from the first sample of the file.

```@example cliref
sigmatau_json("dynamic --synth 1024 --format json"; only = "dynamic") # hide
```

`dev` has a row for each window and a column for each averaging time.

## `drift`

| Option | Meaning |
|:---|:---|
| `--model MODEL` | `linear` (the default), `quadratic`, `log` or `diffusion` |

The fit is that of [`drift_fit`](@ref), which defines the models and their
coefficients. The CSV output is one row:

```text
source,quantity,model,a,a_err,b,b_err,c,c_err,drift_per_day,drift_per_day_err,rms,n
```

`a`, `b` and `c` are the coefficients, with `c` empty for the linear model,
and `n` is the number of frequency samples fitted.

```@example cliref
sigmatau_json("drift --synth 1024 --series 2 --format json"; only = "drift") # hide
```

With `--series`, `series` holds the frequency samples kept and the fit at
their times.

## `holdover`

| Option | Meaning |
|:---|:---|
| `--train SECONDS` | length of the training stretch, on which the frequency is fitted. By default a quarter of the record |
| `--grow K` | in place of `--train`, a training stretch `K` times as long as each holdover time. `K` is 2 or more, and 3 or more with `--fit-drift` |
| `--tau LIST` | holdover times in seconds: `60,3600,14400`. By default the octave multiples of ``τ_0`` up to the training stretch, or to half of what the record has left after it where that is less. With `--grow` they are those that leave half the record for the trials, and with `--kalman` those up to half of what the record has left after the stretch |
| `--fit-drift` | fit a parabola to the stretch, which also carries the frequency drift on. By default a line |
| `--no-anchor` | start each trial from the value of the fit at the end of the stretch. By default from the last reading |
| `--kalman` | predict with a Kalman filter in place of the least-squares fit. It goes with none of `--grow`, `--no-anchor` and `--fit` |
| `--limit SECONDS` | also give the holdover time at which the error reaches this limit |
| `--confidence P` | share of the trials that `quantile` covers. By default 0.95 for this command |
| `--fit` | with JSON, add the error that the power-law model fitted to the record predicts for the same test |

The test is that of [`holdover`](@ref), which [Holdover](holdover.md)
describes. The training stretch and the holdover times must be whole
multiples of the sample interval, and the record needs 16 samples or more.
The CSV output has a row for each holdover time:

```@example cliref
sigmatau("holdover clock.dat --train 512 --tau 16,64,256 --limit 2e-10") # hide
```

`rms`, `quantile`, `worst` and `mean` are the statistics of the time error
over the `n` trials, in seconds. With `--grow`, `train` is `K` times the
holdover time of the row. `limit` is empty without `--limit`.
`time_to_limit` is the holdover time at which `quantile` reaches the limit,
the same on every row: 0 when `quantile` exceeds the limit at the first
holdover time, and empty when it does not reach it at the last. `claimed`
is empty without `--kalman`.

```@example cliref
sigmatau_json("holdover clock.dat --train 512 --tau 16,64,256 --limit 2e-10 --fit --format json"; only = "holdover") # hide
```

`drift` and `anchor` tell which fit and which start the test used. `grow`
is the `K` of `--grow`, with `train` then `null`, and 0 without it. `limit`
is `null` without `--limit`, and `time_to_limit` when the limit is not
reached.

`reference` is ``\sqrt{6}`` times the time deviation at each holdover
time, and ``\sqrt{20}`` times the Hadamard time deviation with
`--fit-drift`. By their definitions these are the RMS errors of predicting
the next ``τ``-average of the phase from the two before it with a line,
and from the three before it with a parabola. The test of `--grow 2`, and
of `--grow 3` with `--fit-drift`, is the counterpart on single readings,
with the line or the parabola fitted by least squares. An entry is `null`
where the record is too short for the deviation at that ``τ``.

`predicted` is there with `--fit`, from the
[`HoldoverPrediction`](@ref) of the model that [`powerlaw_fit`](@ref) fits
to the record: `sigma` and `quantile` are its RMS error and its quantile at
each holdover time, and `budget` holds, for each noise type of `alphas`, the
variance in s² that it adds at each holdover time. When the model cannot be
fitted, `predicted` holds an `error` message.

With `--kalman` the predictor is the Kalman filter of a
[`ClockModel`](@ref), as in `holdover(data, model)`. The filter passes once
over the record, and from every sample past the first training stretch its
estimate is carried ``τ`` ahead without further readings and compared with
the reading there. The model is the [`clock_model`](@ref) of the power-law
model that [`powerlaw_fit`](@ref) fits to the first training stretch alone,
so its noise levels are not learned from the readings that the predictions
are compared with. It holds white FM and random-walk FM, with white PM as
the noise of the readings, and it keeps the flicker FM and the flicker PM
that the fit finds as
[Markov states](../explanation/clock_model.md#Flicker-states), over the
averaging times at which each matters, up to the length of the record.
With `--fit-drift` it also has a drift state.

```@example cliref
sigmatau("holdover clock.dat --train 512 --tau 16,64,256 --kalman") # hide
```

`claimed` is the standard deviation that the filter itself gives the error
at each holdover time. The JSON output has it as an array, with a member
`model` that holds the `R`, `q1`, `q2`, `q3`, `flicker_fm` and
`flicker_pm` of the model, and with `anchor` `false`. `fit_tau_max` is the
longest averaging time of the fit that the model comes from, between an
eighth and a quarter of the training stretch. At a longer holdover time
`claimed` extrapolates the model
([With a Kalman filter as the predictor](@ref)). The two columns agree
where the model fits the clock, and an `rms` above `claimed` says that the
clock has more noise at that holdover time than the model.

## `steer`

| Option | Meaning |
|:---|:---|
| `--time-constant SECONDS` | the gains of a critically damped loop with this time constant, by pole placement |
| `--integral` | with `--time-constant`, a third gain, on the running sum of the phase |
| `--cost P:F:U` | in place of `--time-constant`, the gains that minimize the sum of `P`·phase², `F`·frequency² and `U`·steer². `P` and `U` are above 0 |
| `--interval SECONDS` | time between steers, a whole multiple of the sample interval. By default every sample |
| `--train SECONDS` | length of the first stretch, to which the clock model is fitted. By default a quarter of the record |
| `--fit-drift` | give the clock model a drift state, which the loop steers out |
| `--grid GRID`, `--m LIST` | the averaging times, as for `analyze` |

The record is the clock running free against its reference, and the
command replays it through the loop that [Steering](steering.md)
describes. One of `--time-constant` and `--cost` is required:
[`steering_gains`](@ref) or [`lqg_gains`](@ref) gives the gains, and
[`steer`](@ref) the steered record. The clock model is that of
`holdover --kalman`, fitted to the first stretch alone, so its noise levels
are not learned from the readings that the loop is judged on. The loop runs
over the whole record, and every result is taken from the end of that
stretch on. The record needs 64 samples or more, and the stretch must
leave 16 on either side of its end. The CSV output has a row for each
averaging time:

```@example cliref
sigmatau("steer clock.dat --time-constant 64 --m 1,16,256") # hide
```

`free` and `steered` are the Allan deviations of the record and of the
steered record. `predicted` is that of the steady state that
`steer(model, gains; readings=true)` gives for the clock model: the
readings of the settled loop, with their noise, as the record has it. It
is empty when the steers come less often than every sample or the model
has random-run FM.

```@example cliref
sigmatau_json("steer clock.dat --time-constant 64 --m 1,16,256 --format json"; only = "steer") # hide
```

`gains` holds the phase gain in s⁻¹, the frequency gain and the integral
gain in s⁻¹, as in [`SteeringGains`](@ref). `time_constants`, in seconds,
and `frequencies`, in Hz, are those of the poles of the loop, from
[`steering_poles`](@ref). `drift` tells whether the model has a drift
state, and `model` holds its coefficients as for `holdover --kalman`.
`phase_rms` is the RMS of the steered phase in seconds, and `steers_rms`
that of the steers, in fractional frequency. `predicted` holds the standard
deviations `phase`, `frequency` and `steers` and the `adev` of that
[`SteeringPrediction`](@ref), `phase` and `adev` being those of the
readings, or an `error` message where there is none.
With `--series N`, a member `series` holds the times `t`, the `steered`
phase and the frequency `correction` in force, at up to `N` points over the
whole record.

## `cross`

| Option | Meaning |
|:---|:---|
| `--with FILE2` | the second channel, read with the same options as `FILE` |
| `--grid GRID`, `--m LIST` | averaging factors, as for `analyze` |
| `--nperseg N` | length of the segments of the cross-spectrum, as for `spectrum` |

The two records must have the same length and sample interval, 16 samples
or more and no gaps ([Two-channel measurements](../guide/cross.md)). The
CSV output is the result of [`cross_adev`](@ref), with a row for each
averaging time:

```text
tau,covar,dev,uncertainty,dev_a,dev_b,neff
```

The JSON output holds `with`, the name of the second file, the same result
as `cross_adev`, and the result of [`cross_psd_y`](@ref):

```@example cliref
sigmatau_json("cross a.dat --with b.dat --freq --format json"; only = "cross_spectrum") # hide
```

## `timescale`

| Option | Meaning |
|:---|:---|
| `--value-cols LIST` | columns of the clocks, two or more, each read against the same reference: `2,3,4`. In a table with named columns each is read by its number as `--column` reads one by its name, with the unit in its name |
| `--algorithm ALG` | `at1` (the default), `at2`, `ewfa`, `greenhall`, `jones_tryon` or `stein`, the algorithms of [`timescale`](@ref) |
| `--wmax W` | largest weight that one clock may take, from 0 to 1. By default 1, which sets no limit |
| `--steps MODE` | with `at2`, detect frequency steps of the clocks: `best` or `real` |
| `--pivot K` | with `greenhall`, `jones_tryon` or `stein`, the clock at which the scale is read. By default 1 |
| `--flicker` | with `greenhall` or `jones_tryon`, fit flicker FM as well, and keep the flicker FM and the flicker PM of each clock as [Markov states](../explanation/clock_model.md#Flicker-states) of its model. By default flicker FM is not fitted and flicker PM counts as noise of the readings |

Each column is read with the options of the record. With `--synth` the
clocks are simulated: four of them, or those that `--clocks` or
`--noise-model` describe. `greenhall`, `jones_tryon` and `stein` need three
clocks or more, and fit a model to each with [`clock_models`](@ref). The
record needs 16 samples or more.

The CSV output has a row for each clock:

```text
clock,weight,weight_mean,error_rms
```

`weight` is the weight at the last sample. `weight_mean` and `error_rms`
are the mean weight and the root mean square of the prediction error over
the second half of the record.

```@example cliref
sigmatau_json("timescale --synth 1024 --series 2 --format json"; only = "timescale") # hide
```

`adev.scale` is the Allan deviation of the scale against the reference, and
`adev.clocks` that of each clock against the scale. `series` holds the
scale and, for each clock, its weights at the times kept. With `--flicker`
the object has a member `flicker`, whose arrays `fm` and `pm` give the
number of Markov states in the model of each clock.

## `model`

`model` simulates nothing and reads no file. It writes, as JSON, the
modified Allan, Allan and modified Hadamard deviations that
[`powerlaw_sigma`](@ref) gives for each `--noise-model`, at the octave
averaging times of a record of `--synth N` samples. Besides these two
options it takes `--tau0`, by default 1 s, `--out` and `--quiet`:

```@example cliref
sigmatau_json("model --synth 64 --noise-model 1e-22:1e-24:1e-30:0"; only = "model") # hide
```

With random-run FM in the model (`Q3` above zero), `mdev` and `adev` are
`null`: only the Hadamard deviations converge on that noise.

## `report`

`report` takes the options of the record and those of `analyze` in its
first table, and `--out`. It writes one HTML file: the summary of the
record and its preparation, the settings, a plot of the deviations in
fractional frequency and one of those in seconds, each with its confidence
band, a table for each deviation and the command line that made the
report. Styles and plots are inside the file.

## `live`

| Option | Meaning |
|:---|:---|
| `--phase`, `--freq`, `--hz` | the readings are phase in seconds (the default), fractional frequency, or frequencies in hertz |
| `--nominal HZ` | nominal frequency for `--hz`. By default the first reading, rounded to four significant figures |
| `--tau0 SECONDS` | time between readings. By default 1 s |
| `--devs LIST` | deviations, of `adev`, `mdev`, `hdev` and `mhdev`. By default all four |
| `--grid GRID`, `--m LIST` | averaging factors, as for `analyze` |
| `--query CMD` | query sent before each reading. By default `READ?` for a `tcp://` source and none for a path |
| `--every N` | write a line every `N` readings. By default 1 |
| `--count N` | stop after `N` readings. By default the run ends with the source |
| `--fast` | run the simulator without waiting a sample interval for each reading |
| `--seed S` | seed of the simulator. By default 42 |

`SOURCE` is `sim`, `tcp://HOST:PORT` or a path
([Live measurement](../guide/cli.md#Live-measurement)). Each line of output
is one object:

```@example cliref
sigmatau_json("live sim --count 4 --every 4 --fast --m 1,2 --devs adev") # hide
```

`readings` holds the readings since the line before, `n` their number so
far, and `dev` is `null` at an averaging time that does not have two terms
yet. The last line has `"done": true`. Without `--count`, a grid reaches
32768 sample intervals.

## Exit codes

| Code | Meaning |
|:---|:---|
| 0 | success |
| 1 | usage: an unknown command or option, a value that cannot be used, no file |
| 2 | input: a file that is missing or cannot be read, a column that it does not have, a source that cannot be reached |
| 3 | analysis: a record too short, or with gaps where the analysis takes none |
| 4 | output: the file of `--out` cannot be written |

An error prints one line that begins with `sigmatau:` on standard error. A
reader that closes the pipe early, as in `sigmatau analyze clock.dat | head`,
is not an error.

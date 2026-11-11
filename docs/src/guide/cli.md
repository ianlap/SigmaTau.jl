# Command line

`sigmatau` is the package as one compiled program. It reads a record from a
file, runs one analysis and writes CSV or JSON. It needs no Julia and
starts in a few tens of milliseconds. It is compiled from the code of the
Julia functions, and the tests require its results to equal theirs within
a relative difference of ``10^{-12}``. The
[explorer](../explorer/index.md) runs this program for every analysis it
shows.

This page covers the common tasks. The
[Command-line reference](../reference/cli.md) lists every command and
option, the columns of the CSV output, the JSON schema and the exit codes.

```@setup cli
using SigmaTau
include(joinpath(pkgdir(SigmaTau), "docs", "transcripts.jl"))
```

## Install

A release carries the program for Linux and macOS, as
`sigmatau-<version>-<linux|macos>-<arch>.tar.gz`:

```bash
tar -xzf sigmatau-*.tar.gz
xattr -dr com.apple.quarantine sigmatau-*/     # macOS: clear the download flag
sigmatau-*/bin/sigmatau version
```

The folder holds the program and the Julia runtime libraries it loads, about
150 MB in all. Keep it whole, and put its `bin` directory on the `PATH`.

To build the program from the repository, with Julia 1.13 or later:

```bash
julia -e 'using Pkg; Pkg.Apps.add("JuliaC")'      # installs ~/.julia/bin/juliac
julia --project=cli -e 'using Pkg; Pkg.instantiate()'
juliac --output-exe sigmatau --bundle cli/build --trim=safe --experimental cli
cli/build/bin/sigmatau version
```

Without a build, `julia --project=cli -m SigmaTauCLI analyze clock.dat`
runs the same commands, at the cost of Julia's start-up on each call.

## A first analysis

Every command that reads a file takes `--synth N` in its place: a
simulated record of `N` phase samples at 1 s, white FM with random-walk FM
below it. It serves to try a command before there is a file:

```@example cli
sigmatau("analyze --synth 4096 --devs adev,mdev --m 1,10,100") # hide
```

There is a row for each deviation and averaging time: the deviation, the
limits of its confidence interval, the EDF, the number of terms averaged
and the noise type identified. Without `--devs` the deviations are
`adev,mdev,hdev,mhdev`, and without `--m` the averaging factors are the
octave grid of each one. The program writes its numbers in full. The
transcripts on this page cut them to four significant digits.

## Your own file

```bash
sigmatau analyze clock.dat                   # time in column 1, phase in column 2
sigmatau analyze clock.dat --tau0 0.1        # the sample interval, stated
sigmatau analyze freq.txt --freq             # fractional frequency
sigmatau analyze log.txt --value-col 3 --scaling 1e-9        # column 3, in ns
sigmatau analyze bare.txt --time-col 0 --value-col 1 --tau0 1    # no time column
sigmatau analyze trace.tim                   # a TimeLab file
sigmatau analyze gnss.clk --column G08       # one clock of a RINEX clock file
```

A file of plain columns holds phase in seconds unless `--freq` says
otherwise. The sample interval is `--tau0`, else the one in the file's
header, else the step of the time column. The files of instruments, and
tables whose first line names the columns, are recognized and read as
[`read_data`](@ref) reads them, with the kind of record and the unit taken
from the file.

`check` reads the record the same way and summarizes it, which shows
whether the file was understood before any analysis:

```@example cli
sigmatau("check clock.dat") # hide
```

This file has 4096 samples at 1 s, of which the time column skips ten.

## Preparing the record

The steps of [Gaps, outliers and drift](preparing.md) are options of every
command:

```bash
sigmatau analyze clock.dat --fillgaps                  # fill the missing samples
sigmatau analyze clock.dat --remove-outliers           # replace the outliers
sigmatau analyze clock.dat --remove-drift quadratic    # remove offset and linear drift
sigmatau analyze clock.dat --range 3600:86400          # seconds 3600 to 86400 only
sigmatau analyze clocks.dat --value-col 3 --reference-col 2   # column 3 less column 2
```

The Allan and Hadamard families, MTIE and TIE rms skip gaps in a phase
record. The other deviations need the record filled first, and say so:

```@example cli
sigmatau("analyze clock.dat --devs adev,totdev --m 1,100") # hide
```

```@example cli
sigmatau("analyze clock.dat --devs adev,totdev --m 1,100 --fillgaps") # hide
```

`--remove-drift quadratic` fits a parabola to the phase by least squares
and removes it, which takes out a frequency offset and a linear frequency
drift. `linear` removes the offset alone. The fit leaves out the samples
that `--fillgaps` put in.

A frequency file is integrated to phase as it is read, so one with missing
samples needs `--fillgaps` for every command.

## Deviations, averaging times and intervals

```bash
sigmatau analyze clock.dat --devs adev,mdev,tdev,totdev
sigmatau analyze clock.dat --grid quarter-octave       # a finer grid
sigmatau analyze clock.dat --m 1,2,5,10,20,50,100      # the averaging factors
sigmatau analyze clock.dat --tau-min 10 --tau-max 1000 # part of the grid
sigmatau analyze clock.dat --confidence 0.95           # 95 % intervals
sigmatau analyze clock.dat --noise WHFM                # the noise type, stated
sigmatau analyze clock.dat --bias stable32             # Stable32's bias convention
sigmatau analyze clock.dat --no-ci                     # deviations alone
```

`--noise` and `--bias` are the keywords `noise` and `correct_bias` of the
estimators, with `exact` for `true` and `none` for `false`
([Confidence intervals](../explanation/confidence.md),
[Bias correction](../explanation/bias.md)). The intervals are at 68.3 %
unless `--confidence` says otherwise.

## Other analyses

```bash
sigmatau spectrum clock.dat --nperseg 4096             # S_y(f) and S_x(f)
sigmatau dynamic clock.dat --dev dhdev --window 4096   # deviation in a sliding window
sigmatau drift ocxo.dat --freq --model log             # a model of frequency drift
sigmatau drift clock.dat --method bisection            # a linear drift by the three-point fit
sigmatau holdover gpsdo.dat --train 86400 --limit 1e-6 # time error without a reference
sigmatau steer maser.dat --time-constant 86400         # the record steered to its reference
sigmatau cross a.dat --with b.dat                      # one device, two channels
sigmatau timescale clocks.dat --value-cols 2,3,4       # a time scale from three clocks
sigmatau report clock.dat --out clock.html             # an HTML report
sigmatau phasenoise ocxo.txt --carrier 10e6            # σ(τ) from a table of ℒ(f)
```

`report` takes the options of `analyze` and writes one HTML file with the
record, the settings, the plots, a table for each deviation and the command
line that made it. The file refers to nothing outside itself and prints to
PDF from a browser.

`phasenoise` reads a table of offset frequency and ``ℒ(f)`` in place of a
record, and integrates it to the deviations over the band of the table
([Power-law noise](../explanation/noise.md#A-measured-spectrum)). Each
line of its answer tells how much the noise outside the table could add.

## Live measurement

`sigmatau live SOURCE` keeps the Allan, modified Allan, Hadamard and
modified Hadamard deviations up to date while readings arrive
([`LiveStability`](@ref)). The source is one of:

- `sim`, a simulated clock that gives one phase reading each sample
  interval, or as fast as they compute with `--fast`;
- `tcp://HOST:PORT`, the SCPI socket of an instrument on the network,
  usually port 5025. The query `READ?` is sent before each reading, or the
  one that `--query` gives. This source needs macOS or Linux;
- a path: a serial port, set up with `stty` first (`stty -F` on Linux,
  `stty -f` on macOS), or a file of readings, read to its end. Without
  `--query` the instrument is expected to send readings on its own.

```bash
sigmatau live sim --count 1000 --fast
sigmatau live tcp://192.168.1.20:5025 --hz --nominal 10e6 --tau0 1 --devs adev,mdev
stty -F /dev/ttyUSB0 9600 raw && sigmatau live /dev/ttyUSB0 --freq --tau0 0.1
```

The first number on each line is the reading: phase in seconds, fractional
frequency with `--freq`, or a counter's frequency in hertz with `--hz`,
taken against `--nominal`. `--tau0` is the time between readings, 1 s by
default. The program writes one line of JSON for every reading, or every
`--every N` readings, and stops after `--count N` readings or when the
source ends. The values carry no confidence intervals. Run `analyze` on
the saved readings for those.

A reading of `9.91E+37`, with which a counter reports a gate without a
measurement, stops the run: a running estimate cannot leave a sample out.

## In scripts

`--format json` gives every answer as one JSON object in a versioned schema,
and `--out PATH` writes to a file:

```bash
sigmatau analyze clock.dat --format json --out result.json
sigmatau analyze clock.dat --format json |
    jq -r '.results[0].points[] | [.tau, .dev] | @tsv'      # τ and σ, tab-separated
```

The exit code is 0 on success, 1 for a mistake in the command line, 2 for
a file that cannot be read, 3 for a record that the analysis cannot use and
4 for output that cannot be written. An error prints one line on standard
error, and with `--format json` also an object with `error` and `exit_code`
on standard output, so that a program reading the output always gets JSON.

## Shell completion

`sigmatau completion SHELL` prints a script that completes the commands,
the options and their values, and leaves file names to the shell:

```bash
# bash, for this session
source <(sigmatau completion bash)
# bash, for every session
sigmatau completion bash > /usr/local/etc/bash_completion.d/sigmatau
# zsh
sigmatau completion zsh > "${fpath[1]}/_sigmatau"
# fish
sigmatau completion fish > ~/.config/fish/completions/sigmatau.fish
# PowerShell, in the profile
sigmatau completion powershell | Out-String | Invoke-Expression
```

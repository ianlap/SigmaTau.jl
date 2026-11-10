# Coming from Stable32

Stable32 and SigmaTau compute the same statistics from the same definitions
(NIST SP1065), and on the shipped fixture they agree to Stable32's five
printed figures ([Validation](../explanation/validation.md)). This page maps
a Stable32 session onto SigmaTau and lists the few defaults that differ.

## Reading the file

Stable32's `.phd` / `.frd` files, and text files with its header (`Tau:`
lines closed by `# Header End`, or the older bare `File:` … `Tau:` lines),
read with no keywords:

```julia
using SigmaTau
p = read_data("clock.phd")          # PhaseData, τ₀ from the Tau: line
y = read_data("clock.frd")          # FrequencyData
p = read_phase("export.dat")        # a Stable32 text export
```

| Stable32 | SigmaTau |
|---|---|
| Phase or frequency data type | `PhaseData` or `FrequencyData`; every estimator takes either |
| Tau (sample interval) | `tau0`, from the header or the time column, or `tau0=` |
| Multiplier | `scaling=` on the readers (`1e-9` for ns) |
| MJD time tags | the time column; `time_col=` picks it, τ₀ follows from it and [`start_time`](@ref) gives the first sample's date |
| Converting phase ↔ frequency | not needed for the statistics, which take either; `FrequencyData(p)` and `PhaseData(y)` convert a record, `N` frequency samples integrating to `N + 1` phase points as in Stable32 |

## Running the statistics

A Stable32 *Run* is one SigmaTau call; [`stability`](@ref) runs several at
once, and the [explorer](explorer.md) and the [command line](cli.md) do the
same from a window or a shell.

| Stable32 statistic | SigmaTau |
|---|---|
| Normal (non-overlapping) Allan | `adev(p; overlapping=false)` |
| Overlapping Allan | [`adev`](@ref) |
| Modified Allan | [`mdev`](@ref) |
| Time | [`tdev`](@ref) |
| Hadamard (non-overlapping) | `hdev(p; overlapping=false)` |
| Overlapping Hadamard | [`hdev`](@ref) |
| Total | `totdev(p; correct_bias=:stable32)` |
| Modified Total | `mtotdev(p; correct_bias=:stable32)` |
| Time Total | `ttotdev(p; correct_bias=:stable32)` |
| Hadamard Total | `htotdev(p; correct_bias=:stable32, remove_drift=false)` |
| Thêo1 | `theo1(p; correct_bias=:stable32)` |
| ThêoBR | [`theo1`](@ref) (the default) |
| ThêoH | [`theoh`](@ref) |
| MTIE, TIE rms | [`mtie`](@ref), [`tierms`](@ref) |
| Dynamic stability | [`dadev`](@ref), [`dhdev`](@ref) |
| Power spectrum | [`psd_y`](@ref), [`psd_x`](@ref), [`ssb_phase_noise`](@ref) |
| Drift analysis | [`drift_fit`](@ref) (`:linear`, `:quadratic`, `:log`, `:diffusion`) |
| Outlier check | [`find_outliers`](@ref), [`remove_outliers`](@ref) (median absolute deviation, 5σ) |
| Gap fill | [`fillgaps`](@ref) (Howe's reflect-and-filter method) |
| Three-cornered hat | [`ncornered_hat`](@ref) (any number of clocks ≥ 3) |

The modified Hadamard class ([`mhdev`](@ref), [`htdev`](@ref),
[`mhtotdev`](@ref)) has no Stable32 counterpart.

## Averaging times

Stable32's octave spacing is SigmaTau's `Octave`, the default; `Decade`,
`AllTaus`, `HalfOctave` and the others are [`TauMode`](@ref)s, and an
explicit list of averaging factors works too:

```julia
adev(p)                    # τ = 1, 2, 4, … × τ₀
adev(p, Decade)
adev(p, [1, 10, 100])      # m, so τ = m·τ₀
```

The longest τ follows each estimator's own limit, as in Stable32.

## Confidence intervals and noise identification

Every call identifies the noise type at each τ (lag-1 autocorrelation, with
the B1 and R(n) ratios where those apply) and returns χ² intervals from the
Greenhall–Riley EDF, the method Stable32 uses. The default level is 68.3 %
(Stable32's 1σ); `confidence=0.95` changes it and `ci=false` skips the
intervals for speed.

```julia
r = adev(p; confidence=0.95)
r.tau, r.dev, ci_lower(r), ci_upper(r), r.noise_type, r.edf
```

## Defaults that differ

`correct_bias=:stable32` on any estimator, on [`deviation`](@ref) or on
[`stability`](@ref) applies Stable32's bias convention: the first two rows
and the last three of this table. The command line takes `--bias stable32`,
and the explorer "As Stable32" under Bias correction.

| | Stable32 | SigmaTau default | To match Stable32 |
|---|---|---|---|
| MTOTDEV | raw estimator, which reads low | raised to MDEV's expected value (Howe and Peppler 2001, Table 2) | `correct_bias=:stable32` |
| Thêo1 | ThêoBR over the ladder SP1065 prints, `i = 0…⌊N/6 − 3⌋`; raw with Alpha set to 0 | ThêoBR over `i = 0…⌊N/30 − 3⌋` (Howe 2006) | `correct_bias=:stable32`; `correct_bias=false` for the raw run |
| Noise identification | the subseries as they are | each decimated subseries quadratically detrended first | `identify_noise(x, m; remove_drift=false)` |
| Noise type for HDEV and HTOTDEV | identified on the record as given | identified on the record less its least-squares frequency drift, without the B1 ratio: below 30 averaged points an FM type is the one at the longest τ that has 30 | `hdev(p; remove_drift=false)`, `htotdev(p; remove_drift=false)` |
| Flicker PM | the lag-1 result as it comes | a white-PM result confirmed with R(n) | no switch |
| PM below 30 averaged points | the B1 ratio and R(n) | PM only if it is PM at the longest τ with 30 points and R(n) at the τ itself agrees | `identify_noise(x, m; b1_length=:record)`; no switch on the estimators, and `remove_drift=false` does not change it |
| White or flicker PM, τ ≥ 16τ₀ | the result of the lag-1 or B1 test as it comes | stands only if R(n) agrees; otherwise white FM | no switch |
| Noise type with fewer than 30 averaged points | B1 ratio, labels as with B1(N, μ) at the record length | B1 ratio against B1(N, μ) at the number of averaged points; below 6 points the label measured at 6 | `identify_noise(x, m; b1_length=:record)`; no switch on the estimators |
| Limits between the FM noises, τ ≥ 2τ₀ | lag-1 estimate rounded to the nearest α | limits halfway between what the estimate reads on decimated phase | no switch |
| TOTDEV, HTOTDEV on FM noise | bias-corrected | bias-corrected | the same from τ = 16τ₀ |
| HTOTDEV on FM noise, τ = 2τ₀ to 15τ₀ | the factor published for τ ≥ 16τ₀ | the exact factor at that averaging factor, up to 3 % away in σ (white FM at 3τ₀) | `correct_bias=:stable32` |
| TOTDEV, HTOTDEV on white and flicker PM | raw estimator, which reads high | lowered to the expected value of ADEV, HDEV ([derivation](../explanation/total_pm.md)) | `correct_bias=:stable32` |

The detrending choices change the identified noise type only at long τ, where
few averaged samples remain; when it flips, the EDF and the interval at that τ
follow it, and for HTOTDEV the bias factor too.
[Validation](../explanation/validation.md) shows each case on the shipped
fixture.

The flicker-PM row shows at medium τ. On a pure flicker-PM record of 99 999
frequency samples made by Stable32's own generator, Stable32's HDEV run
labels averaging factors 16, 32 and 128 to 512 white PM and the others up to
4096 flicker PM; SigmaTau labels all of them flicker PM. The deviations agree
to 2e-5 at every τ and the intervals to 3e-4 wherever the labels do. Where
they do not, Stable32's interval is ±0.34 % against a simulated scatter of
±0.5 % (m = 16) to ±1.8 % (m = 512), and SigmaTau's is ±0.55 % to ±1.7 %.

The FM-limits row shows from medium τ. On decimated phase the lag-1 estimate
of flicker FM settles at α = −1.45, next to the −1.5 at which rounding turns it
into random-walk FM: in simulation (N = 4096) pure flicker FM keeps its label in
72 % of records at m = 8 and 48 % at m = 64 with rounding, and in 100 % and
87 % with SigmaTau's limits. Random-walk FM keeps its label in 98 % or more
of records at those m either way.

The B1 row shows at the longest τ, and there it moves σ as well as the
interval for the bias-corrected total estimators. At the record length the
B1 limit for random-walk FM sits far above what a few averaged points can
reach, so that convention almost never reports it: in simulation (N = 1024,
6 to 24 averaged points) random-walk FM is found in 1 to 2 % of records and
flicker FM in 24 to 48 %, against 60 to 75 % and 34 to 79 % at the number of
averaged points. White FM is found less often that way, in 33 to 82 % of
records against 65 to 88 %. Below 6 averaged points the ratio is close to
a coin toss, so SigmaTau reports the label measured at the longest averaging
factor that still has 6. On the 8192-point validation fixture Stable32
reads white FM at τ = 1024 s and beyond, where ADEV is flat and then rising
and SigmaTau reads flicker FM and random-walk FM; TOTDEV is then 3.1 %,
10.9 % and 26.5 % above Stable32's at 1024, 2048 and 4096 s.
[Validation](../explanation/validation.md) has the rows.

## Gaps

Stable32 needs gaps filled before most runs. Its files mark a gap with a
zero (any zero in frequency data, any but the first and last point in phase
data), and SigmaTau's Stable32 reader reads those zeros as gaps. In SigmaTau
`NaN` marks a missing sample; ADEV, MDEV, HDEV, MHDEV, TDEV, HTDEV, MTIE and TIE rms skip
the terms that touch it, and [`fillgaps`](@ref) (or `fillgaps=true` on the
readers) fills the record for the rest ([Choosing an
estimator](choosing.md)).

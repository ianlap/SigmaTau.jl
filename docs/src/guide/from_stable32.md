# Coming from Stable32

Stable32 and SigmaTau compute the same statistics from the same definitions
(NIST SP1065), and on the shipped fixture they agree to Stable32's five
printed figures ([Validation](../explanation/validation.md)). This page maps
a Stable32 session onto SigmaTau and lists the few defaults that differ.

## Reading the file

Stable32's `.phd` / `.frd` files, and text files with its header (`Tau:`
lines closed by `# Header End`), read with no keywords:

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
| MJD time tags | the time column; `time_col=` picks it and τ₀ follows from it |
| Converting phase ↔ frequency | not needed; `FrequencyData` integrates to `N + 1` phase points, as Stable32 does |

## Running the statistics

A Stable32 *Run* is one SigmaTau call; [`stability`](@ref) runs several at
once, and the [explorer](explorer.md) and the [command line](cli.md) do the
same from a window or a shell.

| Stable32 statistic | SigmaTau |
|---|---|
| Normal (non-overlapping) Allan | not provided; the overlapping form has more degrees of freedom |
| Overlapping Allan | [`adev`](@ref) |
| Modified Allan | [`mdev`](@ref) |
| Time | [`tdev`](@ref) |
| Hadamard (non-overlapping) | not provided |
| Overlapping Hadamard | [`hdev`](@ref) |
| Total | [`totdev`](@ref) |
| Modified Total | `mtotdev(p; correct_bias=false)`; see below |
| Time Total | [`ttotdev`](@ref) |
| Hadamard Total | [`htotdev`](@ref) |
| Thêo1 | `theo1(p; correct_bias=false)` |
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

| | Stable32 | SigmaTau default | To match Stable32 |
|---|---|---|---|
| MTOTDEV | raw estimator | divided by √B(α) (SP1065 Table 11), unbiased | `mtotdev(p; correct_bias=false)` |
| Thêo1 | raw Thêo1 run, ThêoBR a separate run | ThêoBR | `theo1(p; correct_bias=false)` |
| Noise identification | the subseries as they are | each decimated subseries quadratically detrended first | `identify_noise(x, m; remove_drift=false)` |
| TOTDEV, HTOTDEV | bias-corrected | bias-corrected | the same |

The detrending choice changes the identified noise type only where a
drift dominates a long τ; when it flips, the EDF and the interval at that τ
follow it. [Validation](../explanation/validation.md) shows each case
on the shipped fixture.

## Gaps

Stable32 needs gaps filled before most runs. In SigmaTau `NaN` marks a
missing sample; ADEV, MDEV, HDEV, MHDEV, TDEV, HTDEV, MTIE and TIE rms skip
the terms that touch it, and [`fillgaps`](@ref) (or `fillgaps=true` on the
readers) fills the record for the rest ([Choosing an
estimator](choosing.md)).

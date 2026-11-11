# Data and files

The readers make a record from a text file, and [`save`](@ref) and
[`load`](@ref) store records, results and suites.
[Reading and saving data](../tutorials/02_reading_your_data.md) works
through them on a file.

| File | Read with |
|:-----|:----------|
| Columns of numbers, with or without a time column | [`read_phase`](@ref), [`read_frequency`](@ref) |
| Stable32 and TimeLab files, counter readings, tables with named columns, BIPM files of UTC − UTC(k), RINEX clock files | [`read_data`](@ref), which finds the format |
| A file that [`save`](@ref) wrote | [`load`](@ref) |

The command line and the explorer read the same files
([The record](cli.md#The-record);
[Opening a file](../explorer/data.md#Opening-a-file)).

## Columns of numbers

```@docs
read_phase
read_frequency
```

## Timestamps and gaps

A time column places the samples. The readers fit its steps to an even
grid of step ``τ_0``, and a reading that is missing from the file leaves a
missing sample, `NaN`, in its place, so the samples on either side of it
keep their distance in time. A file without a time column is taken as
evenly sampled, in file order. [`read_phase`](@ref) states how the grid is
found and when a column is not taken for a time axis.

ISO 8601 timestamps and MJD values also give the record the UTC time of
its first sample, which [`save`](@ref) writes and the readers read back:

```@docs
start_time
```

What a function does with missing samples depends on what it computes:

| Functions | With missing samples |
|:----------|:---------------------|
| [`adev`](@ref), [`mdev`](@ref), [`tdev`](@ref), [`hdev`](@ref), [`mhdev`](@ref), [`htdev`](@ref), [`tierms`](@ref) | Leave out every term that reads a missing sample |
| [`mtie`](@ref) | Takes the peak-to-peak of the samples measured in each window |
| [`totdev`](@ref), [`mtotdev`](@ref), [`ttotdev`](@ref), [`htotdev`](@ref), [`mhtotdev`](@ref), [`pdev`](@ref), [`theo1`](@ref), [`theoh`](@ref), [`dadev`](@ref), [`dhdev`](@ref), [`cross_adev`](@ref), [`groslambert`](@ref), [`cross_psd_y`](@ref), [`cross_psd_x`](@ref), and `adev` and `hdev` with `overlapping=false` | Throw an error |
| [`StreamingStability`](@ref), [`LiveStability`](@ref) | Throw an error at the sample and stay as they were |
| [`psd_y`](@ref), [`psd_x`](@ref), [`ssb_phase_noise`](@ref) | Return `NaN` at every frequency |
| [`detrend`](@ref), [`drift_fit`](@ref), [`find_outliers`](@ref), [`summarize`](@ref) | Work on the samples measured |
| [`fillgaps`](@ref) | Returns a record without missing samples |
| [`remove_outliers`](@ref) | Replaces the outliers and leaves the missing samples missing |

The estimators of the first row follow
[SP1065 §10.4](@cite riley-2008-sp1065): a term counts when every sample
it reads was measured. For ADEV these are three phase samples and for HDEV
four, and MDEV and MHDEV need every sample of the averages they difference.
The noise type is identified on the record with each gap bridged by a
straight line. `neff` is the number of terms summed, and the EDF is that
of a complete record of the length that gives as many terms.

In a frequency record a missing sample loses the phase change across it,
so a term counts only when no missing sample lies within its span, and one
gap in the middle of the record removes most of the terms at long
averaging times. Converting such a record with `PhaseData` throws an
error.

[Gaps](../guide/preparing.md#Gaps) weighs the three ways to handle a
record with gaps: skip them, analyze a complete stretch, or fill them.

## Instrument files

[`read_data`](@ref) reads the files that instruments and their programs
write and returns phase or frequency as the file has it.
[`detect_format`](@ref) names the format it finds.

```@docs
read_data
detect_format
INPUT_FORMATS
```

## A first look

```@docs
summarize
```

## Saving and loading

```@docs
save
load
```

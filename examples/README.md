# Examples

Runnable scripts in [Literate.jl](https://github.com/fredrikekre/Literate.jl)
format. The docs build renders each one into a tutorial page.

| Script | Topic |
|---|---|
| `01_drifting_clock_hadamard.jl` | Drift on a rubidium-class record; `hdev`, `mhdev`, `htdev`, and `mhtotdev` on a short record. |
| `02_reading_your_data.jl` | `.DAT` files and counter CSVs; gaps, outliers, detrending; saving results. |
| `03_three_cornered_hat.jl` | Separating individual clock noise from pairwise comparisons; `nch`. |

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'   # first run only
julia --project=examples examples/01_drifting_clock_hadamard.jl
```

`examples/Project.toml` adds `Plots`; the package itself has no plotting
dependency.

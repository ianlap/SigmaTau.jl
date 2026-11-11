# Examples

Each script is a tutorial of the
[documentation](https://ianlap.github.io/SigmaTau.jl/dev/), written for
[Literate.jl](https://github.com/fredrikekre/Literate.jl): the comments are
the text of the page, and the documentation build runs the code and shows
its results. The scripts also run on their own. Their data are simulated,
so each one runs as it stands.

| Script | Subject |
|:---|:---|
| `01_drifting_clock_hadamard.jl` | A clock with a frequency drift: the drift in the Allan deviation, the Hadamard family, removing the drift, time error, and the total estimators on a short record |
| `02_reading_your_data.jl` | Reading columns, tables and instrument files, missing readings, a first look at a record, and saving records and results |
| `03_three_cornered_hat.jl` | Three clocks separated by the three-cornered hat and by the Groslambert covariance, and where each method stops working |
| `04_clock_ensemble.jl` | One clock from its noise model to a Kalman filter and a holdover prediction, then a time scale formed from four clocks |

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=examples examples/01_drifting_clock_hadamard.jl
```

`examples/Project.toml` adds `Plots`, which the package itself does not
depend on.

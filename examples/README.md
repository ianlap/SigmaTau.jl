# Examples

Each `.jl` file in this directory is a runnable Julia script that
doubles as a docs page on the published site. They're written in
[Literate.jl](https://github.com/fredrikekre/Literate.jl) format —
plain Julia code, with prose comments lifted into Markdown headings
when the docs build. Edit the script, rebuild, and the corresponding
`docs/src/tutorials/*.md` page regenerates.

| Script                              | Topic                                                                |
|-------------------------------------|----------------------------------------------------------------------|
| `00_julia_for_metrologists.jl`      | Stable32-to-Julia primer: install, load, plot, save.                |
| `01_phase_data.jl`                  | Build a `PhaseData` / `FrequencyData` record; phase ↔ frequency.    |
| `02_compute_adev.jl`                | Run `adev`; explore `StabilityResult` (CI, EDF, noise type).         |
| `06_three_cornered_hat.jl`          | Recover individual clock σ_y(τ) from three pairwise differences.     |

For clock state-space tutorials (Kalman filter, PID steering,
holdover budgets) see the sister package
[ClockEnsemble.jl](https://github.com/ianlap/ClockEnsemble.jl).

## Running

From the repo root:

```bash
julia --project=examples examples/01_phase_data.jl
```

`examples/Project.toml` is the env that carries `Plots` + `PGFPlotsX`
(the package's runtime `Project.toml` deliberately doesn't depend on
either — they're only needed for visualisation). On the first run,
`julia --project=examples -e 'using Pkg; Pkg.instantiate()'` will
resolve and download deps; subsequent runs skip straight to compile.
A clean checkout under `--project=.` would fail with
`Package Plots not found` — use `--project=examples`.

## Docs build

`docs/make.jl` walks this directory and runs `Literate.markdown` on
each `.jl` file, dropping the resulting `.md` into
`docs/src/tutorials/`. The generated pages are gitignored — the
`.jl` files are the single source of truth.

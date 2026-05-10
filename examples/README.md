# Examples

Each `.jl` file in this directory is a runnable Julia script that
doubles as a docs page on the published site. They're written in
[Literate.jl](https://github.com/fredrikekre/Literate.jl) format —
plain Julia code, with prose comments lifted into Markdown headings
when the docs build. Edit the script, rebuild, and the corresponding
`docs/src/tutorials/*.md` page regenerates.

| Script                              | Topic                                                                |
|-------------------------------------|----------------------------------------------------------------------|
| `01_phase_data.jl`                  | Build a `PhaseData` / `FrequencyData` record; phase ↔ frequency.    |
| `02_compute_adev.jl`                | Run `adev`; explore `StabilityResult` (CI, EDF, noise type).         |
| `03_kalman_single_clock.jl`         | Track a `ThreeStateClock` with `StandardKalmanFilter` (open-loop).   |
| `04_kalman_pid_steering.jl`         | Closed-loop PID steering with critical-damping gains.                |
| `05_holdover_comparison.jl`         | Holdover budget four ways: TDEV / HTDEV / KF RMS / KF 1σ via `prop!`.|
| `06_three_cornered_hat.jl`          | Recover individual clock σ_y(τ) from three pairwise differences.     |

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

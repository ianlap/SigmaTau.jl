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

## Running

From the repo root:

```bash
julia --project=. examples/01_phase_data.jl
```

The path-based project file (`Project.toml` in the repo root) is the
right environment for the examples. They use `Plots` + `PGFPlotsX` if
present, but they don't crash without them — the plot calls at the
end will simply error and you'll get the computational results in
the REPL.

## Docs build

`docs/make.jl` walks this directory and runs `Literate.markdown` on
each `.jl` file, dropping the resulting `.md` into
`docs/src/tutorials/`. The generated pages are gitignored — the
`.jl` files are the single source of truth.

# # Julia for metrologists — a Stable32 user's guide
#
# If you've been using Stable32 to calculate Allan deviations, this page
# is for you. In 10 minutes you'll load phase data, compute your first
# `adev`, inspect the result, and save it to disk — without clicking a
# single dialog box.

# ## Why Julia?
#
# Three things make Julia a compelling upgrade from Stable32:
#
# - **It's fast without configuration.** Julia compiles to native machine code.
#   Large datasets that bring Python or MATLAB to a crawl run at C speed after
#   the first call.
# - **It's interactive like MATLAB, free like Python.** The REPL gives you an
#   instant feedback loop. No license server, no toolbar to configure.
# - **It reads the same `.DAT` files Stable32 does.** Phase records exported
#   from your counter or measurement system work as-is.

# ## Installation in five minutes
#
# **Step 1** — install Julia via `juliaup` (the recommended version manager):
#
# ```bash
# # macOS / Linux
# curl -fsSL https://install.julialang.org | sh
#
# # Windows (PowerShell)
# winget install julia -s msstore
# ```
#
# **Step 2** — open a Julia REPL and install SigmaTau.jl:
#
# ```julia
# julia> ]   # enter package mode — the prompt changes to pkg>
# pkg> add SigmaTau
# pkg> <Backspace>   # return to normal julia> prompt
# julia> using SigmaTau
# ```
#
# That's it. SigmaTau is now on your path.

# ## Three ways to use Julia
#
# **REPL (interactive)** — type `julia` at a terminal. Results appear
# immediately. This is the closest thing to a Stable32 command window.
# Great for exploration and one-off calculations.
#
# **Script** — put your code in a `.jl` file and run it with
# `julia myscript.jl`. Useful for reproducible, documented workflows
# that you want to run again later.
#
# **Pluto.jl (reactive notebook)** — a browser-based notebook where every
# cell reruns automatically when you change any value. Think Excel but for
# Julia code. Install it once with `] add Pluto` and launch with:
#
# ```julia
# using Pluto
# Pluto.run()
# ```
#
# For getting started, the REPL is the fastest path. The tutorials in this
# documentation are all runnable scripts; everything here can also be pasted
# directly into a Pluto cell.

# ## Typing τ, σ, and friends
#
# Julia identifiers can be Greek letters, and this documentation uses them
# the way the field does: `τ` for averaging time, `σ` for a deviation,
# `τ₀` for the sample interval. You type them with LaTeX-style
# abbreviations followed by Tab. In the REPL and in VS Code (with the
# Julia extension), `\tau` + Tab gives `τ`, `\sigma` + Tab gives `σ`,
# and a subscript follows the same pattern: `\tau` + Tab, then `\_0` +
# Tab gives `τ₀`. Pluto notebooks support the same completion.
#
# None of this is required — every Greek identifier in SigmaTau's API has
# a plain-ASCII spelling (`tau0`, `dev`, …), and your own variables can be
# named however you like. The completion is just how you read and write
# the symbols comfortably when you want them.

# ## Loading phase data from a file
#
# Stable32 exports phase records as plain-text `.DAT` files: a 10-line header
# followed by one phase value per line (in seconds). Here is how to load one:
#
# ```julia
# lines = readlines("my_data.DAT")
# x     = parse.(Float64, strip.(lines[11:end]))  # skip the 10-line header
# pd    = PhaseData(x, 1.0)                        # 1.0 s sample interval
# ```
#
# `PhaseData` takes your phase vector `x` and the sample interval `τ₀` in
# seconds. If you have fractional-frequency data instead of phase, use
# `FrequencyData(y, τ₀)` — the deviation functions accept both.

# ## Your first Allan deviation
#
# Here we generate synthetic white-FM phase data and compute the overlapping
# Allan deviation at four averaging times:

import Random
Random.seed!(20260510)
using SigmaTau

pd = PhaseData(randn(512) .* 1e-9, 1.0)   # 512-point, 1 s sample interval

result = adev(pd, [1, 2, 4, 8, 16, 32])

# The return value is a `StabilityResult`. The two most important fields:
println("τ values (s): ", result.tau)
println("ADEV values:  ", result.dev)

# ## What's inside a StabilityResult
#
# | Field | Stable32 equivalent | Notes |
# |---|---|---|
# | `tau` | τ column | Averaging times in seconds |
# | `dev` | σ_y(τ) column | Allan deviation values |
# | `noise_type` | noise ID column | e.g. `:WHFM`, `:RWFM` |
# | `ci` | error bars | `(lo, hi)` 1σ confidence bounds per τ |
# | `edf` | EDF column | Equivalent degrees of freedom |
# | `neff` | # column | Number of analysis windows averaged |
# | `deviation_type` | window title | `:adev`, `:mdev`, `:totdev`, etc. |
#
# Confidence intervals are computed by default (set `ci=false` to skip).
# The noise type is identified automatically by the lag-1 autocorrelation
# method (NIST SP1065 §4.2). Each `ci` entry is a named tuple — read one
# bound as `result.ci[1].lo`, or get whole vectors with the `ci_lower` /
# `ci_upper` accessor functions.

println("\nNoise types:  ", result.noise_type)
println("CI lower:     ", ci_lower(result))
println("CI upper:     ", ci_upper(result))
println("EDF:          ", result.edf)
println("neff:         ", result.neff)

# ## Overlaying multiple deviations on one plot
#
# This is something Stable32 makes awkward — overlaying ADEV and MDEV from the
# same dataset on a single log-log stability plot. In Julia it's two lines:
#
# ```julia
# using Plots
#
# r_adev  = adev(pd,  [1, 2, 4, 8, 16, 32])
# r_mdev  = mdev(pd,  [1, 2, 4, 8, 16, 32])
# r_tdev  = tdev(pd,  [1, 2, 4, 8, 16, 32])
#
# plot(r_adev, label="ADEV")     # uses the built-in plot recipe
# plot!(r_mdev, label="MDEV")    # ! means "add to the existing plot"
# plot!(r_tdev, label="TDEV")
# ```
#
# `plot(result)` automatically uses a log-log axis and draws shaded error bars
# from the `ci` bounds. Every additional `plot!(...)` call overlays on
# the same figure.

# ## Saving your result to disk
#
# `save_result` writes a tab-separated file you can open in Excel or any
# text editor. `load_result` reads it back into a `StabilityResult`:

tmpdir = mktempdir()
path   = joinpath(tmpdir, "my_adev.tsv")

save_result(path, result)
result2 = load_result(path)

println("\nRound-trip check — devs match: ", result2.dev ≈ result.dev)

# The file format is human-readable plain text: two comment lines encode the
# deviation type and format version, followed by a header row and one data row
# per τ value.

# ## Where to go next
#
# - [Tutorial 01: Phase and frequency data](01_phase_data.md) — deep dive into
#   `PhaseData` / `FrequencyData` construction, `tau0`, and the
#   `_freq_to_phase` conversion that all deviation functions use internally.
#
# - [Tutorial 02: Computing ADEV and reading the result](02_compute_adev.md) —
#   full walkthrough of `adev` options, the `ci` and `confidence`
#   keyword arguments, and how `StabilityResult` fields map to the Stable32
#   output window.
#
# - [Theory: overview](../theory/overview.md) — the five power-law noise types,
#   estimator families, and the mathematical framework behind SigmaTau.

### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ bf66ba3b-78ff-4ec4-9a5e-da844f0e4758
begin
    import Pkg
    Pkg.activate(@__DIR__)      # app/Project.toml: SigmaTau (path ..), Plots, PlutoUI
    Pkg.instantiate()
    using SigmaTau, Plots, PlutoUI, Random
end

# ╔═╡ d38c9b49-c2c7-42b2-b204-1b6155c68d17
md"""
# SigmaTau stability explorer

Load a phase or frequency record (or synthesize one), choose the estimators,
τ grid and confidence level, and read off σ(τ) with χ² error bars, the
identified noise type and the equivalent degrees of freedom. Every control
re-runs the analysis.
"""

# ╔═╡ 4caa5ccc-13cf-4f6e-951c-498dd6903f3b
md"""
## 1. Data

Source: $(@bind source Select(["synth" => "Synthesize power-law noise", "file" => "Load a data file"]; default="synth"))
"""

# ╔═╡ c6d1f2f1-ef48-4b12-a067-cbce4759d81c
md"""
**File** (used when *Load a data file* is selected)

Upload: $(@bind upload FilePicker()) or path: $(@bind path_text TextField(50; default=""))

Record type: $(@bind kind Select(["phase" => "phase x(t), seconds", "freq" => "fractional frequency y(t)"]))
Time column (0 = none): $(@bind time_col NumberField(0:20; default=0))
Value column: $(@bind value_col NumberField(1:20; default=1))
Header rows to skip: $(@bind header NumberField(0:1000; default=0))
τ₀ in seconds (leave blank to infer from the time column): $(@bind tau0_text TextField(8; default="1.0"))
"""

# ╔═╡ fb6d2dc3-9ddf-4a98-a49f-51ed53108a8f
md"""
**Synthetic record** (used when *Synthesize* is selected). Each term is
specified by its deviation at the sample interval, σ_y(τ₀); the slider sets
log₁₀ of that level.

Samples N: $(@bind N_syn Select([4096, 8192, 16384, 32768, 86400]; default=16384))
τ₀ (s): $(@bind tau0_syn NumberField(0.001:0.001:1000.0; default=1.0))
Seed: $(@bind seed NumberField(1:10_000; default=42))

| | on | log₁₀ σ_y(τ₀) |
|---|---|---|
| White PM (α = +2) | $(@bind wpm_on CheckBox(default=true)) | $(@bind wpm_exp Slider(-16:-9; default=-11, show_value=true)) |
| White FM (α = 0) | $(@bind wfm_on CheckBox(default=true)) | $(@bind wfm_exp Slider(-16:-9; default=-12, show_value=true)) |
| Flicker FM (α = −1) | $(@bind ffm_on CheckBox(default=false)) | $(@bind ffm_exp Slider(-16:-9; default=-13, show_value=true)) |
| Random-walk FM (α = −2) | $(@bind rwfm_on CheckBox(default=true)) | $(@bind rwfm_exp Slider(-18:-10; default=-14, show_value=true)) |
"""

# ╔═╡ dd2dae12-0261-4490-895e-d298ea6db75e
md"""
## 2. Preprocessing

Detrend: $(@bind detrend_mode Select([:none, :mean, :linear, :quadratic]; default=:none))
Remove 5σ MAD outliers: $(@bind outliers_on CheckBox(default=false))
Fill gaps (NaN samples): $(@bind fill_on CheckBox(default=false))
"""

# ╔═╡ 51307329-150f-4c9a-bfdf-6246bbe8ec82
md"""
## 3. Analysis

Deviations: $(@bind devs MultiCheckBox(["adev", "mdev", "tdev", "hdev", "mhdev", "htdev", "totdev", "mtotdev", "ttotdev", "htotdev", "mhtotdev", "theo1", "theoh", "pdev", "mtie", "tierms"]; default=["adev", "mdev", "hdev", "mhdev"]))

τ grid: $(@bind taumode_name Select(["Octave", "HalfOctave", "QuarterOctave", "Decade", "HalfDecade", "AllTaus"]; default="Octave"))
Confidence intervals: $(@bind ci_on CheckBox(default=true))
level: $(@bind confidence Select([0.683 => "68.3 % (1σ)", 0.90 => "90 %", 0.95 => "95 %", 0.99 => "99 %"]; default=0.683))
Bias-correct the total family and Thêo1 (ThêoBR): $(@bind correct_bias CheckBox(default=true))
Draw CI as shaded band instead of error bars: $(@bind ci_band_on CheckBox(default=false))
"""

# ╔═╡ d6f381f4-6363-4069-a46e-3dc2d4feeb6f
# The record as loaded or synthesized, before preprocessing.
data = let
    if source == "file"
        bytes = upload === nothing ? UInt8[] : get(upload, "data", UInt8[])
        path = if !isempty(bytes)
            tmp = joinpath(mktempdir(), get(upload, "name", "upload.dat"))
            write(tmp, bytes)
            tmp
        else
            strip(path_text)
        end
        isempty(path) && error("Upload a file or enter a path.")
        t0 = isempty(strip(tau0_text)) ? nothing : parse(Float64, tau0_text)
        reader = kind == "phase" ? read_phase : read_frequency
        reader(path; tau0=t0, time_col=time_col, value_col=value_col, header=header)
    else
        sigma1 = Dict{Int,Float64}()
        wpm_on  && (sigma1[2]  = 10.0^wpm_exp)
        wfm_on  && (sigma1[0]  = 10.0^wfm_exp)
        ffm_on  && (sigma1[-1] = 10.0^ffm_exp)
        rwfm_on && (sigma1[-2] = 10.0^rwfm_exp)
        isempty(sigma1) && error("Switch on at least one noise term.")
        noise_gen(PhaseData, N_syn, tau0_syn; sigma1=sigma1, rng=Xoshiro(seed))
    end
end

# ╔═╡ 6f3523b8-aac6-446e-be6e-67fb675e2315
# Preprocessed record that every estimator below sees.
record = let d = data
    outliers_on && (d = remove_outliers(d))
    fill_on && (d = fillgaps(d))
    detrend_mode != :none && (d = detrend(d; method=detrend_mode))
    d
end

# ╔═╡ 011d48c8-404f-42ed-91ba-b4d748588f11
suite = let
    isempty(devs) && error("Select at least one deviation.")
    taumode = Dict("Octave" => Octave, "HalfOctave" => HalfOctave,
                   "QuarterOctave" => QuarterOctave, "Decade" => Decade,
                   "HalfDecade" => HalfDecade, "AllTaus" => AllTaus)[taumode_name]
    stability(record; devs=Tuple(Symbol.(devs)), taus=taumode,
              ci=ci_on, confidence=confidence, correct_bias=correct_bias)
end

# ╔═╡ 27b74d46-ce84-4c01-9ee8-fa81cfbdfa1e
# σ_y (dimensionless) and σ_x (seconds) results go on separate axes.
let
    is_sigma_x(r) = r.deviation_type in (:tdev, :htdev, :ttotdev, :mtie, :tierms)
    sy = filter(!is_sigma_x, suite.results)
    sx = filter(is_sigma_x, suite.results)
    ttl = string(suite.n, " samples, τ₀ = ", suite.tau0, " s")
    panels = Any[]
    isempty(sy) || push!(panels, plot(sy; ci_band=ci_band_on, legend=:bottomleft,
                                      ylabel="σ_y(τ)", title=ttl))
    isempty(sx) || push!(panels, plot(sx; ci_band=ci_band_on, legend=:topleft,
                                      ylabel="σ_x(τ) (s)", title=isempty(sy) ? ttl : ""))
    plot(panels...; layout=(length(panels), 1), size=(820, 520 * length(panels)))
end

# ╔═╡ f3ad8ee2-4d70-434f-b209-240b76a7976c
md"""
## 4. Table

Show: $(@bind table_dev Select(string.(keys(suite)); default=string(first(keys(suite)))))
"""

# ╔═╡ 879d17ed-9117-401f-be13-cd4a7c8134f1
let r = suite[Symbol(table_dev)]
    has_ci = !isempty(r.edf)
    g(v) = isfinite(v) ? string(round(v; sigdigits=5)) : "NaN"
    io = IOBuffer()
    println(io, has_ci ? "| τ (s) | σ | neff | CI lo | CI hi | EDF | noise |" :
                         "| τ (s) | σ | neff |")
    println(io, has_ci ? "|---:|---:|---:|---:|---:|---:|---|" : "|---:|---:|---:|")
    for i in eachindex(r.tau)
        if has_ci
            println(io, "| ", g(r.tau[i]), " | ", g(r.dev[i]), " | ", r.neff[i], " | ",
                    g(r.ci[i].lo), " | ", g(r.ci[i].hi), " | ", g(r.edf[i]), " | ",
                    r.noise_type[i], " |")
        else
            println(io, "| ", g(r.tau[i]), " | ", g(r.dev[i]), " | ", r.neff[i], " |")
        end
    end
    Markdown.parse(String(take!(io)))
end

# ╔═╡ a5be2956-e65a-49d8-9fd0-4040e4f8aeec
md"""
## 5. Export

The whole suite as a self-describing TSV (`load_suite` reads it back):
$(DownloadButton(let p = joinpath(mktempdir(), "suite.tsv"); save_suite(p, suite; source_file=string(record.source)); read(p) end, "sigmatau_suite.tsv"))
"""

# ╔═╡ f05f370d-0a8a-4988-b259-3bca164cf29a
md"""
## 6. Dynamic stability (optional)

Slide a window along the record and watch σ_y(τ) change with time.

Compute: $(@bind dyn_on CheckBox(default=false))
Estimator: $(@bind dyn_dev Select(["dadev" => "dynamic ADEV", "dhdev" => "dynamic HDEV"]))
Window as a fraction of the record: $(@bind dyn_frac Select([1/4, 1/8, 1/16, 1/32]; default=1/8))
View: $(@bind dyn_view Select(["heatmap" => "heatmap σ(t, τ)", "path3d" => "stacked σ(τ) curves"]))
"""

# ╔═╡ cfb825dc-bfae-4a20-a695-0268c361aadd
if dyn_on
    n = record isa PhaseData ? length(record.x) : length(record.y)
    window = min(n, max(64, round(Int, n * dyn_frac)))   # never longer than the record
    f = dyn_dev == "dadev" ? dadev : dhdev
    dyn = f(record, Octave; window=window)
    dyn_view == "path3d" ? plot(dyn; seriestype=:path3d, size=(820, 560)) :
                           plot(dyn; size=(820, 520))
else
    md"*Dynamic map off.*"
end

# ╔═╡ Cell order:
# ╟─bf66ba3b-78ff-4ec4-9a5e-da844f0e4758
# ╟─d38c9b49-c2c7-42b2-b204-1b6155c68d17
# ╟─4caa5ccc-13cf-4f6e-951c-498dd6903f3b
# ╟─c6d1f2f1-ef48-4b12-a067-cbce4759d81c
# ╟─fb6d2dc3-9ddf-4a98-a49f-51ed53108a8f
# ╟─dd2dae12-0261-4490-895e-d298ea6db75e
# ╟─51307329-150f-4c9a-bfdf-6246bbe8ec82
# ╠═d6f381f4-6363-4069-a46e-3dc2d4feeb6f
# ╠═6f3523b8-aac6-446e-be6e-67fb675e2315
# ╠═011d48c8-404f-42ed-91ba-b4d748588f11
# ╟─27b74d46-ce84-4c01-9ee8-fa81cfbdfa1e
# ╟─f3ad8ee2-4d70-434f-b209-240b76a7976c
# ╟─879d17ed-9117-401f-be13-cd4a7c8134f1
# ╟─a5be2956-e65a-49d8-9fd0-4040e4f8aeec
# ╟─f05f370d-0a8a-4988-b259-3bca164cf29a
# ╟─cfb825dc-bfae-4a20-a695-0268c361aadd

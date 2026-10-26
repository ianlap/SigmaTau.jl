# SigmaTau stability explorer: a desktop window for loading or synthesizing a
# record, choosing estimators and options, and reading the plots and table.
#
#   julia --project=app app/StabilityExplorer.jl            # open the window
#   julia --project=app app/StabilityExplorer.jl --smoke    # headless self-test
#
# Three pages share the plot area: Deviations (σ_y and σ_x curves with χ²
# intervals, and the per-τ table), Spectrum (S_y, S_x and ℒ(f)), and Dynamic
# map (σ(t, τ)). Export renders a clean figure through CairoMakie (PNG, PDF
# or SVG) and writes tables and records as the text files `load_suite` and
# `read_phase` read back. Every number shown comes from the public SigmaTau
# API.

using SigmaTau
using GLMakie
using CairoMakie
using NativeFileDialog
using Random
using Printf

GLMakie.activate!()

const DEVIATIONS = [:adev, :mdev, :tdev, :hdev, :mhdev, :htdev,
                    :totdev, :mtotdev, :ttotdev, :htotdev, :mhtotdev,
                    :theo1, :theoh, :pdev, :mtie, :tierms]
const DEFAULT_ON = Set([:adev, :mdev, :hdev, :mhdev])
const SIGMA_X    = Set([:tdev, :htdev, :ttotdev, :mtie, :tierms])
const TAU_MODES  = ["Octave" => Octave, "Half octave" => HalfOctave,
                    "Quarter octave" => QuarterOctave, "Decade" => Decade,
                    "Half decade" => HalfDecade, "All τ" => AllTaus]
const NOISES     = [(2, "White PM", -11), (1, "Flicker PM", -16),
                    (0, "White FM", -12), (-1, "Flicker FM", -16),
                    (-2, "Random-walk FM", -14)]
const PAGES      = ["Deviations", "Spectrum", "Dynamic map"]
const MAX_TABLE_ROWS = 40

# ── Analysis (no GUI state) ───────────────────────────────────────────────

nsamples(rec::PhaseData) = length(rec.x)
nsamples(rec::FrequencyData) = length(rec.y)

function load_record(path, kind, tau0_text, time_col, value_col, header)
    isempty(strip(path)) && error("Choose a file first (Browse…).")
    t0 = isempty(strip(tau0_text)) ? nothing : parse(Float64, tau0_text)
    reader = kind == "phase" ? read_phase : read_frequency
    return reader(strip(path); tau0=t0, time_col=time_col, value_col=value_col, header=header)
end

function synthesize(N, tau0, seed, levels)
    sigma1 = Dict{Int,Float64}(α => 10.0^e for (α, e) in levels)
    isempty(sigma1) && error("Switch on at least one noise term.")
    return noise_gen(PhaseData, N, tau0; sigma1=sigma1, rng=Xoshiro(seed))
end

function preprocess(rec; outliers=false, gaps=false, detrend_mode=:none)
    outliers && (rec = remove_outliers(rec))
    gaps && (rec = fillgaps(rec))
    detrend_mode != :none && (rec = detrend(rec; method=detrend_mode))
    return rec
end

function analyze(rec, devs, taumode, ci, confidence, correct_bias)
    isempty(devs) && error("Select at least one deviation.")
    return stability(rec; devs=Tuple(devs), taus=taumode, ci=ci,
                     confidence=confidence, correct_bias=correct_bias)
end

function spectra(rec, f_carrier)
    sy = Sy(rec); sx = Sx(rec)
    l = f_carrier > 0 ? L(rec; f_carrier=f_carrier) : nothing
    return (; sy, sx, l)
end

function dynamic_map(rec, dev_sym, frac)
    n = nsamples(rec)
    window = min(n, max(64, round(Int, n * frac)))
    f = dev_sym == :dadev ? dadev : dhdev
    return f(rec, Octave; window=window)
end

# ── Drawing (shared by the window and the export figure) ──────────────────

function draw_results!(ax, results; band=false, legend_position=:lb)
    empty!(ax)
    palette = Makie.wong_colors()
    for (i, r) in enumerate(results)
        color = palette[mod1(i, length(palette))]
        keep = findall(v -> isfinite(v) && v > 0, r.dev)
        isempty(keep) && continue
        τ, σ = r.tau[keep], r.dev[keep]
        if !isempty(r.ci)
            lo, hi = ci_lower(r)[keep], ci_upper(r)[keep]
            ok = isfinite.(lo) .& isfinite.(hi) .& (lo .> 0)
            if band
                band!(ax, τ[ok], lo[ok], hi[ok]; color=(color, 0.22))
            else
                errorbars!(ax, τ[ok], σ[ok], (σ .- lo)[ok], (hi .- σ)[ok];
                           color=color, whiskerwidth=7, linewidth=1)
            end
        end
        scatterlines!(ax, τ, σ; color=color, markersize=9, linewidth=1.6,
                      strokecolor=:white, strokewidth=0.8,
                      label=uppercase(string(r.deviation_type)))
    end
    autolimits!(ax)
    return isempty(results) ? nothing : axislegend(ax; position=legend_position, framevisible=false)
end

function draw_spectrum!(ax, s::SpectralResult; color)
    empty!(ax)
    keep = s.freq .> 0
    lines!(ax, s.freq[keep], s.psd[keep]; color=color, linewidth=1.4)
    autolimits!(ax)
end

# `cb` is a Ref holding the previous Colorbar (or nothing); `cb_cell` the grid
# position it goes in. The colorbar is rebuilt from the heatmap each time.
function draw_dynamic!(ax, dyn; cb=nothing, cb_cell=nothing)
    empty!(ax)
    z = map(v -> (isfinite(v) && v > 0) ? log10(v) : NaN, dyn.dev)
    ly = log10.(dyn.tau)
    hm = heatmap!(ax, dyn.t, ly, z; nan_color=:transparent, colormap=:viridis)
    ks = floor(Int, minimum(ly)):ceil(Int, maximum(ly))
    ax.yticks = (Float64.(ks), ["1e$k" for k in ks])
    ax.title = uppercase(string(dyn.deviation_type)) * ", window = $(dyn.window) samples"
    autolimits!(ax)
    if cb !== nothing
        cb[] === nothing || delete!(cb[])
        cb[] = Colorbar(cb_cell, hm; label="log₁₀ σ")
    end
    return hm
end

function fill_table!(grid, r::StabilityResult)
    foreach(delete!, contents(grid))
    trim!(grid)
    has_ci = !isempty(r.edf)
    header = has_ci ? ["τ (s)", "σ", "neff", "CI lo", "CI hi", "EDF", "noise"] :
                      ["τ (s)", "σ", "neff"]
    g(v) = isfinite(v) ? @sprintf("%.5g", v) : "NaN"
    for (j, h) in enumerate(header)
        Label(grid[1, j], h; font=:bold, halign=:right, tellwidth=true, color=:gray20)
    end
    nrows = min(length(r.tau), MAX_TABLE_ROWS)
    for i in 1:nrows
        cells = has_ci ?
            [g(r.tau[i]), g(r.dev[i]), string(r.neff[i]), g(r.ci[i].lo),
             g(r.ci[i].hi), g(r.edf[i]), string(r.noise_type[i])] :
            [g(r.tau[i]), g(r.dev[i]), string(r.neff[i])]
        for (j, c) in enumerate(cells)
            Label(grid[i + 1, j], c; halign=:right, tellwidth=true,
                  color=isodd(i) ? :gray10 : :gray35)
        end
    end
    if length(r.tau) > nrows
        Label(grid[nrows + 2, 1:length(header)],
              "… $(length(r.tau) - nrows) more rows; export the table for all of them";
              halign=:left, tellwidth=false, color=:gray50)
    end
    colgap!(grid, 14)
    rowgap!(grid, 2)
end

# Publication figure for the current page, rendered through CairoMakie.
function export_figure(page, state, opts)
    fig = Figure(size=(900, page == "Spectrum" ? 900 : 560), fontsize=13)
    if page == "Deviations"
        suite = state.suite[]
        sy = filter(r -> !(r.deviation_type in SIGMA_X), suite.results)
        sx = filter(r ->   r.deviation_type in SIGMA_X,  suite.results)
        row = 1
        if !isempty(sy)
            ax = Axis(fig[row, 1]; xscale=log10, yscale=log10,
                      xlabel="Averaging time τ (s)", ylabel="σ_y(τ)", title=opts.title)
            draw_results!(ax, sy; band=opts.band); row += 1
        end
        if !isempty(sx)
            ax = Axis(fig[row, 1]; xscale=log10, yscale=log10,
                      xlabel="Averaging time τ (s)", ylabel="σ_x(τ) (s)",
                      title=isempty(sy) ? opts.title : "")
            draw_results!(ax, sx; band=opts.band, legend_position=:lt)
        end
    elseif page == "Spectrum"
        sp = state.spectra[]
        ax1 = Axis(fig[1, 1]; xscale=log10, yscale=log10, ylabel="S_y(f) (1/Hz)", title=opts.title)
        ax2 = Axis(fig[2, 1]; xscale=log10, yscale=log10, ylabel="S_x(f) (s²/Hz)")
        draw_spectrum!(ax1, sp.sy; color=Makie.wong_colors()[1])
        draw_spectrum!(ax2, sp.sx; color=Makie.wong_colors()[2])
        if sp.l !== nothing
            ax3 = Axis(fig[3, 1]; xscale=log10, xlabel="Frequency f (Hz)", ylabel="ℒ(f) (dBc/Hz)")
            draw_spectrum!(ax3, sp.l; color=Makie.wong_colors()[3])
        else
            ax2.xlabel = "Frequency f (Hz)"
        end
    else
        ax = Axis(fig[1, 1]; xlabel="Time t (s)", ylabel="τ (s)")
        draw_dynamic!(ax, state.dynamic[]; cb=Ref{Any}(nothing), cb_cell=fig[1, 2])
    end
    return fig
end

# ── Window ────────────────────────────────────────────────────────────────

all_blocks(x) = Makie.Block[]
all_blocks(b::Makie.Block) = Makie.Block[b]
all_blocks(gl::GridLayout) = reduce(vcat, (all_blocks(c) for c in contents(gl)); init=Makie.Block[])
function set_visible!(blocks, v)
    for b in blocks
        b.blockscene.visible[] = v
        b isa Axis && (b.scene.visible[] = v)
    end
end

function build()
    fig = Figure(size=(1600, 1000), fontsize=12, backgroundcolor=:gray98)
    controls = GridLayout(fig[1, 1]; tellheight=false, valign=:top)
    plots    = GridLayout(fig[1, 2])
    colsize!(fig.layout, 1, Fixed(500))
    colgap!(fig.layout, 30)
    rowgap!(controls, 2)
    colgap!(controls, 8)
    const_h = 26
    small(str) = Textbox(fig; stored_string=str, width=58, height=const_h)
    menu(args...; kw...) = Menu(fig, args...; height=const_h, kw...)
    textbox(; kw...) = Textbox(fig; height=const_h, kw...)

    row = Ref(1)
    nextrow() = (r = row[]; row[] += 1; r)
    function section(text)
        r = nextrow()
        Box(controls[r, 1:2]; color=:gray88, strokevisible=false, cornerradius=4)
        Label(controls[r, 1:2], text; font=:bold, halign=:left, tellwidth=false,
              padding=(8, 8, 4, 4), color=:gray15)
    end
    function labelled(text, block)
        r = nextrow()
        Label(controls[r, 1], text; halign=:right, tellwidth=true, color=:gray25)
        controls[r, 2] = block
        return block
    end

    section("Data")
    source = menu(; options=["Synthesize", "File"], default="Synthesize", width=110)
    kind   = menu(; options=["phase", "frequency"], default="phase", width=100)
    let r = nextrow()
        Label(controls[r, 1], "Source, record type"; halign=:right, tellwidth=true, color=:gray25)
        inner = GridLayout(controls[r, 2]); inner[1, 1] = source; inner[1, 2] = kind
    end
    path_box = textbox(; placeholder="choose a file with Browse…", tellwidth=false)
    open_btn = Button(fig; label="Browse…", height=const_h)
    let r = nextrow()
        Label(controls[r, 1], "File"; halign=:right, tellwidth=true, color=:gray25)
        inner = GridLayout(controls[r, 2]); inner[1, 1] = path_box; inner[1, 2] = open_btn
    end
    tau0_box = small("1.0")
    time_col, value_col, header = small("0"), small("1"), small("0")
    let r = nextrow()
        Label(controls[r, 1], "τ₀ (s)  |  time, value, header columns"; halign=:right, tellwidth=true, color=:gray25)
        inner = GridLayout(controls[r, 2])
        inner[1, 1] = tau0_box; inner[1, 2] = time_col; inner[1, 3] = value_col; inner[1, 4] = header
    end

    section("Synthesize  (level = σ_y at τ₀)")
    n_menu   = menu(; options=["4096", "8192", "16384", "32768", "86400"], default="16384", width=90)
    seed_box = small("42")
    let r = nextrow()
        Label(controls[r, 1], "Samples, seed"; halign=:right, tellwidth=true, color=:gray25)
        inner = GridLayout(controls[r, 2]); inner[1, 1] = n_menu; inner[1, 2] = seed_box
    end
    noise_toggles = Dict{Int,Toggle}()
    noise_sliders = Dict{Int,Slider}()
    for (α, name, start) in NOISES
        r = nextrow()
        noise_toggles[α] = Toggle(fig; active=(α in (2, 0, -2)))
        noise_sliders[α] = Slider(fig; range=-18:-9, startvalue=start)
        Label(controls[r, 1], name; halign=:right, color=:gray25)
        inner = GridLayout(controls[r, 2])
        inner[1, 1] = noise_toggles[α]
        inner[1, 2] = noise_sliders[α]
        inner[1, 3] = Label(fig, lift(v -> @sprintf("1e%d", v), noise_sliders[α].value);
                            width=40, halign=:left, color=:gray25)
    end

    section("Preprocessing")
    detrend_menu = labelled("Detrend", menu(; options=["none", "mean", "linear", "quadratic"], default="none"))
    outliers = Toggle(fig; active=false)
    gaps     = Toggle(fig; active=false)
    let r = nextrow()
        Label(controls[r, 1], "Remove 5σ outliers"; halign=:right, tellwidth=true, color=:gray25)
        inner = GridLayout(controls[r, 2]); inner[1, 1] = outliers
        inner[1, 2] = Label(fig, "fill NaN gaps"; halign=:right, color=:gray25); inner[1, 3] = gaps
    end

    section("Deviations")
    dev_toggles = Dict{Symbol,Toggle}()
    devgrid = GridLayout(controls[nextrow(), 1:2])
    for (i, d) in enumerate(DEVIATIONS)
        rr, cc = fldmod1(i, 4)
        dev_toggles[d] = Toggle(fig; active=(d in DEFAULT_ON))
        cell = GridLayout(devgrid[rr, cc])
        cell[1, 1] = dev_toggles[d]
        cell[1, 2] = Label(fig, string(d); halign=:left, color=:gray20)
    end

    section("Analysis")
    tau_menu    = labelled("τ grid", menu(; options=first.(TAU_MODES), default="Octave"))
    ci_toggle   = Toggle(fig; active=true)
    conf_menu   = menu(; options=["0.683", "0.90", "0.95", "0.99"], default="0.683", width=80)
    let r = nextrow()
        Label(controls[r, 1], "Confidence intervals, level"; halign=:right, tellwidth=true, color=:gray25)
        inner = GridLayout(controls[r, 2]); inner[1, 1] = ci_toggle; inner[1, 2] = conf_menu
    end
    bias_toggle = Toggle(fig; active=true)
    band_toggle = Toggle(fig; active=false)
    let r = nextrow()
        Label(controls[r, 1], "Bias-correct totals, Thêo1"; halign=:right, tellwidth=true, color=:gray25)
        inner = GridLayout(controls[r, 2]); inner[1, 1] = bias_toggle
        inner[1, 2] = Label(fig, "CI as band"; halign=:right, color=:gray25); inner[1, 3] = band_toggle
    end
    carrier_box = textbox(; stored_string="10e6", width=90)
    dyn_menu  = menu(; options=["dadev", "dhdev"], default="dadev", width=80)
    frac_menu = menu(; options=["1/4", "1/8", "1/16", "1/32"], default="1/8", width=70)
    run_btn = Button(fig; label="Run analysis", buttoncolor=RGBf(0.33, 0.47, 0.72),
                     labelcolor=:white, buttoncolor_hover=RGBf(0.4, 0.55, 0.8),
                     font=:bold, height=30)
    controls[nextrow(), 1:2] = run_btn

    status = Label(fig[2, 1:2], "Ready. Press Run analysis."; halign=:left, tellwidth=false,
                   color=:gray30, padding=(6, 6, 4, 4))

    # Plot area: a header bar (page selector, title, export buttons) and
    # three page layouts in collapsible rows.
    bar = GridLayout(plots[1, 1])
    page_menu = Menu(bar[1, 1]; options=PAGES, default="Deviations", width=150, tellwidth=true)
    Label(bar[1, 2], ""; tellwidth=false)
    export_fig_btn = Button(bar[1, 3]; label="Export figure…", tellwidth=true)
    export_tsv_btn = Button(bar[1, 4]; label="Export table…", tellwidth=true)
    export_rec_btn = Button(bar[1, 5]; label="Export record…", tellwidth=true)
    title_lbl = Label(bar[2, 1:5], "No analysis yet."; halign=:left, tellwidth=false, color=:gray30)

    dev_page = GridLayout(plots[2, 1])
    spec_page = GridLayout(plots[3, 1])
    dyn_page = GridLayout(plots[4, 1])
    ax_y = Axis(dev_page[1, 1]; xscale=log10, yscale=log10,
                xlabel="Averaging time τ (s)", ylabel="σ_y(τ)")
    ax_x = Axis(dev_page[2, 1]; xscale=log10, yscale=log10,
                xlabel="Averaging time τ (s)", ylabel="σ_x(τ) (s)")
    table_menu = Menu(dev_page[3, 1][1, 1]; options=["adev"], default="adev", width=140,
                      tellwidth=true, halign=:left)
    table_grid = GridLayout(dev_page[3, 1][2, 1]; halign=:left, valign=:top)
    rowsize!(dev_page, 3, Auto(0.8))

    let opts = GridLayout(spec_page[1, 1]; halign=:left)
        opts[1, 1] = Label(fig, "Carrier frequency for ℒ(f) (Hz)"; color=:gray25)
        opts[1, 2] = carrier_box
        opts[1, 3] = Label(fig, "(re-run to apply)"; color=:gray50)
    end
    ax_sy = Axis(spec_page[2, 1]; xscale=log10, yscale=log10, ylabel="S_y(f) (1/Hz)")
    ax_sx = Axis(spec_page[3, 1]; xscale=log10, yscale=log10, ylabel="S_x(f) (s²/Hz)")
    ax_l  = Axis(spec_page[4, 1]; xscale=log10, xlabel="Frequency f (Hz)", ylabel="ℒ(f) (dBc/Hz)")

    let opts = GridLayout(dyn_page[1, 1:2]; halign=:left)
        opts[1, 1] = Label(fig, "Estimator"; color=:gray25); opts[1, 2] = dyn_menu
        opts[1, 3] = Label(fig, "window / record length"; color=:gray25); opts[1, 4] = frac_menu
        opts[1, 5] = Label(fig, "(re-run to apply)"; color=:gray50)
    end
    ax_d = Axis(dyn_page[2, 1]; xlabel="Time t (s)", ylabel="τ (s)")
    cb_d = Ref{Any}(nothing)

    pages = Dict("Deviations" => (dev_page, 2), "Spectrum" => (spec_page, 3), "Dynamic map" => (dyn_page, 4))
    legends = Dict{String,Vector{Makie.Block}}(p => Makie.Block[] for p in PAGES)
    function show_page(name)
        for (p, (gl, r)) in pages
            on_page = p == name
            set_visible!(all_blocks(gl), on_page)
            set_visible!(legends[p], on_page)
            rowsize!(plots, r, on_page ? Auto(1) : Fixed(0))
        end
        cb_d[] === nothing || set_visible!([cb_d[]], name == "Dynamic map")
    end
    on(page_menu.selection) do name; name === nothing || show_page(name); end

    state = (record=Ref{Any}(nothing), suite=Ref{Any}(nothing),
             spectra=Ref{Any}(nothing), dynamic=Ref{Any}(nothing))

    say(msg; error=false) = (status.text = msg; status.color = error ? :firebrick : :gray30)
    text(tb) = something(tb.stored_string[], "")

    function current_record()
        rec = if source.selection[] == "File"
            load_record(text(path_box), kind.selection[], text(tau0_box),
                        parse(Int, text(time_col)), parse(Int, text(value_col)), parse(Int, text(header)))
        else
            levels = [(α, noise_sliders[α].value[]) for (α, _, _) in NOISES if noise_toggles[α].active[]]
            synthesize(parse(Int, n_menu.selection[]), parse(Float64, text(tau0_box)),
                       parse(Int, text(seed_box)), levels)
        end
        return preprocess(rec; outliers=outliers.active[], gaps=gaps.active[],
                          detrend_mode=Symbol(detrend_menu.selection[]))
    end

    function run()
        try
            rec = current_record()
            devs = [d for d in DEVIATIONS if dev_toggles[d].active[]]
            taumode = Dict(TAU_MODES)[tau_menu.selection[]]
            suite = analyze(rec, devs, taumode, ci_toggle.active[],
                            parse(Float64, conf_menu.selection[]), bias_toggle.active[])
            sp  = spectra(rec, parse(Float64, text(carrier_box)))
            dyn = dynamic_map(rec, Symbol(dyn_menu.selection[]),
                              1 / parse(Int, split(frac_menu.selection[], "/")[2]))
            state.record[] = rec; state.suite[] = suite; state.spectra[] = sp; state.dynamic[] = dyn

            foreach(delete!, legends["Deviations"]); empty!(legends["Deviations"])
            sy = filter(r -> !(r.deviation_type in SIGMA_X), suite.results)
            sx = filter(r ->   r.deviation_type in SIGMA_X,  suite.results)
            for (ax, rs, pos) in ((ax_y, sy, :lb), (ax_x, sx, :lt))
                leg = draw_results!(ax, rs; band=band_toggle.active[], legend_position=pos)
                leg === nothing || push!(legends["Deviations"], leg)
            end
            rowsize!(dev_page, 2, isempty(sx) ? Fixed(0) : Auto(1))
            set_visible!([ax_x], !isempty(sx))
            names = string.(keys(suite))
            table_menu.options = names
            table_menu.i_selected = 1
            fill_table!(table_grid, suite[1])

            draw_spectrum!(ax_sy, sp.sy; color=Makie.wong_colors()[1])
            draw_spectrum!(ax_sx, sp.sx; color=Makie.wong_colors()[2])
            sp.l === nothing ? empty!(ax_l) : draw_spectrum!(ax_l, sp.l; color=Makie.wong_colors()[3])
            draw_dynamic!(ax_d, dyn; cb=cb_d, cb_cell=dyn_page[2, 2])

            src = rec.source == "user" ? "synthetic record" : basename(rec.source)
            title_lbl.text = "$src: $(suite.n) samples, τ₀ = $(suite.tau0) s, CI $(suite.confidence)"
            show_page(page_menu.selection[])
            say("Computed " * join(names, ", ") * " on $(suite.n) samples.")
        catch err
            say("Error: " * sprint(showerror, err); error=true)
        end
    end

    on(run_btn.clicks) do _; run(); end
    on(table_menu.selection) do name
        s = state.suite[]
        (s === nothing || name === nothing) && return
        fill_table!(table_grid, s[Symbol(name)])
    end
    on(open_btn.clicks) do _
        p = pick_file(; filterlist="dat,txt,csv,tsv")
        isempty(p) && return
        path_box.displayed_string[] = p; path_box.stored_string[] = p
        source.i_selected = 2
    end
    on(export_fig_btn.clicks) do _
        state.suite[] === nothing && (say("Run an analysis first."; error=true); return)
        p = save_file("sigmatau_$(lowercase(replace(page_menu.selection[], " " => "_"))).png";
                      filterlist="png,pdf,svg")
        isempty(p) && return
        try
            exfig = export_figure(page_menu.selection[], state,
                                  (; title=title_lbl.text[], band=band_toggle.active[]))
            Makie.save(p, exfig; backend=CairoMakie, px_per_unit=2)
            say("Wrote $p.")
        catch err
            say("Error: " * sprint(showerror, err); error=true)
        end
    end
    on(export_tsv_btn.clicks) do _
        state.suite[] === nothing && (say("Run an analysis first."; error=true); return)
        p = save_file("sigmatau_suite.tsv"; filterlist="tsv,txt")
        isempty(p) || (SigmaTau.save(p, state.suite[]); say("Wrote $p."))
    end
    on(export_rec_btn.clicks) do _
        state.record[] === nothing && (say("Run an analysis first."; error=true); return)
        p = save_file("sigmatau_record.txt"; filterlist="txt,csv,dat")
        isempty(p) || (SigmaTau.save(p, state.record[]); say("Wrote $p."))
    end

    DataInspector(fig)
    show_page("Deviations")
    return (; fig, run, state, show_page, export_figure=(page) -> export_figure(page, state,
            (; title=title_lbl.text[], band=band_toggle.active[])),
            controls=(; dev_toggles, band_toggle, source, page_menu))
end

# ── Entry point ───────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    app = build()
    if "--smoke" in ARGS
        app.controls.dev_toggles[:tdev].active[] = true
        app.controls.dev_toggles[:mtie].active[] = true
        app.run()
        s = app.state.suite[]
        s === nothing && error("smoke test: no suite produced")
        dir = get(ENV, "SIGMATAU_SMOKE_DIR", tempdir())
        for page in PAGES
            app.controls.page_menu.i_selected = findfirst(==(page), PAGES)
            Makie.save(joinpath(dir, "explorer_" * replace(lowercase(page), " " => "_") * ".png"),
                 app.fig; px_per_unit=1)
            Makie.save(joinpath(dir, "export_" * replace(lowercase(page), " " => "_") * ".pdf"),
                 app.export_figure(page); backend=CairoMakie)
        end
        println("smoke test OK: ", join(string.(keys(s)), ", "), " → ", dir)
    else
        screen = display(app.fig)
        wait(screen)
    end
end

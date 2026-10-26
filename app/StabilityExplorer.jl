# SigmaTau stability explorer: a desktop window for loading or synthesizing a
# record, choosing estimators and options, and reading the plots and tables.
#
#   julia --project=app app/StabilityExplorer.jl            # open the window
#   julia --project=app app/StabilityExplorer.jl --smoke    # headless self-test
#
# Pages: Record (data source, synthesis and the raw time series), Deviations
# (σ_y and σ_x curves with χ² intervals and the per-τ table), Spectrum (S_y,
# S_x, ℒ(f)), Dynamic deviation (a waterfall of σ(τ) per window) and Help. Export renders a clean figure
# through CairoMakie (PNG, PDF, SVG) and writes tables and records as the
# text files `load_suite` and `read_phase` read back. Every number shown comes
# from the public SigmaTau API.

using SigmaTau
using GLMakie
using CairoMakie
using NativeFileDialog
using Random
using Printf
using Statistics

GLMakie.activate!()

const VERSION_STRING = string(pkgversion(SigmaTau))
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
const PAGES      = ["Record", "Deviations", "Spectrum", "Dynamic deviation", "Help"]
const MAX_TABLE_ROWS = 40
const MAX_SERIES_POINTS = 20_000

const HELP_SECTIONS = [
    ("Workflow",
     "1. Record page: synthesize a power-law record or load a file (Browse…, then Load file). " *
     "The time series and its summary appear below the controls.\n" *
     "2. Choose deviations, τ grid, confidence level and preprocessing in the left column.\n" *
     "3. Press Run analysis, or Enter in a field (on the Record page Enter also reloads the record). The Deviations, Spectrum and Dynamic deviation pages " *
     "fill in and the status line at the bottom reports what was computed.\n" *
     "4. Export figure… writes the current page as PNG, PDF or SVG. Export table… writes the suite " *
     "as a tab-separated file. Export record… writes the preprocessed record as text.\n" *
     "Keyboard: Tab and Shift+Tab move between the fields of a page, Enter in a field or F5 runs the " *
     "analysis, Page Up and Page Down switch pages, F1 opens this page. Start with --fontsize=14 or 16 " *
     "for larger text."),
    ("Deviations",
     "Frequency deviations are dimensionless; time deviations (tdev, htdev, ttotdev, mtie, tierms) are in " *
     "seconds and plot on their own axes. Hadamard estimators remove linear frequency drift; modified " *
     "estimators separate white PM from flicker PM."),
    ("Confidence intervals",
     "χ² bounds from the Greenhall–Riley equivalent degrees of freedom, with the noise type identified " *
     "at each τ by the lag-1 autocorrelation test. Bias correction applies the NIST SP1065 / Howe " *
     "unbias factors to the total family and ThêoBR removal to Thêo1; switch it off to reproduce " *
     "allantools (totals) or raw Thêo1. mtie and tierms have no interval model."),
    ("File formats",
     "Stable32 .DAT exports: header rows 10, time column 0, value column 1, τ₀ from the header's Tau " *
     "line. Two-column time/value text or CSV: τ₀ is inferred from the time column. Files written by " *
     "Export record… need no settings."),
    ("Synthetic records",
     "Each term is a power-law noise specified by its deviation at the sample interval, σ(τ₀); the " *
     "slider sets the exponent of 10. Terms are independent and the seed makes the record reproducible."),
    ("Links",
     "Documentation: https://ianlap.github.io/SigmaTau.jl\nSource: https://github.com/ianlap/SigmaTau.jl"),
]
const GLOSSARY = [
    ("adev", "overlapping Allan"), ("tdev", "time deviation, τ/√3 · MDEV"),
    ("mdev", "modified Allan"), ("htdev", "Hadamard time deviation"),
    ("hdev", "overlapping Hadamard"), ("ttotdev", "time total deviation"),
    ("mhdev", "modified Hadamard"), ("mtie", "maximum time interval error"),
    ("totdev", "total"), ("tierms", "RMS time interval error"),
    ("mtotdev", "modified total"), ("theo1", "Thêo1 (ThêoBR when bias-corrected)"),
    ("htotdev", "Hadamard total"), ("theoh", "ThêoH composite of ADEV and ThêoBR"),
    ("mhtotdev", "modified Hadamard total"), ("pdev", "parabolic deviation"),
]

# ── Analysis (no GUI state) ───────────────────────────────────────────────

nsamples(rec::PhaseData) = length(rec.x)
nsamples(rec::FrequencyData) = length(rec.y)
samples(rec::PhaseData) = rec.x
samples(rec::FrequencyData) = rec.y

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

# Summary lines for the Record page.
function record_summary(rec)
    v = samples(rec)
    n = length(v)
    T = (n - 1) * rec.tau0
    y = rec isa PhaseData ? diff(v) ./ rec.tau0 : v
    finite = filter(isfinite, y)
    drift = length(finite) > 2 ? begin
        t = (0:length(finite)-1) .* rec.tau0
        cov(t, finite) / var(t)
    end : NaN
    kind = rec isa PhaseData ? "phase x(t), seconds" : "fractional frequency y(t)"
    span = T >= 86400 ? @sprintf("%.2f d", T / 86400) : T >= 3600 ? @sprintf("%.2f h", T / 3600) : @sprintf("%.0f s", T)
    return [
        "Record: $kind" * (rec.source == "user" ? ", synthetic" : ", " * basename(rec.source)),
        @sprintf("N = %d samples, τ₀ = %g s, span T = %s", n, rec.tau0, span),
        @sprintf("mean y = %.3e, std y = %.3e, NaN samples = %d", mean(finite), std(finite), count(isnan, v)),
        @sprintf("least-squares frequency drift = %.3e per s (%.3e per day)", drift, drift * 86400),
    ]
end

# ── Drawing (shared by the window and the export figure) ──────────────────

# Axis and table labels with real subscripts.
const LBL_SIGMA_Y   = rich("σ", subscript("y"), "(τ)")
const LBL_SIGMA_X = rich("σ", subscript("x"), "(τ)  (s)")
const LBL_SIGMA_MIN = rich("σ", subscript("min"))
const LBL_SIGMA_MAX = rich("σ", subscript("max"))
const N_EFF     = rich("n", subscript("eff"))
const S_Y       = rich("S", subscript("y"), "(f)  (1/Hz)")
const S_X       = rich("S", subscript("x"), "(f)  (s²/Hz)")

# Integer-decade ticks for log axes (10⁰, 10¹, …); Makie's default picker
# skips decades on wide ranges. Falls back to the linear picker inside a decade.
struct DecadeTicks end
function Makie.get_tickvalues(::DecadeTicks, vmin, vmax)
    ks = ceil(Int, vmin):floor(Int, vmax)
    length(ks) >= 2 && return Float64.(ks)
    return Makie.get_tickvalues(WilkinsonTicks(5, k_min=3), vmin, vmax)
end
const LOG_TICKS = LogTicks(DecadeTicks())
log_axis_kw() = (xticks=LOG_TICKS, yticks=LOG_TICKS,
                 xminorticksvisible=true, yminorticksvisible=true,
                 xminorticks=IntervalsBetween(9), yminorticks=IntervalsBetween(9),
                 xminorgridvisible=true, yminorgridvisible=true,
                 xminorgridstyle=:dot, yminorgridstyle=:dot)

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

# Dynamic deviation as the SP1065 waterfall: one σ(τ) curve per window time,
# stacked along the time axis, coloured by window centre. Axis3 has no log
# scales, so log10(τ) and log10(σ) are plotted and the decade ticks relabelled.
# `cb` is a Ref holding the previous Colorbar (or nothing) and `cb_cell` its
# grid position; the colorbar is rebuilt on every draw.
decade_labels(ks) = [Makie.rich("10", Makie.superscript(string(k))) for k in ks]
dynamic_axis(gp) = Axis3(gp; xlabel="τ (s)", ylabel="Time t (s)", zlabel="σ",
                         azimuth=1.6π, elevation=0.2π, perspectiveness=0.3,
                         ytickformat=vs -> [@sprintf("%g", v) for v in vs],
                         xlabeloffset=36, ylabeloffset=48, zlabeloffset=56)
function draw_dynamic!(ax, dyn; cb=nothing, cb_cell=nothing)
    empty!(ax)
    z = map(v -> (isfinite(v) && v > 0) ? log10(v) : NaN, dyn.dev)
    lx = log10.(dyn.tau)
    tmin, tmax = Float64.(extrema(dyn.t))
    cmap = cgrad(:viridis)
    for (i, tc) in enumerate(dyn.t)
        c = tmax > tmin ? cmap[(tc - tmin) / (tmax - tmin)] : cmap[0.5]
        lines!(ax, lx, fill(Float64(tc), length(lx)), z[i, :]; color=c, linewidth=1.6)
    end
    kx = floor(Int, minimum(lx)):ceil(Int, maximum(lx))
    zfin = filter(isfinite, z)
    kz = isempty(zfin) ? (0:0) : (floor(Int, minimum(zfin)):ceil(Int, maximum(zfin)))
    ax.xticks = (Float64.(kx), decade_labels(kx))
    ax.zticks = (Float64.(kz), decade_labels(kz))
    ax.title = uppercase(string(dyn.deviation_type)) * ", window = $(dyn.window) samples, $(length(dyn.t)) windows"
    autolimits!(ax)
    if cb !== nothing
        cb[] === nothing || delete!(cb[])
        cb[] = Colorbar(cb_cell; colormap=:viridis, limits=(tmin, max(tmax, tmin + 1)),
                        label="window centre t (s)", height=Relative(0.6))
    end
    return ax
end

function draw_record!(ax, rec)
    empty!(ax)
    v = samples(rec)
    n = length(v)
    step = max(1, cld(n, MAX_SERIES_POINTS))
    idx = 1:step:n
    t = (idx .- 1) .* rec.tau0
    lines!(ax, t, v[idx]; color=Makie.wong_colors()[1], linewidth=1)
    ax.ylabel = rec isa PhaseData ? "x(t) (s)" : "y(t)"
    ax.title = step > 1 ? "every $(step)th sample shown" : ""
    autolimits!(ax)
end

function fill_table!(grid, r::StabilityResult)
    foreach(delete!, contents(grid))
    trim!(grid)
    has_ci = !isempty(r.edf)
    header = has_ci ? ["τ (s)", "σ", N_EFF, LBL_SIGMA_MIN, LBL_SIGMA_MAX, "EDF", "noise"] :
                      ["τ (s)", "σ", N_EFF]
    ncol = length(header)
    g(v) = isfinite(v) ? @sprintf("%.5g", v) : "NaN"
    Box(grid[1, 1:ncol]; color=:gray85, strokevisible=false)
    for (j, h) in enumerate(header)
        Label(grid[1, j], h; font=:bold, halign=:right, tellwidth=true, color=:gray15,
              padding=(6, 6, 3, 3))
    end
    nrows = min(length(r.tau), MAX_TABLE_ROWS)
    for i in 1:nrows
        iseven(i) && Box(grid[i + 1, 1:ncol]; color=:gray94, strokevisible=false)
        cells = has_ci ?
            [g(r.tau[i]), g(r.dev[i]), string(r.neff[i]), g(r.ci[i].lo),
             g(r.ci[i].hi), g(r.edf[i]), string(r.noise_type[i])] :
            [g(r.tau[i]), g(r.dev[i]), string(r.neff[i])]
        for (j, c) in enumerate(cells)
            Label(grid[i + 1, j], c; halign=:right, tellwidth=true, color=:gray10,
                  padding=(6, 6, 1, 1))
        end
    end
    if length(r.tau) > nrows
        Label(grid[nrows + 2, 1:ncol],
              "$(length(r.tau) - nrows) more rows not shown; Export table… writes all of them";
              halign=:left, tellwidth=false, color=:gray35)
    end
    colgap!(grid, 0)
    rowgap!(grid, 0)
end

# Publication figure for a page, rendered through CairoMakie.
function export_figure(page, state, opts)
    fig = Figure(size=(900, page == "Spectrum" ? 900 : page == "Dynamic deviation" ? 640 : 560), fontsize=13)
    if page == "Deviations"
        suite = state.suite[]
        sy = filter(r -> !(r.deviation_type in SIGMA_X), suite.results)
        sx = filter(r ->   r.deviation_type in SIGMA_X,  suite.results)
        row = 1
        if !isempty(sy)
            ax = Axis(fig[row, 1]; xscale=log10, yscale=log10, log_axis_kw()...,
                      xlabel="Averaging time τ (s)", ylabel=LBL_SIGMA_Y, title=opts.title)
            draw_results!(ax, sy; band=opts.band); row += 1
        end
        if !isempty(sx)
            ax = Axis(fig[row, 1]; xscale=log10, yscale=log10, log_axis_kw()...,
                      xlabel="Averaging time τ (s)", ylabel=LBL_SIGMA_X,
                      title=isempty(sy) ? opts.title : "")
            draw_results!(ax, sx; band=opts.band, legend_position=:lt)
        end
    elseif page == "Spectrum"
        sp = state.spectra[]
        ax1 = Axis(fig[1, 1]; xscale=log10, yscale=log10, log_axis_kw()..., ylabel=S_Y, title=opts.title)
        ax2 = Axis(fig[2, 1]; xscale=log10, yscale=log10, log_axis_kw()..., ylabel=S_X)
        draw_spectrum!(ax1, sp.sy; color=Makie.wong_colors()[1])
        draw_spectrum!(ax2, sp.sx; color=Makie.wong_colors()[2])
        if sp.l !== nothing
            ax3 = Axis(fig[3, 1]; xscale=log10, xticks=LOG_TICKS, xminorticksvisible=true,
                       xminorticks=IntervalsBetween(9), xlabel="Frequency f (Hz)", ylabel="ℒ(f) (dBc/Hz)")
            draw_spectrum!(ax3, sp.l; color=Makie.wong_colors()[3])
        else
            ax2.xlabel = "Frequency f (Hz)"
        end
    elseif page == "Dynamic deviation"
        ax = dynamic_axis(fig[1, 1]); ax.title = opts.title
        draw_dynamic!(ax, state.dynamic[]; cb=Ref{Any}(nothing), cb_cell=fig[1, 2])
    else
        ax = Axis(fig[1, 1]; xlabel="Time t (s)", title=opts.title)
        draw_record!(ax, state.record[])
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
        b isa Union{Axis,Axis3} && (b.scene.visible[] = v)
    end
end

function build(; size=(1500, 940), fontsize=12)
    fig = Figure(; size=size, fontsize=fontsize, backgroundcolor=:gray98)
    controls = GridLayout(fig[1, 1]; tellheight=false, valign=:top)
    plots    = GridLayout(fig[1, 2]; tellwidth=false, tellheight=false)
    colsize!(fig.layout, 1, Fixed(round(Int, 430 * fontsize / 12)))
    colgap!(fig.layout, 24)
    rowgap!(controls, 2)
    colgap!(controls, 8)

    const_h = round(Int, 2.2 * fontsize)
    hint = fontsize - 1
    small(str) = Textbox(fig; stored_string=str, width=62, height=const_h, tellwidth=false, halign=:left)
    menu(; kw...) = Menu(fig; height=const_h, kw...)
    textbox(; kw...) = Textbox(fig; height=const_h, kw...)
    caption(text) = Label(fig, text; halign=:right, tellwidth=true, color=:gray25)

    row = Ref(1)
    nextrow() = (r = row[]; row[] += 1; r)
    function section(text)
        r = nextrow()
        Box(controls[r, 1:2]; color=:gray88, strokevisible=false, cornerradius=4)
        Label(controls[r, 1:2], text; font=:bold, halign=:left, tellwidth=false,
              padding=(8, 8, 4, 4), color=:gray15)
    end
    labelled(text, block) = (r = nextrow(); controls[r, 1] = caption(text); controls[r, 2] = block; block)
    function labelled_row(text, blocks...)
        r = nextrow()
        controls[r, 1] = caption(text)
        inner = GridLayout(controls[r, 2]; halign=:left)
        for (j, b) in enumerate(blocks); inner[1, j] = b; end
    end

    # ── Analysis column ──
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
    Label(controls[nextrow(), 1:2], rich("σ", subscript("y"), ": adev … theoh, pdev.   σ", subscript("x"), " (seconds): tdev, htdev, ttotdev, mtie, tierms.");
          halign=:left, tellwidth=false, color=:gray35, fontsize=hint)

    section("Analysis")
    tau_menu    = labelled("τ grid", menu(; options=first.(TAU_MODES), default="Octave", width=170))
    ci_toggle   = Toggle(fig; active=true)
    conf_menu   = menu(; options=["0.683", "0.90", "0.95", "0.99"], default="0.683", width=84)
    labelled_row("Confidence intervals, level", ci_toggle, conf_menu)
    bias_toggle = Toggle(fig; active=true)
    band_toggle = Toggle(fig; active=false)
    labelled_row("Bias-correct totals, Thêo1", bias_toggle, Label(fig, "CI as band"; color=:gray25), band_toggle)

    section("Preprocessing")
    detrend_menu = labelled("Detrend", menu(; options=["none", "mean", "linear", "quadratic"], default="none", width=170))
    outliers = Toggle(fig; active=false)
    gaps     = Toggle(fig; active=false)
    labelled_row("Remove 5σ outliers", outliers, Label(fig, "fill NaN gaps"; color=:gray25), gaps)

    run_btn = Button(fig; label="Run analysis", buttoncolor=RGBf(0.33, 0.47, 0.72),
                     labelcolor=:white, buttoncolor_hover=RGBf(0.4, 0.55, 0.8),
                     font=:bold, height=const_h + 6)
    controls[nextrow(), 1:2] = run_btn
    Label(controls[nextrow(), 1:2], "Enter in a field also runs the analysis.";
          halign=:center, tellwidth=false, color=:gray35, fontsize=hint)

    status = Label(fig[2, 1:2], "Ready."; halign=:left, tellwidth=false, color=:gray30,
                   padding=(6, 6, 4, 4))

    # ── Plot area: header bar and pages in collapsible rows ──
    bar = GridLayout(plots[1, 1])
    rowgap!(plots, 0)
    page_menu = Menu(bar[1, 1]; options=PAGES, default="Record", width=150, tellwidth=true)
    Label(bar[1, 2], ""; tellwidth=false)
    export_fig_btn = Button(bar[1, 3]; label="Export figure…", tellwidth=true)
    export_tsv_btn = Button(bar[1, 4]; label="Export table…", tellwidth=true)
    export_rec_btn = Button(bar[1, 5]; label="Export record…", tellwidth=true)
    title_lbl = Label(bar[2, 1:5], "No record loaded."; halign=:left, tellwidth=false, color=:gray30)

    # Record page: data source and synthesis controls, then the time series.
    rec_page = GridLayout(plots[2, 1]; alignmode=Outside())
    src = GridLayout(rec_page[1, 1]; valign=:top)
    Box(src[1, 1:3]; color=:gray88, strokevisible=false, cornerradius=4)
    Label(src[1, 1:3], "Load a file"; font=:bold, halign=:left, tellwidth=false, padding=(8, 8, 4, 4), color=:gray15)
    path_box = textbox(; placeholder="no file chosen", tellwidth=false, halign=:left)
    open_btn = Button(fig; label="Browse…", height=const_h)
    src[2, 1] = Label(fig, "File"; halign=:right, color=:gray25)
    src[2, 2:3] = path_box
    kind = menu(; options=["phase", "frequency"], default="phase", width=110, tellwidth=false, halign=:left)
    src[3, 1] = Label(fig, "Record type"; halign=:right, color=:gray25)
    src[3, 2] = kind
    tau0_box = Textbox(fig; stored_string="1.0", width=62, height=const_h, tellwidth=false, halign=:left)
    time_col = small("0"); value_col = small("1"); header = small("0")
    src[4, 1] = Label(fig, "τ₀ (s)"; halign=:right, color=:gray25)
    src[4, 2] = tau0_box
    src[5, 1] = Label(fig, "Time, value columns"; halign=:right, color=:gray25)
    let g = GridLayout(src[5, 2]; tellwidth=false, halign=:left); g[1, 1] = time_col; g[1, 2] = value_col; end
    src[6, 1] = Label(fig, "Header rows"; halign=:right, color=:gray25)
    src[6, 2] = header
    Label(src[7, 1:3], "Stable32 .DAT: header rows 10, time column 0, value column 1.";
          halign=:left, tellwidth=false, color=:gray35, fontsize=hint)
    load_btn = Button(fig; label="Load file", height=const_h)
    let g = GridLayout(src[8, 2:3]; halign=:right, tellwidth=false); g[1, 1] = open_btn; g[1, 2] = load_btn; end

    syn = GridLayout(rec_page[1, 2]; valign=:top)
    Box(syn[1, 1:2]; color=:gray88, strokevisible=false, cornerradius=4)
    Label(syn[1, 1:2], rich("Synthesize a record  (level = σ", subscript("y"), " at τ₀)"); font=:bold, halign=:left,
          tellwidth=false, padding=(8, 8, 4, 4), color=:gray15)
    n_menu   = menu(; options=["4096", "8192", "16384", "32768", "86400"], default="16384", width=90, tellwidth=false, halign=:left)
    seed_box = Textbox(fig; stored_string="42", width=62, height=const_h, tellwidth=false, halign=:left)
    syn[2, 1] = Label(fig, "Samples, seed"; halign=:right, color=:gray25)
    let g = GridLayout(syn[2, 2]; tellwidth=false, halign=:left); g[1, 1] = n_menu; g[1, 2] = seed_box; colsize!(g, 1, Fixed(90)); colsize!(g, 2, Fixed(62)); end
    noise_toggles = Dict{Int,Toggle}()
    noise_sliders = Dict{Int,Slider}()
    for (k, (α, name, start)) in enumerate(NOISES)
        noise_toggles[α] = Toggle(fig; active=(α in (2, 0, -2)))
        noise_sliders[α] = Slider(fig; range=-18:-9, startvalue=start, tellwidth=false)
        syn[2 + k, 1] = Label(fig, name; halign=:right, color=:gray25)
        g = GridLayout(syn[2 + k, 2]; tellwidth=false)
        g[1, 1] = noise_toggles[α]
        g[1, 2] = noise_sliders[α]
        g[1, 3] = Label(fig, lift(v -> @sprintf("1e%d", v), noise_sliders[α].value); width=40,
                        halign=:left, color=:gray25)
    end
    synth_btn = Button(fig; label="Synthesize", height=const_h, tellwidth=false, halign=:right)
    syn[3 + length(NOISES), 2] = synth_btn
    colsize!(rec_page, 1, Relative(0.5))
    colgap!(rec_page, 24)

    ax_rec = Axis(rec_page[2, 1:2]; xlabel="Time t (s)", xtickformat=vs -> [@sprintf("%g", v) for v in vs])
    summary_lbls = [Label(rec_page[2 + i, 1:2], ""; halign=:left, tellwidth=false, color=:gray20)
                    for i in 1:4]
    rowsize!(rec_page, 2, Auto(1.0))
    rowgap!(rec_page, 6)

    # Deviations page.
    dev_page = GridLayout(plots[3, 1]; alignmode=Outside())
    ax_y = Axis(dev_page[1, 1]; xscale=log10, yscale=log10, log_axis_kw()...,
                xlabel="Averaging time τ (s)", ylabel=LBL_SIGMA_Y)
    ax_x = Axis(dev_page[2, 1]; xscale=log10, yscale=log10, log_axis_kw()...,
                xlabel="Averaging time τ (s)", ylabel=LBL_SIGMA_X)
    tbl = GridLayout(dev_page[1:2, 2]; halign=:left, valign=:top)
    tbl_head = GridLayout(tbl[1, 1]; halign=:left)
    table_menu = Menu(tbl_head[1, 1]; options=["adev"], default="adev", width=140, tellwidth=true, halign=:left)
    tbl_head[1, 2] = Label(fig, "per-τ table"; halign=:left, color=:gray35)
    table_grid = GridLayout(tbl[2, 1]; halign=:left, valign=:top)
    colgap!(dev_page, 20)

    # Spectrum page.
    spec_page = GridLayout(plots[4, 1]; alignmode=Outside())
    carrier_box = textbox(; stored_string="10e6", width=90)
    let opts = GridLayout(spec_page[1, 1]; halign=:left, tellwidth=false)
        opts[1, 1] = Label(fig, "Carrier frequency for ℒ(f) (Hz)"; color=:gray25)
        opts[1, 2] = carrier_box
        opts[1, 3] = Label(fig, "(takes effect at the next run)"; color=:gray35)
    end
    ax_sy = Axis(spec_page[2, 1]; xscale=log10, yscale=log10, log_axis_kw()..., ylabel=S_Y)
    ax_sx = Axis(spec_page[3, 1]; xscale=log10, yscale=log10, log_axis_kw()..., ylabel=S_X)
    ax_l  = Axis(spec_page[4, 1]; xscale=log10, xticks=LOG_TICKS, xminorticksvisible=true,
                 xminorticks=IntervalsBetween(9), xlabel="Frequency f (Hz)", ylabel="ℒ(f) (dBc/Hz)")

    # Dynamic deviation page.
    dyn_page = GridLayout(plots[5, 1]; alignmode=Outside())
    dyn_menu  = menu(; options=["dadev", "dhdev"], default="dadev", width=84)
    frac_menu = menu(; options=["1/4", "1/8", "1/16", "1/32"], default="1/8", width=72)
    let opts = GridLayout(dyn_page[1, 1:2]; halign=:left, tellwidth=false)
        opts[1, 1] = Label(fig, "Estimator"; color=:gray25); opts[1, 2] = dyn_menu
        opts[1, 3] = Label(fig, "window as a fraction of the record"; color=:gray25); opts[1, 4] = frac_menu
        opts[1, 5] = Label(fig, "(takes effect at the next run; drag the plot to rotate)"; color=:gray35)
    end
    ax_d = dynamic_axis(dyn_page[2, 1])
    cb_d = Ref{Any}(nothing)

    # Help page: two columns of sections, the glossary under "Deviations".
    help_page = GridLayout(plots[6, 1]; tellwidth=false, tellheight=false, valign=:top, alignmode=Outside())
    Label(help_page[1, 1:2], "SigmaTau explorer " * VERSION_STRING; font=:bold, fontsize=16,
          halign=:left, tellwidth=false, color=:gray15, padding=(0, 0, 10, 0))
    help_cols = (GridLayout(help_page[2, 1]; valign=:top, tellwidth=false, tellheight=false),
                 GridLayout(help_page[2, 2]; valign=:top, tellwidth=false, tellheight=false))
    colgap!(help_page, 40)
    let rows = [1, 1]
        for (title, body) in HELP_SECTIONS
            c = title in ("Workflow", "Deviations") ? 1 : 2
            col = help_cols[c]
            Label(col[rows[c], 1], title; font=:bold, halign=:left, tellwidth=false, color=:gray15,
                  padding=(0, 0, 2, 12)); rows[c] += 1
            Label(col[rows[c], 1], body; halign=:left, justification=:left, tellwidth=false,
                  word_wrap=true, color=:gray15, lineheight=1.15); rows[c] += 1
            if title == "Deviations"
                g = GridLayout(col[rows[c], 1]; halign=:left, tellwidth=false); rows[c] += 1
                for (i, (name, desc)) in enumerate(GLOSSARY)
                    Label(g[i, 1], name; font=:bold, halign=:left, color=:gray15, padding=(0, 10, 1, 1))
                    Label(g[i, 2], desc; halign=:left, color=:gray25, padding=(0, 0, 1, 1))
                end
                rowgap!(g, 0)
            end
        end
        foreach(col -> rowgap!(col, 4), help_cols)
    end

    pages = Dict("Record" => (rec_page, 2), "Deviations" => (dev_page, 3), "Spectrum" => (spec_page, 4),
                 "Dynamic deviation" => (dyn_page, 5), "Help" => (help_page, 6))
    legends = Dict{String,Vector{Makie.Block}}(p => Makie.Block[] for p in PAGES)
    sx_shown = Ref(true)          # whether the σ_x axis currently has content
    function show_page(name)
        for (p, (gl, r)) in pages
            on_page = p == name
            set_visible!(all_blocks(gl), on_page)
            set_visible!(legends[p], on_page)
            rowsize!(plots, r, on_page ? Auto(1) : Fixed(0))
            rowgap!(plots, r - 1, on_page ? 16 : 0)
        end
        cb_d[] === nothing || set_visible!([cb_d[]], name == "Dynamic deviation")
        sx_shown[] || set_visible!([ax_x], false)     # keep the empty σ_x axis collapsed
    end
    on(page_menu.selection) do name; name === nothing || show_page(name); end

    state = (record=Ref{Any}(nothing), suite=Ref{Any}(nothing),
             spectra=Ref{Any}(nothing), dynamic=Ref{Any}(nothing), last_dir=Ref(pwd()))

    say(msg; error=false) = (status.text = msg; status.color = error ? :firebrick : :gray30)
    # What the field shows counts, Enter or not (Makie only commits on Enter).
    text(tb) = something(tb.displayed_string[], "")
    parse_int(tb, what) = (v = tryparse(Int, text(tb)); v === nothing ? error("$what must be an integer, got \"$(text(tb))\".") : v)
    parse_float(tb, what) = (v = tryparse(Float64, text(tb)); v === nothing ? error("$what must be a number, got \"$(text(tb))\".") : v)

    function show_record!(rec)
        state.record[] = rec
        draw_record!(ax_rec, rec)
        for (lbl, line) in zip(summary_lbls, record_summary(rec))
            lbl.text = line
        end
        src_name = rec.source == "user" ? "synthetic record" : basename(rec.source)
        title_lbl.text = "$src_name: $(nsamples(rec)) samples, τ₀ = $(rec.tau0) s"
    end

    function load_data(from_file)
        rec = if from_file
            load_record(text(path_box), kind.selection[], text(tau0_box),
                        parse_int(time_col, "time column"), parse_int(value_col, "value column"),
                        parse_int(header, "header rows"))
        else
            levels = [(α, noise_sliders[α].value[]) for (α, _, _) in NOISES if noise_toggles[α].active[]]
            synthesize(parse(Int, n_menu.selection[]), parse_float(tau0_box, "τ₀"),
                       parse_int(seed_box, "seed"), levels)
        end
        rec = preprocess(rec; outliers=outliers.active[], gaps=gaps.active[],
                         detrend_mode=Symbol(detrend_menu.selection[]))
        show_record!(rec)
        return rec
    end

    function run()
        try
            say("Computing…"); yield()
            rec = state.record[]
            rec === nothing && (rec = load_data(false))
            devs = [d for d in DEVIATIONS if dev_toggles[d].active[]]
            taumode = Dict(TAU_MODES)[tau_menu.selection[]]
            suite = analyze(rec, devs, taumode, ci_toggle.active[],
                            parse(Float64, conf_menu.selection[]), bias_toggle.active[])
            sp  = spectra(rec, parse_float(carrier_box, "carrier frequency"))
            dyn = dynamic_map(rec, Symbol(dyn_menu.selection[]),
                              1 / parse(Int, split(frac_menu.selection[], "/")[2]))
            state.suite[] = suite; state.spectra[] = sp; state.dynamic[] = dyn

            foreach(delete!, legends["Deviations"]); empty!(legends["Deviations"])
            sy = filter(r -> !(r.deviation_type in SIGMA_X), suite.results)
            sx = filter(r ->   r.deviation_type in SIGMA_X,  suite.results)
            for (ax, rs, pos) in ((ax_y, sy, :lb), (ax_x, sx, :lt))
                leg = draw_results!(ax, rs; band=band_toggle.active[], legend_position=pos)
                leg === nothing || push!(legends["Deviations"], leg)
            end
            sx_shown[] = !isempty(sx)
            rowsize!(dev_page, 2, sx_shown[] ? Auto(1) : Fixed(0))
            rowgap!(dev_page, 1, sx_shown[] ? 16 : 0)
            names = string.(keys(suite))
            table_menu.options = names
            table_menu.i_selected = 1
            fill_table!(table_grid, suite[1])

            draw_spectrum!(ax_sy, sp.sy; color=Makie.wong_colors()[1])
            draw_spectrum!(ax_sx, sp.sx; color=Makie.wong_colors()[2])
            sp.l === nothing ? empty!(ax_l) : draw_spectrum!(ax_l, sp.l; color=Makie.wong_colors()[3])
            draw_dynamic!(ax_d, dyn; cb=cb_d, cb_cell=dyn_page[2, 2])

            src_name = rec.source == "user" ? "synthetic record" : basename(rec.source)
            title_lbl.text = "$src_name: $(suite.n) samples, τ₀ = $(suite.tau0) s, " *
                             (ci_toggle.active[] ? "CI $(suite.confidence)" : "no CI")
            show_page(page_menu.selection[])
            say("Computed " * join(names, ", ") * " on $(suite.n) samples.")
        catch err
            say("Error: " * sprint(showerror, err); error=true)
        end
    end

    function reload(from_file)
        try
            load_data(from_file)
            say("Loaded $(nsamples(state.record[])) samples. Press Run analysis.")
        catch err
            say("Error: " * sprint(showerror, err); error=true)
        end
    end

    on(run_btn.clicks) do _; run(); end
    from_file = Ref(false)        # how the current record was produced
    # Keyboard: Tab / Shift+Tab move between the fields of the current page,
    # F5 runs, Page Up / Page Down switch pages, F1 opens Help.
    page_fields = Dict("Record" => [path_box, tau0_box, time_col, value_col, header, seed_box],
                       "Spectrum" => [carrier_box])
    on(events(fig).keyboardbutton; priority=100) do ev
        ev.action == Keyboard.press || return Consume(false)
        page = page_menu.selection[]
        if ev.key == Keyboard.tab
            fields = get(page_fields, page, Textbox[])
            isempty(fields) && return Consume(false)   # no custom order here; leave Tab to Makie
            i = findfirst(tb -> tb.focused[], fields)
            back = ispressed(fig, Keyboard.left_shift | Keyboard.right_shift)
            nxt = i === nothing ? (back ? length(fields) : 1) : mod1(i + (back ? -1 : 1), length(fields))
            foreach(Makie.defocus!, fields)
            Makie.focus!(fields[nxt])
            return Consume(true)
        elseif ev.key == Keyboard.f5
            run(); return Consume(true)
        elseif ev.key == Keyboard.f1
            page_menu.i_selected = findfirst(==("Help"), PAGES); return Consume(true)
        elseif ev.key in (Keyboard.page_down, Keyboard.page_up)
            step = ev.key == Keyboard.page_down ? 1 : -1
            page_menu.i_selected = mod1(something(page_menu.i_selected[], 1) + step, length(PAGES))
            return Consume(true)
        end
        return Consume(false)
    end
    on(load_btn.clicks) do _; from_file[] = true; reload(true); end
    on(synth_btn.clicks) do _; from_file[] = false; reload(false); end
    for tb in (tau0_box, time_col, value_col, header, seed_box)
        on(tb.stored_string) do _; reload(from_file[]); state.record[] === nothing || run(); end
    end
    on(carrier_box.stored_string) do _; run(); end
    on(table_menu.selection) do name
        s = state.suite[]
        (s === nothing || name === nothing) && return
        fill_table!(table_grid, s[Symbol(name)])
    end
    on(open_btn.clicks) do _
        p = pick_file(state.last_dir[]; filterlist="dat,txt,csv,tsv")
        isempty(p) && return
        state.last_dir[] = dirname(p)
        path_box.displayed_string[] = p; path_box.stored_string[] = p
        from_file[] = true
        reload(true)
    end
    function export_to(default_name, filterlist, writer)
        p = save_file(joinpath(state.last_dir[], default_name); filterlist=filterlist)
        isempty(p) && return
        try
            writer(p); state.last_dir[] = dirname(p); say("Wrote $p.")
        catch err
            say("Error: " * sprint(showerror, err); error=true)
        end
    end
    write_figure(p, page) = Makie.save(p, export_figure(page, state, (; title=title_lbl.text[], band=band_toggle.active[]));
                                       backend=CairoMakie, px_per_unit=2)
    write_table(p) = SigmaTau.save(p, state.suite[])
    write_record(p) = SigmaTau.save(p, state.record[])
    on(export_fig_btn.clicks) do _
        page = page_menu.selection[]
        page == "Help" && (say("Nothing to export on the Help page."); return)
        (page == "Record" ? state.record[] : state.suite[]) === nothing &&
            (say("Run an analysis first."; error=true); return)
        export_to("sigmatau_$(lowercase(replace(page, " " => "_"))).png", "png,pdf,svg", p -> write_figure(p, page))
    end
    on(export_tsv_btn.clicks) do _
        state.suite[] === nothing && (say("Run an analysis first."; error=true); return)
        export_to("sigmatau_suite.tsv", "tsv,txt", write_table)
    end
    on(export_rec_btn.clicks) do _
        state.record[] === nothing && (say("Load or synthesize a record first."; error=true); return)
        export_to("sigmatau_record.txt", "txt,csv,dat", write_record)
    end

    DataInspector(fig)
    show_page("Record")
    return (; fig, run, reload, state, show_page, write_figure, write_table, write_record,
            status=() -> status.text[],
            controls=(; dev_toggles, band_toggle, page_menu, path_box, kind, tau0_box, time_col, value_col, header))
end

# ── Entry point ───────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    # Options: --smoke (headless self-test), --fontsize=N (default 12; 14 or 16
    # for projection), --size=WxH (window size). SIGMATAU_FONTSIZE also works.
    option(name, default) = (i = findfirst(startswith("--$name="), ARGS); i === nothing ? default : split(ARGS[i], "=")[2])
    smoke = "--smoke" in ARGS
    fontsize = parse(Int, option("fontsize", get(ENV, "SIGMATAU_FONTSIZE", "12")))
    # `win_size`, not `size`: a top-level `size` would shadow Base.size.
    win_size = Tuple(parse.(Int, split(option("size", "1500x940"), "x")))
    app = build(; size=win_size, fontsize=fontsize)
    app.reload(false)       # open with a synthetic record on screen
    app.run()
    if smoke
        dir = get(ENV, "SIGMATAU_SMOKE_DIR", tempdir())
        c = app.controls
        # Load path: the Stable32 fixture through the Record-page fields.
        fixture = joinpath(pkgdir(SigmaTau), "test", "fixtures", "validation", "stable32gen.DAT")
        c.path_box.displayed_string[] = fixture; c.header.displayed_string[] = "10"   # typed, no Enter
        app.reload(true)
        n_fixture = nsamples(app.state.record[])
        n_fixture == 8192 || error("smoke test: Stable32 fixture read $n_fixture samples")
        c.dev_toggles[:tdev].active[] = true
        c.dev_toggles[:mtie].active[] = true
        app.run()
        s = app.state.suite[]
        (s === nothing || s.n != 8192) && error("smoke test: no suite from the loaded file")
        # Export paths and the record round trip.
        app.write_table(joinpath(dir, "suite.tsv"))
        length(load_suite(joinpath(dir, "suite.tsv")).results) == length(s.results) || error("smoke test: suite round trip")
        app.write_record(joinpath(dir, "record.txt"))
        c.path_box.displayed_string[] = joinpath(dir, "record.txt"); c.header.displayed_string[] = "0"
        app.reload(true)
        nsamples(app.state.record[]) == 8192 || error("smoke test: record round trip")
        # Error path: a bad field reports in the status line and keeps the record.
        c.tau0_box.displayed_string[] = "abc"
        app.reload(false)
        startswith(app.status(), "Error:") || error("smoke test: bad τ₀ not reported")
        c.tau0_box.displayed_string[] = "1.0"
        app.reload(false); app.run()
        s = app.state.suite[]
        for page in PAGES
            c.page_menu.i_selected = findfirst(==(page), PAGES)
            Makie.save(joinpath(dir, "explorer_" * replace(lowercase(page), " " => "_") * ".png"),
                       app.fig; px_per_unit=1)
            page == "Help" || app.write_figure(joinpath(dir, "export_" * replace(lowercase(page), " " => "_") * ".pdf"), page)
        end
        println("smoke test OK: ", join(string.(keys(s)), ", "), " → ", dir)
    else
        screen = display(app.fig; title="SigmaTau explorer " * VERSION_STRING)
        wait(screen)
    end
end

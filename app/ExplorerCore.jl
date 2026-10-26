# Shared core of the SigmaTau explorers: constants, analysis helpers and the
# Makie drawing functions used by both the GLMakie window
# (StabilityExplorer.jl) and the browser app (StabilityWeb.jl). Backend-free:
# the including file activates GLMakie or WGLMakie.

using SigmaTau
using Makie
using CairoMakie
using Random
using Printf
using Statistics

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


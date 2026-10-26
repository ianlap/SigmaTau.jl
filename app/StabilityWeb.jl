# SigmaTau stability explorer, browser edition. Same analysis and plots as
# the GLMakie window (shared through ExplorerCore.jl), with native HTML
# controls: every input has a label, the table is a real <table>, status
# messages are announced through a live region, and Tab moves through the
# controls in reading order. Plots are WGLMakie canvases.
#
#   julia --project=app app/StabilityWeb.jl             # serve on http://127.0.0.1:9384
#   julia --project=app app/StabilityWeb.jl --port=8000
#   julia --project=app app/StabilityWeb.jl --smoke     # start, wait for SIGMATAU_SMOKE_DONE, exit
#
# This is a local application: it reads files from the machine it runs on
# (type a path, or pick a file in the browser, which uploads it to the
# server). Exports download through the browser.

using Bonito
using WGLMakie
using Base64
include("ExplorerCore.jl")

# Bonito and Makie both export these widget names; the browser edition
# wants Bonito's.
const Checkbox  = Bonito.Checkbox
const Button    = Bonito.Button
const TextField = Bonito.TextField
const Dropdown  = Bonito.Dropdown

WGLMakie.activate!()

const CSS_TEXT = """
:root { color-scheme: light; }
body { font-family: system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif; font-size: 15px;
       color: #1d1d1f; background: #f7f7f8; margin: 0; }
header { background: #2f4e8f; color: white; padding: 10px 24px; display: flex; gap: 24px; align-items: baseline; }
header h1 { font-size: 18px; margin: 0; } header nav a { color: #dfe7f7; margin-right: 16px; }
main { display: grid; grid-template-columns: 360px 1fr; gap: 24px; padding: 16px 24px 72px; align-items: start; }
@media (max-width: 1000px) { main { grid-template-columns: 1fr; } }
fieldset { border: 1px solid #c9ccd3; border-radius: 6px; margin: 0 0 14px; padding: 8px 12px 10px; }
legend { font-weight: 600; padding: 0 6px; }
label { display: inline-flex; align-items: center; gap: 6px; margin: 3px 12px 3px 0; }
.grid4 { display: grid; grid-template-columns: repeat(4, 1fr); } .grid4 label { margin-right: 2px; }
.row { display: flex; flex-wrap: wrap; gap: 10px 18px; align-items: center; margin: 4px 0; }
input[type=text], input[type=number], select { font: inherit; padding: 4px 6px; border: 1px solid #9aa0ab; border-radius: 4px; background: white; }
input[type=text] { min-width: 6em; } .wide { width: 100%; box-sizing: border-box; }
button { font: inherit; padding: 6px 14px; border: 1px solid #7f8794; border-radius: 5px; background: #eef0f3; cursor: pointer; }
button.primary { background: #2f4e8f; color: white; border-color: #2f4e8f; font-weight: 600; }
button:focus-visible, input:focus-visible, select:focus-visible, a:focus-visible { outline: 3px solid #f0a500; outline-offset: 2px; }
section { background: white; border: 1px solid #e1e3e8; border-radius: 8px; padding: 12px 18px 18px; margin-bottom: 20px; }
section h2 { font-size: 17px; margin: 0 0 10px; }
table { border-collapse: collapse; font-variant-numeric: tabular-nums; margin-top: 8px; }
caption { text-align: left; font-weight: 600; padding: 4px 0; }
th, td { padding: 3px 10px; text-align: right; } th { background: #e8eaef; }
tbody tr:nth-child(even) { background: #f2f3f6; } th:last-child, td:last-child { text-align: left; }
#status { padding: 8px 24px; border-top: 1px solid #e1e3e8; background: white; position: sticky; bottom: 0; z-index: 10; }
#status.error { color: #b3261e; } .hint { color: #595d66; font-size: 14px; }
.summary p { margin: 2px 0; } .plot { margin-top: 8px; }
"""

# ── Widgets with labels ───────────────────────────────────────────────────

# Checkboxes read "[x] name"; everything else reads "name: [control]".
labelled(text, w::Checkbox; kw...) = DOM.label(w, text; kw...)
labelled(text, w; kw...) = DOM.label(text, w; kw...)
hint(text) = DOM.p(text; class="hint")

# Browser-side download of bytes produced on the server.
function download!(session, name, mime, bytes::Vector{UInt8})
    data = base64encode(bytes)
    evaljs(session, js"""
        const a = document.createElement('a');
        a.href = 'data:' + $(mime) + ';base64,' + $(data);
        a.download = $(name);
        document.body.appendChild(a); a.click(); a.remove();
    """)
end

function figure_bytes(fig, ext)
    path = tempname() * "." * ext
    Makie.save(path, fig; backend=CairoMakie, px_per_unit=2)
    bytes = read(path); rm(path; force=true)
    return bytes
end

# ── The app ───────────────────────────────────────────────────────────────

function explorer_app()
    App(; title="SigmaTau explorer " * VERSION_STRING) do session
        # State
        state = (record=Ref{Any}(nothing), suite=Ref{Any}(nothing),
                 spectra=Ref{Any}(nothing), dynamic=Ref{Any}(nothing))
        status = Observable("Ready.")
        status_class = Observable("")
        say(msg; error=false) = (status_class[] = error ? "error" : ""; status[] = msg)

        # Deviation toggles
        dev_boxes = Dict(d => Checkbox(d in DEFAULT_ON; style=Styles()) for d in DEVIATIONS)
        dev_grid = DOM.div([labelled(string(d), dev_boxes[d]) for d in DEVIATIONS]...; class="grid4")

        # Analysis options
        tau_menu  = Dropdown(first.(TAU_MODES); style=Styles())
        ci_box    = Checkbox(true; style=Styles())
        conf_menu = Dropdown(["0.683", "0.90", "0.95", "0.99"]; style=Styles())
        bias_box  = Checkbox(true; style=Styles())
        band_box  = Checkbox(false; style=Styles())
        detrend_menu = Dropdown(["none", "mean", "linear", "quadratic"]; style=Styles())
        outliers_box = Checkbox(false; style=Styles())
        gaps_box     = Checkbox(false; style=Styles())
        run_btn = Button("Run analysis"; style=Styles(), class="primary")

        # Record: file
        path_tf   = TextField(""; style=Styles(), class="wide", placeholder="path on the server, or pick a file below")
        upload    = Observable{Any}(nothing)
        file_in   = DOM.input(; type="file", id="file-pick", accept=".dat,.txt,.csv,.tsv",
                              onchange=js"""e => { const f = e.target.files[0]; if (!f) return;
                                             f.text().then(t => $(upload).notify([f.name, t])); }""")
        kind_menu = Dropdown(["phase", "frequency"]; style=Styles())
        tau0_tf   = TextField("1.0"; style=Styles())
        tcol_tf   = TextField("0"; style=Styles())
        vcol_tf   = TextField("1"; style=Styles())
        hdr_tf    = TextField("0"; style=Styles())
        load_btn  = Button("Load file"; style=Styles())
        # Record: synthesis
        n_menu    = Dropdown(["4096", "8192", "16384", "32768", "86400"]; index=3, style=Styles())
        seed_tf   = TextField("42"; style=Styles())
        noise_box = Dict(α => Checkbox(α in (2, 0, -2); style=Styles()) for (α, _, _) in NOISES)
        noise_exp = Dict(α => Dropdown(["1e$k" for k in -18:-9]; index=start + 19, style=Styles()) for (α, _, start) in NOISES)
        synth_btn = Button("Synthesize"; style=Styles())
        # Spectrum / dynamic options
        carrier_tf = TextField("10e6"; style=Styles())
        dyn_menu   = Dropdown(["dadev", "dhdev"]; style=Styles())
        frac_menu  = Dropdown(["1/4", "1/8", "1/16", "1/32"]; index=2, style=Styles())
        # Exports
        fig_fmt    = Dropdown(["png", "pdf", "svg"]; style=Styles())
        exp_fig_btn = Button("Download figure"; style=Styles())
        exp_tsv_btn = Button("Download table"; style=Styles())
        exp_rec_btn = Button("Download record"; style=Styles())
        table_menu = Dropdown(["adev"]; style=Styles())

        # Dynamic content
        title_obs   = Observable("No record loaded.")
        summary_obs = Observable{Any}(DOM.div())
        record_plot = Observable{Any}(DOM.div())
        dev_plot    = Observable{Any}(DOM.div("Run the analysis to see the deviations."))
        table_obs   = Observable{Any}(DOM.div())
        spec_plot   = Observable{Any}(DOM.div())
        dyn_plot    = Observable{Any}(DOM.div())
        last_title  = Ref("")

        text(tf) = String(something(tf.value[], ""))
        parse_int(tf, what) = (v = tryparse(Int, text(tf)); v === nothing ? error("$what must be an integer, got \"$(text(tf))\".") : v)
        parse_float(tf, what) = (v = tryparse(Float64, text(tf)); v === nothing ? error("$what must be a number, got \"$(text(tf))\".") : v)

        function show_record!(rec)
            state.record[] = rec
            fig = Figure(size=(880, 300), fontsize=13)
            ax = Axis(fig[1, 1]; xlabel="Time t (s)", xtickformat=vs -> [@sprintf("%g", v) for v in vs])
            draw_record!(ax, rec)
            record_plot[] = DOM.div(fig; class="plot")
            summary_obs[] = DOM.div([DOM.p(l) for l in record_summary(rec)]...; class="summary")
            name = rec.source == "user" ? "synthetic record" : basename(rec.source)
            title_obs[] = "$name: $(nsamples(rec)) samples, τ₀ = $(rec.tau0) s"
        end

        function load_data(from_file)
            rec = if from_file
                load_record(text(path_tf), kind_menu.value[], text(tau0_tf),
                            parse_int(tcol_tf, "time column"), parse_int(vcol_tf, "value column"),
                            parse_int(hdr_tf, "header rows"))
            else
                levels = [(α, parse(Int, replace(noise_exp[α].value[], "1e" => ""))) for (α, _, _) in NOISES if noise_box[α].value[]]
                synthesize(parse(Int, n_menu.value[]), parse_float(tau0_tf, "τ₀"), parse_int(seed_tf, "seed"), levels)
            end
            rec = preprocess(rec; outliers=outliers_box.value[], gaps=gaps_box.value[],
                             detrend_mode=Symbol(detrend_menu.value[]))
            show_record!(rec)
            return rec
        end

        function table_dom(r::StabilityResult)
            has_ci = !isempty(r.edf)
            g(v) = isfinite(v) ? @sprintf("%.5g", v) : "NaN"
            head = has_ci ? ["τ (s)", "σ", "n_eff", "σ_min", "σ_max", "EDF", "noise"] : ["τ (s)", "σ", "n_eff"]
            rows = [DOM.tr([DOM.td(c) for c in (has_ci ?
                        [g(r.tau[i]), g(r.dev[i]), string(r.neff[i]), g(r.ci[i].lo), g(r.ci[i].hi), g(r.edf[i]), string(r.noise_type[i])] :
                        [g(r.tau[i]), g(r.dev[i]), string(r.neff[i])])]...) for i in eachindex(r.tau)]
            DOM.table(DOM.caption("Per-τ values of " * uppercase(string(r.deviation_type))),
                      DOM.thead(DOM.tr([DOM.th(h; scope="col") for h in head]...)),
                      DOM.tbody(rows...))
        end

        function run()
            try
                say("Computing…")
                rec = state.record[]
                rec === nothing && (rec = load_data(false))
                devs = [d for d in DEVIATIONS if dev_boxes[d].value[]]
                taumode = Dict(TAU_MODES)[tau_menu.value[]]
                suite = analyze(rec, devs, taumode, ci_box.value[], parse(Float64, conf_menu.value[]), bias_box.value[])
                sp  = spectra(rec, parse_float(carrier_tf, "carrier frequency"))
                dyn = dynamic_map(rec, Symbol(dyn_menu.value[]), 1 / parse(Int, split(frac_menu.value[], "/")[2]))
                state.suite[] = suite; state.spectra[] = sp; state.dynamic[] = dyn
                name = rec.source == "user" ? "synthetic record" : basename(rec.source)
                last_title[] = "$name: $(suite.n) samples, τ₀ = $(suite.tau0) s, " * (ci_box.value[] ? "CI $(suite.confidence)" : "no CI")
                title_obs[] = last_title[]
                opts = (; title=last_title[], band=band_box.value[])
                dev_plot[]  = DOM.div(export_figure("Deviations", state, opts); class="plot")
                spec_plot[] = DOM.div(export_figure("Spectrum", state, opts); class="plot")
                dyn_plot[]  = DOM.div(export_figure("Dynamic deviation", state, opts); class="plot")
                names = string.(keys(suite))
                table_menu.options[] = names
                table_menu.option_index[] = 1
                table_obs[] = table_dom(suite[1])
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

        on(run_btn.value) do _; run(); end
        on(load_btn.value) do _; reload(true); end
        on(synth_btn.value) do _; reload(false); end
        on(upload) do v
            v === nothing && return
            name, content = String(v[1]), String(v[2])
            path = joinpath(mktempdir(), name)
            write(path, content)
            path_tf.value[] = path
            reload(true)
        end
        on(table_menu.value) do name
            s = state.suite[]
            (s === nothing || name === nothing) && return
            table_obs[] = table_dom(s[Symbol(name)])
        end
        on(exp_fig_btn.value) do _
            state.suite[] === nothing && (say("Run an analysis first."; error=true); return)
            ext = fig_fmt.value[]
            fig = export_figure("Deviations", state, (; title=last_title[], band=band_box.value[]))
            mime = Dict("png" => "image/png", "pdf" => "application/pdf", "svg" => "image/svg+xml")[ext]
            download!(session, "sigmatau_deviations.$ext", mime, figure_bytes(fig, ext))
            say("Figure sent to the browser as sigmatau_deviations.$ext.")
        end
        on(exp_tsv_btn.value) do _
            state.suite[] === nothing && (say("Run an analysis first."; error=true); return)
            path = tempname() * ".tsv"; SigmaTau.save(path, state.suite[])
            download!(session, "sigmatau_suite.tsv", "text/tab-separated-values", read(path)); rm(path)
            say("Table sent to the browser as sigmatau_suite.tsv.")
        end
        on(exp_rec_btn.value) do _
            state.record[] === nothing && (say("Load or synthesize a record first."; error=true); return)
            path = tempname() * ".txt"; SigmaTau.save(path, state.record[])
            download!(session, "sigmatau_record.txt", "text/plain", read(path)); rm(path)
            say("Record sent to the browser as sigmatau_record.txt.")
        end

        # Open with a synthetic record analysed, like the window.
        load_data(false); run()

        controls = DOM.div(
            DOM.fieldset(DOM.legend("Deviations"), dev_grid,
                         hint("σ_y: adev … theoh, pdev. σ_x (seconds): tdev, htdev, ttotdev, mtie, tierms.")),
            DOM.fieldset(DOM.legend("Analysis"),
                         DOM.div(labelled("τ grid", tau_menu); class="row"),
                         DOM.div(labelled("Confidence intervals", ci_box), labelled("Level", conf_menu); class="row"),
                         DOM.div(labelled("Bias-correct totals and Thêo1", bias_box), labelled("CI as band", band_box); class="row")),
            DOM.fieldset(DOM.legend("Preprocessing"),
                         DOM.div(labelled("Detrend", detrend_menu); class="row"),
                         DOM.div(labelled("Remove 5σ outliers", outliers_box), labelled("Fill NaN gaps", gaps_box); class="row")),
            DOM.div(run_btn; class="row"),
            DOM.fieldset(DOM.legend("Export"),
                         DOM.div(labelled("Figure format", fig_fmt), exp_fig_btn; class="row"),
                         DOM.div(exp_tsv_btn, exp_rec_btn; class="row")),
        )

        record_section = DOM.section(DOM.h2("Record"; id="record"),
            DOM.div(
                DOM.fieldset(DOM.legend("Load a file"),
                    DOM.div(labelled("File path", path_tf); class="row"),
                    DOM.div(labelled("Or pick a file", file_in); class="row"),
                    DOM.div(labelled("Record type", kind_menu), labelled("τ₀ (s)", tau0_tf); class="row"),
                    DOM.div(labelled("Time column", tcol_tf), labelled("Value column", vcol_tf), labelled("Header rows", hdr_tf); class="row"),
                    hint("Stable32 .DAT: header rows 10, time column 0, value column 1."),
                    DOM.div(load_btn; class="row")),
                DOM.fieldset(DOM.legend("Synthesize a record (level = σ_y at τ₀)"),
                    DOM.div(labelled("Samples", n_menu), labelled("Seed", seed_tf); class="row"),
                    [DOM.div(labelled(name, noise_box[α]), labelled("level", noise_exp[α]); class="row") for (α, name, _) in NOISES]...,
                    DOM.div(synth_btn; class="row"));
                style="display:grid;grid-template-columns:1fr 1fr;gap:16px"),
            record_plot, summary_obs)

        dev_section = DOM.section(DOM.h2("Deviations"; id="deviations"), dev_plot,
                                  DOM.div(labelled("Per-τ table for", table_menu); class="row"), table_obs)
        spec_section = DOM.section(DOM.h2("Spectrum"; id="spectrum"),
                                   DOM.div(labelled("Carrier frequency for ℒ(f) (Hz)", carrier_tf), hint("takes effect at the next run"); class="row"),
                                   spec_plot)
        dyn_section = DOM.section(DOM.h2("Dynamic deviation"; id="dynamic"),
                                  DOM.div(labelled("Estimator", dyn_menu), labelled("Window as a fraction of the record", frac_menu),
                                          hint("takes effect at the next run; drag the plot to rotate"); class="row"),
                                  dyn_plot)
        help_section = DOM.section(DOM.h2("Help"; id="help"),
                                   [DOM.div(DOM.h3(t), [DOM.p(l) for l in split(b, "\n")]...) for (t, b) in HELP_SECTIONS]...,
                                   DOM.h3("Deviations"),
                                   DOM.dl([DOM.div(DOM.dt(n), DOM.dd(d)) for (n, d) in GLOSSARY]...))

        DOM.div(
            DOM.style(CSS_TEXT),
            DOM.header(DOM.h1("SigmaTau explorer " * VERSION_STRING),
                       DOM.nav(DOM.a("Record"; href="#record"), DOM.a("Deviations"; href="#deviations"),
                               DOM.a("Spectrum"; href="#spectrum"), DOM.a("Dynamic deviation"; href="#dynamic"),
                               DOM.a("Help"; href="#help"); var"aria-label"="Pages")),
            DOM.main(DOM.div(controls), DOM.div(DOM.p(title_obs; class="hint"), record_section, dev_section, spec_section, dyn_section, help_section)),
            DOM.div(status; id="status", class=status_class, role="status", var"aria-live"="polite"),
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    option(name, default) = (i = findfirst(startswith("--$name="), ARGS); i === nothing ? default : split(ARGS[i], "=")[2])
    port = parse(Int, option("port", "9384"))
    server = Bonito.Server(explorer_app(), "127.0.0.1", port)
    println("SigmaTau explorer serving on http://127.0.0.1:$port")
    if "--smoke" in ARGS
        done = get(ENV, "SIGMATAU_SMOKE_DONE", "")
        while isempty(done) || !isfile(done); sleep(1); end
        close(server)
    else
        wait(server.server.serve_task)
    end
end

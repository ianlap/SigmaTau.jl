# sigmatau report: one self-contained HTML file (styles and SVG plots inline)
# with the record, the settings, σ(τ) plots and per-τ tables. A browser
# prints it to PDF. Written by hand like the JSON: no templating package.

const PALETTE = ["#0b6fb0", "#e69f00", "#009e73", "#cc79a7", "#9e1b32", "#56b4e9", "#d55e00", "#555555"]

function html_escape(io::IO, s::AbstractString)
    for c in s
        c == '&' ? print(io, "&amp;") : c == '<' ? print(io, "&lt;") : c == '>' ? print(io, "&gt;") :
        c == '"' ? print(io, "&quot;") : print(io, c)
    end
end

# A POSIX shell word for `a`: as is when it is plainly safe, else single-quoted.
function shell_quote(a::String)
    !isempty(a) && all(c -> isletter(c) || isdigit(c) || c in "-_.,:/=+@%", a) && return a
    return "'" * replace(a, "'" => "'\\''") * "'"
end

# x to `d` significant figures in e-notation (1.234e-12), "–" if not finite.
function sci(io::IO, x::Float64, d::Int=4)
    isfinite(x) || return print(io, "–")
    x == 0 && return print(io, "0")
    e = floor(Int, log10(abs(x)))
    m = round(x / exp10(e); digits=d - 1)
    abs(m) >= 10 && (m /= 10; e += 1)
    print(io, round(m; digits=d - 1)); print(io, 'e'); print(io, e)
end

# Plain to `d` significant figures for moderate values, else e-notation.
function plain(io::IO, x::Float64, d::Int=5)
    (isfinite(x) && x != 0 && (abs(x) >= 1e5 || abs(x) < 1e-3)) && return sci(io, x, d)
    isfinite(x) || return print(io, "–")
    v = round(x; sigdigits=d)
    v == round(v) && abs(v) < 1e15 ? print(io, round(Int, v)) : print(io, v)
end

function duration_text(io::IO, s::Float64)
    s < 120 ? (plain(io, s, 3); print(io, " s")) :
    s < 7200 ? (plain(io, s / 60, 3); print(io, " min")) :
    s < 3 * 86400 ? (plain(io, s / 3600, 3); print(io, " h")) :
    (plain(io, s / 86400, 3); print(io, " d"))
end

# ── SVG log-log plot ──────────────────────────────────────────────────────

const PW, PH, ML, MR, MT, MB = 760, 420, 70, 16, 14, 46

# Decade lo at pixel a, decade hi at pixel b.
struct LogScale
    lo::Int
    hi::Int
    a::Float64
    b::Float64
end
(s::LogScale)(v::Float64) = s.a + (log10(v) - s.lo) / (s.hi - s.lo) * (s.b - s.a)

function svg_plot(io::IO, results::Vector{StabilityResult}, labels::Vector{String}, ylabel::String)
    xs = Float64[]; ys = Float64[]
    for r in results, i in eachindex(r.tau)
        isfinite(r.dev[i]) && r.dev[i] > 0 || continue
        push!(xs, r.tau[i]); push!(ys, r.dev[i])
        if !isempty(r.ci)
            isfinite(r.ci[i].lo) && r.ci[i].lo > 0 && push!(ys, r.ci[i].lo)
            isfinite(r.ci[i].hi) && push!(ys, r.ci[i].hi)
        end
    end
    isempty(xs) && return
    x0 = floor(Int, log10(minimum(xs))); x1 = max(ceil(Int, log10(maximum(xs))), x0 + 1)
    y0 = floor(Int, log10(minimum(ys))); y1 = max(ceil(Int, log10(maximum(ys))), y0 + 1)
    px = LogScale(x0, x1, Float64(ML), Float64(PW - MR))
    py = LogScale(y1, y0, Float64(MT), Float64(PH - MB))   # decades run upwards
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 ", PW, ' ', PH, "\" role=\"img\" font-family=\"sans-serif\" font-size=\"12\">")
    print(io, "<rect x=\"", ML, "\" y=\"", MT, "\" width=\"", PW - ML - MR, "\" height=\"", PH - MT - MB, "\" fill=\"#fff\" stroke=\"#999\"/>")
    for k in x0:x1
        x = px(exp10(k))
        print(io, "<line x1=\"", x, "\" x2=\"", x, "\" y1=\"", MT, "\" y2=\"", PH - MB, "\" stroke=\"#e4e4e4\"/>")
        print(io, "<text x=\"", x, "\" y=\"", PH - MB + 16, "\" text-anchor=\"middle\">10<tspan dy=\"-5\" font-size=\"9\">", k, "</tspan></text>")
    end
    for k in y0:y1
        y = py(exp10(k))
        print(io, "<line x1=\"", ML, "\" x2=\"", PW - MR, "\" y1=\"", y, "\" y2=\"", y, "\" stroke=\"#e4e4e4\"/>")
        print(io, "<text x=\"", ML - 6, "\" y=\"", y + 4, "\" text-anchor=\"end\">10<tspan dy=\"-5\" font-size=\"9\">", k, "</tspan></text>")
    end
    print(io, "<text x=\"", (ML + PW - MR) / 2, "\" y=\"", PH - 8, "\" text-anchor=\"middle\">Averaging time τ (s)</text>")
    print(io, "<text transform=\"translate(16 ", (MT + PH - MB) / 2, ") rotate(-90)\" text-anchor=\"middle\">")
    html_escape(io, ylabel); print(io, "</text>")
    legend = IOBuffer()
    for (j, r) in enumerate(results)
        col = PALETTE[mod1(j, length(PALETTE))]
        ok = [i for i in eachindex(r.tau) if isfinite(r.dev[i]) && r.dev[i] > 0]
        isempty(ok) && continue
        if !isempty(r.ci)
            band = [i for i in ok if isfinite(r.ci[i].lo) && r.ci[i].lo > 0 && isfinite(r.ci[i].hi)]
            if length(band) > 1
                print(io, "<polygon fill=\"", col, "\" fill-opacity=\"0.15\" stroke=\"none\" points=\"")
                for i in band; print(io, px(r.tau[i]), ',', py(r.ci[i].hi), ' '); end
                for i in reverse(band); print(io, px(r.tau[i]), ',', py(r.ci[i].lo), ' '); end
                print(io, "\"/>")
            end
        end
        print(io, "<polyline fill=\"none\" stroke=\"", col, "\" stroke-width=\"2\" points=\"")
        for i in ok; print(io, px(r.tau[i]), ',', py(r.dev[i]), ' '); end
        print(io, "\"/>")
        for i in ok
            print(io, "<circle r=\"2.5\" fill=\"", col, "\" cx=\"", px(r.tau[i]), "\" cy=\"", py(r.dev[i]), "\"/>")
        end
        ly = MT + 8 + 18j
        print(legend, "<line x1=\"", PW - MR - 110, "\" x2=\"", PW - MR - 88, "\" y1=\"", ly, "\" y2=\"", ly, "\" stroke=\"", col, "\" stroke-width=\"2\"/>")
        print(legend, "<text x=\"", PW - MR - 82, "\" y=\"", ly + 4, "\">"); html_escape(legend, labels[j]); print(legend, "</text>")
    end
    print(io, "<rect x=\"", PW - MR - 118, "\" y=\"", MT + 10, "\" width=\"110\" height=\"", 18 * length(results) + 6,
          "\" fill=\"#fff\" fill-opacity=\"0.9\" stroke=\"#ccc\"/>")
    write(io, take!(legend))
    print(io, "</svg>")
end

# ── Tables ────────────────────────────────────────────────────────────────

function html_result_table(io::IO, r::StabilityResult, lab::String)
    has_ci, has_edf, has_nt = !isempty(r.ci), !isempty(r.edf), !isempty(r.noise_type)
    print(io, "<table class=\"results\"><caption>"); html_escape(io, lab)
    print(io, "</caption><thead><tr><th>τ (s)</th><th>σ</th>")
    has_ci && print(io, "<th>CI low</th><th>CI high</th>")
    has_edf && print(io, "<th>EDF</th>")
    has_nt && print(io, "<th>Noise</th>")
    print(io, "<th>Terms</th></tr></thead><tbody>")
    for i in eachindex(r.tau)
        print(io, "<tr><td class=\"n\">"); plain(io, r.tau[i])
        print(io, "</td><td class=\"n\">"); sci(io, r.dev[i], 4)
        if has_ci
            print(io, "</td><td class=\"n\">"); sci(io, r.ci[i].lo, 3)
            print(io, "</td><td class=\"n\">"); sci(io, r.ci[i].hi, 3)
        end
        has_edf && (print(io, "</td><td class=\"n\">"); plain(io, r.edf[i], 4))
        has_nt && (print(io, "</td><td>"); print(io, noise_short(r.noise_type[i])))
        print(io, "</td><td class=\"n\">"); print(io, r.neff[i]); print(io, "</td></tr>")
    end
    print(io, "</tbody></table>")
end

noise_short(s::Symbol) = s === :WHPM ? "WPM" : s === :FLPM ? "FPM" : s === :WHFM ? "WFM" :
                         s === :FLFM ? "FFM" : s === :RWFM ? "RWFM" : String(s)

const REPORT_STYLE = """
body{font:13px/1.5 system-ui,sans-serif;color:#1d1d1f;background:#fff;margin:32px auto;max-width:900px;padding:0 16px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 8px;border-bottom:1px solid #ddd;padding-bottom:4px}
.meta{color:#666;margin:0 0 16px}table{border-collapse:collapse;margin:8px 0 16px;font-size:12px}
th,td{text-align:left;padding:3px 8px;border-bottom:1px solid #eee;vertical-align:top}thead th{border-bottom:1px solid #bbb}
td.n{text-align:right;font-family:ui-monospace,monospace;white-space:nowrap}caption{text-align:left;font-weight:600;padding:4px 0}
dl{display:grid;grid-template-columns:max-content 1fr;gap:2px 16px;margin:0}dt{color:#666}dd{margin:0}
pre{font:11.5px ui-monospace,monospace;background:#f6f6f6;padding:8px 10px;white-space:pre-wrap;word-break:break-all}
svg{width:100%;height:auto}.results{break-inside:avoid}
@media print{body{margin:0;max-width:none}h2{break-after:avoid}svg{break-inside:avoid}}
"""

function cmd_report(io::IO, o::Options, data::PhaseData, args::Vector{String})
    results = analysis_results(o, data)
    isempty(results) && throw(CLIError(EXIT_ANALYSIS, "the record is too short for any τ"))
    s = summarize(data)
    src = source_name(o, data)
    print(io, "<!doctype html>\n<html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">")
    print(io, "<title>Stability report</title><style>", REPORT_STYLE, "</style></head><body>\n")
    print(io, "<h1>Frequency stability report</h1><p class=\"meta\">"); html_escape(io, src)
    print(io, " · SigmaTau ", SIGMATAU_VERSION, "</p>\n<h2>Record</h2><dl>")
    print(io, "<dt>Source</dt><dd>"); html_escape(io, src); print(io, "</dd>")
    print(io, "<dt>Type</dt><dd>", quantity(o), "</dd>")
    print(io, "<dt>Samples</dt><dd>", s.n, "</dd>")
    print(io, "<dt>τ<sub>0</sub></dt><dd>"); plain(io, Float64(data.tau0), 6); print(io, " s</dd>")
    print(io, "<dt>Span</dt><dd>"); duration_text(io, Float64(s.span)); print(io, "</dd>")
    print(io, "<dt>Gaps</dt><dd>", s.gaps, "</dd><dt>Outliers (5σ)</dt><dd>", length(s.outliers), "</dd>")
    print(io, "<dt>Drift / day</dt><dd>"); sci(io, Float64(s.drift_per_day), 3); print(io, "</dd>")
    print(io, "<dt>Preprocessing</dt><dd>detrend ", o.detrend, ", remove drift ", o.remove_drift,
          o.fillgaps ? ", gaps filled" : "", o.remove_outliers ? ", outliers removed" : "",
          o.reference_col > 0 ? ", minus column " * string(o.reference_col) : "")
    if !(isnan(o.range_start) && isnan(o.range_end))
        print(io, ", seconds ")
        isnan(o.range_start) ? print(io, "0") : plain(io, o.range_start, 6)
        print(io, " to ")
        isnan(o.range_end) ? print(io, "end") : plain(io, o.range_end, 6)
    end
    print(io, "</dd></dl>\n")
    print(io, "<h2>Analysis</h2><dl><dt>Estimators</dt><dd>")
    for (k, r) in enumerate(results)
        k > 1 && print(io, ", ")
        print(io, registry_entry(r.deviation_type).label)
    end
    print(io, "</dd><dt>τ</dt><dd>", isempty(o.m) ? o.grid * " grid" : "m = " * join(string.(o.m), ", "))
    isnan(o.tau_min) || (print(io, ", from "); plain(io, o.tau_min, 6); print(io, " s"))
    isnan(o.tau_max) || (print(io, ", up to "); plain(io, o.tau_max, 6); print(io, " s"))
    print(io, "</dd>")
    print(io, "<dt>Confidence</dt><dd>")
    o.ci ? (plain(io, 100 * o.confidence, 3); print(io, " % (Greenhall EDF, χ² intervals)")) : print(io, "off")
    print(io, "</dd><dt>Bias correction</dt><dd>", o.bias ? "on" : "off", "</dd></dl>\n<h2>Stability</h2>\n")
    for (axis, ylabel) in ((:sigma_y, "σy(τ)"), (:sigma_x, "σx(τ) (s)"))
        sel = StabilityResult[r for r in results if registry_entry(r.deviation_type).axis === axis]
        isempty(sel) && continue
        svg_plot(io, sel, String[registry_entry(r.deviation_type).label for r in sel], ylabel)
        print(io, '\n')
    end
    print(io, "<h2>Results</h2>\n")
    for r in results
        html_result_table(io, r, registry_entry(r.deviation_type).label)
        print(io, '\n')
    end
    print(io, "<h2>Reproduce</h2><pre>sigmatau")
    for i in eachindex(args)
        (args[i] == "--out" || (i > 1 && args[i-1] == "--out")) && continue
        print(io, ' '); html_escape(io, shell_quote(args[i]))
    end
    print(io, "</pre>\n</body></html>\n")
    return nothing
end

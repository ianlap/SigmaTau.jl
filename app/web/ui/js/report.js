// HTML report of the current analysis: one self-contained file (styles and
// SVG figures inline) that prints to PDF from any browser. It records what
// was analysed and how, so the numbers can be traced and reproduced.
import { state, visibleClocks, referenceClock } from "./state.js";
import { byId, label, gridName, noiseLabel } from "./estimators.js";
import { loadPlotly } from "./plot.js";
import { analyzeCommand } from "./core.js";
import { figure } from "./pages/stability.js";
import { sci, num, thousands, duration, escapeHtml as h } from "./format.js";

const AXIS_TITLE = { sigma_y: "Frequency stability σ<sub>y</sub>(τ)", sigma_x: "Time stability σ<sub>x</sub>(τ)" };
const DRIFT_TEXT = { none: "kept", linear: "offset removed", quadratic: "linear drift removed" };

async function sha256(text) {
  if (!crypto?.subtle) return "";
  const buf = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text));
  return [...new Uint8Array(buf)].map((b) => b.toString(16).padStart(2, "0")).join("");
}

// The Stability page's figure for one axis, drawn off screen as SVG.
async function figureSvg(clocks, axis) {
  const devs = state.devs.filter((d) => byId[d].axis === axis && clocks.some((c) => c.result?.results.some((r) => r.dev === d)));
  if (!devs.length) return "";
  const P = await loadPlotly();
  const { data, layout } = figure(clocks, devs, axis);
  layout.legend = { ...layout.legend, traceorder: "normal" };
  Object.assign(layout, { showlegend: true, width: 820, height: 460, paper_bgcolor: "#ffffff", plot_bgcolor: "#ffffff" });
  layout.font = { ...layout.font, color: "#1d1d1f" };
  const el = document.createElement("div");
  el.style.cssText = "position:fixed;left:-10000px;top:0;width:820px;height:460px";
  document.body.appendChild(el);
  try {
    await P.newPlot(el, data, layout, { staticPlot: true });
    const url = await P.toImage(el, { format: "svg", width: 820, height: 460 });
    return `<figure><img alt="${h(axis)} plot" src="${url}"><figcaption>${AXIS_TITLE[axis]}</figcaption></figure>`;
  } finally {
    P.purge(el);
    el.remove();
  }
}

function sourceText(c) {
  const s = c.source;
  if (s.type === "synth") return `synthetic, ${thousands(s.n)} samples, seed ${s.seed}`;
  return `${h(s.name)}${c.format ? ` (${h(c.format)})` : ` (${c.kind}, column ${c.value_col})`}`;
}

async function recordRows(clocks) {
  const ref = referenceClock();
  const rows = [];
  for (const c of clocks) {
    const r = c.record, k = r.check;
    const hash = c.source.type === "file" ? await sha256(c.source.content) : "";
    const prep = [`drift ${DRIFT_TEXT[c.drift] ?? c.drift}`, c.fillgaps ? "gaps filled" : "", c.remove_outliers ? "outliers removed" : "",
      ref && ref !== c && ref.group === c.group ? `minus ${h(ref.name)}` : ""].filter(Boolean).join(", ");
    rows.push(`<tr><th>${h(c.name)}</th><td>${sourceText(c)}${hash ? `<br><code class="hash">SHA-256 ${hash}</code>` : ""}</td>
      <td>${h(r.kind)}</td><td class="n">${thousands(r.n)}</td><td class="n">${num(r.tau0, 6)} s</td><td class="n">${duration(r.n * r.tau0)}</td>
      <td class="n">${k.gaps ?? 0}</td><td class="n">${k.outliers ?? 0}</td><td class="n">${sci(k.drift_per_day, 3)}</td><td>${prep}</td></tr>`);
  }
  return rows.join("");
}

function resultTable(c, r) {
  const ci = r.ci_lo.length === r.tau.length;
  const rows = r.tau.map((t, i) => `<tr><td class="n">${num(t, 5)}</td><td class="n">${sci(r.sigma[i], 4)}</td>${ci ? `<td class="n">${sci(r.ci_lo[i], 3)}</td><td class="n">${sci(r.ci_hi[i], 3)}</td>` : ""}${r.edf.length ? `<td class="n">${num(r.edf[i], 4)}</td>` : ""}${r.noise.length ? `<td>${h(noiseLabel(r.noise[i] ?? ""))}</td>` : ""}<td class="n">${r.neff[i] ?? ""}</td></tr>`).join("");
  return `<table class="results"><caption>${h(c.name)} · ${h(label(r.dev))}</caption>
    <thead><tr><th>τ (s)</th><th>σ</th>${ci ? "<th>CI low</th><th>CI high</th>" : ""}${r.edf.length ? "<th>EDF</th>" : ""}${r.noise.length ? "<th>Noise</th>" : ""}<th>Terms</th></tr></thead>
    <tbody>${rows}</tbody></table>`;
}

const STYLE = `
  body { font: 13px/1.5 "IBM Plex Sans", system-ui, sans-serif; color: #1d1d1f; background: #fff; margin: 32px auto; max-width: 900px; padding: 0 16px; }
  h1 { font-size: 22px; margin: 0 0 4px; } h2 { font-size: 16px; margin: 28px 0 8px; border-bottom: 1px solid #ddd; padding-bottom: 4px; }
  .meta { color: #666; margin: 0 0 16px; }
  table { border-collapse: collapse; width: 100%; margin: 8px 0 16px; font-size: 12px; }
  th, td { text-align: left; padding: 3px 8px; border-bottom: 1px solid #eee; vertical-align: top; }
  thead th { border-bottom: 1px solid #bbb; font-weight: 600; }
  td.n { text-align: right; font-family: "IBM Plex Mono", ui-monospace, monospace; white-space: nowrap; }
  caption { text-align: left; font-weight: 600; padding: 4px 0; }
  dl { display: grid; grid-template-columns: max-content 1fr; gap: 2px 16px; margin: 0; } dt { color: #666; } dd { margin: 0; }
  code, pre { font-family: "IBM Plex Mono", ui-monospace, monospace; font-size: 11.5px; }
  .hash { color: #666; word-break: break-all; }
  pre { background: #f6f6f6; padding: 8px 10px; white-space: pre-wrap; word-break: break-all; }
  figure { margin: 12px 0; } figure img { width: 100%; height: auto; border: 1px solid #eee; } figcaption { color: #666; font-size: 12px; }
  .results { width: auto; min-width: 60%; break-inside: avoid; }
  @media print { body { margin: 0; max-width: none; } h2 { break-after: avoid; } figure { break-inside: avoid; } }
`;

// Build the report from the clocks with results. Returns the HTML text.
export async function buildReport() {
  const clocks = visibleClocks().filter((c) => c.result && c.record);
  if (!clocks.length) throw new Error("Run the analysis first; there are no results to report.");
  const doc = clocks[0].result.doc;
  const when = new Date();
  const conf = state.ci ? `${num(state.confidence * 100, 3)} %` : "off";
  const ref = referenceClock();
  const figures = (await figureSvg(clocks, "sigma_y")) + (await figureSvg(clocks, "sigma_x"));
  const tables = clocks.flatMap((c) => c.result.results.map((r) => resultTable(c, r))).join("");
  const commands = clocks.map((c) => analyzeCommand(c)).join("\n");
  return `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Stability report</title><style>${STYLE}</style></head><body>
<h1>Frequency stability report</h1>
<p class="meta">${h(when.toISOString().replace("T", " ").slice(0, 19))} UTC · SigmaTau ${h(doc.sigmatau ?? "")} · ${clocks.length === 1 ? h(clocks[0].name) : `${clocks.length} clocks`}</p>
<h2>Analysis</h2>
<dl>
  <dt>Estimators</dt><dd>${state.devs.map((d) => h(label(d))).join(", ")}</dd>
  <dt>τ grid</dt><dd>${h(gridName(state.grid))}</dd>
  <dt>Confidence</dt><dd>${conf}${state.ci ? " (Greenhall EDF, χ² intervals)" : ""}</dd>
  <dt>Bias correction</dt><dd>${state.bias ? "on" : "off"}</dd>
  <dt>Reference clock</dt><dd>${ref ? h(ref.name) : "none (as measured)"}</dd>
</dl>
<h2>Records</h2>
<table><thead><tr><th>Clock</th><th>Source</th><th>Type</th><th>Samples</th><th>τ<sub>0</sub></th><th>Span</th><th>Gaps</th><th>Outliers</th><th>Drift / day</th><th>Preprocessing</th></tr></thead>
<tbody>${await recordRows(clocks)}</tbody></table>
<h2>Stability</h2>
${figures}
<h2>Results</h2>
${tables}
<h2>Reproduce</h2>
<pre>${h(commands)}</pre>
</body></html>
`;
}

// Ensemble page: a timescale (AT1 or EWFA) formed from several clocks read
// against one reference. The clocks are the value columns of one loaded file;
// without such a file the core simulates four. The scale is computed when this
// page needs it and recomputed when the settings or the clocks change.
import { state, update, isReference } from "../state.js";
import { timescale, synthArgs } from "../core.js";
import { loadPlotly, baseLayout, logAxis, linAxis, CONFIG, seriesColor, css } from "../plot.js";
import { sci, num, escapeHtml as h } from "../format.js";

let root, ctx, Plotly;
let inflight = null;           // key of the run in progress
let inflightCtl = null;        // stops it (the core process too) when the settings change
const ALGS = [["at1", "AT1"], ["at2", "AT2"], ["ewfa", "EWFA"]];

export const ownsPanel = true;

// The clocks the scale is formed from: the first loaded file with two or
// more phase or frequency columns (without the reference clock), else the
// first ensemble simulated on the Data page. [] = the core simulates four.
function fileClocks() {
  const groups = new Map();
  for (const c of state.clocks) {
    if (c.source.type !== "file" || c.group === null || c.format || !c.record || isReference(c)) continue;
    groups.set(c.group, [...(groups.get(c.group) ?? []), c]);
  }
  const file = [...groups.values()].find((g) => g.length >= 2);
  if (file) return file;
  const sims = new Map();
  for (const c of state.clocks) {
    const src = c.source;
    if (src.type !== "synth" || !(src.clocks > 1) || !c.record) continue;
    const key = JSON.stringify([src.n, src.seed, src.tau0, src.clocks, src.models ?? null]);
    sims.set(key, [...(sims.get(key) ?? []), c]);
  }
  const ens = [...sims.values()].find((g) => g.length === g[0].source.clocks);
  return ens ? ens.sort((a, b) => a.source.index - b.source.index) : [];
}

function request() {
  const E = state.ens, clocks = fileClocks();
  const a = ["timescale"];
  let file = null;
  if (clocks.length && clocks[0].source.type === "synth") {
    const src = clocks[0].source;
    a.push(...synthArgs(src), "--tau0", src.tau0);
  } else if (clocks.length) {
    const c = clocks[0];
    file = { name: c.source.name, content: c.source.content };
    a.push("--value-cols", clocks.map((k) => k.value_col).join(","), "--time-col", c.time_col, "--header", c.header, "--scaling", c.scaling);
    if (c.tau0) a.push("--tau0", c.tau0);
    if (c.kind === "frequency") a.push("--freq");
  } else a.push("--synth", 16384, "--seed", 42);
  a.push("--algorithm", E.algorithm, "--wmax", E.wmax, ...(E.algorithm === "ewfa" && E.tauY ? ["--tau-y", E.tauY] : []), "--format", "json", "--series", 2000);
  return { clocks, body: { args: a.map(String), file } };
}

export function controls() {
  const E = state.ens, clocks = fileClocks();
  return `<h2>Timescale</h2>
    <div class="ctl-block"><label class="lbl">Algorithm</label><div class="seg" role="group" aria-label="Algorithm">${ALGS.map(([id, label]) =>
      `<button type="button" data-ctl="alg" data-v="${id}" class="${E.algorithm === id ? "on" : ""}">${label}</button>`).join("")}</div></div>
    ${E.algorithm === "ewfa" ? `<label class="opt"><span>Frequency averaging τ<sub>y</sub> (s)</span><input type="text" data-ctl="tauY" value="${h(E.tauY)}" placeholder="100 τ₀" inputmode="decimal"></label>` : ""}
    <label class="opt"><span>Largest weight</span><input type="text" data-ctl="wmax" value="${h(String(E.wmax))}" inputmode="decimal"></label>
    <div class="rule"></div>
    <p class="muted small">${clocks.length
      ? (clocks[0].source.type === "synth"
        ? `${clocks.length} clocks simulated on the Data page, each against a perfect reference.`
        : `${clocks.length} clocks from ${h(clocks[0].source.name)}, each read against the same reference.`)
      : "No file with several clocks is loaded and none simulated, so four simulated clocks are used. Load a file with one column per clock, or simulate several clocks, on the Data page."}</p>
    <p class="muted small">AT1 weighs each clock by the inverse of its prediction-error variance, learned as it runs. EWFA weighs every clock the same.</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl;
  if (k === "alg" && e.type === "click") update((s) => { s.ens.algorithm = t.dataset.v; });
  if (e.type !== "change") return;
  const v = t.value.trim();
  if (k === "tauY") {
    if (v === "" || Number(v) > 0) update((s) => { s.ens.tauY = v; });
    else { t.value = state.ens.tauY; ctx.status("error", "τ_y must be a time above 0 s, or empty for 100 τ₀."); }
  }
  if (k === "wmax") {
    if (Number(v) > 0 && Number(v) <= 1) update((s) => { s.ens.wmax = Number(v); });
    else { t.value = state.ens.wmax; ctx.status("error", "The largest weight is a fraction between 0 and 1."); }
  }
}

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="card plot-card"><div class="card-row"><h2>ADEV of the scale and of each clock against it</h2><span class="muted small" data-ref="sub"></span></div>
        <div class="plot" data-ref="adev"></div><div class="empty" data-ref="empty" hidden></div></div>
      <div class="card plot-card ens-weights"><div class="card-row"><h2>Weights</h2></div><div class="plot" data-ref="weights"></div></div>
    </div>
    <div class="card side-card ens-table" data-ref="table"></div>`;
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

// Ask the core for the scale when it is missing or out of date.
function ensure() {
  const { body } = request();
  const key = JSON.stringify(body);
  if (state.ens.result?.key === key || inflight === key || !state.core.ok) return;
  inflightCtl?.abort();
  const ctl = inflightCtl = new AbortController();
  inflight = key;
  state.ens.error = null;
  timescale(body, ctl.signal).then((doc) => {
    if (inflight !== key) return;
    inflight = null;
    update((s) => { s.ens.result = { key, ts: doc.timescale }; });
  }).catch((err) => {
    if (err.name === "AbortError" || inflight !== key) return;
    inflight = null;
    update((s) => { s.ens.error = err.message; s.ens.result = null; });
  });
}

function showEmpty(text) {
  const e = root.querySelector("[data-ref=empty]");
  e.innerHTML = text ? `<p>${text}</p>` : "";
  e.hidden = !text;
  for (const p of root.querySelectorAll(".plot")) p.style.visibility = text ? "hidden" : "visible";
}

export function refresh() {
  if (!root) return;
  ensure();
  const { clocks } = request();
  const r = state.ens.result?.ts;
  const table = root.querySelector("[data-ref=table]");
  if (state.ens.error) { table.innerHTML = ""; return showEmpty(h(state.ens.error)); }
  if (!r) { table.innerHTML = ""; return showEmpty(state.core.ok ? "Forming the timescale…" : "No analysis core. Start the explorer with serve.py."); }
  showEmpty("");
  const names = r.clocks.map((n, i) => clocks[i]?.name ?? n);
  const color = (i) => seriesColor(clocks[i]?.slot ?? i);
  root.querySelector("[data-ref=sub]").textContent = `${r.algorithm.toUpperCase()} · ${r.n} samples · τ₀ = ${num(r.tau0, 4)} s · against ${r.reference || "the file's reference"}`;
  table.innerHTML = `<div class="card-row"><h2>Clocks</h2></div>
    <div class="tbl"><div class="tbl-row tbl-head"><span>Clock</span><span>Weight</span><span>Mean</span><span>RMS error</span></div>
    ${names.map((n, i) => `<div class="tbl-row tbl-body"><span class="ens-name"><span class="clock-dot" style="background:${color(i)}"></span>${h(n)}</span>
      <span>${num(r.weight[i], 3)}</span><span>${num(r.weight_mean[i], 3)}</span><span>${sci(r.error_rms[i], 2)} s</span></div>`).join("")}</div>
    <p class="muted small ens-note">Mean weight and RMS prediction error cover the second half of the record.</p>`;
  if (!Plotly) return;
  const scale = css("--accent");
  Plotly.react(root.querySelector("[data-ref=adev]"), [
    { x: r.adev.tau, y: r.adev.scale, name: "Timescale", mode: "lines+markers", line: { color: scale, width: 3 }, marker: { color: scale, size: 6 },
      hovertemplate: "<b>Timescale</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>" },
    ...r.adev.clocks.map((dev, i) => ({ x: r.adev.tau, y: dev, name: names[i], mode: "lines", line: { color: color(i), width: 1.5, dash: "dot" },
      hovertemplate: `<b>${h(names[i])}</b> − scale<br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>` })),
  ], { ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis("σ<sub>y</sub>(τ)"), showlegend: true, uirevision: "ens" }, CONFIG);
  const s = r.series;
  Plotly.react(root.querySelector("[data-ref=weights]"), s ? s.weights.map((w, i) => ({
    x: s.t, y: w, name: names[i], mode: "lines", line: { color: color(i), width: 1.5 },
    hovertemplate: `<b>${h(names[i])}</b><br>t = %{x:.6g} s<br>w = %{y:.3f}<extra></extra>` })) : [],
  { ...baseLayout(), margin: { l: 62, r: 14, t: 8, b: 44 }, xaxis: linAxis("Time (s)"), yaxis: { ...linAxis("Weight"), range: [0, 1] }, showlegend: false, uirevision: "ens-w" }, CONFIG);
}

export function plotElement() { return root?.querySelector("[data-ref=adev]") ?? null; }
export function hint() { return "Solid = the scale against the reference · dotted = each clock against the scale"; }

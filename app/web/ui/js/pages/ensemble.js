// Ensemble page: a timescale (AT1, AT2, EWFA or a Kalman scale) formed from
// several clocks read against one reference. The clocks are the value columns
// of one loaded file; without such a file the core simulates four. The scale
// is computed when this page needs it and recomputed when the settings or the
// clocks change.
import { state, update, isReference } from "../state.js";
import { timescale, synthArgs } from "../core.js";
import { loadPlotly, baseLayout, logAxis, linAxis, CONFIG, seriesColor, css } from "../plot.js";
import { sci, num, escapeHtml as h } from "../format.js";

let root, ctx, Plotly;
let inflight = null;           // key of the run in progress
let inflightCtl = null;        // stops it (the core process too) when the settings change
const ALGS = [["at1", "AT1"], ["at2", "AT2"], ["ewfa", "EWFA"], ["greenhall", "Greenhall"], ["jones_tryon", "Jones–Tryon"], ["stein", "Stein"]];
const KALMAN = new Set(["greenhall", "jones_tryon", "stein"]);
const capped = (alg) => !KALMAN.has(alg) || alg === "stein";   // takes --wmax
const ABOUT = {
  at1: "AT1 weighs each clock by the inverse of its prediction-error variance, learned as it runs.",
  at2: "AT2 keeps AT1's weights and follows each clock's frequency with a Kalman update. Frequency steps: best reruns the scale from a detected step; real keeps the values already issued and corrects from the step on.",
  ewfa: "EWFA weighs every clock the same.",
  greenhall: "Greenhall's Kalman scale fits a model to each clock and reads the scale at the pivot clock. It needs three or more clocks.",
  jones_tryon: "The Jones–Tryon Kalman scale: Greenhall's filter without the covariance reduction, so it has no per-clock weights. It needs three or more clocks.",
  stein: "Stein's Kalman scale (KAS-2) fits a model to each clock, filters each clock's difference from the pivot, and weighs each clock's phase, frequency and aging by the inverse of its noise. The weights shown are the phase weights. It needs three or more clocks.",
};
// The pivot as the core takes it: 1-based, back to the first clock when out of range.
const pivotOf = (n) => (state.ens.pivot <= n ? state.ens.pivot : 1);

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
  a.push("--algorithm", E.algorithm, ...(KALMAN.has(E.algorithm) ? ["--pivot", pivotOf(clocks.length || 4)] : []), ...(capped(E.algorithm) ? ["--wmax", E.wmax] : []),
    ...(E.algorithm === "at2" && E.steps ? ["--steps", E.steps] : []), "--format", "json", "--series", 2000);
  return { clocks, body: { args: a.map(String), file } };
}

export function controls() {
  const E = state.ens, clocks = fileClocks();
  return `<h2>Timescale</h2>
    <label class="opt"><span>Algorithm</span><select data-ctl="alg">${ALGS.map(([id, label]) =>
      `<option value="${id}" ${E.algorithm === id ? "selected" : ""}>${label}</option>`).join("")}</select></label>
    ${KALMAN.has(E.algorithm)
      ? `<label class="opt"><span>Pivot clock</span><select data-ctl="pivot">${(clocks.length ? clocks.map((c) => c.name) : [1, 2, 3, 4].map((k) => `clock ${k}`)).map((n, i) =>
        `<option value="${i + 1}" ${pivotOf(clocks.length || 4) === i + 1 ? "selected" : ""}>${h(n)}</option>`).join("")}</select></label>` : ""}
    ${capped(E.algorithm) ? `<label class="opt"><span>Largest weight</span><input type="text" data-ctl="wmax" value="${h(String(E.wmax))}" inputmode="decimal"></label>` : ""}
    ${E.algorithm === "at2" ? `<label class="opt"><span>Frequency steps</span><select data-ctl="steps">${[["", "Off"], ["best", "Best"], ["real", "Real"]].map(([id, label]) =>
      `<option value="${id}" ${E.steps === id ? "selected" : ""}>${label}</option>`).join("")}</select></label>` : ""}
    <div class="rule"></div>
    ${KALMAN.has(E.algorithm) && (clocks.length || 4) < 3 ? `<p class="small ens-warn">${ALGS.find(([id]) => id === E.algorithm)[1]} needs three or more clocks.</p>` : ""}
    <p class="muted small">${clocks.length
      ? (clocks[0].source.type === "synth"
        ? `${clocks.length} clocks simulated on the Data page, each against a perfect reference.`
        : `${clocks.length} clocks from ${h(clocks[0].source.name)}, each read against the same reference.`)
      : "No file with several clocks is loaded and none simulated, so four simulated clocks are used. Load a file with one column per clock, or simulate several clocks, on the Data page."}</p>
    <p class="muted small">${ABOUT[E.algorithm]}</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl;
  if (e.type !== "change") return;
  const v = t.value.trim();
  if (k === "alg") update((s) => { s.ens.algorithm = v; });
  if (k === "steps") update((s) => { s.ens.steps = v; });
  if (k === "pivot") update((s) => { s.ens.pivot = Number(v); });
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
  if (state.ens.result?.key === key || state.ens.error?.key === key || inflight === key || !state.core.ok) return;
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
    update((s) => { s.ens.error = { key, message: err.message }; s.ens.result = null; });
  });
}

const weight = (w) => (w === null ? "–" : num(w, 3));

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
  if (state.ens.error) { table.innerHTML = ""; return showEmpty(h(state.ens.error.message)); }
  if (!r) { table.innerHTML = ""; return showEmpty(state.core.ok ? "Forming the timescale…" : "No analysis core. Start the explorer with serve.py."); }
  showEmpty("");
  const names = r.clocks.map((n, i) => clocks[i]?.name ?? n);
  const color = (i) => seriesColor(clocks[i]?.slot ?? i);
  // Jones–Tryon has no weights: its plot goes and the ADEV plot takes the room.
  const weighted = r.weight.some((w) => w !== null);
  const wcard = root.querySelector(".ens-weights"), relaid = wcard.hidden === weighted;
  wcard.hidden = !weighted;
  const alg = ALGS.find(([id]) => id === r.algorithm)?.[1] ?? r.algorithm;
  root.querySelector("[data-ref=sub]").textContent = `${alg}${r.pivot ? ` at ${names[r.pivot - 1]}` : ""}${r.steps ? `, ${r.steps} steps` : ""} · ${r.n} samples · τ₀ = ${num(r.tau0, 4)} s · against ${r.reference || "the file's reference"}`;
  table.innerHTML = `<div class="card-row"><h2>Clocks</h2></div>
    <div class="tbl"><div class="tbl-row tbl-head"><span>Clock</span><span>Weight</span><span>Mean</span><span>RMS error</span></div>
    ${names.map((n, i) => `<div class="tbl-row tbl-body"><span class="ens-name"><span class="clock-dot" style="background:${color(i)}"></span>${h(n)}</span>
      <span>${weight(r.weight[i])}</span><span>${weight(r.weight_mean[i])}</span><span>${sci(r.error_rms[i], 2)} s</span></div>`).join("")}</div>
    <p class="muted small ens-note">${weighted ? "" : "The Jones–Tryon scale has no per-clock weights. "}Mean weight and RMS prediction error cover the second half of the record.</p>`;
  if (!Plotly) return;
  const scale = css("--accent");
  Plotly.react(root.querySelector("[data-ref=adev]"), [
    { x: r.adev.tau, y: r.adev.scale, name: "Timescale", mode: "lines+markers", line: { color: scale, width: 3 }, marker: { color: scale, size: 6 },
      hovertemplate: "<b>Timescale</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>" },
    ...r.adev.clocks.map((dev, i) => ({ x: r.adev.tau, y: dev, name: names[i], mode: "lines", line: { color: color(i), width: 1.5, dash: "dot" },
      hovertemplate: `<b>${h(names[i])}</b> − scale<br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>` })),
  ], { ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis("σ<sub>y</sub>(τ)"), showlegend: true, uirevision: "ens" }, CONFIG);
  if (relaid) Plotly.Plots.resize(root.querySelector("[data-ref=adev]"));
  const s = r.series;
  Plotly.react(root.querySelector("[data-ref=weights]"), s && weighted ? s.weights.map((w, i) => ({
    x: s.t, y: w, name: names[i], mode: "lines", line: { color: color(i), width: 1.5 },
    hovertemplate: `<b>${h(names[i])}</b><br>t = %{x:.6g} s<br>w = %{y:.3f}<extra></extra>` })) : [],
  { ...baseLayout(), margin: { l: 62, r: 14, t: 8, b: 44 }, xaxis: linAxis("Time (s)"), yaxis: { ...linAxis("Weight"), range: [0, 1] }, showlegend: false, uirevision: "ens-w" }, CONFIG);
}

export function plotElement() { return root?.querySelector("[data-ref=adev]") ?? null; }
export function hint() { return "Solid = the scale against the reference · dotted = each clock against the scale"; }

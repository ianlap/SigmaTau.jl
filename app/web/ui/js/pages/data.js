// Data page: load a file (with the import dialog) or synthesize a record,
// the record check, and the phase plot with outliers in orange.
import { state, update, activeClock, multi, markStale } from "../state.js";
import { loadPlotly, baseLayout, linAxis, CONFIG, css, seriesColor, timeUnit } from "../plot.js";
import { thousands, duration, sci, escapeHtml as h } from "../format.js";

const FILE_ICON = `<svg width="22" height="26" viewBox="0 0 22 26" fill="none" stroke="currentColor" stroke-width="1.8" aria-hidden="true"><path d="M3 2h10l6 6v16H3z"/><path d="M13 2v6h6"/></svg>`;
// Fitted to the phase and removed before every analysis. On phase a line is a
// frequency offset and a parabola adds a linear frequency drift.
const DRIFTS = [
  ["none", "None", "Analyse the record as measured"],
  ["linear", "Linear", "Remove a frequency offset (a straight line in phase)"],
  ["quadratic", "Quadratic", "Remove a frequency offset and a linear frequency drift (a parabola in phase), as SP1065 asks before the total estimators"],
];

let root, plotEl, Plotly, ctx, relayoutBound = false;

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="toolbar">
        <div class="seg" role="group" aria-label="Source" data-ref="mode">
          <button type="button" data-mode="file">Load a file</button><button type="button" data-mode="synth">Synthesize</button></div>
        <span class="caption" data-ref="caption"></span>
      </div>
      <div class="data-top" data-ref="top"></div>
      <div class="card record-card">
        <div class="card-row"><h2 data-ref="plot-title">Phase x(t)</h2><span class="range-bar" data-ref="range"></span></div>
        <div class="plot" data-ref="plot"></div><div class="empty" data-ref="empty" hidden></div>
      </div>
    </div>`;
  plotEl = root.querySelector("[data-ref=plot]");
  relayoutBound = false;
  root.addEventListener("click", onClick);
  root.addEventListener("change", onChange);
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

function fileCard(c) {
  const isFile = c?.source.type === "file";
  const nameField = multi()
    ? `<input class="clock-name" data-field="name" value="${h(c.name)}" aria-label="Clock name" title="Rename this clock">`
    : `<b>${h(isFile ? c.source.name : c?.name ?? "")}</b>`;
  return `<div class="card file-card">
    ${isFile ? `<div class="file-title">${FILE_ICON}<div style="flex:1;min-width:0">${nameField}
        <span class="muted small">${multi() ? `${h(c.source.name)} · ` : ""}${h(c.sourceInfo)}</span></div></div>
      <div class="file-actions"><button class="btn" data-act="import">Import settings…</button><button class="btn" data-act="open">Open another file…</button>${addClockBtn()}</div>`
    : `<div class="file-actions"><button class="btn btn-primary" data-act="open">Open a file…</button>${addClockBtn()}</div>`}
    <div class="dropzone" data-ref="drop">Drop a data file anywhere in the window</div>
    ${c?.error ? `<p class="error small">${h(c.error)}</p>` : ""}
  </div>`;
}

const addClockBtn = () => state.clocks.some((c) => c.record) ? `<button class="btn" data-act="add-clock">+ Add clock</button>` : "";

function synthCard(c) {
  const s = c?.source.type === "synth" ? c.source : { n: 16384, seed: 42, tau0: 1 };
  return `<div class="card synth-card">
    <div><b>Synthetic record</b><div class="muted small">White plus random-walk frequency noise, made in the core. Same seed, same record.</div></div>
    <div class="synth-grid">
      <label>Samples<input type="number" data-synth="n" min="64" max="10000000" step="1" value="${s.n}"></label>
      <label>τ<sub>0</sub> (s)<input type="number" data-synth="tau0" min="0" step="any" value="${s.tau0}"></label>
      <label>Seed<input type="number" data-synth="seed" min="0" step="1" value="${s.seed}"></label>
    </div>
    <div class="file-actions"><button class="btn btn-primary" data-act="synth">${c?.source.type === "synth" ? "Make record" : "Use a synthetic record"}</button>${addClockBtn()}</div>
  </div>`;
}

// NaN samples from the core plus samples missing from the time column, which
// the readers skip unless gaps are filled.
function gapRow(c, k) {
  const missing = c.fillgaps ? 0 : (c.missing ?? 0);
  const n = k.gaps + missing;
  const tip = `${thousands(k.gaps)} NaN sample${k.gaps === 1 ? "" : "s"} · ${thousands(missing)} missing from the time column${missing ? " (turn on Fill gaps to put them back)" : ""}`;
  return `<div class="kv" title="${tip}"><span>Gaps</span><span class="mono ${n ? "warnv" : ""}">${thousands(n)}</span></div>`;
}

function checkCard(c) {
  const k = c?.record?.check;
  if (!k) return `<div class="card check-card"><h2>Record check</h2><p class="muted small">${c?.error ? "The core could not read this record." : c ? "Reading…" : "No record yet."}</p></div>`;
  const drift = Number.isFinite(k.drift_per_day) ? `${sci(k.drift_per_day, 2)} /day` : "–";
  const out = c.remove_outliers
    ? `<span class="mono">removed · <button class="link-btn" data-act="keep-outliers">undo</button></span>`
    : k.outliers ? `<span class="mono warnv">${k.outliers} · <button class="link-btn" data-act="remove-outliers">remove</button></span>`
    : `<span class="mono">0</span>`;
  return `<div class="card check-card">
    <h2>Record check</h2>
    <div class="kv"><span>Samples</span><span class="mono">${thousands(k.samples)}</span></div>
    <div class="kv"><span>Span</span><span class="mono">${duration(k.span_s)}</span></div>
    ${gapRow(c, k)}
    <div class="kv"><span>Points past 5σ</span>${out}</div>
    <div class="kv" title="Least-squares linear frequency drift of the record as analysed"><span>${c.drift === "quadratic" ? "Drift left" : "Frequency drift"}</span><span class="mono">${drift}</span></div>
    <div class="rule"></div>
    ${driftBlock(c, k)}
    <label class="check" title="${c.hasTimeCol ? "Put missing samples back on an even grid and fill them (Howe's method)" : "Needs a time column in the file"}">
      <input type="checkbox" data-field="fillgaps" ${c.fillgaps ? "checked" : ""} ${c.hasTimeCol ? "" : "disabled"}>Fill gaps by interpolation</label>
  </div>`;
}

function driftBlock(c, k) {
  const mode = c.drift ?? "none";
  const removed = [];
  if (Number.isFinite(k.removed_offset)) removed.push(`offset ${sci(k.removed_offset, 3)}`);
  if (Number.isFinite(k.removed_drift)) removed.push(`drift ${sci(k.removed_drift, 3)} /day`);
  return `<div class="drift-block">
    <div class="opt"><span title="Fitted to the phase by least squares and removed before every analysis">Drift removal</span>
      <div class="seg seg-sm" role="group" aria-label="Drift removal">${DRIFTS.map(([v, l, tip]) =>
        `<button type="button" data-drift="${v}" class="${v === mode ? "on" : ""}" title="${tip}">${l}</button>`).join("")}</div></div>
    <div class="muted small drift-note">${mode === "none" ? "Frequency offset and drift stay in the record." : removed.length ? `Removed ${removed.join(" · ")}` : "Fitting…"}</div>
  </div>`;
}

export function refresh() {
  if (!root) return;
  const c = activeClock();
  const mode = state.dataMode;
  root.querySelectorAll("[data-mode]").forEach((b) => b.classList.toggle("on", b.dataset.mode === mode));
  root.querySelector("[data-ref=caption]").innerHTML = mode === "file"
    ? "Phase or frequency records: text columns, Stable32, TimeLab .tim, counter logs"
    : "A test record with known noise, for trying estimators";
  const focus = document.activeElement?.dataset?.field === "name";
  if (!focus) root.querySelector("[data-ref=top]").innerHTML = (mode === "file" ? fileCard(c) : synthCard(c)) + checkCard(c);
  else root.querySelector(".check-card").outerHTML = checkCard(c);
  const title = c?.drift && c.drift !== "none" ? `Phase residual x(t), ${c.drift} drift removed` : "Phase x(t)";
  root.querySelector("[data-ref=plot-title]").textContent = multi() && c ? `${title} · ${c.name}` : title;
  root.querySelector("[data-ref=range]").innerHTML = rangeBar(c);
  drawRecord(c);
}

// The analysis range: drag across the plot in "Select range" mode. The check
// above and this plot keep the whole record; the runs use the range.
function rangeBar(c) {
  if (!c?.record) return "";
  const r = c.range;
  return `${r ? `<span class="small">Analysing <span class="mono">${duration(r[0])} – ${duration(r[1])}</span> · <button class="link-btn" data-act="clear-range">whole record</button></span>` : `<span class="muted small">Outliers in orange</span>`}
    <span class="seg seg-sm" role="group" aria-label="Drag on the plot to">${[["zoom", "Zoom", "Drag to zoom"], ["range", "Select range", "Drag across the plot to analyse only that part of the record"]].map(([v, l, tip]) =>
      `<button type="button" data-drag="${v}" class="${state.dataDrag === v ? "on" : ""}" title="${tip}">${l}</button>`).join("")}</span>`;
}

function setRange(c, r) {
  update((s) => { c.range = r; c.result = null; c.dynamic = null; markStale(s, "record"); });
}

function drawRecord(c) {
  const empty = root.querySelector("[data-ref=empty]");
  const s = c?.record?.series;
  if (!s || !s.t.length) {
    empty.hidden = false;
    empty.innerHTML = `<p>${c?.error ? h(c.error) : c ? "Reading the record…" : "Open a file or drop one anywhere in the window."}</p>`;
    plotEl.style.visibility = "hidden";
    return;
  }
  empty.hidden = true;
  plotEl.style.visibility = "visible";
  if (!Plotly) return;
  const k = c.record.check;
  const finite = s.x.filter(Number.isFinite);
  const maxAbs = Math.max(...finite.map(Math.abs), ...(k?.outlier_x ?? []).map(Math.abs));
  const [f, u] = timeUnit(maxAbs || 1);
  const traces = [{
    x: s.t, y: s.x.map((v) => v / f), mode: "lines", line: { color: multi() ? seriesColor(c.slot) : css("--s1"), width: 1.4 },
    hovertemplate: `t = %{x:.6g} s<br>x = %{y:.5g} ${u}<extra></extra>`, name: "phase",
  }];
  if (k?.outlier_t?.length && !c.remove_outliers) {
    traces.push({
      x: k.outlier_t, y: k.outlier_x.map((v) => v / f), mode: "markers", name: "outliers",
      marker: { color: css("--outlier") || css("--s2"), size: 8 },
      hovertemplate: `outlier<br>t = %{x:.6g} s<br>x = %{y:.5g} ${u}<extra></extra>`,
    });
  }
  const picking = state.dataDrag === "range";
  const yaxis = { ...linAxis(""), ticksuffix: ` ${u}`, fixedrange: picking };
  const shapes = c.range ? [{ type: "rect", xref: "x", yref: "paper", x0: c.range[0], x1: c.range[1], y0: 0, y1: 1,
    fillcolor: css("--accent-soft") || "rgba(158,27,50,0.12)", line: { width: 0 }, layer: "below" }] : [];
  Plotly.react(plotEl, traces, {
    ...baseLayout(), margin: { l: 70, r: 18, t: 10, b: 46 }, xaxis: linAxis("Time t (s)"), yaxis, shapes,
    dragmode: "zoom", uirevision: `${c.id}-${picking}`,
  }, { ...CONFIG, showTips: false });
  if (!relayoutBound) {
    relayoutBound = true;
    // In "Select range" mode a horizontal drag is a zoom we turn into the range.
    plotEl.on("plotly_relayout", (ev) => {
      const k = activeClock();
      const a = ev["xaxis.range[0]"], b = ev["xaxis.range[1]"];
      if (state.dataDrag !== "range" || !k?.record || a === undefined || b === undefined) return;
      const end = (k.record.n - 1) * k.record.tau0;
      const t0 = Math.max(0, Math.min(a, b)), t1 = Math.min(end, Math.max(a, b));
      const snap = (t) => Math.round(t / k.record.tau0) * k.record.tau0;
      Plotly.relayout(plotEl, { "xaxis.autorange": true });
      if (snap(t1) - snap(t0) >= 16 * k.record.tau0) setRange(k, [snap(t0), snap(t1)]);
      else ctx.status("warn", "The range needs at least 16 samples.");
    });
  }
}

function onClick(e) {
  const dr = e.target.closest("[data-drift]")?.dataset.drift;
  if (dr) { const c = activeClock(); if (c && c.drift !== dr) ctx.setRecord(c, { drift: dr }); return; }
  const dg = e.target.closest("[data-drag]")?.dataset.drag;
  if (dg) { update((s) => { s.dataDrag = dg; }); return; }
  const m = e.target.closest("[data-mode]")?.dataset.mode;
  if (m) { update((s) => { s.dataMode = m; }); return; }
  const act = e.target.closest("[data-act]")?.dataset.act;
  const c = activeClock();
  if (act === "open") ctx.openFile(c);
  else if (act === "import") ctx.reimport(c);
  else if (act === "add-clock") ctx.addClock();
  else if (act === "remove-outliers") ctx.setRecord(c, { remove_outliers: true });
  else if (act === "keep-outliers") ctx.setRecord(c, { remove_outliers: false });
  else if (act === "clear-range") setRange(c, null);
  else if (act === "synth") {
    const v = (k) => Number(root.querySelector(`[data-synth=${k}]`).value);
    const n = Math.round(v("n")), tau0 = v("tau0"), seed = Math.round(v("seed"));
    if (!(n >= 64) || !(tau0 > 0) || !(seed >= 0)) { ctx.status("error", "Samples must be at least 64 and τ₀ above zero."); return; }
    ctx.useSynth(c, { n, tau0, seed });
  }
}

function onChange(e) {
  const f = e.target.dataset.field;
  const c = activeClock();
  if (!c || !f) return;
  if (f === "name") { const v = e.target.value.trim() || c.name; update(() => { c.name = v; }); }
  else if (f === "fillgaps") ctx.setRecord(c, { fillgaps: e.target.checked });
}

export function plotElement() { return plotEl; }
export function hint() {
  return state.dataDrag === "range"
    ? "Drag across the plot to choose the part of the record to analyse"
    : "Drop a file anywhere in the window · drag on the plot to zoom · double-click to reset";
}

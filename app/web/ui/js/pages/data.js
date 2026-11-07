// Data page: load a file (with the import dialog) or synthesize a record,
// the record check, and the phase plot with outliers in orange.
import { state, update, activeClock, multi } from "../state.js";
import { loadPlotly, baseLayout, linAxis, CONFIG, css, seriesColor, timeUnit } from "../plot.js";
import { thousands, duration, sci, escapeHtml as h } from "../format.js";

const FILE_ICON = `<svg width="22" height="26" viewBox="0 0 22 26" fill="none" stroke="currentColor" stroke-width="1.8" aria-hidden="true"><path d="M3 2h10l6 6v16H3z"/><path d="M13 2v6h6"/></svg>`;
const DETRENDS = [["none", "none"], ["mean", "mean"], ["endpoint", "endpoint"], ["linear", "linear"]];

let root, plotEl, Plotly, ctx;

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
        <div class="card-row"><h2 data-ref="plot-title">Phase x(t)</h2><span class="muted small">Outliers marked in orange · drag to zoom</span></div>
        <div class="plot" data-ref="plot"></div><div class="empty" data-ref="empty" hidden></div>
      </div>
    </div>`;
  plotEl = root.querySelector("[data-ref=plot]");
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
    <div class="dropzone" data-ref="drop">Drop a .dat, .txt, .csv or .tsv file anywhere in the window</div>
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
    <div class="kv"><span>Gaps (NaN)</span><span class="mono ${k.gaps ? "warnv" : ""}">${thousands(k.gaps)}</span></div>
    <div class="kv"><span>Points past 5σ</span>${out}</div>
    <div class="kv"><span>Frequency drift</span><span class="mono">${drift}</span></div>
    <div class="rule"></div>
    <label class="opt"><span>Detrend</span><select data-field="detrend">${DETRENDS.map(([v, l]) => `<option value="${v}" ${v === c.detrend ? "selected" : ""}>${l}</option>`).join("")}</select></label>
    <label class="check" title="${c.hasTimeCol ? "Put missing samples back on an even grid and fill them (Howe's method)" : "Needs a time column in the file"}">
      <input type="checkbox" data-field="fillgaps" ${c.fillgaps ? "checked" : ""} ${c.hasTimeCol ? "" : "disabled"}>Fill gaps by interpolation</label>
  </div>`;
}

export function refresh() {
  if (!root) return;
  const c = activeClock();
  const mode = state.dataMode;
  root.querySelectorAll("[data-mode]").forEach((b) => b.classList.toggle("on", b.dataset.mode === mode));
  root.querySelector("[data-ref=caption]").innerHTML = mode === "file"
    ? "Phase or frequency records as text: .dat, .txt, .csv, .tsv"
    : "A test record with known noise, for trying estimators";
  const focus = document.activeElement?.dataset?.field === "name";
  if (!focus) root.querySelector("[data-ref=top]").innerHTML = (mode === "file" ? fileCard(c) : synthCard(c)) + checkCard(c);
  else root.querySelector(".check-card").outerHTML = checkCard(c);
  root.querySelector("[data-ref=plot-title]").textContent = multi() && c ? `Phase x(t) · ${c.name}` : "Phase x(t)";
  drawRecord(c);
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
  const yaxis = { ...linAxis(""), ticksuffix: ` ${u}` };
  Plotly.react(plotEl, traces, { ...baseLayout(), margin: { l: 70, r: 18, t: 10, b: 46 }, xaxis: linAxis("Time t (s)"), yaxis, uirevision: c.id }, CONFIG);
}

function onClick(e) {
  const m = e.target.closest("[data-mode]")?.dataset.mode;
  if (m) { update((s) => { s.dataMode = m; }); return; }
  const act = e.target.closest("[data-act]")?.dataset.act;
  const c = activeClock();
  if (act === "open") ctx.openFile(c);
  else if (act === "import") ctx.reimport(c);
  else if (act === "add-clock") ctx.addClock();
  else if (act === "remove-outliers") ctx.setRecord(c, { remove_outliers: true });
  else if (act === "keep-outliers") ctx.setRecord(c, { remove_outliers: false });
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
  else if (f === "detrend") ctx.setRecord(c, { detrend: e.target.value });
  else if (f === "fillgaps") ctx.setRecord(c, { fillgaps: e.target.checked });
}

export function plotElement() { return plotEl; }
export function hint() {
  return "Drop a file anywhere in the window · drag on the plot to zoom · double-click to reset";
}

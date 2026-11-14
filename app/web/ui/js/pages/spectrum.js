// Spectrum page: S_y(f) on top and ℒ(f) (or S_x(f)) below, all clocks
// overlaid, both plots on one frequency axis. The core returns S_y and S_x
// with each run; the carrier, phase-noise view, averaging and slope guides
// are display settings applied here without a run.
import { state, update, visibleClocks, multi } from "../state.js";
import { loadPlotly, baseLayout, logAxis, linAxis, CONFIG, seriesColor, css } from "../plot.js";
import { num, escapeHtml as h } from "../format.js";

let root, ctx, Plotly, top, bottom, syncing = false;

const SLOPES = { "-2": "random-walk FM", "-1": "flicker FM", "0": "white FM", "1": "flicker PM", "2": "white PM" };
const SUP = { "-2": "⁻²", "-1": "⁻¹", "0": "⁰", "1": "¹", "2": "²" };

export function carrierLabel(hz) {
  if (!(hz > 0)) return "–";
  const units = [[1e9, "GHz"], [1e6, "MHz"], [1e3, "kHz"], [1, "Hz"]];
  const [f, u] = units.find(([f]) => hz >= f) ?? [1, "Hz"];
  return `${num(hz / f, 4)} ${u}`;
}

// The settings, as panel content (one clock) or a toolbar (several clocks).
export function controls(inline = false) {
  const sp = state.spec;
  const carrier = `<input type="text" data-ctl="carrier" value="${h(sp.carrierText)}" inputmode="decimal" aria-label="Carrier frequency in Hz" class="carrier">`;
  const phase = `<div class="seg" role="group" aria-label="Phase noise shown as">
      <button type="button" data-ctl="phase" data-v="L" class="${sp.phase === "L" ? "on" : ""}">ℒ(f) dBc/Hz</button>
      <button type="button" data-ctl="phase" data-v="Sx" class="${sp.phase === "Sx" ? "on" : ""}">S<sub>x</sub>(f)</button></div>`;
  const avg = `<select data-ctl="avg" aria-label="Averaging"><option value="log" ${sp.avg === "log" ? "selected" : ""}>Log-binned</option><option value="raw" ${sp.avg === "raw" ? "selected" : ""}>Welch, every bin</option></select>`;
  const guides = `<label class="check"><input type="checkbox" data-ctl="guides" ${sp.guides ? "checked" : ""}>Show power-law slope guides</label>`;
  if (inline) {
    return `<label class="tool-field">Carrier for ℒ(f) ${carrier} Hz</label>${phase}${avg}`;
  }
  return `<h2>Spectrum</h2>
    <div class="ctl-block"><label class="lbl">Carrier frequency for ℒ(f)</label><div class="tau-line">${carrier}<span class="muted">Hz</span></div></div>
    <div class="ctl-block"><label class="lbl">Phase noise shown as</label>${phase}</div>
    <label class="opt"><span>Averaging</span>${avg}</label>
    ${guides}
    <div class="rule"></div>
    <p class="muted small">Spectrum settings apply as soon as you change them. The estimator choice on the Stability page doesn't affect this page.</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl;
  if (k === "phase" && e.type === "click") update((s) => { s.spec.phase = t.dataset.v; });
  if (k === "avg" && e.type === "change") update((s) => { s.spec.avg = t.value; });
  if (k === "guides" && e.type === "change") update((s) => { s.spec.guides = t.checked; });
  if (k === "carrier" && e.type === "change") {
    const v = Number(t.value.replace(/\s/g, ""));
    if (v > 0) update((s) => { s.spec.carrier = v; s.spec.carrierText = t.value.trim(); });
    else { t.value = state.spec.carrierText; ctx.status("error", "The carrier must be a frequency above 0 Hz, e.g. 10e6."); }
  }
}

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="toolbar" data-ref="toolbar"></div>
      <p class="spec-note muted small" data-ref="note"></p>
      <div class="card plot-card spec-card"><div class="card-row"><h2>Fractional frequency S<sub>y</sub>(f)</h2><span class="muted small">1/Hz</span></div>
        <div class="plot" data-ref="top"></div><div class="empty" data-ref="empty" hidden></div></div>
      <div class="card plot-card spec-card"><div class="card-row"><h2 data-ref="bottom-title"></h2><span class="muted small" data-ref="bottom-units"></span></div>
        <div class="plot" data-ref="bottom"></div></div>
    </div>`;
  top = root.querySelector("[data-ref=top]");
  bottom = root.querySelector("[data-ref=bottom]");
  root.addEventListener("click", (e) => {
    if (e.target.closest("[data-act=load]")) ctx.goData();
    if (e.target.closest("[data-act=run]")) ctx.run();
  });
  loadPlotly().then((P) => { Plotly = P; refresh(); bindSync(); }).catch(() => {});
}

// Zooming either plot moves the frequency axis of both.
function bindSync() {
  const link = (from, to) => from.on?.("plotly_relayout", (ev) => {
    if (syncing) return;
    const upd = {};
    if ("xaxis.range[0]" in ev) { upd["xaxis.range[0]"] = ev["xaxis.range[0]"]; upd["xaxis.range[1]"] = ev["xaxis.range[1]"]; }
    if (ev["xaxis.autorange"]) upd["xaxis.autorange"] = true;
    if (!Object.keys(upd).length) return;
    syncing = true;
    Promise.resolve(Plotly.relayout(to, upd)).finally(() => { syncing = false; });
  });
  link(top, bottom);
  link(bottom, top);
}

// Mean in log-spaced bins (10 per decade): smooths the Welch estimate where
// it is densest, at high frequency.
function logBin(f, p) {
  const out = { f: [], p: [] };
  let i = 0;
  while (i < f.length) {
    const edge = f[i] * Math.pow(10, 0.1);
    let sf = 0, sp = 0, n = 0;
    while (i < f.length && f[i] < edge) {
      if (Number.isFinite(p[i]) && p[i] > 0) { sf += Math.log(f[i]); sp += p[i]; n++; }
      i++;
    }
    if (n) { out.f.push(Math.exp(sf / n)); out.p.push(sp / n); }
  }
  return out;
}

const shown = (f, p) => (state.spec.avg === "log" ? logBin(f, p) : { f, p });
const toL = (sx) => sx.map((v) => 10 * Math.log10(0.5 * Math.pow(2 * Math.PI * state.spec.carrier, 2) * v));

// Label stretches whose log-log slope is near a whole power law.
function slopeGuides(f, p) {
  const n = f.length;
  if (n < 9) return [];
  const notes = [];
  for (let k = 0; k < 3; k++) {
    const a = Math.floor((k * n) / 3), b = Math.floor(((k + 1) * n) / 3);
    const xs = [], ys = [];
    for (let i = a; i < b; i++) if (p[i] > 0) { xs.push(Math.log10(f[i])); ys.push(Math.log10(p[i])); }
    if (xs.length < 3) continue;
    const mx = xs.reduce((s, v) => s + v, 0) / xs.length, my = ys.reduce((s, v) => s + v, 0) / ys.length;
    let sxy = 0, sxx = 0;
    xs.forEach((x, i) => { sxy += (x - mx) * (ys[i] - my); sxx += (x - mx) ** 2; });
    const slope = sxy / sxx, r = Math.round(slope);
    if (Math.abs(slope - r) > 0.35 || !SLOPES[r]) continue;
    const mid = Math.floor((a + b) / 2);
    if (notes.length && notes[notes.length - 1].r === r) continue;
    notes.push({ r, x: f[mid], y: p[mid] });
  }
  const accent = css("--accent");
  return notes.map((g) => ({
    x: Math.log10(g.x), y: Math.log10(g.y), xref: "x", yref: "y", showarrow: false, yshift: 16,
    text: `<b>f${SUP[g.r]} ${SLOPES[g.r]}</b>`, font: { color: accent, size: 13 },
  }));
}

function showEmpty(html) {
  const e = root.querySelector("[data-ref=empty]");
  e.innerHTML = html;
  e.hidden = !html;
  for (const p of [top, bottom]) p.style.visibility = html ? "hidden" : "visible";
}

export function refresh() {
  if (!root) return;
  const several = multi();
  const tb = root.querySelector("[data-ref=toolbar]");
  tb.hidden = !several;
  if (several && !tb.querySelector("input:focus, select:focus")) {   // a field in use keeps its text
    const at = tb.contains(document.activeElement) ? document.activeElement.dataset : null;
    tb.innerHTML = controls(true);
    if (at) tb.querySelector(`[data-ctl=${at.ctl}][data-v="${at.v}"]`)?.focus();
  }
  const isL = state.spec.phase === "L";
  root.querySelector("[data-ref=bottom-title]").innerHTML = isL ? `Phase noise ℒ(f) at ${carrierLabel(state.spec.carrier)}` : "Phase S<sub>x</sub>(f)";
  root.querySelector("[data-ref=bottom-units]").textContent = isL ? "dBc/Hz" : "s²/Hz";

  // The core analyzes a record with gaps without its spectrum.
  const analyzed = visibleClocks().filter((c) => c.result?.spectrum);
  const clocks = analyzed.filter((c) => !c.result.spectrum.error);
  const gapped = analyzed.filter((c) => c.result.spectrum.error).map((c) => h(c.name));
  const why = "A spectrum needs every sample: turn on <i>Fill gaps by reflection</i> on the Data page, or select a range without gaps.";
  root.querySelector("[data-ref=note]").innerHTML = clocks.length && gapped.length ? `No spectrum for ${gapped.join(", ")}, which ${gapped.length > 1 ? "have" : "has"} gaps. ${why}` : "";
  if (!state.clocks.some((c) => c.record)) return showEmpty(`<h2>No record yet</h2><p>Load or synthesize a record on the Data page.</p><button class="btn" data-act="load">Go to Data</button>`);
  if (!analyzed.length) return showEmpty(`<h2>Not analyzed yet</h2><p>Run the analysis to compute the spectra.</p><button class="btn btn-primary" data-act="run">Run analysis</button>`);
  if (!clocks.length) return showEmpty(`<h2>No spectrum: the record has gaps</h2><p>${why}</p><button class="btn" data-act="load">Go to Data</button>`);
  showEmpty("");
  if (!Plotly) return;

  const topTraces = [], botTraces = [];
  let notes = [];
  clocks.forEach((c, i) => {
    const sp = c.result.spectrum;
    const color = several ? seriesColor(c.slot) : css("--s1");
    const sy = shown(sp.f, sp.sy), sx = shown(sp.f, sp.sx);
    const name = several ? c.name : "S_y";
    topTraces.push({ x: sy.f, y: sy.p, mode: "lines", name, line: { color, width: 2 },
      hovertemplate: `${several ? `<b>${h(c.name)}</b><br>` : ""}f = %{x:.4g} Hz<br>S<sub>y</sub> = %{y:.3e} /Hz<extra></extra>` });
    botTraces.push({ x: sx.f, y: isL ? toL(sx.p) : sx.p, mode: "lines", name,
      line: { color: several ? color : css("--s3"), width: 2 },
      hovertemplate: `${several ? `<b>${h(c.name)}</b><br>` : ""}f = %{x:.4g} Hz<br>${isL ? "ℒ = %{y:.1f} dBc/Hz" : "S<sub>x</sub> = %{y:.3e} s²/Hz"}<extra></extra>` });
    if (!several && i === 0 && state.spec.guides) notes = slopeGuides(logBin(sp.f, sp.sy).f, logBin(sp.f, sp.sy).p);
  });
  const base = { ...baseLayout(), margin: { l: 70, r: 18, t: 8, b: 30 }, showlegend: false };
  const range = top.layout?.xaxis?.range;
  const xa = (title) => ({ ...logAxis(title), ...(range && !top.layout.xaxis.autorange ? { range, autorange: false } : {}) });
  Plotly.react(top, topTraces, { ...base, xaxis: xa(""), yaxis: logAxis(""), annotations: notes, uirevision: "spec" }, CONFIG);
  Plotly.react(bottom, botTraces, {
    ...base, margin: { l: 70, r: 18, t: 8, b: 48 }, xaxis: xa("Fourier frequency f (Hz)"),
    yaxis: isL ? { ...linAxis(""), ticksuffix: "" } : logAxis(""), uirevision: `spec-${state.spec.phase}`,
  }, CONFIG);
}

export function hint() {
  return multi() ? "Color = clock · both plots share the frequency axis" : "Both plots share the frequency axis · drag either to zoom both";
}

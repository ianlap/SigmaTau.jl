// Steering page: each clock steered to its reference, in replay. The core
// fits a clock model to the first stretch of the record, runs a Kalman
// filter and the steering loop over the whole record and reports from the
// end of that stretch on. One clock: the Allan deviation free-running,
// steered and as the model predicts for the settled loop, or the steered
// phase with the frequency correction in force. Several clocks: one colour
// per clock. The gains come from a time constant (pole placement) or from a
// cost. Computed when this page needs it; no Run is needed.
import { state, update, visibleClocks, multi } from "../state.js";
import { steer as fetchSteer, steerKey, steerDefault } from "../core.js";
import { loadPlotly, baseLayout, logAxis, linAxis, CONFIG, seriesColor, css, timeUnit } from "../plot.js";
import { sci, num, duration, modelText, escapeHtml as h } from "../format.js";

let root, ctx, Plotly;
const inflight = new Map();     // clock id → key being fetched
const FRACS = [8, 4, 2];
const EVERY = [1, 2, 4, 8, 16, 32, 64];
const WEIGHTS = ["phase²", "frequency²", "steer²"];

export function controls(inline = false) {
  const o = state.steer, c = visibleClocks()[0], rec = c?.record;
  const from = `<div class="seg" role="group" aria-label="Gains from">
      <button type="button" data-ctl="from" data-v="time" class="${o.from === "time" ? "on" : ""}" title="A critically damped loop with this time constant: gains by pole placement">Time constant</button>
      <button type="button" data-ctl="from" data-v="cost" class="${o.from === "cost" ? "on" : ""}" title="The gains that minimize a weighted sum of phase², frequency² and steer² (LQG)">Cost</button></div>`;
  const T = `<input type="text" data-ctl="T" value="${h(o.TText)}" placeholder="${rec ? num(steerDefault(c), 4) : "auto"}" inputmode="decimal" aria-label="Loop time constant in seconds" class="carrier">`;
  const integral = `<label class="check" title="A third gain, on the running sum of the phase: it removes the phase offset a frequency drift leaves, at the cost of more noise"><input type="checkbox" data-ctl="integral" ${o.integral ? "checked" : ""}>Integral${inline ? "" : " gain"}</label>`;
  const cost = WEIGHTS.map((w, i) => `<label class="${inline ? "tool-field" : "opt"}"><span>${w}</span><input type="text" data-ctl="cost" data-i="${i}" value="${h(o.cost[i])}" inputmode="decimal" aria-label="Weight on ${w}" class="carrier steer-w"></label>`).join("");
  const every = `<select data-ctl="every" aria-label="Samples between steers">${EVERY.map((k) => `<option value="${k}" ${k === o.every ? "selected" : ""}>${k === 1 ? "sample" : `${k} samples`}</option>`).join("")}</select>`;
  const frac = inline
    ? `<select data-ctl="frac" aria-label="Share of the record the clock model is fitted to">${FRACS.map((f) => `<option value="${f}" ${f === o.frac ? "selected" : ""}>1/${f}</option>`).join("")}</select>`
    : `<div class="pills">${FRACS.map((f) => `<button type="button" class="pill ${f === o.frac ? "on" : ""}" data-ctl="frac" data-v="${f}">1/${f}</button>`).join("")}</div>`;
  const drift = `<label class="check" title="The model carries a frequency drift, and the loop steers it out ahead"><input type="checkbox" data-ctl="drift" ${o.drift ? "checked" : ""}>Drift${inline ? "" : " state"}</label>`;
  const view = `<div class="seg" role="group" aria-label="Plot">
      <button type="button" data-ctl="view" data-v="sigma" class="${o.view === "sigma" ? "on" : ""}">σ(τ)</button>
      <button type="button" data-ctl="view" data-v="record" class="${o.view === "record" ? "on" : ""}">Record</button></div>`;
  const gains = o.from === "time" ? `<label class="tool-field">T ${T} s</label>${integral}` : cost;
  if (inline) return `${from}${gains}<label class="tool-field">Steer&nbsp;every ${every}</label><label class="tool-field">Model&nbsp;from ${frac}</label>${drift}${view}`;
  return `<h2>Steering</h2>
    <div class="ctl-block"><label class="lbl">Gains from</label>${from}</div>
    ${o.from === "time"
      ? `<div class="ctl-block"><label class="lbl">Time constant</label><div class="tau-line">${T}<span class="muted">s</span></div></div>${integral}`
      : `<div class="ctl-block"><label class="lbl">Weights on</label>${cost}</div>`}
    <label class="opt"><span>Steer every</span>${every}</label>
    <div class="ctl-block"><label class="lbl">Clock model from the first</label>${frac}
      ${rec ? `<div class="muted small">${duration(Math.floor(rec.n / o.frac) * rec.tau0)} of the record; the results count from its end</div>` : ""}</div>
    ${drift}
    <div class="ctl-block"><label class="lbl">Plot</label>${view}</div>
    <div class="rule"></div>
    <p class="muted small">Changes here recompute this page right away. The record is the clock running free against its reference; the loop is replayed on it.</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl, click = e.type === "click", change = e.type === "change";
  if ((k === "from" || k === "view") && click) update((s) => { s.steer[k] = t.dataset.v; });
  if (k === "frac" && (click && t.dataset.v || change && t.tagName === "SELECT")) update((s) => { s.steer.frac = Number(t.dataset.v ?? t.value); });
  if (k === "every" && change) update((s) => { s.steer.every = Number(t.value); });
  if ((k === "integral" || k === "drift") && change) update((s) => { s.steer[k] = t.checked; });
  if (k === "T" && change) {
    const text = t.value.trim(), v = Number(text);
    if (!text) update((s) => { s.steer.T = null; s.steer.TText = ""; });
    else if (v > 0) update((s) => { s.steer.T = v; s.steer.TText = text; });
    else { t.value = state.steer.TText; ctx.status("error", "The time constant is in seconds, above 0."); }
  }
  if (k === "cost" && change) {
    const i = Number(t.dataset.i), text = t.value.trim(), v = Number(text);
    if (text && (i === 1 ? v >= 0 : v > 0)) update((s) => { s.steer.cost[i] = text; });
    else { t.value = state.steer.cost[i]; ctx.status("error", "The weights on phase² and steer² are above 0; the weight on frequency² may be 0."); }
  }
}

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="toolbar hold-bar steer-bar" data-ref="toolbar"></div>
      <div class="card plot-card"><div class="card-row"><h2 data-ref="title"></h2><span class="muted small" data-ref="sub"></span></div>
        <div class="plot" data-ref="plot"></div><div class="empty" data-ref="empty" hidden></div></div>
    </div>
    <div class="card side-card steer-card" data-ref="side" hidden></div>`;
  root.addEventListener("click", (e) => { if (e.target.closest("[data-act=load]")) ctx.goData(); });
  // Plotly only follows window resizes; follow the card too (the side card comes and goes).
  const plot = root.querySelector("[data-ref=plot]");
  new ResizeObserver(() => { if (Plotly && plot.layout && plot.clientWidth) Plotly.Plots.resize(plot); }).observe(plot);
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

// Ask the core for any replay that is missing or out of date.
function ensure(clocks) {
  for (const c of clocks) {
    const key = steerKey(c);
    if (c.steer?.key === key || inflight.get(c.id) === key) continue;
    inflight.set(c.id, key);
    const done = (fields) => {
      if (inflight.get(c.id) !== key) return;
      inflight.delete(c.id);
      update(() => { c.steer = { key, ...fields }; });
    };
    fetchSteer(c).then(done, (err) => done({ error: err.message }));
  }
}

function showEmpty(html) {
  const e = root.querySelector("[data-ref=empty]");
  e.innerHTML = html;
  e.hidden = !html;
  root.querySelector("[data-ref=plot]").style.visibility = html ? "hidden" : "visible";
}

// The loop's time constants as shown, each once: a critically damped loop has one.
const timeConstants = (st) => [...new Set(st.time_constants.map(duration))];
const settled = (st) => (st.predicted.error ? null : st.predicted);

function renderSide(box, clocks) {
  box.hidden = !clocks.length;
  if (!clocks.length) return;
  if (clocks.length > 1) {
    box.className = "card side-card hold-cmp";
    box.innerHTML = `<div class="card-head"><h2>Steered clocks</h2></div>
      <div class="tbl"><div class="tbl-row tbl-head"><span>Clock</span><span>Phase RMS</span><span>Steer RMS</span></div>
      ${clocks.map((c) => `<div class="tbl-row tbl-body"><span class="name"><span class="clock-dot" style="background:${seriesColor(c.slot)}"></span>${h(c.name)}</span>
        <span>${sci(c.steer.phase_rms, 3)} s</span><span>${sci(c.steer.steers_rms, 3)}</span></div>`).join("")}</div>`;
    return;
  }
  const c = clocks[0], st = c.steer, g = st.gains, p = settled(st);
  // Poles that are complex beyond rounding: a critically damped loop has a double pole.
  const osc = Math.max(...st.frequencies), rings = 2 * Math.PI * osc * Math.max(...st.time_constants) > 0.05;
  const kv = (name, value, title = "") => `<div class="kv" ${title ? `title="${title}"` : ""}><span>${name}</span><span class="mono">${value}</span></div>`;
  const both = (measured, predicted, unit) => `${sci(measured, 3)}${unit}${p ? ` <span class="muted">(${sci(predicted, 3)})</span>` : ""}`;
  box.className = "card side-card steer-card";
  box.innerHTML = `<div class="card-head"><h2>Steering loop</h2></div>
    <div class="cell-kv">
      ${kv("Phase gain", `${sci(g.phase, 4)} s⁻¹`, "Frequency steer per second of estimated phase")}
      ${kv("Frequency gain", num(g.frequency, 4), "Share of the estimated frequency taken out at each steer")}
      ${g.integral ? kv("Integral gain", `${sci(g.integral, 4)} s⁻¹`, "Frequency steer per second of the running sum of the estimated phase") : ""}
      ${kv(`Time constant${timeConstants(st).length > 1 ? "s" : ""}`, timeConstants(st).join(", "), "How fast the loop pulls an offset in: from its poles")}
      ${rings ? kv("Rings with period", duration(1 / osc), "The poles are complex: the loop overshoots") : ""}
      ${kv("Steer", st.interval > c.record.tau0 ? `every ${duration(st.interval)}` : "every sample")}
      <div class="rule"></div>
      ${kv("Steered phase, RMS", both(st.phase_rms, p?.phase, " s"), "From the end of the model's stretch on; in brackets, the settled loop as the model predicts it")}
      ${kv("Steers, RMS", both(st.steers_rms, p?.steers, ""), "Fractional frequency per steer")}
    </div>
    <p class="muted small hold-note">Model fitted to the first ${duration(st.train)}: ${modelText(st.model)}${st.drift ? ", with a drift state" : ""}.</p>
    ${p ? "" : `<p class="muted small hold-note">No predicted curve: ${h(st.predicted.error)}</p>`}
    <div class="tbl"><div class="tbl-row tbl-head"><span>τ (s)</span><span>free</span><span>steered</span><span>predicted</span></div>
    ${st.tau.map((t, i) => `<div class="tbl-row tbl-body"><span>${num(t, 5)}</span><span>${sci(st.free[i], 3)}</span><span>${sci(st.steered[i], 3)}</span><span>${sci(p?.adev[i], 3)}</span></div>`).join("")}</div>`;
}

// σ(τ) from the end of the model's stretch on.
function sigmaPlot(ready, several) {
  const traces = [];
  for (const c of ready) {
    const st = c.steer, who = several ? `<b>${h(c.name)}</b><br>` : "";
    const line = (name, y, color, dash) => traces.push({
      x: st.tau, y, name: several ? `${c.name}, ${name}` : name, mode: "lines+markers", line: { color, dash, width: 2 }, marker: { color, size: 5 },
      hovertemplate: `${who}${name}<br>τ = %{x:.4g} s<br>σ<sub>y</sub> = %{y:.3e}<extra></extra>`,
    });
    line("Free-running", st.free, several ? seriesColor(c.slot) : css("--s1"), several ? "dash" : "solid");
    line("Steered", st.steered, several ? seriesColor(c.slot) : css("--s2"), "solid");
    const p = settled(st);
    if (!several && p) traces.push({
      x: st.tau, y: p.adev, name: "Predicted (clock model)", mode: "lines", line: { color: css("--muted"), width: 1.5, dash: "dashdot" },
      hovertemplate: "The settled loop, from the clock model<br>τ = %{x:.4g} s<br>σ<sub>y</sub> = %{y:.3e}<extra></extra>",
    });
  }
  // The loop holds the clock to its reference beyond its time constant.
  const T = several ? null : Math.max(...ready[0].steer.time_constants), accent = css("--accent");
  return [traces, {
    ...baseLayout(), margin: { l: 70, r: 18, t: 10, b: 48 },
    xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis("Allan deviation σ<sub>y</sub>(τ)"),
    showlegend: !several, legend: { ...baseLayout().legend, x: 0.02, xanchor: "left", y: 0.02, yanchor: "bottom" },
    shapes: T ? [{ type: "line", xref: "x", yref: "paper", x0: T, x1: T, y0: 0, y1: 1, line: { color: accent, width: 1, dash: "dot" } }] : [],
    annotations: T ? [{ xref: "x", x: Math.log10(T), xanchor: "left", yref: "paper", y: 0.99, yanchor: "top", showarrow: false, text: " time constant", font: { color: accent, size: 11.5 } }] : [],
    uirevision: "steer-sigma",
  }];
}

// The steered phase over the whole record, with the correction in force for one clock.
function recordPlot(ready, several) {
  const all = ready.flatMap((c) => c.steer.series.steered).filter(Number.isFinite).map(Math.abs);
  const [f, u] = timeUnit(Math.max(...all) || 1);
  const traces = ready.map((c) => ({
    x: c.steer.series.t, y: c.steer.series.steered.map((v) => v / f), name: several ? c.name : "Steered phase", mode: "lines",
    line: { color: several ? seriesColor(c.slot) : css("--s2"), width: 1.4 },
    hovertemplate: `${several ? `<b>${h(c.name)}</b><br>` : ""}t = %{x:.6g} s<br>x = %{y:.5g} ${u}<extra></extra>`,
  }));
  const st = ready[0].steer;
  if (!several) traces.push({
    x: st.series.t, y: st.series.correction, name: "Frequency correction", yaxis: "y2", mode: "lines", line: { color: css("--s4"), width: 1.2 },
    hovertemplate: "t = %{x:.6g} s<br>correction = %{y:.4e}<extra></extra>",
  });
  const muted = css("--muted");
  return [traces, {
    ...baseLayout(), margin: { l: 70, r: several ? 18 : 76, t: 10, b: 46 },
    xaxis: linAxis("Time t (s)"), yaxis: { ...linAxis("Steered phase"), ticksuffix: ` ${u}` },
    ...(several ? {} : { yaxis2: { ...linAxis("Frequency correction"), overlaying: "y", side: "right", showgrid: false, exponentformat: "e" } }),
    showlegend: !several, legend: { ...baseLayout().legend, x: 0.98, xanchor: "right", y: 0.02, yanchor: "bottom" },
    // Left of the line the model was fitted; the results count from it on.
    shapes: several ? [] : [{ type: "line", xref: "x", yref: "paper", x0: st.train, x1: st.train, y0: 0, y1: 1, line: { color: muted, width: 1, dash: "dot" } }],
    annotations: several ? [] : [{ xref: "x", x: st.train, xanchor: "left", yref: "paper", y: 0.99, yanchor: "top", showarrow: false, text: " results from here", font: { color: muted, size: 11.5 } }],
    uirevision: "steer-record",
  }];
}

export function refresh() {
  if (!root) return;
  const several = multi(), record = state.steer.view === "record";
  const tb = root.querySelector("[data-ref=toolbar]");
  tb.hidden = !several;
  if (several && !tb.querySelector("input:focus, select:focus")) {   // a field in use keeps its text
    const at = tb.contains(document.activeElement) ? document.activeElement.dataset : null;
    tb.innerHTML = controls(true);
    if (at) tb.querySelector(`[data-ctl=${at.ctl}][data-v="${at.v}"]`)?.focus();
  }
  root.querySelector("[data-ref=title]").textContent = record ? "Steered phase" : "Stability, free-running and steered";
  const side = root.querySelector("[data-ref=side]");
  const clocks = visibleClocks().filter((c) => c.record);
  if (!clocks.length) {
    side.hidden = true;
    return showEmpty(`<h2>No record yet</h2><p>Load or synthesize a record on the Data page.</p><button class="btn" data-act="load">Go to Data</button>`);
  }
  ensure(clocks);
  const current = clocks.filter((c) => c.steer?.key === steerKey(c));
  const ready = current.filter((c) => !c.steer.error);
  const failed = current.filter((c) => c.steer.error);
  renderSide(side, ready);
  if (!ready.length) return showEmpty(`<p>${failed.length ? h(failed.map((c) => `${several ? `${c.name}: ` : ""}${c.steer.error}`).join(" · ")) : "Replaying the steering loop…"}</p>`);
  showEmpty("");
  if (!Plotly) return;
  root.querySelector("[data-ref=sub]").textContent = record ? "the whole record, as the loop would have left it"
    : several ? "solid = steered · dashed = free-running" : "from the end of the stretch the clock model was fitted to";
  Plotly.react(root.querySelector("[data-ref=plot]"), ...(record ? recordPlot : sigmaPlot)(ready, several), CONFIG);
}

export function plotElement() { return root?.querySelector("[data-ref=plot]") ?? null; }
export function hint() {
  if (state.steer.view === "record") return multi() ? "Colour = clock · the steered phase over the whole record" : "Steered phase, with the frequency correction in force on the right axis";
  return multi() ? "Colour = clock · solid = steered · dashed = free-running" : "Free-running and steered σ(τ) · grey = the settled loop as the clock model predicts it";
}

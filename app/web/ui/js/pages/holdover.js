// Holdover page: how far each clock wanders once its reference is gone. The
// core learns the frequency on a training stretch (a least-squares fit, or a
// Kalman filter that reads on up to each start), lets the clock run and
// compares with the record, from every start the record allows. One clock:
// RMS, quantile and worst case against holdover time, the curve its noise
// model predicts or the filter claims, √6·TDEV (√20·HTDEV with drift) for reference, and the time
// to a limit. Several clocks: one colour per clock. Computed when this page
// needs it; no Run is needed.
import { state, update, visibleClocks, multi } from "../state.js";
import { holdover as fetchHoldover, holdoverKey } from "../core.js";
import { loadPlotly, baseLayout, logAxis, CONFIG, seriesColor, css } from "../plot.js";
import { sci, num, thousands, duration, modelText, escapeHtml as h } from "../format.js";

let root, ctx, Plotly;
const inflight = new Map();     // clock id → key being fetched
const FRACS = [8, 4, 2, 0];
const LEVELS = [[0.683, "68.3 %"], [0.9, "90 %"], [0.95, "95 %"], [0.99, "99 %"]];
const NOISE = { WHPM: "white PM", FLPM: "flicker PM", WHFM: "white FM", FLFM: "flicker FM", RWFM: "random-walk FM" };
const pct = (p) => `${num(100 * p, 3)} %`;
// The stretch as the control names it: a share of the record, or (0) as many
// holdover times as the fit has terms.
const stretch = (f, drift) => (f ? `1/${f}` : `${drift ? 3 : 2}τ`);
const reference = (drift) => (drift ? "√20·HTDEV" : "√6·TDEV");

export function controls(inline = false) {
  const o = state.hold;
  const train = inline
    ? `<select data-ctl="frac" aria-label="Training stretch, as a share of the record or of the holdover time">${FRACS.map((f) => `<option value="${f}" ${f === o.frac ? "selected" : ""} ${!f && o.kalman ? "disabled" : ""}>${stretch(f, o.drift)}</option>`).join("")}</select>`
    : `<div class="pills">${FRACS.map((f) => `<button type="button" class="pill ${f === o.frac ? "on" : ""}" data-ctl="frac" data-v="${f}" ${!f && o.kalman ? `disabled title="A Kalman filter trains on a share of the record"` : ""}>${stretch(f, o.drift)}</button>`).join("")}</div>`;
  const by = `<div class="seg" role="group" aria-label="Predicted by">
      <button type="button" data-ctl="kalman" data-v="0" class="${o.kalman ? "" : "on"}" title="A least-squares line or parabola through the training stretch">${inline ? "Fit" : "Least-squares fit"}</button>
      <button type="button" data-ctl="kalman" data-v="1" class="${o.kalman ? "on" : ""}" title="A Kalman filter that has read the record up to each start; its clock model is fitted to the first training stretch alone">${inline ? "Kalman" : "Kalman filter"}</button></div>`;
  const fit = `<div class="seg" role="group" aria-label="Learned from the stretch">
      <button type="button" data-ctl="drift" data-v="0" class="${o.drift ? "" : "on"}">Frequency</button>
      <button type="button" data-ctl="drift" data-v="1" class="${o.drift ? "on" : ""}">${inline ? "+ drift" : "Frequency and drift"}</button></div>`;
  const limit = `<input type="text" data-ctl="limit" value="${h(o.limitText)}" placeholder="none" inputmode="decimal" aria-label="Time error limit in seconds" class="carrier">`;
  const level = `<select data-ctl="confidence" aria-label="Share of trials">${LEVELS.map(([v, l]) => `<option value="${v}" ${v === o.confidence ? "selected" : ""}>${l}</option>`).join("")}</select>`;
  if (inline) {
    return `<label class="tool-field">Train&nbsp;on ${train}</label>${by}${fit}<label class="tool-field">Limit ${limit} s</label><label class="tool-field">Trials ${level}</label>`;
  }
  const c = visibleClocks()[0];
  const n = c?.record?.n ?? 0;
  return `<h2>Holdover</h2>
    <div class="ctl-block"><label class="lbl">Training stretch</label>${train}
      ${n ? `<div class="muted small">${o.frac ? `${duration(Math.floor(n / o.frac) * c.record.tau0)} of the record` : `${o.drift ? "three times" : "twice"} each holdover time, the test behind ${o.drift ? "HTDEV" : "TDEV"}`}, moved along the record one sample at a time</div>` : ""}</div>
    <div class="ctl-block"><label class="lbl">Predicted by</label>${by}</div>
    <div class="ctl-block"><label class="lbl">Learned from the stretch</label>${fit}</div>
    <div class="ctl-block"><label class="lbl">Time error limit</label><div class="tau-line">${limit}<span class="muted">s</span></div></div>
    <label class="opt"><span>Share of trials</span>${level}</label>
    <label class="check" ${o.kalman ? `title="The filter starts from its own estimate"` : ""}><input type="checkbox" data-ctl="anchor" ${o.anchor && !o.kalman ? "checked" : ""} ${o.kalman ? "disabled" : ""}>Start from the last reading</label>
    <div class="rule"></div>
    <p class="muted small">Changes here recompute this page right away. The dashed curve is the error that this share of the trials stayed within; the limit is reached where it crosses.</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl;
  if (k === "frac" && (e.type === "click" && t.dataset.v || e.type === "change" && t.tagName === "SELECT")) update((s) => { s.hold.frac = Number(t.dataset.v ?? t.value); });
  if (k === "kalman" && e.type === "click") update((s) => { s.hold.kalman = t.dataset.v === "1"; if (s.hold.kalman && !s.hold.frac) s.hold.frac = 4; });
  if (k === "drift" && e.type === "click") update((s) => { s.hold.drift = t.dataset.v === "1"; });
  if (k === "anchor" && e.type === "change") update((s) => { s.hold.anchor = t.checked; });
  if (k === "confidence" && e.type === "change") update((s) => { s.hold.confidence = Number(t.value); });
  if (k === "limit" && e.type === "change") {
    const text = t.value.trim(), v = Number(text);
    if (!text) update((s) => { s.hold.limit = null; s.hold.limitText = ""; });
    else if (v > 0) update((s) => { s.hold.limit = v; s.hold.limitText = text; });
    else { t.value = state.hold.limitText; ctx.status("error", "The limit is a time error in seconds above 0, e.g. 1e-6 for a microsecond."); }
  }
}

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="toolbar hold-bar" data-ref="toolbar"></div>
      <div class="card plot-card"><div class="card-row"><h2>Time error in holdover</h2><span class="muted small" data-ref="sub"></span></div>
        <div class="plot" data-ref="plot"></div><div class="empty" data-ref="empty" hidden></div></div>
    </div>
    <div class="card side-card hold-card" data-ref="side" hidden></div>`;
  root.addEventListener("click", (e) => { if (e.target.closest("[data-act=load]")) ctx.goData(); });
  // Plotly only follows window resizes; follow the card too (the side card comes and goes).
  const plot = root.querySelector("[data-ref=plot]");
  new ResizeObserver(() => { if (Plotly && plot.layout && plot.clientWidth) Plotly.Plots.resize(plot); }).observe(plot);
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

// Ask the core for any test that is missing or out of date.
function ensure(clocks) {
  for (const c of clocks) {
    const key = holdoverKey(c);
    if (c.holdover?.key === key || inflight.get(c.id) === key) continue;
    inflight.set(c.id, key);
    const done = (fields) => {
      if (inflight.get(c.id) !== key) return;
      inflight.delete(c.id);
      update(() => { c.holdover = { key, ...fields }; });
    };
    fetchHoldover(c).then((doc) => done(doc.holdover), (err) => done({ error: err.message }));
  }
}

// "within 1e-06 s for 2.3 h", or why there is no such time.
function limitText(hd) {
  if (hd.limit === null) return "–";
  if (hd.time_to_limit === null) return `longer than ${duration(hd.tau[hd.tau.length - 1])}`;
  return hd.time_to_limit === 0 ? `under ${duration(hd.tau[0])}` : duration(hd.time_to_limit);
}

function showEmpty(html) {
  const e = root.querySelector("[data-ref=empty]");
  e.innerHTML = html;
  e.hidden = !html;
  root.querySelector("[data-ref=plot]").style.visibility = html ? "hidden" : "visible";
}

// The RMS over the reference at the shortest, a middle and the longest τ that has one.
function ratios(hd) {
  const at = hd.tau.map((_, i) => i).filter((i) => hd.reference[i] > 0);
  if (!at.length) return "";
  const pick = [...new Set([at[0], at[at.length >> 1], at[at.length - 1]])];
  const name = reference(hd.drift);
  return `<p class="muted small hold-note" title="${name} is the RMS error of extrapolating the next τ-average of the phase from the ${hd.drift ? "three" : "two"} before it. With a stretch of ${hd.drift ? 3 : 2}τ the test equals it at τ₀ and stays a constant 1.2 to 1.4 above on FM noise; white PM parts them.">RMS over ${name}: ${pick.map((i) => `${(hd.rms[i] / hd.reference[i]).toFixed(2)} at ${duration(hd.tau[i])}`).join(", ")}.</p>`;
}

function renderSide(box, clocks) {
  box.hidden = !clocks.length;
  if (!clocks.length) return;
  if (clocks.length > 1) {
    box.className = "card side-card hold-card hold-cmp";
    box.innerHTML = `<div class="card-head"><h2>Time within the limit</h2></div>
      <div class="tbl"><div class="tbl-row tbl-head"><span>Clock</span><span>${pct(state.hold.confidence)} of trials</span><span>RMS at longest τ</span></div>
      ${clocks.map((c) => `<div class="tbl-row tbl-body"><span class="name"><span class="clock-dot" style="background:${seriesColor(c.slot)}"></span>${h(c.name)}</span>
        <span>${limitText(c.holdover)}</span><span>${sci(c.holdover.rms[c.holdover.rms.length - 1], 3)} s</span></div>`).join("")}</div>
      ${state.hold.limit ? "" : `<p class="muted small hold-note">Set a time error limit to compare how long each clock stays within it.</p>`}`;
    return;
  }
  const hd = clocks[0].holdover, p = hd.predicted, last = hd.tau.length - 1;
  const model = p && !p.error;
  const shares = model ? p.alphas.map((_, j) => [p.noise_types[j], p.budget[j][last] / (p.sigma[last] ** 2)]).filter(([, s]) => s >= 0.005) : [];
  box.className = "card side-card hold-card";
  box.innerHTML = `<div class="card-head"><h2>Holdover</h2></div>
    <div class="cell-kv">
      <div class="kv"><span>Training stretch</span><span class="mono">${hd.grow ? `${hd.grow}τ` : duration(hd.train)}</span></div>
      <div class="kv" title="Holdover time at which the error that ${pct(hd.confidence)} of the trials stayed within reaches the limit"><span>Within ${hd.limit === null ? "the limit" : `${sci(hd.limit, 3)} s`} for</span><span class="mono">${limitText(hd)}</span></div>
      ${model && hd.limit !== null ? `<div class="kv" title="The same time from the power-law noise model fitted to the record"><span>Predicted from the noise model</span><span class="mono">${limitText({ ...hd, time_to_limit: p.time_to_limit })}</span></div>` : ""}
    </div>
    ${shares.length ? `<p class="muted small hold-note" title="Share of the predicted error variance at the longest holdover time">Predicted error at ${duration(hd.tau[last])}: ${shares.map(([n, s]) => `${NOISE[n] ?? n} ${Math.round(100 * s)} %`).join(", ")}.</p>` : ""}
    ${ratios(hd)}
    <div class="tbl"><div class="tbl-row tbl-head"><span>τ (s)</span><span>RMS (s)</span><span>${pct(hd.confidence)}</span><span>worst</span><span>trials</span></div>
    ${hd.tau.map((t, i) => `<div class="tbl-row tbl-body"><span>${num(t, 5)}</span><span>${sci(hd.rms[i], 3)}</span><span>${sci(hd.quantile[i], 3)}</span><span>${sci(hd.worst[i], 3)}</span><span>${thousands(hd.n[i])}</span></div>`).join("")}</div>
    ${p?.error ? `<p class="muted small hold-note">No predicted curve: ${h(p.error)}</p>` : ""}
    ${hd.claimed ? `<p class="muted small hold-note" title="Above 1, the clock did worse than the filter's model said it would">RMS over the filter's own claim: ${(hd.rms[0] / hd.claimed[0]).toFixed(2)} at ${duration(hd.tau[0])}, ${(hd.rms[last] / hd.claimed[last]).toFixed(2)} at ${duration(hd.tau[last])}. Model fitted to the first ${duration(hd.train)}: ${modelText(hd.model)}.</p>` : ""}`;
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
  const side = root.querySelector("[data-ref=side]");
  const clocks = visibleClocks().filter((c) => c.record);
  if (!clocks.length) {
    side.hidden = true;
    return showEmpty(`<h2>No record yet</h2><p>Load or synthesize a record on the Data page.</p><button class="btn" data-act="load">Go to Data</button>`);
  }
  ensure(clocks);
  const current = clocks.filter((c) => c.holdover?.key === holdoverKey(c));
  const ready = current.filter((c) => !c.holdover.error);
  const failed = current.filter((c) => c.holdover.error);
  renderSide(side, ready);
  if (!ready.length) return showEmpty(`<p>${failed.length ? h(failed.map((c) => `${several ? `${c.name}: ` : ""}${c.holdover.error}`).join(" · ")) : "Running the holdover test…"}</p>`);
  showEmpty("");
  if (!Plotly) return;

  const level = pct(state.hold.confidence);
  const traces = [];
  for (const c of ready) {
    const hd = c.holdover, who = several ? `<b>${h(c.name)}</b><br>` : "";
    const line = (name, y, color, dash, width = 2) => traces.push({
      x: hd.tau, y, name: several ? `${c.name}, ${name}` : name, mode: "lines+markers", line: { color, dash, width }, marker: { color, size: 5 },
      hovertemplate: `${who}${name}<br>τ = %{x:.4g} s<br>%{y:.3e} s<extra></extra>`,
    });
    const col = (k) => (several ? seriesColor(c.slot) : css(`--s${k}`));
    line("RMS", hd.rms, col(1), "solid");
    line(`${level} of trials`, hd.quantile, col(2), "dash");
    if (!several) {
      line("Worst trial", hd.worst, col(3), "dot");
      const p = hd.predicted, grey = { mode: "lines", line: { color: css("--muted"), width: 1.5, dash: "dashdot" } };
      if (p && !p.error) traces.push({
        x: hd.tau, y: p.sigma, name: "Predicted RMS (noise model)", ...grey,
        hovertemplate: "Predicted from the noise model (h<sub>α</sub> fitted to MHDEV)<br>τ = %{x:.4g} s<br>%{y:.3e} s<extra></extra>",
      });
      if (hd.claimed) traces.push({
        x: hd.tau, y: hd.claimed, name: "Claimed by the filter", ...grey,
        hovertemplate: "The filter's own standard deviation<br>τ = %{x:.4g} s<br>%{y:.3e} s<extra></extra>",
      });
      traces.push({
        x: hd.tau, y: hd.reference, name: reference(hd.drift), mode: "lines", line: { color: css("--s4"), width: 1.5 },
        hovertemplate: `${reference(hd.drift)}: the same prediction from τ-averages<br>τ = %{x:.4g} s<br>%{y:.3e} s<extra></extra>`,
      });
    }
  }
  const limit = state.hold.limit;
  const accent = css("--accent");
  const shapes = limit ? [{ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: limit, y1: limit, line: { color: accent, width: 1.5 } }] : [];
  if (limit) for (const c of ready) {
    const t = c.holdover.time_to_limit;
    if (t > 0) shapes.push({ type: "line", xref: "x", yref: "paper", x0: t, x1: t, y0: 0, y1: 1, line: { color: several ? seriesColor(c.slot) : accent, width: 1, dash: "dot" } });
  }
  root.querySelector("[data-ref=sub]").textContent = several ? `solid = RMS · dashed = ${level} of trials` : state.hold.kalman ? "from the estimate of a Kalman filter that has read up to each start" : `from the last ${state.hold.anchor ? "reading" : "value of the fit"} of each training stretch`;
  Plotly.react(root.querySelector("[data-ref=plot]"), traces, {
    ...baseLayout(), margin: { l: 70, r: 18, t: 10, b: 48 },
    xaxis: logAxis("Holdover time τ (s)"), yaxis: logAxis("Time error (s)"),
    showlegend: !several, legend: { ...baseLayout().legend, x: 0.98, xanchor: "right", y: 0.02, yanchor: "bottom" },
    shapes, annotations: limit ? [{ xref: "paper", x: 0.01, xanchor: "left", yref: "y", y: Math.log10(limit), yanchor: "bottom", showarrow: false, text: `limit ${sci(limit, 3)} s`, font: { color: accent, size: 11.5 } }] : [],
    uirevision: "hold",
  }, CONFIG);
}

export function plotElement() { return root?.querySelector("[data-ref=plot]") ?? null; }
export function hint() {
  return multi() ? "Colour = clock · solid = RMS · dashed = the chosen share of trials" : `Solid = RMS · dashed = the chosen share of trials · dotted = worst trial · grey = predicted from the noise model, or claimed by the filter · thin = ${reference(state.hold.drift)}`;
}

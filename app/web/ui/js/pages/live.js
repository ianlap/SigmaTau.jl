// Live page: readings from an instrument (or the simulator) streamed through
// `sigmatau live`, with σ(τ) and the fractional frequency updating as they
// arrive. The settings live in the side panel. The readings are kept here,
// not in the app state, so polling redraws only this page; "Analyse" hands
// them to the other pages as a new clock, with confidence intervals and
// every estimator.
import { state, update } from "../state.js";
import { liveStart, livePoll, liveStop } from "../core.js";
import { loadPlotly, baseLayout, logAxis, linAxis, CONFIG, seriesColor } from "../plot.js";
import { num, sci, thousands, duration, escapeHtml as h } from "../format.js";

let root, ctx, Plotly;
const DEVS = [["adev", "ADEV"], ["mdev", "MDEV"], ["hdev", "HDEV"], ["mhdev", "MHDEV"]];
const SOURCES = [["sim", "Simulator"], ["tcp", "Network"], ["device", "Serial / file"]];
const READINGS = [["hz", "Counter frequency (Hz)"], ["frequency", "Fractional frequency"], ["phase", "Phase (s)"]];
const STRIP = 5000;          // points kept on the frequency strip
const POLL_MS = 250;

// The measurement in progress (or the last one).
const run = { id: null, since: 0, running: false, error: null, readings: [], last: null, quantity: "phase", tau0: 1, started: 0, timer: null };

export const ownsPanel = true;     // the settings stay in the panel with several clocks too

export function controls() {
  const L = state.live;
  const busy = run.running;
  const dis = busy ? "disabled" : "";
  const seg = (k, items, v) => `<div class="seg" role="group">${items.map(([id, label]) =>
    `<button type="button" data-ctl="${k}" data-v="${id}" class="${v === id ? "on" : ""}" ${dis}>${label}</button>`).join("")}</div>`;
  const field = (k, label, v, hint = "") => `<label class="opt live-field"><span>${label}</span><input type="text" data-ctl="${k}" value="${h(v)}" ${dis} spellcheck="false" ${hint ? `placeholder="${h(hint)}"` : ""}></label>`;
  const where = L.source === "tcp" ? `${field("address", "Address", L.address, "host:port")}${field("query", "Query", L.query, "none")}
      <p class="muted small">An LXI instrument's SCPI socket, usually port 5025. The query is sent before every reading.</p>`
    : L.source === "device" ? `${field("device", "Path", L.device)}${field("query", "Query", L.source === "device" && L.query === "READ?" ? "" : L.query, "none")}
      <p class="muted small">A serial port set up with <code>stty</code> first, or a file being written. One reading per line.</p>`
    : `<p class="muted small">A simulated clock: white PM, then white FM, then random-walk FM, one reading per τ₀.</p>`;
  return `<h2>Live measurement</h2>
    <div class="ctl-block"><label class="lbl">Source</label>${seg("source", SOURCES, L.source)}</div>
    ${where}
    ${L.source === "sim" ? "" : `<label class="opt"><span>Readings</span><select class="wide" data-ctl="reading" ${dis}>${READINGS.map(([id, label]) =>
      `<option value="${id}" ${L.reading === id ? "selected" : ""}>${label}</option>`).join("")}</select></label>
      ${L.reading === "hz" ? field("nominal", "Nominal (Hz)", L.nominal, "first reading") : ""}`}
    ${field("tau0", "τ₀ (s)", String(L.tau0))}
    <div class="ctl-block"><label class="lbl">Estimators</label><div class="live-devs">${DEVS.map(([id, label]) =>
      `<label class="check"><input type="checkbox" data-ctl="dev" data-v="${id}" ${L.devs.includes(id) ? "checked" : ""} ${dis}>${label}</label>`).join("")}</div></div>
    <button class="btn ${busy ? "" : "btn-primary"} live-go" type="button" data-ctl="go">${busy ? "Stop" : "Start"}</button>
    <div class="rule"></div>
    <p class="muted small">σ(τ) updates with every reading, without intervals. Analyse the readings afterwards for confidence intervals, noise identification and every other estimator.</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl;
  if (k === "go" && e.type === "click") return run.running ? stop() : start();
  if (run.running) return;
  if (k === "source" && e.type === "click") update((s) => { s.live.source = t.dataset.v; });
  if (k === "dev" && e.type === "change") update((s) => {
    s.live.devs = DEVS.map(([id]) => id).filter((id) => (id === t.dataset.v ? t.checked : s.live.devs.includes(id)));
  });
  if (e.type !== "change") return;
  if (k === "reading") update((s) => { s.live.reading = t.value; });
  if (["address", "device", "query", "nominal"].includes(k)) update((s) => { s.live[k] = t.value.trim(); });
  if (k === "tau0") {
    const v = Number(t.value);
    if (v > 0) update((s) => { s.live.tau0 = v; });
    else { t.value = state.live.tau0; ctx.status("error", "τ₀ must be a time above 0 s, e.g. 1 or 0.1."); }
  }
}

function args() {
  const L = state.live;
  const a = ["live"];
  if (L.source === "sim") a.push("sim");
  else if (L.source === "tcp") a.push(`tcp://${L.address}`, "--query", L.query);
  else a.push(L.device, ...(L.query ? ["--query", L.query] : []));
  if (L.source !== "sim") {
    if (L.reading === "hz") a.push("--hz", ...(L.nominal ? ["--nominal", L.nominal] : []));
    if (L.reading === "frequency") a.push("--freq");
  }
  // About four updates a second, whatever the reading rate.
  a.push("--tau0", String(L.tau0), "--devs", L.devs.join(","), "--every", String(Math.max(1, Math.round(0.25 / L.tau0))));
  return a;
}

async function start() {
  if (!state.live.devs.length) return ctx.status("warn", "Pick at least one estimator to stream.");
  if (state.live.source === "device" && !state.live.device) return ctx.status("warn", "Give the device or file path.");
  Object.assign(run, { since: 0, error: null, readings: [], last: null, started: Date.now(), running: true });
  update(() => {});
  try {
    run.id = await liveStart(args());
  } catch (err) {
    Object.assign(run, { running: false, error: err.message });
    update(() => {});
    return;
  }
  poll();
}

async function stop() {
  clearTimeout(run.timer);
  await liveStop().catch(() => {});
  run.running = false;
  update(() => {});
}

async function poll() {
  const id = run.id;
  let out;
  try {
    out = await livePoll(id, run.since);
  } catch (err) {
    out = { lines: [], running: false, error: err.message };
  }
  if (id !== run.id) return;
  for (const line of out.lines) {
    run.readings.push(...line.readings);
    Object.assign(run, { last: line, quantity: line.quantity, tau0: line.tau0 });
  }
  run.since = out.next ?? run.since;
  if (out.lines.length && state.page === "live") draw();
  if (out.running && run.running) {
    run.timer = setTimeout(poll, POLL_MS);
  } else if (run.running) {
    Object.assign(run, { running: false, error: out.error });
    update(() => {});
  }
}

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="card live-bar" data-ref="bar"></div>
      <div class="card plot-card live-card"><div class="card-row"><h2>σ(τ) so far</h2><span class="muted small" data-ref="sub"></span></div>
        <div class="plot" data-ref="sigma"></div><div class="empty" data-ref="empty"></div></div>
      <div class="card plot-card live-card live-strip"><div class="card-row"><h2>Fractional frequency</h2><span class="muted small" data-ref="strip-sub"></span></div>
        <div class="plot" data-ref="strip"></div></div>
    </div>`;
  root.addEventListener("click", (e) => {
    const act = e.target.closest("[data-act]")?.dataset.act;
    if (act === "analyse") analyse();
    if (act === "save") save();
  });
  loadPlotly().then((P) => { Plotly = P; draw(); }).catch(() => {});
}

export function refresh() {
  if (root) draw();
}

// Fractional frequency of each reading: frequency as is, counter readings
// against their nominal, phase as first differences.
function frequencies() {
  const r = run.readings, q = run.quantity;
  if (q === "frequency") return r;
  if (q === "hz") { const f0 = run.last?.nominal; return f0 ? r.map((f) => f / f0 - 1) : []; }
  return r.slice(1).map((x, i) => (x - r[i]) / run.tau0);
}

function draw() {
  const bar = root.querySelector("[data-ref=bar]");
  const n = run.readings.length;
  const state_ = run.running ? `<span class="live-dot on"></span>Measuring` : run.error ? `<span class="live-dot err"></span>Stopped` : n ? "Stopped" : "Not started";
  bar.innerHTML = `<div class="live-stats">
      <span class="live-state">${state_}</span>
      <span><b>${thousands(n)}</b> readings</span>
      <span>${n ? duration(n * run.tau0) : "–"} of data</span>
      <span>last ${n ? h(num(run.readings[n - 1], 12)) : "–"}</span>
    </div>
    <span class="grow"></span>
    <button class="btn" type="button" data-act="save" ${n ? "" : "disabled"}>Save readings</button>
    <button class="btn btn-primary" type="button" data-act="analyse" ${n >= 16 && !run.running ? "" : "disabled"} title="Load the readings as a clock and run the full analysis">Analyse readings</button>
    ${run.error ? `<div class="live-error">${h(run.error)}</div>` : ""}`;
  const empty = root.querySelector("[data-ref=empty]");
  const sigma = root.querySelector("[data-ref=sigma]");
  const ready = run.last && run.last.results.some((r) => r.dev.some((v) => v !== null));
  empty.hidden = !!ready;
  empty.innerHTML = ready ? "" : `<p>${run.running ? "Waiting for enough readings…" : "Set up the source in the panel and press Start."}</p>`;
  sigma.style.visibility = ready ? "visible" : "hidden";
  if (!Plotly || !ready) return;
  const L = run.last;
  root.querySelector("[data-ref=sub]").textContent = `${thousands(L.n)} readings · τ₀ = ${num(L.tau0, 4)} s`;
  Plotly.react(sigma, L.results.map((r, i) => ({
    x: L.tau, y: r.dev, name: r.deviation.toUpperCase(), mode: "lines+markers",
    line: { color: seriesColor(i), width: 2 }, marker: { color: seriesColor(i), size: 6 },
    customdata: r.neff, hovertemplate: `<b>${r.deviation.toUpperCase()}</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}<br>%{customdata} terms<extra></extra>`,
  })), {
    ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis("σ(τ)"),
    showlegend: true, uirevision: "live",
  }, CONFIG);
  const y = frequencies();
  const k0 = Math.max(0, y.length - STRIP);
  const off = n - y.length;
  root.querySelector("[data-ref=strip-sub]").textContent = y.length > STRIP ? `last ${thousands(STRIP)} readings` : "";
  Plotly.react(root.querySelector("[data-ref=strip]"), [{
    x: y.slice(k0).map((_, i) => (k0 + i + off) * L.tau0), y: y.slice(k0), mode: "lines",
    line: { color: seriesColor(0), width: 1 }, hovertemplate: "t = %{x:.6g} s<br>y = %{y:.4e}<extra></extra>",
  }], {
    ...baseLayout(), margin: { l: 70, r: 14, t: 8, b: 44 },
    xaxis: linAxis("Time (s)"), yaxis: { ...linAxis("y"), exponentformat: "e" }, uirevision: "strip",
  }, CONFIG);
}

// The readings as a one-column record: phase or fractional frequency.
function record() {
  const phase = run.quantity === "phase";
  return { phase, text: (phase ? run.readings : frequencies()).map((v) => String(v)).join("\n") + "\n" };
}

function stamp() {
  return new Date(run.started).toISOString().slice(0, 19).replace(/[-:]/g, "").replace("T", "-");
}

async function analyse() {
  const { phase, text } = record();
  await ctx.addRecord({
    name: `Live ${new Date(run.started).toLocaleTimeString()}`,
    source: { type: "file", name: `live-${stamp()}.txt`, content: text },
    sourceInfo: `live: ${run.last?.source ?? ""}`,
    kind: phase ? "phase" : "frequency", tau0: run.tau0, time_col: 0, value_col: 1, hasTimeCol: false, missing: 0,
  });
}

function save() {
  const { phase, text } = record();
  ctx.save(`live-${stamp()}.txt`, new Blob([`# ${phase ? "phase (s)" : "fractional frequency"}, tau0 = ${run.tau0} s, source ${run.last?.source ?? ""}\n${text}`], { type: "text/plain" }));
}

export function plotElement() {
  return root?.querySelector("[data-ref=sigma]") ?? null;
}
export function hint() {
  return run.running ? "Stop the measurement to analyse the readings" : "Start reads from the source in the panel";
}

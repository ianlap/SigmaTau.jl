// Entry point: routing, loading records, the Run action, the header chip,
// status bar, export menu, drag and drop and keyboard shortcuts.
import { state, update, subscribe, newClock, activeClock, multi, referenceClock, isReference } from "./state.js";
import { health, load, analyze } from "./core.js";
import { initPanel, renderPanel } from "./panel.js";
import { openImport } from "./importer.js";
import * as stability from "./pages/stability.js";
import * as data from "./pages/data.js";
import * as spectrum from "./pages/spectrum.js";
import * as dynamicPage from "./pages/dynamic.js";
import * as holdover from "./pages/holdover.js";
import * as live from "./pages/live.js";
import * as ensemble from "./pages/ensemble.js";
import * as synthesis from "./pages/synthesis.js";
import * as help from "./pages/help.js";
import { label, gridName } from "./estimators.js";
import { openFile, saveFile, readFile } from "./platform.js";
import { buildReport } from "./report.js";
import { sessionText, parseSession, applySession, isSession } from "./session.js";
import { loadPlotly, css } from "./plot.js";
import { thousands, num, duration, escapeHtml as h } from "./format.js";

const PAGES = { data, synthesis, stability, spectrum, dynamic: dynamicPage, holdover, ensemble, live, help };
const ORDER = ["data", "synthesis", "stability", "spectrum", "dynamic", "holdover", "ensemble", "live"];
let pageEl = document.getElementById("page");
let page = null;
let transient = null;      // { kind, text } shown in the status bar until the next change

const ctx = {
  run, addClock, openFile: openInto, reimport, useSynth, setRecord, addRecord, save: saveFile,
  goData: () => go("data"),
  go,
  status: (kind, text) => { transient = { kind, text, until: kind === "error" ? Infinity : Date.now() + 8000 }; renderStatus(); },
};

// ── Routing ──
function go(name) {
  if (!PAGES[name]) name = "stability";
  if (location.hash !== `#${name}`) history.replaceState(null, "", `#${name}`);
  if (page === PAGES[name]) return;
  state.page = name;
  page = PAGES[name];
  // Each visit mounts on a fresh element, so the listeners of earlier mounts go with the old one.
  const fresh = pageEl.cloneNode(false);
  pageEl.replaceWith(fresh);
  pageEl = fresh;
  pageEl.className = `page page-${name}`;
  page.mount(pageEl, ctx);
  render();
}
window.addEventListener("hashchange", () => go(location.hash.slice(1)));

// ── Records ──
const controllers = new Map();

async function loadClock(c) {
  controllers.get(c.id)?.abort();
  const ctl = new AbortController();
  controllers.set(c.id, ctl);
  update((s) => { c.loading = true; c.error = null; markRecordStale(s); });
  try {
    const out = await load(c, ctl.signal);
    update(() => { c.record = out.record; c.loading = false; });
  } catch (err) {
    if (err.name === "AbortError") return;
    update(() => { c.record = null; c.result = null; c.loading = false; c.error = err.message; });
  }
}

function markRecordStale(s) { s.run.stale = "record"; s.run.gen = (s.run.gen ?? 0) + 1; }

const DRIFT_FACT = { none: "drift kept", linear: "offset removed", quadratic: "drift removed" };
const baseName = (n) => n.replace(/\.[^.]+$/, "");

// A frequency record cannot be read with samples missing, so a file with gaps is filled from the start.
const fillNeeded = (settings) => settings.kind === "frequency" && (settings.gaps > 0 || settings.blanks > 0);

let nextGroup = 1;

// Load a file through the import dialog. The first value column goes into
// `clock` (or a new clock); each further value column becomes its own clock,
// all in one group so a reference clock can be chosen among them.
async function importFile(file, clock, prev) {
  const settings = await openImport(file, prev);
  if (!settings) return;
  const group = prev?.group ?? nextGroup++;
  const touched = [];
  update((s) => {
    settings.value_cols.forEach((col, k) => {
      let c = k === 0 ? clock : null;
      const given = settings.labels?.[k] ?? null;     // a name typed in the dialog
      const colName = given ?? settings.names[k] ?? (settings.value_cols.length > 1 ? `${baseName(file.name)} col ${col}` : baseName(file.name));
      if (!c) {
        c = newClock({ name: colName });
        s.clocks.push(c);
      } else if (given) {
        c.name = given;
      } else if (!multi() || c.source.type !== "file" || c.name === baseName(c.source.name ?? "")) {
        c.name = settings.value_cols.length > 1 ? colName : baseName(file.name);
      }
      Object.assign(c, {
        source: { type: "file", name: file.name, content: file.text, sha256: file.sha256 ?? "" }, group,
        kind: settings.kind, format: settings.format, tau0: settings.tau0, time_col: settings.time_col, timeUnit: settings.timeUnit,
        value_col: col, column_name: settings.names[k] ?? null, against: settings.against, header: settings.header, scaling: settings.scaling,
        sourceInfo: settings.sourceInfo,
        keep_zeros: settings.keep_zeros, fillgaps: fillNeeded(settings), remove_outliers: false, range: null, result: null,
      });
      touched.push(c);
    });
    s.activeClock = touched[0].id;
    s.dataMode = "file";
  });
  go("data");
  if (settings.remembered) ctx.status("good", `Read ${file.name} with the settings remembered for files like this. Import settings… changes them.`);
  await Promise.all(touched.map(loadClock));
}

// Choose the reference clock (null = as measured). Clocks from the same file
// are re-read against it; the reference itself leaves the plots.
async function setReference(id) {
  const old = referenceClock();
  update((s) => { s.reference = id; markRecordStale(s); });
  const r = referenceClock();
  const groups = new Set([old?.group, r?.group].filter((g) => g !== null && g !== undefined));
  const hadResults = state.clocks.some((c) => c.result);
  await Promise.all(state.clocks.filter((c) => groups.has(c.group)).map((c) => { c.result = null; return loadClock(c); }));
  if (hadResults) await run();
}

async function openInto(clock) {
  const f = await openFile();
  if (f) await importFile(f, clock ?? null);
}

async function addClock() {
  const f = await openFile();
  if (f) await importFile(f, null);
}

// Reopen the import dialog for a file's clocks. Columns dropped there drop
// their clocks; new ones are added.
async function reimport(c) {
  if (c?.source.type !== "file") return;
  const group = state.clocks.filter((k) => k.group === c.group);
  const prev = { ...c, value_cols: group.map((k) => k.value_col), clockNames: Object.fromEntries(group.map((k) => [k.value_col, k.name])) };
  const settings = await openImport({ name: c.source.name, text: c.source.content }, prev);
  if (!settings) return;
  const keep = new Map(group.map((k) => [k.value_col, k]));
  update((s) => {
    s.clocks = s.clocks.filter((k) => k.group !== c.group || settings.value_cols.includes(k.value_col));
    if (s.reference !== null && !s.clocks.some((k) => k.id === s.reference)) s.reference = null;
  });
  const touched = [];
  update((s) => {
    settings.value_cols.forEach((col, i) => {
      let k = keep.get(col);
      const given = settings.labels?.[i] ?? null;
      if (!k) { k = newClock({ name: given ?? settings.names[i] ?? `${baseName(c.source.name)} col ${col}` }); s.clocks.push(k); }
      else if (given) k.name = given;
      Object.assign(k, {
        group: c.group, kind: settings.kind, format: settings.format, tau0: settings.tau0, time_col: settings.time_col, timeUnit: settings.timeUnit,
        value_col: col, column_name: settings.names[i] ?? null, against: settings.against, header: settings.header, scaling: settings.scaling,
        sourceInfo: settings.sourceInfo, source: c.source,
        keep_zeros: settings.keep_zeros, fillgaps: k.fillgaps || fillNeeded(settings), result: null,
      });
      touched.push(k);
    });
    if (!s.clocks.some((k) => k.id === s.activeClock)) s.activeClock = touched[0]?.id ?? null;
  });
  await Promise.all(touched.map(loadClock));
}

// Simulated clocks that belong together share the numbers they were made from.
const ensembleKey = (src) => (src.type === "synth" && src.clocks > 1 ? JSON.stringify([src.n, src.seed, src.tau0, src.clocks, src.models ?? null]) : null);
const DEFAULT_NAMES = /^(Clock|Synthetic record|Sim \d+)$/;

// Make `clocks` synthetic clocks (one, or an ensemble of simulated clocks
// against a perfect reference), reusing the clock you are on and the other
// members of its ensemble.
async function useSynth(c, { n, tau0, seed, clocks = 1, models = null, names = null }) {
  const k = Math.max(1, Math.round(clocks));
  const made = [];
  update((s) => {
    const key = c ? ensembleKey(c.source) : null;
    const sibs = key ? s.clocks.filter((x) => ensembleKey(x.source) === key).sort((a, b) => a.source.index - b.source.index) : [];
    const pool = sibs.length ? sibs : c ? [c] : [];
    s.clocks = s.clocks.filter((x) => !pool.slice(k).includes(x));        // fewer clocks than before
    for (let i = 0; i < k; i++) {
      let x = pool[i];
      const wasFile = x?.source.type === "file";
      if (!x) { x = newClock(); s.clocks.push(x); }
      Object.assign(x, {
        source: { type: "synth", n, seed, tau0, ...(models ? { models } : {}), ...(k > 1 ? { clocks: k, index: i + 1 } : {}) },
        sourceInfo: k > 1 ? `simulated clock ${i + 1} of ${k}` : models ? "simulated clock" : "white + random-walk FM", kind: "phase", group: null,
        tau0, fillgaps: false, remove_outliers: false, drift: "none", range: null, result: null,
      });
      if (names) x.name = names[i];
      else if (wasFile || DEFAULT_NAMES.test(x.name)) x.name = k > 1 ? `Sim ${i + 1}` : "Synthetic record";
      if (s.reference === x.id) s.reference = null;
      made.push(x);
    }
    s.activeClock = made[0].id;
  });
  await Promise.all(made.map(loadClock));
}

async function setRecord(c, fields) {
  if (!c) return;
  update(() => Object.assign(c, fields));
  await loadClock(c);
}

// A record made on another page (the Live page's readings): a new clock,
// analysed straight away.
async function addRecord(fields) {
  let c;
  update((s) => { c = newClock(fields); s.clocks.push(c); s.activeClock = c.id; s.dataMode = "file"; });
  go("stability");
  await loadClock(c);
  if (c.record) await run();
}

// Replace everything with a saved session, read its records and rerun.
async function openSession(file) {
  let doc, clocks;
  try { doc = parseSession(file.text); } catch (err) { return ctx.status("error", `${file.name}: ${err.message}`); }
  for (const c of state.clocks) controllers.get(c.id)?.abort();
  runCtl?.abort();     // a run in flight belongs to the clocks being replaced
  update((s) => { clocks = applySession(s, doc); markRecordStale(s); s.run.state = "idle"; });
  nextGroup = Math.max(nextGroup, ...clocks.map((c) => (c.group ?? 0) + 1));
  go("stability");
  await Promise.all(clocks.map(loadClock));
  await run();
  const failed = clocks.find((c) => c.error);
  if (failed) ctx.status("error", `Opened ${file.name}, but ${failed.name} failed: ${failed.error}`);
  else ctx.status("good", `Opened the session ${file.name}.`);
}

// Take a clock out of the explorer. Clocks simulated together share one random
// stream, so the rest are simulated again without it, each with its own model.
async function removeClock(id) {
  const c = state.clocks.find((k) => k.id === id);
  if (!c) return;
  const key = ensembleKey(c.source);
  const sibs = key ? state.clocks.filter((k) => ensembleKey(k.source) === key).sort((a, b) => a.source.index - b.source.index) : [];
  if (sibs.length > 1) {
    const { n, tau0, seed } = c.source;
    const models = c.source.models ?? sibs.map((_, j) => [0, 1e-24 * 4 ** j, 1e-30 * 4 ** j, 0]);   // the core's default clocks
    const rest = (xs) => xs.filter((_, j) => j !== c.source.index - 1);
    return useSynth(sibs[0], { n, tau0, seed, clocks: sibs.length - 1, models: rest(models), names: rest(sibs).map((k) => k.name) });
  }
  if (state.reference === id) await setReference(null);
  update((s) => {
    s.clocks = s.clocks.filter((k) => k.id !== id);
    if (!s.clocks.some((k) => k.id === s.activeClock)) s.activeClock = s.clocks[0]?.id ?? null;
  });
}

function selectClock(id) {
  update((s) => { s.activeClock = id; });
  go("data");
}

// ── Run ──
let runCtl = null;
let runGen = 0;      // the input generation the run in flight started from

// Stop the run in flight. Results from before it stay as they were.
function cancelRun() {
  if (state.run.state !== "running") return;
  runCtl?.abort();
  update((s) => { s.run.state = "idle"; s.run.note = "Run cancelled"; });
  ctx.status("warn", "Run cancelled. The results are from before it.");
}

// Run analysis (the button, F5). While a run is going it cancels it, or, when
// the inputs changed since it started, restarts it on the new ones.
async function run() {
  const clocks = state.clocks.filter((c) => c.record && !isReference(c));
  if (state.run.state === "running") {
    const outdated = (state.run.gen ?? 0) !== runGen;
    cancelRun();
    if (!outdated) return;
  }
  if (!clocks.length || !state.devs.length) return;
  const ctl = runCtl = new AbortController();
  transient = null;
  update((s) => { s.run.state = "running"; });
  const gen = runGen = state.run.gen ?? 0;
  const t0 = performance.now();
  let failed = 0;
  await Promise.all(clocks.map(async (c) => {
    try {
      const out = await analyze(c, ctl.signal);
      update(() => { c.result = out; c.error = null; });
    } catch (err) {
      if (err.name === "AbortError") return;
      failed++;
      update(() => { c.result = null; c.error = err.message; });
    }
  }));
  if (ctl.signal.aborted) return;       // cancelled or replaced: whoever did that has the state
  const s = (performance.now() - t0) / 1000;
  const nd = state.devs.length, ok = clocks.length - failed;
  const r = clocks[0].range, n = r ? Math.floor((r[1] - r[0]) / clocks[0].record.tau0) + 1 : clocks[0].record.n;
  const what = clocks.length > 1 ? `${ok === clocks.length ? "" : `${ok} of `}${clocks.length} clocks` : `${thousands(n)} samples`;
  update((st) => {
    st.run.state = failed ? "error" : "idle";
    // Inputs changed while this ran: the results are already out of date.
    if (!failed && (st.run.gen ?? 0) === gen) st.run.stale = null;
    st.run.note = ok ? `Ran ${nd} estimator${nd === 1 ? "" : "s"} on ${what} in ${num(s, 2)} s` : "";
  });
  if (failed) ctx.status("error", `${failed} of ${clocks.length} clocks failed: ${clocks.find((c) => c.error)?.error}`);
}

// ── Header chip and status bar ──
function renderChip() {
  const dot = document.getElementById("record-dot");
  const text = document.getElementById("record-chip-text");
  const withRec = state.clocks.filter((c) => c.record);
  const fresh = !state.run.stale && state.clocks.some((c) => c.result);
  dot.className = `dot ${withRec.length ? (fresh ? "good" : "warn") : ""}`;
  if (!withRec.length) { text.innerHTML = "No record"; return; }
  if (withRec.length === 1 && !multi()) {
    const c = withRec[0], r = c.record;
    const facts = [`${thousands(r.n)} samples`, `τ<sub>0</sub> ${num(r.tau0, 4)} s`, duration(r.n * r.tau0),
      fresh ? DRIFT_FACT[c.drift] ?? "drift kept" : "not analysed yet"];
    text.innerHTML = `${h(c.source.type === "file" ? c.source.name : "synthetic record")} <span class="facts">${facts.join(" · ")}</span>`;
    return;
  }
  const n = withRec.reduce((a, c) => a + c.record.n, 0);
  text.innerHTML = `${state.clocks.length} clocks <span class="facts">${thousands(n)} samples in all${fresh ? "" : " · not analysed yet"}</span>`;
}

function renderStatus() {
  const dot = document.getElementById("status-dot");
  const text = document.getElementById("status-text");
  const hintEl = document.getElementById("status-hint");
  let kind, msg;
  const anyRecord = state.clocks.some((c) => c.record);
  if (transient) [kind, msg] = [transient.kind, transient.text];
  else if (!state.core.ok) [kind, msg] = ["error", "No analysis core. Start the explorer with serve.py."];
  else if (state.run.state === "running") [kind, msg] = ["busy", `Running ${state.devs.map(label).join(", ")}…`];
  else if (state.clocks.some((c) => c.loading)) [kind, msg] = ["busy", "Reading the record…"];
  else if (!anyRecord) [kind, msg] = ["", "No record. Open a file on the Data page."];
  else if (state.run.stale) [kind, msg] = ["warn", state.run.stale === "record" ? "New record" : "Inputs changed"];
  else [kind, msg] = ["good", "Up to date"];
  dot.className = `status-dot ${kind}`;
  text.textContent = msg;
  hintEl.textContent = kind === "warn" ? "Run analysis (F5) to update Stability, Spectrum and Dynamic" : (page?.hint() ?? "");
}

function render() {
  document.querySelectorAll(".tabs a").forEach((a) => a.classList.toggle("on", a.dataset.page === state.page));
  document.getElementById("help-btn").classList.toggle("on", state.page === "help");
  // One clock: Spectrum, Dynamic and Holdover put their settings in the panel. Several
  // clocks: the panel keeps the clocks and the page shows a toolbar.
  renderPanel(page?.controls && (!multi() || page.ownsPanel) ? page.controls() : null);
  renderChip();
  page?.refresh();
  renderStatus();
}
subscribe(() => { if (transient && Date.now() > transient.until) transient = null; render(); });
window.addEventListener("themechange", render);

// ── File menu: exports, report, session ──
const exportBtn = document.getElementById("export-btn");
const exportMenu = document.getElementById("export-menu");
exportBtn.addEventListener("click", () => {
  exportMenu.hidden = !exportMenu.hidden;
  exportBtn.setAttribute("aria-expanded", String(!exportMenu.hidden));
});
document.addEventListener("click", (e) => {
  if (!e.target.closest(".menu-wrap")) { exportMenu.hidden = true; exportBtn.setAttribute("aria-expanded", "false"); }
});
exportMenu.addEventListener("click", async (e) => {
  const kind = e.target.closest("[data-export]")?.dataset.export;
  if (!kind) return;
  exportMenu.hidden = true;
  const analysed = state.clocks.filter((c) => c.result);
  if (kind === "save-session") {
    if (!state.clocks.length) return ctx.status("warn", "There is nothing to save yet.");
    return saveFile("sigmatau-session.json", new Blob([sessionText()], { type: "application/json" }));
  }
  if (kind === "open-session") {
    const f = await openFile(".json,application/json");
    return f && openSession(f);
  }
  if (kind === "report") {
    try {
      const html = await buildReport();
      return saveFile("sigmatau-report.html", new Blob([html], { type: "text/html" }));
    } catch (err) {
      return ctx.status("warn", err.message);
    }
  }
  if (kind === "png" || kind === "svg") {
    const blob = await plotImage(kind);
    if (!blob) return ctx.status("warn", "Nothing to export on this page yet.");
    saveFile(`sigmatau-${state.page}.${kind}`, blob);
  } else if (!analysed.length) {
    ctx.status("warn", "Run the analysis first; there are no results to export.");
  } else if (kind === "json") {
    const out = analysed.map((c) => ({ clock: c.name, ...c.result.doc }));
    saveFile("sigmatau-results.json", new Blob([JSON.stringify(analysed.length === 1 ? out[0] : out, null, 2)], { type: "application/json" }));
  } else {
    const lines = ["clock\testimator\ttau_s\tsigma\tci_lo\tci_hi\tedf\tnoise\tneff"];
    for (const c of analysed) for (const r of c.result.results) r.tau.forEach((t, i) => {
      lines.push([c.name, r.dev, t, r.sigma[i], r.ci_lo[i] ?? "", r.ci_hi[i] ?? "", r.edf[i] ?? "", r.noise[i] ?? "", r.neff[i] ?? ""].join("\t"));
    });
    saveFile("sigmatau-results.tsv", new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" }));
  }
});

// Every plot shown on the page (both spectra, each clock's map), placed as on
// screen, in one PNG (at twice the size) or SVG.
async function plotImage(kind) {
  const els = [...pageEl.querySelectorAll(".js-plotly-plot")]
    .filter((el) => el.data?.length && el.offsetParent && getComputedStyle(el).visibility !== "hidden");
  if (!els.length) return null;
  const P = await loadPlotly();
  // A dark plot keeps its card behind it; light plots export on clear paper.
  const paper = document.documentElement.dataset.theme === "dark" ? css("--bg-card") : null;
  const boxes = els.map((el) => el.getBoundingClientRect());
  const x0 = Math.min(...boxes.map((b) => b.left)), y0 = Math.min(...boxes.map((b) => b.top));
  const W = Math.round(Math.max(...boxes.map((b) => b.right)) - x0), H = Math.round(Math.max(...boxes.map((b) => b.bottom)) - y0);
  const k = kind === "png" ? 2 : 1;
  const parts = await Promise.all(els.map(async (el, i) => {
    const b = boxes[i], fig = paper ? { data: el.data, layout: { ...el.layout, paper_bgcolor: paper } } : el;
    const url = await P.toImage(fig, { format: kind, scale: k, width: Math.round(b.width), height: Math.round(b.height) });
    return { x: Math.round(b.left - x0), y: Math.round(b.top - y0), w: Math.round(b.width), h: Math.round(b.height), url };
  }));
  if (kind === "svg") {
    const inner = parts.map((p) => decodeURIComponent(p.url.slice(p.url.indexOf(",") + 1)).replace("<svg ", `<svg x="${p.x}" y="${p.y}" `));
    const bg = paper ? `<rect width="100%" height="100%" fill="${paper}"/>` : "";
    return new Blob([`<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${W}" height="${H}">${bg}${inner.join("")}</svg>`], { type: "image/svg+xml" });
  }
  const canvas = Object.assign(document.createElement("canvas"), { width: W * k, height: H * k });
  const g = canvas.getContext("2d");
  if (paper) { g.fillStyle = paper; g.fillRect(0, 0, canvas.width, canvas.height); }
  for (const p of parts) {
    const img = new Image();
    img.src = p.url;
    await img.decode();
    g.drawImage(img, p.x * k, p.y * k, p.w * k, p.h * k);
  }
  return new Promise((resolve) => canvas.toBlob(resolve, "image/png"));
}

// ── Header buttons and keys ──
document.getElementById("record-chip").addEventListener("click", () => go("data"));
document.getElementById("help-btn").addEventListener("click", () => go("help"));
// Page settings (data-ctl) live in the panel or in a page toolbar.
for (const type of ["click", "change"]) {
  document.addEventListener(type, (e) => { if (e.target.closest?.("[data-ctl]")) page?.onControl?.(e); });
}
document.addEventListener("keydown", (e) => {
  // F5 is Run here, so it never reloads the page and loses the loaded clocks (Ctrl+F5 still does).
  if (e.key === "F5" && !e.ctrlKey && !e.metaKey && !e.shiftKey) e.preventDefault();
  if (document.querySelector("dialog[open]")) return;
  if (e.key === "F5" || ((e.ctrlKey || e.metaKey) && e.key === "Enter")) { e.preventDefault(); run(); return; }
  if (e.key === "Escape" && state.run.state === "running") { e.preventDefault(); cancelRun(); return; }
  if (e.target.closest("input, select, textarea") || e.ctrlKey || e.metaKey || e.altKey) return;
  const i = "12345678".indexOf(e.key);
  if (i >= 0) go(ORDER[i]);
  if (e.key === "?") go("help");
});

// ── Drag and drop anywhere ──
let dragDepth = 0;
const isFileDrag = (e) => [...(e.dataTransfer?.types ?? [])].includes("Files");
document.addEventListener("dragenter", (e) => { if (isFileDrag(e)) { dragDepth++; document.body.classList.add("dragging"); } });
document.addEventListener("dragleave", () => { if (--dragDepth <= 0) { dragDepth = 0; document.body.classList.remove("dragging"); } });
document.addEventListener("dragover", (e) => { if (isFileDrag(e)) e.preventDefault(); });
document.addEventListener("drop", async (e) => {
  if (!isFileDrag(e)) return;
  e.preventDefault();
  dragDepth = 0;
  document.body.classList.remove("dragging");
  const f = e.dataTransfer.files[0];
  if (!f) return;
  const file = await readFile(f);
  await (isSession(file.text) ? openSession(file) : importFile(file, activeClock()));
});

// ── Start ──
initPanel({ run, cancel: cancelRun, addClock, selectClock, setReference, removeClock });
go(location.hash.slice(1) || "stability");
(async () => {
  const hc = await health();
  update((s) => { s.core = { ok: !!hc.ok, name: hc.core }; });
  if (!hc.ok) return;
  // Start on a synthetic record so the plots show something straight away.
  await useSynth(null, { n: 16384, tau0: 1, seed: 42 });
  await run();
})();

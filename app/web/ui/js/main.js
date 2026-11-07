// Entry point: routing, loading records, the Run action, the header chip,
// status bar, export menu, drag and drop and keyboard shortcuts.
import { state, update, subscribe, newClock, activeClock, multi } from "./state.js";
import { health, load, analyze } from "./core.js";
import { initPanel, renderPanel } from "./panel.js";
import { openImport } from "./importer.js";
import * as stability from "./pages/stability.js";
import * as data from "./pages/data.js";
import { placeholder } from "./pages/placeholder.js";
import { label, gridName } from "./estimators.js";
import { openFile, saveFile } from "./platform.js";
import { loadPlotly } from "./plot.js";
import { thousands, num, duration, escapeHtml as h } from "./format.js";

const PAGES = { data, stability, spectrum: placeholder("spectrum"), dynamic: placeholder("dynamic") };
const ORDER = ["data", "stability", "spectrum", "dynamic"];
const pageEl = document.getElementById("page");
let page = null;
let transient = null;      // { kind, text } shown in the status bar until the next change

const ctx = {
  run, addClock, openFile: openInto, reimport, useSynth, setRecord,
  goData: () => go("data"),
  status: (kind, text) => { transient = { kind, text }; renderStatus(); },
};

// ── Routing ──
function go(name) {
  if (!PAGES[name]) name = "stability";
  if (location.hash !== `#${name}`) history.replaceState(null, "", `#${name}`);
  if (page === PAGES[name]) return;
  state.page = name;
  page = PAGES[name];
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

function markRecordStale(s) { s.run.stale = "record"; }

const baseName = (n) => n.replace(/\.[^.]+$/, "");

async function importFile(file, clock, prev) {
  const settings = await openImport(file, prev);
  if (!settings) return;
  let c = clock;
  update((s) => {
    if (!c) {
      c = newClock({ name: baseName(file.name) });
      s.clocks.push(c);
    } else if (!multi() || c.source.type !== "file" || c.name === baseName(c.source.name ?? "")) {
      c.name = baseName(file.name);
    }
    Object.assign(c, {
      source: { type: "file", name: file.name, content: file.text },
      kind: settings.kind, tau0: settings.tau0, time_col: settings.time_col, timeUnit: settings.timeUnit,
      value_col: settings.value_col, header: settings.header, scaling: settings.scaling,
      hasTimeCol: settings.hasTimeCol, sourceInfo: settings.sourceInfo,
      fillgaps: settings.hasTimeCol ? c.fillgaps : false, remove_outliers: false, result: null,
    });
    s.activeClock = c.id;
    s.dataMode = "file";
  });
  go("data");
  await loadClock(c);
}

async function openInto(clock) {
  const f = await openFile();
  if (f) await importFile(f, clock ?? null);
}

async function addClock() {
  const f = await openFile();
  if (f) await importFile(f, null);
}

async function reimport(c) {
  if (c?.source.type !== "file") return;
  await importFile({ name: c.source.name, text: c.source.content }, c, c);
}

async function useSynth(c, { n, tau0, seed }) {
  update((s) => {
    if (!c) { c = newClock(); s.clocks.push(c); s.activeClock = c.id; }
    Object.assign(c, {
      source: { type: "synth", n, seed, tau0 }, sourceInfo: "white + random-walk FM", kind: "phase",
      tau0, hasTimeCol: false, fillgaps: false, remove_outliers: false, detrend: "none", result: null,
    });
    if (c.name === "Clock" || !multi()) c.name = "Synthetic record";
  });
  await loadClock(c);
}

async function setRecord(c, fields) {
  if (!c) return;
  update(() => Object.assign(c, fields));
  await loadClock(c);
}

function selectClock(id) {
  update((s) => { s.activeClock = id; });
  go("data");
}

// ── Run ──
let runCtl = null;
async function run() {
  const clocks = state.clocks.filter((c) => c.record);
  if (!clocks.length || !state.devs.length || state.run.state === "running") return;
  runCtl?.abort();
  runCtl = new AbortController();
  transient = null;
  update((s) => { s.run.state = "running"; });
  const t0 = performance.now();
  let failed = 0;
  await Promise.all(clocks.map(async (c) => {
    try {
      const out = await analyze(c, runCtl.signal);
      update(() => { c.result = out; c.error = null; });
    } catch (err) {
      if (err.name === "AbortError") return;
      failed++;
      update(() => { c.result = null; c.error = err.message; });
    }
  }));
  const s = (performance.now() - t0) / 1000;
  const nd = state.devs.length;
  const what = clocks.length > 1 ? `${clocks.length} clocks` : `${thousands(clocks[0].record.n)} samples`;
  update((st) => {
    st.run.state = failed ? "error" : "idle";
    st.run.stale = failed ? st.run.stale : null;
    st.run.note = `Ran ${nd} estimator${nd === 1 ? "" : "s"} on ${what} in ${num(s, 2)} s`;
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
      fresh ? (c.detrend === "none" ? "no detrend" : `${c.detrend} detrend`) : "not analysed yet"];
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
  renderPanel();
  renderChip();
  page?.refresh();
  renderStatus();
}
subscribe(() => { if (state.run.state !== "running") transient = transient?.kind === "error" ? transient : null; render(); });

// ── Export ──
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
  if (kind === "png" || kind === "svg") {
    const el = page?.plotElement();
    if (!el?.data) return ctx.status("warn", "Nothing to export on this page yet.");
    const P = await loadPlotly();
    const url = await P.toImage(el, { format: kind, scale: kind === "png" ? 2 : 1, width: el.clientWidth, height: el.clientHeight });
    saveFile(`sigmatau-${state.page}.${kind}`, await (await fetch(url)).blob());
  } else if (!analysed.length) {
    ctx.status("warn", "Run the analysis first; there are no results to export.");
  } else if (kind === "json") {
    const out = analysed.map((c) => ({ clock: c.name, ...c.result }));
    saveFile("sigmatau-results.json", new Blob([JSON.stringify(analysed.length === 1 ? out[0] : out, null, 2)], { type: "application/json" }));
  } else {
    const lines = ["clock\testimator\ttau_s\tsigma\tci_lo\tci_hi\tedf\tnoise\tneff"];
    for (const c of analysed) for (const r of c.result.results) r.tau.forEach((t, i) => {
      lines.push([c.name, r.dev, t, r.sigma[i], r.ci_lo[i] ?? "", r.ci_hi[i] ?? "", r.edf[i] ?? "", r.noise[i] ?? "", r.neff[i] ?? ""].join("\t"));
    });
    saveFile("sigmatau-results.tsv", new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" }));
  }
});

// ── Header buttons and keys ──
document.getElementById("record-chip").addEventListener("click", () => go("data"));
document.getElementById("help-btn").addEventListener("click", () => document.getElementById("help-dialog").showModal());
document.addEventListener("keydown", (e) => {
  if (document.querySelector("dialog[open]")) return;
  if (e.key === "F5" || ((e.ctrlKey || e.metaKey) && e.key === "Enter")) { e.preventDefault(); run(); return; }
  if (e.target.closest("input, select, textarea") || e.ctrlKey || e.metaKey || e.altKey) return;
  const i = "1234".indexOf(e.key);
  if (i >= 0) go(ORDER[i]);
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
  if (f) await importFile({ name: f.name, text: await f.text() }, activeClock());
});

// ── Start ──
initPanel({ run, addClock, selectClock });
go(location.hash.slice(1) || "stability");
(async () => {
  const hc = await health();
  update((s) => { s.core = { ok: !!hc.ok, name: hc.core }; });
  if (!hc.ok) return;
  // Start on a synthetic record so the plots show something straight away.
  await useSynth(null, { n: 16384, tau0: 1, seed: 42 });
  await run();
})();

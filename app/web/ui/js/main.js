// Entry point: page routing, the analysis runner, header, status bar,
// export menu and keyboard shortcuts.
import { state, update, subscribe, newClock, activeClocks } from "./state.js";
import { health, analyze } from "./api.js";
import { loadPlotly } from "./plot.js";
import { mountStability } from "./pages/stability.js";
import { mountData } from "./pages/data.js";
import { mountPlaceholder } from "./pages/placeholder.js";
import { label } from "./estimators.js";
import { saveFile } from "./platform.js";
import { thousands, num, duration } from "./format.js";

const PAGES = ["data", "stability", "spectrum", "dynamic"];
const pageEl = document.getElementById("page");
let current = null;     // { name, refresh }

// ── Runner ──
let controller = null;
let autoTimer = null;

async function run() {
  clearTimeout(autoTimer);
  const clocks = activeClocks();
  if (!clocks.length || !state.devs.length || !state.core.ok) return;
  controller?.abort();
  controller = new AbortController();
  const signal = controller.signal;
  update((s) => { s.run.state = "running"; });
  const t0 = performance.now();
  await Promise.all(clocks.map(async (c) => {
    try {
      const res = await analyze(c, signal);
      update((s) => { const k = s.clocks.find((x) => x.id === c.id); if (k) { k.result = res; k.error = null; } });
    } catch (err) {
      if (err.name === "AbortError") return;
      update((s) => { const k = s.clocks.find((x) => x.id === c.id); if (k) { k.error = err.message; k.result = null; } });
    }
  }));
  if (signal.aborted) return;
  const ms = performance.now() - t0;
  const failed = clocks.filter((c) => state.clocks.find((x) => x.id === c.id)?.error).length;
  const samples = clocks.reduce((n, c) => n + (state.clocks.find((x) => x.id === c.id)?.result?.record.n ?? 0), 0);
  update((s) => {
    s.run.state = failed ? "error" : "ok";
    s.run.lastMs = ms;
    s.run.summary = `Ran ${s.devs.length} estimator${s.devs.length === 1 ? "" : "s"} on ${thousands(samples)} samples in ${num(ms / 1000, 2)} s`;
  });
}

// Re-run automatically when inputs change and the core is fast enough.
function scheduleAuto() {
  update((s) => { if (s.run.state !== "running") s.run.state = "stale"; });
  clearTimeout(autoTimer);
  const fast = state.run.lastMs === null || state.run.lastMs < 1500;
  if (state.run.autoRun && fast) autoTimer = setTimeout(run, 250);
}

// ── Routing ──
function mount(name) {
  if (!PAGES.includes(name)) name = "stability";
  state.page = name;
  document.querySelectorAll(".tabs [data-page]").forEach((a) => {
    const on = a.dataset.page === name;
    a.classList.toggle("on", on);
    a.setAttribute("aria-selected", String(on));
  });
  current = name === "stability" ? mountStability(pageEl, { run })
    : name === "data" ? mountData(pageEl)
    : mountPlaceholder(pageEl, name);
  current.refresh();
}
window.addEventListener("hashchange", () => mount(location.hash.slice(1)));

// ── Header and status ──
function renderChrome() {
  const pill = document.getElementById("record-pill-text");
  const clocks = state.clocks;
  if (!clocks.length) pill.textContent = "No record";
  else if (clocks.length > 1) pill.textContent = `${clocks.length} clocks · ${activeClocks().length} shown`;
  else {
    const c = clocks[0], rec = c.result?.record;
    pill.innerHTML = rec
      ? `<b>${c.name}</b> ${thousands(rec.n)} samples · τ₀ ${num(rec.tau0, 3)} s · ${duration(rec.n * rec.tau0)} · ${c.detrend === "none" ? "no detrend" : c.detrend}`
      : `<b>${c.name}</b> not run yet`;
  }
  const dot = document.getElementById("status-dot");
  const text = document.getElementById("status-text");
  const hint = document.getElementById("status-hint");
  const r = state.run;
  let cls = "idle", msg = "Ready";
  if (!state.core.ok) { cls = "error"; msg = "No analysis core found"; }
  else if (r.state === "running") { cls = "busy"; msg = "Running…"; }
  else if (r.state === "error") { cls = "error"; msg = "The last run failed; see the plot area"; }
  else if (r.state === "stale") { cls = "warn"; msg = "Inputs changed; press Run analysis (F5)"; }
  else if (r.state === "ok") { cls = "ok"; msg = "Up to date"; }
  dot.className = `status-dot ${cls}`;
  text.textContent = msg;
  hint.textContent = state.page === "stability" || state.page === "data"
    ? "Hover a point for its value · drag to zoom · double-click to reset" : "";
}

subscribe((what) => {
  if (what === "analysis") {
    scheduleAuto();
    return;            // scheduleAuto's own update re-renders
  }
  current?.refresh();
  renderChrome();
});

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
exportMenu.addEventListener("click", (e) => {
  const kind = e.target.closest("[data-export]")?.dataset.export;
  if (!kind) return;
  exportMenu.hidden = true;
  const plot = document.querySelector(".plot:not([hidden]) .js-plotly-plot, .plot.js-plotly-plot:not([hidden])");
  if ((kind === "png" || kind === "svg") && plot) {
    window.Plotly.downloadImage(plot, { format: kind, scale: 2, width: 1200, height: 800, filename: `sigmatau-${state.page}` });
  } else if (kind === "tsv") {
    const lines = ["clock\tdev\ttau_s\tsigma\tci_lo\tci_hi\tedf\tnoise\tneff"];
    for (const c of state.clocks) for (const r of c.result?.results ?? []) {
      r.tau.forEach((t, i) => lines.push([c.name, r.dev, t, r.sigma[i], r.ci_lo[i] ?? "", r.ci_hi[i] ?? "", r.edf[i] ?? "", r.noise[i] ?? "", r.neff[i] ?? ""].join("\t")));
    }
    saveFile("sigmatau-results.tsv", new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" }));
  } else if (kind === "json") {
    const out = state.clocks.map((c) => ({ name: c.name, result: c.result }));
    saveFile("sigmatau-results.json", new Blob([JSON.stringify(out, null, 1)], { type: "application/json" }));
  }
});

document.getElementById("help-btn").addEventListener("click", () => document.getElementById("help-dialog").showModal());
document.getElementById("record-pill").addEventListener("click", () => { location.hash = "#data"; });

document.addEventListener("keydown", (e) => {
  if (e.key === "F5") { e.preventDefault(); run(); return; }
  const typing = e.target.closest("input, select, textarea, [contenteditable]") || document.querySelector("dialog[open]");
  if (!typing && /^[1-4]$/.test(e.key)) location.hash = "#" + PAGES[Number(e.key) - 1];
});

// Empty-state button: start with a synthetic record.
pageEl.addEventListener("click", (e) => {
  if (e.target.closest('[data-action="synth"]')) {
    update((s) => { s.clocks.push(newClock({ name: "Synthetic" })); }, "analysis");
  }
});

// ── Start ──
(async () => {
  const [h] = await Promise.all([health(), loadPlotly().catch(() => null)]);
  state.core = { ok: !!h.ok, name: h.core };
  mount(location.hash.slice(1) || "stability");
  renderChrome();
})();

// The side panel: clocks (when there are several), estimators, options and
// the Run button. One clock: colour = estimator. Several: colour = clock and
// line style = estimator.
import { state, update, markStale, multi, isReference, referenceColFor } from "./state.js";
import { byId, GRIDS, CONFIDENCES, gridName } from "./estimators.js";
import { seriesColor, strokeSvg } from "./plot.js";
import { escapeHtml as h, thousands } from "./format.js";
import { openPicker } from "./picker.js";

let handlers = { run() {}, addClock() {}, selectClock() {}, setReference() {} };
export function initPanel(h) {
  handlers = { ...handlers, ...h };
  const side = document.getElementById("side");
  side.addEventListener("click", onClick);
  side.addEventListener("change", onChange);
}

const confLabel = (c) => (CONFIDENCES.find(([v]) => v === c) ?? [c, `${c * 100} %`])[1];

function clockRows() {
  return state.clocks.map((c) => {
    const n = c.record?.n;
    const ref = isReference(c);
    const against = !ref && state.reference !== null && referenceColFor(c) ? "vs reference" : "";
    const sub = ref ? "reference" : [c.source.type === "file" ? c.source.name : "synthetic", n ? `${thousands(n)} pts` : c.error ? "error" : "…", against].filter(Boolean).join(" · ");
    return `<div class="row clock-row ${c.id === state.activeClock ? "active" : ""} ${c.visible ? "" : "hidden-clock"}" data-clock="${c.id}" role="button" tabindex="0">
      <span class="row-main"><span class="clock-dot" style="background:${seriesColor(c.slot)}"></span>
        <span class="clock-text"><b>${h(c.name)}</b><span>${h(sub)}</span></span></span>
      ${ref ? `<span class="ref-tag" title="Reference clock">ref</span>` : `<button class="eye" type="button" data-eye="${c.id}" title="${c.visible ? "Hide" : "Show"} ${h(c.name)} on the plots" aria-pressed="${c.visible}">${c.visible ? "◉" : "○"}</button>`}
    </div>`;
  }).join("");
}

// "As measured" or any clock that shares its original reference and time
// stamps with others: today, clocks read from the same file.
function referenceSelect() {
  const against = [...new Set(state.clocks.map((c) => c.against).filter(Boolean))];
  const groupSize = (g) => state.clocks.filter((c) => c.group !== null && c.group === g).length;
  const opts = state.clocks.map((c) => {
    const ok = c.source.type === "file" && groupSize(c.group) > 1;
    return `<option value="${c.id}" ${c.id === state.reference ? "selected" : ""} ${ok ? "" : "disabled"}>${h(c.name)}${ok ? "" : " (no shared reference)"}</option>`;
  }).join("");
  const r = state.clocks.find((c) => c.id === state.reference);
  const others = r ? state.clocks.filter((c) => c !== r && c.group !== r.group).length : 0;
  return `<div class="ref-block">
    <label class="opt" title="Only clocks measured together, in one file, share a reference and time stamps"><span>Reference</span>
      <select class="wide" data-ref-select><option value="">As measured${against.length === 1 ? ` (vs ${h(against[0])})` : ""}</option>${opts}</select></label>
    ${others ? `<div class="muted small">${others} clock${others > 1 ? "s" : ""} from other files stay as measured.</div>` : ""}
  </div>`;
}

function estimatorRows(several) {
  return state.devs.map((id, i) => {
    const e = byId[id];
    const mark = several ? strokeSvg(i, "var(--ink)", 26) : `<span class="swatch" style="background:${seriesColor(i)}"></span>`;
    return `<div class="row"><span class="row-main">${mark}<span class="row-name">${e.label}</span>${several ? "" : `<span class="row-family">${e.family}</span>`}</span>
      <button class="row-x" type="button" data-remove="${id}" title="Remove ${e.label}">×</button></div>`;
  }).join("");
}

function options() {
  return `<div class="opts">
    <label class="opt"><span>τ grid</span><select class="wide" data-opt="grid">${GRIDS.map(([v, l]) => `<option value="${v}" ${v === state.grid ? "selected" : ""}>${l}</option>`).join("")}</select></label>
    <div class="opt"><label class="check"><input type="checkbox" data-opt="ci" ${state.ci ? "checked" : ""}>Confidence intervals</label>
      <select class="narrow" data-opt="confidence" ${state.ci ? "" : "disabled"}>${CONFIDENCES.map(([v, l]) => `<option value="${v}" ${v === state.confidence ? "selected" : ""}>${l}</option>`).join("")}</select></div>
    <label class="check"><input type="checkbox" data-opt="band" ${state.band ? "checked" : ""} ${state.ci ? "" : "disabled"}>Show intervals as a band</label>
    <label class="check"><input type="checkbox" data-opt="bias" ${state.bias ? "checked" : ""}>Bias-correct totals and Thêo1</label>
  </div>`;
}

export function renderPanel() {
  const side = document.getElementById("side");
  const several = multi();
  const nDev = state.devs.length;
  const run = state.run;
  const stale = run.stale && state.clocks.some((c) => c.record);
  const runLabel = run.state === "running" ? "Running…"
    : run.stale === "record" && state.clocks.some((c) => c.result) ? "New record, run analysis" : "Run analysis";

  const scroll = side.querySelector(".side-scroll");
  const top = scroll ? scroll.scrollTop : 0;
  side.innerHTML = `
    <div class="side-scroll">
      ${several ? `
      <section><div class="side-head"><h2>Clocks</h2><span class="note">colour = clock</span></div></section>
      <div class="rows">${clockRows()}</div>
      ${referenceSelect()}
      <button class="btn-dashed" type="button" data-act="add-clock">+ Add clock</button>
      <div class="rule"></div>` : ""}
      <div class="side-head"><h2>Estimators</h2><span class="note">${several ? "line style = estimator" : `${nDev} selected`}</span></div>
      <div class="rows">${estimatorRows(several)}</div>
      <button class="btn-dashed" type="button" data-act="pick">+ Add estimator</button>
      <div class="rule"></div>
      ${several ? `
        <button class="opts-toggle" type="button" data-act="opts" aria-expanded="${state.optionsOpen}"><h2>Options</h2>
          <span class="note">${gridName(state.grid).toLowerCase()} τ · ${state.ci ? `${confLabel(state.confidence)} CI` : "no CI"} ${state.optionsOpen ? "▾" : "▸"}</span></button>
        ${state.optionsOpen ? options() : ""}`
      : `<h2>Options</h2>${options()}`}
    </div>
    <div class="side-foot">
      <button class="run-btn" type="button" data-act="run" ${run.state === "running" || !state.clocks.length || !nDev ? "disabled" : ""}>
        ${stale ? `<span class="stale-dot" aria-hidden="true"></span>` : ""}${runLabel}<kbd>F5</kbd></button>
      <div class="run-note">${h(run.note || "")}</div>
    </div>`;
  const ns = side.querySelector(".side-scroll");
  if (ns) ns.scrollTop = top;
}

function onClick(e) {
  const t = e.target;
  const eye = t.closest("[data-eye]");
  if (eye) {
    const id = Number(eye.dataset.eye);
    update((s) => { const c = s.clocks.find((k) => k.id === id); c.visible = !c.visible; });
    return;
  }
  const rm = t.closest("[data-remove]");
  if (rm) {
    update((s) => { s.devs = s.devs.filter((d) => d !== rm.dataset.remove); markStale(s, "inputs"); });
    return;
  }
  const row = t.closest("[data-clock]");
  if (row) { handlers.selectClock(Number(row.dataset.clock)); return; }
  const act = t.closest("[data-act]")?.dataset.act;
  if (act === "pick") openPicker();
  else if (act === "run") handlers.run();
  else if (act === "add-clock") handlers.addClock();
  else if (act === "opts") update((s) => { s.optionsOpen = !s.optionsOpen; });
}

function onChange(e) {
  if ("refSelect" in e.target.dataset) { handlers.setReference(e.target.value ? Number(e.target.value) : null); return; }
  const opt = e.target.dataset.opt;
  if (!opt) return;
  const v = e.target.type === "checkbox" ? e.target.checked : e.target.value;
  update((s) => {
    if (opt === "grid") { s.grid = v; markStale(s, "inputs"); }
    if (opt === "confidence") { s.confidence = Number(v); markStale(s, "inputs"); }
    if (opt === "ci") { s.ci = v; markStale(s, "inputs"); }
    if (opt === "bias") { s.bias = v; markStale(s, "inputs"); }
    if (opt === "band") s.band = v;       // display only
  });
}

// Synthesis page: define the noise of one or more simulated clocks, see the
// stability it predicts (closed form, no simulation), and generate the clocks
// on the Data page. Each noise is set as the clock model's coefficient or as
// the power-law h_α; the two follow each other.
import { state, update } from "../state.js";
import { modelPreview } from "../core.js";
import { loadPlotly, baseLayout, logAxis, CONFIG, seriesColor } from "../plot.js";
import { escapeHtml as h } from "../format.js";

const MAX_CLOCKS = 8;
const DEVS = [["mdev", "MDEV"], ["adev", "ADEV"], ["mhdev", "MHDEV"]];
// Per noise: the model coefficient, its power-law h_α as a function of it and of τ₀, and the inverse.
const NOISES = [
  { id: "R", name: "White PM", q: "R (s²)", h: "h₂", toH: (q, t) => 8 * Math.PI ** 2 * t * q, toQ: (v, t) => v / (8 * Math.PI ** 2 * t) },
  { id: "q1", name: "White FM", q: "q₁ (s)", h: "h₀", toH: (q) => 2 * q, toQ: (v) => v / 2 },
  { id: "q2", name: "Random-walk FM", q: "q₂ (1/s)", h: "h₋₂", toH: (q) => q / (2 * Math.PI ** 2), toQ: (v) => v * 2 * Math.PI ** 2 },
  { id: "q3", name: "Random-run FM", q: "q₃ (1/s³)", h: "h₋₄", toH: (q) => q / (8 * Math.PI ** 4), toQ: (v) => v * 8 * Math.PI ** 4 },
];

let root, ctx, Plotly;
let inflight = null, inflightCtl = null, shownKeys = 0;

const fmt = (x) => (x === 0 ? "0" : x.toExponential(3).replace(/\.?0+e/, "e"));
const model = (c) => NOISES.map((n) => c[n.id]);

export const ownsPanel = true;

export function controls() {
  const S = state.synth;
  return `<h2>Record</h2>
    <label class="opt"><span>Samples</span><input type="text" data-ctl="n" value="${S.n}" inputmode="numeric"></label>
    <label class="opt"><span>τ<sub>0</sub> (s)</span><input type="text" data-ctl="tau0" value="${S.tau0}" inputmode="decimal"></label>
    <label class="opt"><span>Seed</span><input type="text" data-ctl="seed" value="${S.seed}" inputmode="numeric"></label>
    <div class="ctl-block"><label class="lbl">Preview</label><div class="seg" role="group" aria-label="Deviation">${DEVS.map(([id, label]) =>
      `<button type="button" data-ctl="dev" data-v="${id}" class="${S.dev === id ? "on" : ""}">${label}</button>`).join("")}</div></div>
    <div class="rule"></div>
    <button class="btn btn-primary run-btn" type="button" data-ctl="generate">Generate ${S.clocks.length > 1 ? `${S.clocks.length} clocks` : "the clock"}</button>
    <p class="muted small">The curves are what the models predict, not a simulation. Generate simulates the clocks on the Data page${S.clocks.length > 1 ? " and forms their time scale on the Ensemble page" : ""}.</p>`;
}

function setPositive(t, key, v, whole) {
  const x = Number(v);
  if (v !== "" && Number.isFinite(x) && x > 0 && (!whole || Number.isInteger(x))) { update((s) => { s.synth[key] = x; }); return; }
  t.value = state.synth[key];
  ctx.status("error", whole ? `${key === "n" ? "Samples" : "Seed"} must be a whole number${key === "n" ? " of at least 64" : ""}.` : "τ₀ must be a time above 0 s.");
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl, S = state.synth;
  if (k === "dev" && e.type === "click") update((s) => { s.synth.dev = t.dataset.v; });
  if (k === "generate" && e.type === "click") generate();
  if (e.type !== "change") return;
  const v = t.value.trim();
  if (k === "n") { if (Number(v) >= 64) setPositive(t, "n", v, true); else { t.value = S.n; ctx.status("error", "Samples must be a whole number of at least 64."); } }
  if (k === "tau0") setPositive(t, "tau0", v, false);
  if (k === "seed") { if (v !== "" && Number.isInteger(Number(v)) && Number(v) >= 0) update((s) => { s.synth.seed = Number(v); }); else { t.value = S.seed; ctx.status("error", "The seed is a whole number, 0 or more."); } }
}

function generate() {
  const S = state.synth;
  if (!S.clocks.some((c) => model(c).some((x) => x > 0))) { ctx.status("error", "Give at least one clock some noise."); return; }
  const k = S.clocks.length;
  ctx.useSynth(null, { n: S.n, tau0: S.tau0, seed: S.seed, clocks: k, models: S.clocks.map(model), names: S.clocks.map((c) => c.name) })
    .then(() => ctx.go(k > 1 ? "ensemble" : "data"));
}

export function mount(el, c) {
  root = el; ctx = c; shownKeys = 0;
  root.innerHTML = `
    <div class="col">
      <div class="card plot-card"><div class="card-row"><h2>Predicted stability</h2><span class="muted small" data-ref="sub"></span></div>
        <div class="plot" data-ref="plot"></div><div class="empty" data-ref="empty" hidden></div></div>
    </div>
    <div class="card side-card synth-clocks" data-ref="clocks"></div>`;
  root.addEventListener("change", onEdit);
  root.addEventListener("click", onClick);
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

// ── the clocks being defined ──
function clocksHtml() {
  const S = state.synth;
  return `<div class="card-row"><h2>Clocks</h2><button class="btn" type="button" data-act="add" ${S.clocks.length >= MAX_CLOCKS ? "disabled" : ""}>Add a clock</button></div>
    ${S.clocks.map((c, i) => `<div class="synth-clock" data-i="${i}">
      <div class="name-row"><span class="clock-dot" style="background:${seriesColor(i)}"></span>
        <input class="clock-name" data-name aria-label="Clock name" value="${h(c.name)}">
        <button class="btn" type="button" data-act="remove" ${S.clocks.length < 2 ? "disabled" : ""} aria-label="Remove ${h(c.name)}">Remove</button></div>
      <div class="noise-grid"><span class="muted small">Noise</span><span class="muted small">Model</span><span class="muted small">Power law</span>
        ${NOISES.map((n) => `<span class="noise-name">${n.name}</span>
          <label><span class="sr-only">${n.q}</span><input type="text" data-q="${n.id}" inputmode="decimal" placeholder="${n.q}" title="${n.q}"></label>
          <label><span class="sr-only">${n.h}</span><input type="text" data-h="${n.id}" inputmode="decimal" placeholder="${n.h}" title="${n.h}"></label>`).join("")}
      </div></div>`).join("")}
    <p class="muted small synth-note">Each noise is given as the clock model's coefficient (S<sub>y</sub> to <i>q</i> as in the ClockModel docs) or as its power-law h<sub>α</sub>; typing one fills in the other. Flicker noise is not offered here.</p>`;
}

// Put the numbers in the inputs, except the one being typed in.
function fill() {
  const S = state.synth;
  root.querySelectorAll(".synth-clock").forEach((row) => {
    const c = S.clocks[Number(row.dataset.i)];
    for (const n of NOISES) {
      const q = row.querySelector(`[data-q=${n.id}]`), hh = row.querySelector(`[data-h=${n.id}]`);
      if (q !== document.activeElement) q.value = fmt(c[n.id]);
      if (hh !== document.activeElement) hh.value = fmt(n.toH(c[n.id], S.tau0));
    }
  });
}

function onEdit(e) {
  const row = e.target.closest(".synth-clock");
  if (!row) return;
  const c = state.synth.clocks[Number(row.dataset.i)];
  if (e.target.dataset.name !== undefined) { update(() => { c.name = e.target.value.trim() || c.name; }); return; }
  const id = e.target.dataset.q ?? e.target.dataset.h;
  const n = NOISES.find((x) => x.id === id);
  const v = Number(e.target.value);
  if (e.target.value.trim() === "" || !Number.isFinite(v) || v < 0) { ctx.status("error", `${n.name} must be a number, 0 or more.`); fill(); return; }
  update(() => { c[id] = e.target.dataset.q !== undefined ? v : n.toQ(v, state.synth.tau0); });
}

function onClick(e) {
  const act = e.target.closest("[data-act]")?.dataset.act;
  if (!act) return;
  const S = state.synth, i = Number(e.target.closest(".synth-clock")?.dataset.i);
  if (act === "add" && S.clocks.length < MAX_CLOCKS) {
    const last = S.clocks[S.clocks.length - 1];
    update((s) => { s.synth.clocks.push({ ...last, name: `Sim ${S.clocks.length + 1}` }); });
  }
  if (act === "remove" && S.clocks.length > 1) update((s) => { s.synth.clocks.splice(i, 1); });
}

// ── the preview ──
function ensure() {
  const S = state.synth;
  const body = { n: S.n, tau0: S.tau0, models: S.clocks.map(model) };
  const key = JSON.stringify(body);
  if (S.preview?.key === key || inflight === key || !state.core.ok) return;
  inflightCtl?.abort();
  const ctl = inflightCtl = new AbortController();
  inflight = key;
  S.error = null;
  modelPreview(body, ctl.signal).then((doc) => {
    if (inflight !== key) return;
    inflight = null;
    update((s) => { s.synth.preview = { key, doc }; });
  }).catch((err) => {
    if (err.name === "AbortError" || inflight !== key) return;
    inflight = null;
    update((s) => { s.synth.error = err.message; s.synth.preview = null; });
  });
}

function showEmpty(text) {
  const e = root.querySelector("[data-ref=empty]");
  e.innerHTML = text ? `<p>${text}</p>` : "";
  e.hidden = !text;
  root.querySelector("[data-ref=plot]").style.visibility = text ? "hidden" : "visible";
}

export function refresh() {
  if (!root) return;
  const S = state.synth;
  // The clock editor is rebuilt only when clocks come or go, so typing keeps its focus.
  const box = root.querySelector("[data-ref=clocks]");
  if (shownKeys !== S.clocks.length && document.activeElement?.tagName !== "INPUT") { box.innerHTML = clocksHtml(); shownKeys = S.clocks.length; }
  fill();
  ensure();
  const doc = S.preview?.doc;
  if (S.error) return showEmpty(h(S.error));
  if (!doc) return showEmpty(state.core.ok ? "Working out the stability…" : "No analysis core. Start the explorer with serve.py.");
  const lines = doc.clocks.map((c, i) => ({ i, y: c[S.dev] })).filter((l) => l.y && l.y.some((v) => v > 0));
  const label = DEVS.find(([id]) => id === S.dev)[1];
  const missing = doc.clocks.some((c) => !c[S.dev]);
  if (!lines.length) return showEmpty(missing ? `Random-run FM has no ${label}. Choose MHDEV.` : "Give a clock some noise to see its stability.");
  showEmpty("");
  root.querySelector("[data-ref=sub]").textContent = `${label} · ${doc.n} samples · τ₀ = ${doc.tau0} s${missing ? ` · random-run FM has no ${label}` : ""}`;
  if (!Plotly) return;
  Plotly.react(root.querySelector("[data-ref=plot]"), lines.map(({ i, y }) => ({
    x: doc.tau, y, name: S.clocks[i].name, mode: "lines+markers", line: { color: seriesColor(i), width: 2 }, marker: { color: seriesColor(i), size: 5 },
    hovertemplate: `<b>${h(S.clocks[i].name)}</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>` })),
  { ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis(`${label}, σ<sub>y</sub>(τ)`), showlegend: true, uirevision: "synth" }, CONFIG);
}

export function hint() { return "Type a noise as the model coefficient or as h_α, and the other follows"; }

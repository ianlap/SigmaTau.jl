// Ensemble page: several clocks read against one reference. One view forms a
// time scale from them (AT1, AT2, EWFA or a Kalman scale); another separates
// the stability of each clock by an N-cornered hat or the Groslambert
// covariance; the third takes two columns as two channels measuring one
// device, for the cross Allan deviation and the cross-spectrum. The clocks
// are the value columns of one loaded file; without such a file the core
// simulates four (not for two channels). Each is computed when this page
// needs it and recomputed when the settings or the clocks change.
import { state, update, isReference } from "../state.js";
import { ensemble, synthArgs } from "../core.js";
import { loadPlotly, baseLayout, logAxis, linAxis, CONFIG, seriesColor, css, LINE_W } from "../plot.js";
import { byId } from "../estimators.js";
import { gridField, confidenceField, onSharedControl } from "../panel.js";
import { sci, num, interval, timeAxisUnit, inWords as plain, escapeHtml as h } from "../format.js";

let root, ctx, Plotly;
let inflight = null;           // key of the run in progress
let inflightCtl = null;        // stops it (the core process too) when the settings change
const ALGS = [["at1", "AT1"], ["at2", "AT2"], ["ewfa", "EWFA"], ["greenhall", "Greenhall"], ["jones_tryon", "Jones–Tryon"], ["stein", "Stein"]];
const KALMAN = new Set(["greenhall", "jones_tryon", "stein"]);
const CAPPED = new Set(["at1", "at2", "stein"]);                // take --wmax
const FLICKER = new Set(["greenhall", "jones_tryon"]);          // take --flicker
const ABOUT_FLICKER = " Flicker states: on keeps each clock's flicker FM and flicker PM in its model as Markov states, across the τ where its fit finds them; off folds flicker PM into the reading noise and leaves flicker FM out.";
const MODELS = [["", "As simulated"], ["fitted", "Fitted"]];
const ABOUT_MODELS = { true: " The models are fitted from the comparisons, as for a file.", false: " The models are those the clocks were simulated from." };
const ABOUT = {
  at1: "AT1 weighs each clock by the inverse of its prediction-error variance, learned as it runs.",
  at2: "AT2 keeps AT1's weights and follows each clock's frequency with a Kalman update. Frequency steps: best reruns the scale from a detected step; real keeps the values already issued and corrects from the step on.",
  ewfa: "EWFA weighs every clock the same.",
  greenhall: "Greenhall's Kalman scale runs one filter over a model of each clock and reads the scale at the pivot clock. It needs three or more clocks.",
  jones_tryon: "The Jones–Tryon Kalman scale: Greenhall's filter without the covariance reduction, so it has no per-clock weights. It needs three or more clocks.",
  stein: "Stein's Kalman scale (KAS-2) takes a model of each clock, filters each clock's difference from the pivot, and weighs each clock's phase, frequency and aging by the inverse of its noise. The weights shown are the phase weights. It needs three or more clocks.",
};
const HAT_DEVS = ["adev", "mdev", "hdev", "mhdev", "pdev", "tdev", "htdev"];
const HAT_METHODS = [["hat", "Hat"], ["groslambert", "Groslambert covariance"], ["klts", "Groslambert, KLTS interval"]];
const ABOUT_HAT = {
  hat: "The hat takes the deviation of every pair of clocks, which holds the noise of both, and solves for each clock.",
  groslambert: "The Groslambert covariance multiplies the differences of two records that share a clock. Noise that each record of a pair has on its own, such as its counter's, stays out. For three clocks read against one reference it equals the hat on ADEV.",
  klts: "The Groslambert covariance with the interval of Lantz et al. (2019), for three clocks. Where a clock is not resolved, ▽ marks the upper limit.",
};
// The pivot as the core takes it: 1-based, back to the first clock when out of range.
const pivotOf = (n) => (state.ens.pivot <= n ? state.ens.pivot : 1);
// What each clock's curve is taken against. Only the reference compares a
// clock with the scale: a clock against the scale has its own share of the
// scale taken out. Simulated clocks have an ideal reference, so that is the
// default for them; a file's reference may be noisier than its clocks.
const AGAINST = [["reference", "Reference"], ["scale", "Scale"]];
const againstOf = (clocks) => state.ens.against ?? (clocks[0]?.source.type === "file" ? "scale" : "reference");
// Whether a Kalman scale fits its models: always for a file, and for simulated
// clocks when asked, in place of the models they were made from.
const simulated = (clocks) => clocks[0]?.source.type !== "file";
const fittedOf = (clocks) => !simulated(clocks) || state.ens.fitted;
// The file goes to the core as read, without the preparation of the Data page,
// so the core's remedy for gaps (--fill-gaps) has no counterpart on this page.
const GAPS = [
  [/^clock_models: .* usable τ points/, "Greenhall, Jones–Tryon and Stein fit a model to each clock, and this file leaves one of them too few samples at which it is read with two others. AT1, AT2 and EWFA form a scale from it."],
  [/^the KLTS interval needs/, "The KLTS interval needs a file without missing samples. The hat and the Groslambert covariance separate the clocks of this one."],
  [/^the frequency record is missing/, "Frequency columns with missing samples cannot be read."],
];
const FEW = [/^hat needs three or more simulated/, "Simulated clocks have a perfect reference, which counts as no clock. Simulate three or more, on the Data or the Synthesis page."];
// Simulated clocks have no gaps: a fit that stops there met noise it has no term for.
const UNFITTED = "The fit leaves one of these clocks too few averaging times. It has terms for white and flicker PM, white FM and random-walk FM, and none for random-run FM. With Models at As simulated the scale takes the models the clocks were made from.";
// Two channels are read with the Data page's gap filling, which names the remedy.
const CROSS_GAPS = [/^cross analysis needs records without gaps/, "The two channels need every sample: turn on Fill gaps by reflection on the Data page."];
const inWords = (message, sim) => (sim && GAPS[0][0].test(message) ? UNFITTED : null)
  ?? (CROSS_GAPS[0].test(message) ? CROSS_GAPS[1] : null)
  ?? GAPS.find(([re]) => re.test(message))?.[1].concat(" The gap filling of the Data page does not apply here.")
  ?? (FEW[0].test(message) ? FEW[1] : plain(message));

export const ownsPanel = true;

// The clocks of the page: the first loaded file with two or more phase or
// frequency columns, else the first ensemble simulated on the Data page.
// [] = the core simulates four. A scale leaves the reference clock out; a hat
// (`all`) keeps it as one of the clocks.
function fileClocks(all = false) {
  const groups = new Map();
  for (const c of state.clocks) {
    if (c.source.type !== "file" || c.group === null || c.format || !c.record || (!all && isReference(c))) continue;
    groups.set(c.group, [...(groups.get(c.group) ?? []), c]);
  }
  const file = [...groups.values()].find((g) => g.length >= 2);
  if (file) return file;
  const sims = new Map();
  for (const c of state.clocks) {
    const src = c.source;
    if (src.type !== "synth" || !(src.clocks > 1) || !c.record) continue;
    const key = JSON.stringify([src.n, src.seed, src.tau0, src.clocks, src.models ?? null]);
    sims.set(key, [...(sims.get(key) ?? []), c]);
  }
  const ens = [...sims.values()].find((g) => g.length === g[0].source.clocks);
  return ens ? ens.sort((a, b) => a.source.index - b.source.index) : [];
}

// The command, the flags that give it its clocks, and their file.
function clockArgs(command, clocks) {
  const a = [command];
  let file = null;
  if (clocks.length && clocks[0].source.type === "synth") {
    const src = clocks[0].source;
    a.push(...synthArgs(src), "--tau0", src.tau0);
  } else if (clocks.length) {
    const c = clocks[0];
    file = { name: c.source.name, content: c.source.content };
    a.push("--column", clocks.map((k) => k.value_col).join(","), "--time-column", c.time_col, "--header", c.header, "--scaling", c.scaling);
    if (c.tau0) a.push("--tau0", c.tau0);
    if (c.kind === "frequency") a.push("--freq");
  } else a.push("--synth", 16384, "--seed", 42);
  return [a, file];
}

// The hat of the page's clocks. Clocks of a file are read against its
// reference, which is one more clock unless the box says otherwise; three
// columns may hold the pairs instead. Simulated clocks have a perfect
// reference, which is none. `corners` is the number of clocks separated.
function hatRequest() {
  const H = state.ens.hat, clocks = fileClocks(true);
  const [a, file] = clockArgs("hat", clocks);
  const pairs = !!file && clocks.length === 3 && H.pairs;
  const reference = !!file && !pairs && (H.reference || clocks.length < 3);
  const corners = (clocks.length || 4) + reference;
  const method = H.method === "klts" && corners !== 3 ? "groslambert" : H.method;
  if (pairs) a.push("--pairs");
  else if (file && !reference) a.push("--no-reference");
  a.push(...(method === "hat" ? ["--devs", H.dev] : [method === "klts" ? "--klts" : "--covariance"]), "--grid", state.grid, "--confidence", state.confidence, "--format", "json");
  return { clocks, file: !!file, pairs, reference, corners, method, body: { args: a.map(String), file } };
}

// Two channels of one device: two columns of the page's file, read as the
// Data page reads them, gaps filled where it fills them. Chosen in the
// panel; by default the first two. `body` is null without a file.
function crossRequest() {
  const clocks = fileClocks(true).filter((c) => c.source.type === "file");
  const pick = (id, k) => clocks.find((c) => c.id === id) ?? clocks[k];
  const a = pick(state.ens.cross.a, 0);
  const b = clocks.length < 2 ? null : pick(state.ens.cross.b === a.id ? null : state.ens.cross.b, a === clocks[1] ? 0 : 1);
  if (!b) return { clocks, a, b, body: null };
  const [args, file] = clockArgs("cross", [a, b]);
  if (a.fillgaps || b.fillgaps) args.push("--fill-gaps");
  args.push("--grid", state.grid, "--format", "json");
  return { clocks, a, b, body: { args: args.map(String), file } };
}

function request() {
  if (state.ens.view === "hat") return hatRequest();
  if (state.ens.view === "cross") return crossRequest();
  const E = state.ens, clocks = fileClocks();
  const [a, file] = clockArgs("timescale", clocks);
  a.push("--algorithm", E.algorithm, ...(KALMAN.has(E.algorithm) ? ["--pivot", pivotOf(clocks.length || 4)] : []), ...(CAPPED.has(E.algorithm) ? ["--wmax", E.wmax] : []),
    ...(E.algorithm === "at2" && E.steps ? ["--steps", E.steps] : []), ...(KALMAN.has(E.algorithm) && simulated(clocks) && E.fitted ? ["--fitted-models"] : []),
    ...(FLICKER.has(E.algorithm) && E.flicker && fittedOf(clocks) ? ["--flicker"] : []), "--format", "json", "--series", 2000);
  return { clocks, body: { args: a.map(String), file } };
}

const views = () => `<div class="seg ens-view" role="group" aria-label="Analysis">${[["scale", "Time scale"], ["hat", "N-cornered hat"], ["cross", "Two channels"]].map(([id, label]) =>
  `<button type="button" data-ctl="view" data-v="${id}" class="${state.ens.view === id ? "on" : ""}">${label}</button>`).join("")}</div>`;

function hatControls() {
  const H = state.ens.hat, { clocks, file, pairs, corners, method } = hatRequest();
  const select = (ctl, options, value, off = () => false) => `<select data-ctl="${ctl}">${options.map(([id, label]) =>
    `<option value="${id}" ${id === value ? "selected" : ""} ${off(id) ? "disabled" : ""}>${label}</option>`).join("")}</select>`;
  return `${views()}
    <label class="opt"><span>Method</span>${select("method", HAT_METHODS, method, (id) => id === "klts" && corners !== 3)}</label>
    <label class="opt"><span>Deviation</span>${method === "hat" ? select("hat-dev", HAT_DEVS.map((d) => [d, byId[d].label]), H.dev) : `<select disabled><option>ADEV</option></select>`}</label>
    ${file && clocks.length === 3 ? `<label class="opt"><span>Columns hold</span>${select("pairs", [["", "Clocks"], ["on", "Pairs"]], pairs ? "on" : "")}</label>` : ""}
    ${file && clocks.length >= 3 && !pairs ? `<label class="check" title="The clocks are read against the file's reference, which is then one more clock of the hat"><input type="checkbox" data-ctl="reference" ${H.reference ? "checked" : ""}>The reference counts as a clock</label>` : ""}
    ${gridField("opt")}${confidenceField("opt")}
    <div class="rule"></div>
    <p class="muted small">${!clocks.length
      ? "No file with several clocks is loaded and none simulated, so four simulated clocks are used. Load a file with one column per clock, or simulate several clocks, on the Data page."
      : !file ? `${clocks.length} clocks simulated on the Data page. Their reference is perfect and counts as no clock.`
      : pairs ? `The three columns of ${h(clocks[0].source.name)} as the differences A−B, B−C and C−A of clocks A, B and C, each measured on its own.`
      : `${clocks.length} clocks from ${h(clocks[0].source.name)}, each read against the same reference${corners > clocks.length ? ", which is one more clock" : ""}.`}</p>
    <p class="muted small">${ABOUT_HAT[method]} The τ grid and confidence level are shared with the Stability page.</p>`;
}

function crossControls() {
  const { clocks, a, b } = crossRequest();
  const channel = (ctl, pick, other) => `<select data-ctl="${ctl}">${clocks.map((c) =>
    `<option value="${c.id}" ${c === pick ? "selected" : ""} ${c === other ? "disabled" : ""}>${h(c.name)}</option>`).join("")}</select>`;
  return `${views()}
    ${b ? `<label class="opt"><span>Channel A</span>${channel("cross-a", a, b)}</label>
    <label class="opt"><span>Channel B</span>${channel("cross-b", b, a)}</label>` : ""}
    ${gridField("opt")}
    <div class="rule"></div>
    <p class="muted small">${b ? `Two columns of ${h(a.source.name)}, each a measurement of one device: through two channels of an instrument, or against two references.`
      : "Load a file that holds the two channels as two columns on the Data page."}</p>
    <p class="muted small">The cross Allan deviation keeps the noise the two channels share, the device's, and averages toward zero the noise each has on its own. Its uncertainty follows the covariance model of Vernotte and Lantz (2019), with each channel's Allan EDF. The τ grid is shared with the Stability page.</p>`;
}

export function controls() {
  if (state.ens.view === "hat") return hatControls();
  if (state.ens.view === "cross") return crossControls();
  const E = state.ens, clocks = fileClocks();
  return `${views()}
    <label class="opt"><span>Algorithm</span><select data-ctl="alg">${ALGS.map(([id, label]) =>
      `<option value="${id}" ${E.algorithm === id ? "selected" : ""}>${label}</option>`).join("")}</select></label>
    ${KALMAN.has(E.algorithm)
      ? `<label class="opt"><span>Pivot clock</span><select data-ctl="pivot">${(clocks.length ? clocks.map((c) => c.name) : [1, 2, 3, 4].map((k) => `clock ${k}`)).map((n, i) =>
        `<option value="${i + 1}" ${pivotOf(clocks.length || 4) === i + 1 ? "selected" : ""}>${h(n)}</option>`).join("")}</select></label>` : ""}
    ${CAPPED.has(E.algorithm) ? `<label class="opt"><span>Largest weight</span><input type="text" data-ctl="wmax" value="${h(String(E.wmax))}" inputmode="decimal"></label>` : ""}
    ${E.algorithm === "at2" ? `<label class="opt"><span>Frequency steps</span><select data-ctl="steps">${[["", "Off"], ["best", "Best"], ["real", "Real"]].map(([id, label]) =>
      `<option value="${id}" ${E.steps === id ? "selected" : ""}>${label}</option>`).join("")}</select></label>` : ""}
    <label class="opt" title="Each clock's curve: the clock as read against the reference, which compares with the scale, or the clock against the scale"><span>Clocks against</span><select data-ctl="against">${AGAINST.map(([id, label]) =>
      `<option value="${id}" ${againstOf(clocks) === id ? "selected" : ""}>${label}</option>`).join("")}</select></label>
    ${KALMAN.has(E.algorithm) && simulated(clocks) ? `<label class="opt"><span>Models</span><select data-ctl="models">${MODELS.map(([id, label]) =>
      `<option value="${id}" ${E.fitted === (id === "fitted") ? "selected" : ""}>${label}</option>`).join("")}</select></label>` : ""}
    ${FLICKER.has(E.algorithm) && fittedOf(clocks) ? `<label class="opt"><span>Flicker states</span><select data-ctl="flicker">${[["", "Off"], ["on", "On"]].map(([id, label]) =>
      `<option value="${id}" ${E.flicker === (id === "on") ? "selected" : ""}>${label}</option>`).join("")}</select></label>` : ""}
    <div class="rule"></div>
    ${KALMAN.has(E.algorithm) && (clocks.length || 4) < 3 ? `<p class="small ens-warn">${ALGS.find(([id]) => id === E.algorithm)[1]} needs three or more clocks.</p>` : ""}
    <p class="muted small">${clocks.length
      ? (clocks[0].source.type === "synth"
        ? `${clocks.length} clocks simulated on the Data page, each against a perfect reference.`
        : `${clocks.length} clocks from ${h(clocks[0].source.name)}, each read against the same reference.`)
      : "No file with several clocks is loaded and none simulated, so four simulated clocks are used. Load a file with one column per clock, or simulate several clocks, on the Data page."}</p>
    <p class="muted small">${ABOUT[E.algorithm]}${KALMAN.has(E.algorithm) && simulated(clocks) ? ABOUT_MODELS[E.fitted] : ""}${FLICKER.has(E.algorithm) && fittedOf(clocks) ? ABOUT_FLICKER : ""}</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl;
  if (k === "view" && e.type === "click") update((s) => { s.ens.view = t.dataset.v; });
  if (e.type !== "change") return;
  onSharedControl(t);
  if (k === "reference") update((s) => { s.ens.hat.reference = t.checked; });
  const v = t.value.trim();
  if (k === "method") update((s) => { s.ens.hat.method = v; });
  if (k === "hat-dev") update((s) => { s.ens.hat.dev = v; });
  if (k === "pairs") update((s) => { s.ens.hat.pairs = v === "on"; });
  if (k === "cross-a" || k === "cross-b") update((s) => { s.ens.cross[k.slice(-1)] = Number(v); });
  if (k === "alg") update((s) => { s.ens.algorithm = v; });
  if (k === "steps") update((s) => { s.ens.steps = v; });
  if (k === "flicker") update((s) => { s.ens.flicker = v === "on"; });
  if (k === "models") update((s) => { s.ens.fitted = v === "fitted"; });
  if (k === "pivot") update((s) => { s.ens.pivot = Number(v); });
  if (k === "against") update((s) => { s.ens.against = v; });
  if (k === "wmax") {
    if (Number(v) > 0 && Number(v) <= 1) update((s) => { s.ens.wmax = Number(v); });
    else { t.value = state.ens.wmax; ctx.status("error", "The largest weight is a fraction between 0 and 1."); }
  }
}

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="card plot-card"><div class="card-row"><h2 data-ref="title"></h2><span class="muted small" data-ref="sub"></span></div>
        <div class="plot" data-ref="adev"></div><div class="empty" data-ref="empty" hidden></div></div>
      <div class="card plot-card ens-weights"><div class="card-row"><h2 data-ref="title2">Weights</h2></div><div class="plot" data-ref="weights"></div></div>
    </div>
    <div class="card side-card ens-table" data-ref="table"></div>`;
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

// Ask the core for the scale or the hat when it is missing or out of date.
function ensure() {
  const { body, clocks } = request();
  if (!body) return;
  const key = JSON.stringify(body);
  if (state.ens.result?.key === key || state.ens.error?.key === key || inflight === key || !state.core.ok) return;
  inflightCtl?.abort();
  const ctl = inflightCtl = new AbortController();
  inflight = key;
  state.ens.error = null;
  ensemble(body, ctl.signal).then((doc) => {
    if (inflight !== key) return;
    inflight = null;
    update((s) => { s.ens.result = { key, ts: doc.timescale, hat: doc.hat, cross: doc.cross_adev && { adev: doc.cross_adev, spec: doc.cross_spectrum, record: doc.record } }; });
  }).catch((err) => {
    if (err.name === "AbortError" || inflight !== key) return;
    inflight = null;
    update((s) => { s.ens.error = { key, message: inWords(err.message, simulated(clocks)) }; s.ens.result = null; });
  });
}

const weight = (w) => (w === null ? "–" : num(w, 3));

// Against an ideal reference: the scale beside its most stable clock at the
// shortest τ, and what the best fixed weights would give independent clocks,
// σ⁻² = Σ σᵢ⁻² (Weiss and Weissert (1991), Eq. (4) and Appendix A); 1/√N for N equal clocks.
function gain(r) {
  const own = r.adev.readings.map((dev) => dev[0]).filter((v) => v > 0), best = Math.min(...own);
  if (!own.length || !(r.adev.scale[0] > 0)) return "";
  const bound = 1 / Math.sqrt(own.reduce((sum, v) => sum + (best / v) ** 2, 0));
  return `At τ = ${num(r.adev.tau[0], 4)} s the scale is ${num(r.adev.scale[0] / best, 2)} of its most stable clock; the best fixed weights for these clocks give ${num(bound, 2)}. `;
}

// Without a scale the plots are hidden, and its table and the line that describes it go.
function showEmpty(text) {
  const e = root.querySelector("[data-ref=empty]");
  e.innerHTML = text ? `<p>${text}</p>` : "";
  e.hidden = !text;
  for (const p of root.querySelectorAll(".plot")) p.style.visibility = text ? "hidden" : "visible";
  if (text) for (const ref of ["table", "sub"]) root.querySelector(`[data-ref=${ref}]`).innerHTML = "";
}

// The core names the file's clocks of a hat by column and its reference
// "reference"; the page names them as the clock list does.
const hatNames = (r, { clocks, pairs }) => r.clocks.map((n, i) => (pairs ? n : clocks[i]?.name ?? (n === "reference" ? clocks[0]?.against || "Reference" : n)));
const scaleNames = (r, { clocks }) => r.clocks.map((n, i) => clocks[i]?.name ?? n);

// The hat: each clock's σ(τ) with its interval, over the pairs as measured.
function refreshHat(r, req) {
  const { clocks, pairs } = req;
  const e = byId[r.results[0].deviation], res = r.results[0];
  const col = (c, k) => c.points.map((p) => p[k]);
  const names = hatNames(r, req);
  const color = (i) => (pairs || clocks[i] ? seriesColor(pairs ? i : clocks[i].slot) : clocks.length ? css("--ink-2") : seriesColor(i));
  const sigma = e.axis === "sigma_x" ? "σ<sub>x</sub>(τ) (s)" : "σ<sub>y</sub>(τ)";
  root.querySelector("[data-ref=title]").textContent = `${e.label} of each clock`;
  root.querySelector("[data-ref=sub]").textContent = `${HAT_METHODS.find(([id]) => id === (r.interval === "klts" ? "klts" : r.method))[1]} · ${r.clocks.length} clocks · ${r.n} samples · τ₀ = ${interval(r.tau0)} · ${num(100 * r.confidence, 3)} % intervals`;
  root.querySelector(".ens-weights").hidden = true;
  const tau = col(res.clocks[0], "tau"), resolved = (c) => c.points.filter((p) => p.dev !== null).length;
  const first = (p) => (p.dev === null && p.ci_upper !== null ? `&lt; ${sci(p.ci_upper, 3)}` : sci(p.dev, 3));
  root.querySelector("[data-ref=table]").innerHTML = `<div class="card-row"><h2>Clocks</h2></div>
    <div class="tbl hat-tbl"><div class="tbl-row tbl-head"><span>Clock</span><span>${e.label} at ${num(tau[0], 4)} s</span><span>Resolved</span></div>
    ${res.clocks.map((c, i) => `<div class="tbl-row tbl-body"><span class="ens-name"><span class="clock-dot" style="background:${color(i)}"></span>${h(names[i])}</span>
      <span>${first(c.points[0])}</span><span>${resolved(c)} of ${tau.length} τ</span></div>`).join("")}</div>
    <p class="muted small ens-note">Resolved counts the averaging times at which the estimate of the clock's variance is positive. At the others the clock is not resolved under the noise of the rest, and it has no point${r.interval === "klts" ? ": ▽ and &lt; give its upper limit" : ""}.</p>`;
  if (!Plotly) return;
  const bar = (c) => ({ type: "data", symmetric: false, thickness: 1, width: 3,
    array: c.points.map((p) => (p.dev !== null && p.ci_upper !== null ? p.ci_upper - p.dev : 0)),
    arrayminus: c.points.map((p) => (p.dev !== null && p.ci_lower ? p.dev - p.ci_lower : 0)) });
  const span = (p) => (p.ci_upper === null ? "" : `<br>${sci(p.ci_lower, 3)} to ${sci(p.ci_upper, 3)}`);
  Plotly.react(root.querySelector("[data-ref=adev]"), [
    ...res.pairs.map((c, k) => ({ x: tau, y: col(c, "dev"), name: "Pairs, as measured", legendgroup: "pairs", showlegend: k === 0, mode: "lines",
      line: { color: css("--muted"), width: 1, dash: "dot" },
      hovertemplate: `<b>${h(r.pairs[k].map((i) => names[i - 1]).join(" − "))}</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>` })),
    ...res.clocks.map((c, i) => ({ x: tau, y: col(c, "dev"), name: names[i], mode: "lines+markers", line: { color: color(i), width: LINE_W }, marker: { color: color(i), size: 6 },
      error_y: bar(c), text: c.points.map(span), hovertemplate: `<b>${h(names[i])}</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}%{text}<extra></extra>` })),
    // KLTS gives an upper limit where the clock is not resolved.
    ...res.clocks.map((c, i) => ({ x: tau, y: c.points.map((p) => (p.dev === null ? p.ci_upper : null)), name: `${names[i]}, upper limit`, showlegend: false, mode: "markers",
      marker: { color: color(i), size: 9, symbol: "triangle-down-open" }, hovertemplate: `<b>${h(names[i])}</b>, not resolved<br>τ = %{x:.4g} s<br>σ < %{y:.3e}<extra></extra>` })),
  ], { ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis(sigma), showlegend: true,
    legend: { ...baseLayout().legend, orientation: "h", x: 0, xanchor: "left", y: 1.02, yanchor: "bottom", bgcolor: "rgba(0,0,0,0)", borderwidth: 0 },   // above the curves, not over them
    margin: { l: 62, r: 14, t: 40, b: 48 }, uirevision: `hat-${res.deviation}` }, CONFIG);
  Plotly.Plots.resize(root.querySelector("[data-ref=adev]"));
}

// Two channels: the cross Allan deviation with its uncertainty over each
// channel's own ADEV, and the cross-spectrum over its floor and each
// channel's spectrum. Where the covariance is not positive, or not above its
// uncertainty, the device is not resolved: ▽ marks √(covariance + uncertainty).
function refreshCross(r, { a, b }) {
  const x = r.adev, sp = r.spec, ca = seriesColor(a.slot), cb = seriesColor(b.slot), dev = css("--accent"), ink = css("--ink-2");
  const resolved = x.covar.filter((v, i) => v > x.uncertainty[i]).length;
  root.querySelector("[data-ref=title]").textContent = "Cross Allan deviation of the two channels";
  root.querySelector("[data-ref=sub]").textContent = `${a.name} × ${b.name} · ${r.record.n} samples · τ₀ = ${num(r.record.tau0, 4)} s`;
  const wcard = root.querySelector(".ens-weights"), relaid = wcard.hidden;
  wcard.hidden = false;
  root.querySelector("[data-ref=title2]").innerHTML = "Cross-spectrum S<sub>y</sub>(f)";
  const row = (i) => `<div class="tbl-row tbl-body"><span>${num(x.tau[i], 4)}</span><span>${x.covar[i] > 0 ? sci(x.dev[i], 3) : "–"}</span><span>${sci(x.dev_a[i], 3)}</span><span>${sci(x.dev_b[i], 3)}</span></div>`;
  root.querySelector("[data-ref=table]").innerHTML = `<div class="card-row"><h2>Per-τ results</h2></div>
    <div class="tbl cross-tbl"><div class="tbl-row tbl-head"><span>τ (s)</span><span>Cross</span><span>A</span><span>B</span></div>${x.tau.map((_, i) => row(i)).join("")}</div>
    <p class="muted small ens-note">The device is resolved at ${resolved} of ${x.tau.length} τ, where the covariance is above its uncertainty. Elsewhere the channels' own noise dominates. The cross-spectrum averages ${sp.nseg} segments of ${sp.nperseg} samples; its floor is the residual of the noise the channels do not share.</p>`;
  if (!Plotly) return;
  const sq = (v) => (v > 0 ? Math.sqrt(v) : null);
  const ok = x.covar.map((v, i) => v > x.uncertainty[i]);
  const hover = (name) => `<b>${h(name)}</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>`;
  Plotly.react(root.querySelector("[data-ref=adev]"), [
    { x: x.tau, y: x.dev_a, name: a.name, mode: "lines", line: { color: ca, width: 1.5, dash: "dot" }, hovertemplate: hover(`${a.name} alone`) },
    { x: x.tau, y: x.dev_b, name: b.name, mode: "lines", line: { color: cb, width: 1.5, dash: "dot" }, hovertemplate: hover(`${b.name} alone`) },
    { x: x.tau, y: x.dev.map((v, i) => (ok[i] ? v : null)), name: "Cross", mode: "lines+markers", line: { color: dev, width: 2.5 }, marker: { color: dev, size: 6 },
      error_y: { type: "data", symmetric: false, thickness: 1, width: 3,
        array: x.covar.map((v, i) => (ok[i] ? sq(v + x.uncertainty[i]) - x.dev[i] : 0)),
        arrayminus: x.covar.map((v, i) => (ok[i] ? x.dev[i] - sq(v - x.uncertainty[i]) : 0)) },
      hovertemplate: hover("Cross") },
    { x: x.tau, y: x.covar.map((v, i) => (ok[i] ? null : sq(v + x.uncertainty[i]))), name: "Cross, upper limit", showlegend: false, mode: "markers",
      marker: { color: dev, size: 9, symbol: "triangle-down-open" }, hovertemplate: "<b>Cross</b>, not resolved<br>τ = %{x:.4g} s<br>σ < %{y:.3e}<extra></extra>" },
  ], { ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis("σ<sub>y</sub>(τ)"), showlegend: true,
    legend: { ...baseLayout().legend, orientation: "h", x: 0, xanchor: "left", y: 1.02, yanchor: "bottom", bgcolor: "rgba(0,0,0,0)", borderwidth: 0 },
    margin: { l: 62, r: 14, t: 40, b: 48 }, uirevision: "cross" }, CONFIG);
  const shover = (name) => `<b>${h(name)}</b><br>f = %{x:.4g} Hz<br>S = %{y:.3e} /Hz<extra></extra>`;
  Plotly.react(root.querySelector("[data-ref=weights]"), [
    { x: sp.f, y: sp.sy_a, name: a.name, mode: "lines", line: { color: ca, width: 1 }, opacity: 0.6, hovertemplate: shover(`${a.name} alone`) },
    { x: sp.f, y: sp.sy_b, name: b.name, mode: "lines", line: { color: cb, width: 1 }, opacity: 0.6, hovertemplate: shover(`${b.name} alone`) },
    { x: sp.f, y: sp.floor, name: "Floor", mode: "lines", line: { color: ink, width: 1, dash: "dash" }, hovertemplate: shover("Floor") },
    // The real part of the cross-spectrum; below zero it has no point on a log axis.
    { x: sp.f, y: sp.sy.map((v) => (v > 0 ? v : null)), name: "Cross", mode: "lines", line: { color: dev, width: 1.5 }, hovertemplate: shover("Cross") },
  ], { ...baseLayout(), margin: { l: 62, r: 14, t: 8, b: 44 }, xaxis: logAxis("Fourier frequency f (Hz)"), yaxis: logAxis("S<sub>y</sub>(f) (1/Hz)"), showlegend: false, uirevision: "cross-s" }, CONFIG);
  if (relaid) for (const p of root.querySelectorAll(".plot")) Plotly.Plots.resize(p);
}

export function refresh() {
  if (!root) return;
  ensure();
  const req = request(), { clocks } = req;
  const hat = state.ens.view === "hat", cross = state.ens.view === "cross";
  const r = hat ? state.ens.result?.hat : cross ? state.ens.result?.cross : state.ens.result?.ts;
  const table = root.querySelector("[data-ref=table]");
  if (cross && !req.body) return showEmpty("Two channels needs a file that holds them as two columns, such as the readings of one device against two references. Load one on the Data page.");
  if (state.ens.error) return showEmpty(h(state.ens.error.message));
  if (!r) return showEmpty(!state.core.ok ? "No analysis core. Start the explorer with serve.py." : hat ? "Separating the clocks…" : cross ? "Computing the cross estimates…" : "Forming the time scale…");
  showEmpty("");
  if (hat) return refreshHat(r, req);
  if (cross) return refreshCross(r, req);
  root.querySelector("[data-ref=title2]").textContent = "Weights";
  const against = againstOf(clocks), ideal = clocks[0]?.source.type !== "file";
  root.querySelector("[data-ref=title]").textContent = `ADEV of the scale and of each clock against ${against === "scale" ? "it" : "the reference"}`;
  const names = scaleNames(r, req);
  const color = (i) => seriesColor(clocks[i]?.slot ?? i);
  // Jones–Tryon has no weights: its plot goes and the ADEV plot takes the room.
  const weighted = r.weight.some((w) => w !== null);
  const wcard = root.querySelector(".ens-weights"), relaid = wcard.hidden === weighted;
  wcard.hidden = !weighted;
  const alg = ALGS.find(([id]) => id === r.algorithm)?.[1] ?? r.algorithm;
  root.querySelector("[data-ref=sub]").textContent = `${alg}${r.pivot ? ` at ${names[r.pivot - 1]}` : ""}${r.steps ? `, ${r.steps} steps` : ""}${r.flicker ? ", flicker states" : ""} · ${r.n} samples · τ₀ = ${interval(r.tau0)} · against ${r.reference || clocks[0]?.against || "the file's reference"}`;
  // With flicker states on: the Markov states each clock's fitted model carries.
  const states = (i) => [[r.flicker.fm[i], "FM"], [r.flicker.pm[i], "PM"]].filter(([k]) => k).map(([k, kind]) => `${k} ${kind}`).join(" + ");
  const carried = r.flicker ? names.map((n, i) => [n, states(i)]).filter(([, s]) => s) : [];
  const flicker = !r.flicker ? "" : carried.length ? `Flicker states: ${carried.map(([n, s]) => `${h(n)} (${s})`).join(", ")}. ` : "No clock's fit has flicker that needs states. ";
  // A Kalman scale shows the white FM of each clock's model; a fit without any takes nearly all of the weight.
  const fitted = r.models, bare = (r.fitted ? fitted : []).flatMap((m, i) => (m.q1 === 0 && m.flicker_fm === 0 ? [names[i]] : []));
  const unresolved = bare.length ? `The fit finds no white FM in ${bare.map(h).join(", ")}: the record does not resolve ${bare.length > 1 ? "them" : "it"}, and the scale follows ${bare.length > 1 ? "them" : "it"}. AT1 takes no model. ` : "";
  table.innerHTML = `<div class="card-row"><h2>Clocks</h2></div>
    <div class="tbl${fitted ? " fitted" : ""}"><div class="tbl-row tbl-head"><span>Clock</span><span>Weight</span><span>Mean</span><span>RMS error</span>${fitted ? "<span>White FM</span>" : ""}</div>
    ${names.map((n, i) => `<div class="tbl-row tbl-body"><span class="ens-name"><span class="clock-dot" style="background:${color(i)}"></span>${h(n)}</span>
      <span>${weight(r.weight[i])}</span><span>${weight(r.weight_mean[i])}</span><span>${sci(r.error_rms[i], 2)} s</span>${fitted ? `<span>${sci(fitted[i].q1, 2)}</span>` : ""}</div>`).join("")}</div>
    <p class="muted small ens-note">${ideal ? gain(r) : ""}${unresolved}${flicker}${weighted ? "" : "The Jones–Tryon scale has no per-clock weights. "}Mean weight and RMS prediction error cover the second half of the record.</p>`;
  if (!Plotly) return;
  const scale = css("--accent");
  Plotly.react(root.querySelector("[data-ref=adev]"), [
    { x: r.adev.tau, y: r.adev.scale, name: "Time scale", mode: "lines+markers", line: { color: scale, width: 3 }, marker: { color: scale, size: 6 },
      hovertemplate: "<b>Time scale</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>" },
    ...r.adev[against === "scale" ? "clocks" : "readings"].map((dev, i) => ({ x: r.adev.tau, y: dev, name: names[i], mode: "lines", line: { color: color(i), width: 1.5, dash: "dot" },
      hovertemplate: `<b>${h(names[i])}</b>${against === "scale" ? " − scale" : ""}<br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>` })),
  ], { ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis("σ<sub>y</sub>(τ)"), showlegend: false, uirevision: "ens" }, CONFIG);
  if (relaid) Plotly.Plots.resize(root.querySelector("[data-ref=adev]"));
  const s = r.series;
  const [tf, tu] = timeAxisUnit(s?.t.at(-1));
  Plotly.react(root.querySelector("[data-ref=weights]"), s && weighted ? s.weights.map((w, i) => ({
    x: s.t.map((t) => t / tf), y: w, name: names[i], mode: "lines", line: { color: color(i), width: 1.5 },
    hovertemplate: `<b>${h(names[i])}</b><br>t = %{x:.6g} ${tu}<br>w = %{y:.3f}<extra></extra>` })) : [],
  { ...baseLayout(), margin: { l: 62, r: 14, t: 8, b: 44 }, xaxis: linAxis(`Time (${tu})`), yaxis: { ...linAxis("Weight"), range: [0, 1] }, showlegend: false, uirevision: "ens-w" }, CONFIG);
}

// The File menu's results, of the view shown: the hat's σ(τ) of each clock
// separated, or the scale's ADEV and each clock's against the scale or the
// reference as plotted, with the clock's weights. The JSON is the core's, with
// the clocks named as on the page.
export function results() {
  const req = request(), r = state.ens.result, hat = state.ens.view === "hat";
  const doc = r?.key === JSON.stringify(req.body) ? (hat ? r.hat : r.ts) : null;
  const base = hat
    ? { what: "Hat results", file: "hat", columns: ["clock", "estimator", "tau_s", "sigma", "ci_lo", "ci_hi", "edf", "neff"] }
    : { what: "Time scale table", file: "timescale", columns: ["clock", "against", "tau_s", "adev", "weight", "weight_mean", "error_rms_s"] };
  if (!doc) return { ...base, rows: [], json: null };
  if (hat) {
    const names = hatNames(doc, req);
    return { ...base, json: { ...doc, clocks: names },
      rows: doc.results.flatMap((res) => res.clocks.flatMap((c, i) => c.points.map((p) => [names[i], res.deviation, p.tau, p.dev, p.ci_lower, p.ci_upper, p.edf, p.neff]))) };
  }
  const names = scaleNames(doc, req), a = doc.adev, against = againstOf(req.clocks), dev = against === "scale" ? a.clocks : a.readings;
  return { ...base, json: { ...doc, clocks: names },
    rows: [...a.tau.map((t, j) => ["time scale", "reference", t, a.scale[j]]),
      ...names.flatMap((n, i) => a.tau.map((t, j) => [n, against, t, dev[i][j], doc.weight[i], doc.weight_mean[i], doc.error_rms[i]]))] };
}

export function hint() {
  if (state.ens.view === "cross") return "Solid = the cross deviation, bars = its uncertainty, ▽ = upper limit · dotted = each channel alone";
  return state.ens.view === "hat" ? "Solid = each clock, separated · dotted = the pairs as measured · bars = the confidence interval"
    : `Solid = the scale against the reference · dotted = each clock against the ${againstOf(fileClocks())}`;
}

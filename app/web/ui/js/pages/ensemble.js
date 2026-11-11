// Ensemble page: several clocks read against one reference. One view forms a
// time scale from them (AT1, AT2, EWFA or a Kalman scale); the other separates
// the stability of each clock by an N-cornered hat or the Groslambert
// covariance. The clocks are the value columns of one loaded file; without
// such a file the core simulates four. Either is computed when this page
// needs it and recomputed when the settings or the clocks change.
import { state, update, isReference } from "../state.js";
import { ensemble, synthArgs } from "../core.js";
import { loadPlotly, baseLayout, logAxis, linAxis, CONFIG, seriesColor, css } from "../plot.js";
import { byId } from "../estimators.js";
import { sci, num, escapeHtml as h } from "../format.js";

let root, ctx, Plotly;
let inflight = null;           // key of the run in progress
let inflightCtl = null;        // stops it (the core process too) when the settings change
const ALGS = [["at1", "AT1"], ["at2", "AT2"], ["ewfa", "EWFA"], ["greenhall", "Greenhall"], ["jones_tryon", "Jones–Tryon"], ["stein", "Stein"]];
const KALMAN = new Set(["greenhall", "jones_tryon", "stein"]);
const CAPPED = new Set(["at1", "at2", "stein"]);                // take --wmax
const FLICKER = new Set(["greenhall", "jones_tryon"]);          // take --flicker
const ABOUT_FLICKER = " Flicker states: on keeps each clock's flicker FM and flicker PM in its model as Markov states, across the τ where its fit finds them; off folds flicker PM into the reading noise and leaves flicker FM out.";
const ABOUT = {
  at1: "AT1 weighs each clock by the inverse of its prediction-error variance, learned as it runs.",
  at2: "AT2 keeps AT1's weights and follows each clock's frequency with a Kalman update. Frequency steps: best reruns the scale from a detected step; real keeps the values already issued and corrects from the step on.",
  ewfa: "EWFA weighs every clock the same.",
  greenhall: "Greenhall's Kalman scale fits a model to each clock and reads the scale at the pivot clock. It needs three or more clocks.",
  jones_tryon: "The Jones–Tryon Kalman scale: Greenhall's filter without the covariance reduction, so it has no per-clock weights. It needs three or more clocks.",
  stein: "Stein's Kalman scale (KAS-2) fits a model to each clock, filters each clock's difference from the pivot, and weighs each clock's phase, frequency and aging by the inverse of its noise. The weights shown are the phase weights. It needs three or more clocks.",
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
// The file goes to the core as read, without the preparation of the Data page,
// so the core's remedy for gaps (--fillgaps) has no counterpart on this page.
const GAPS = [
  [/^clock_models: .* usable τ points/, "Greenhall, Jones–Tryon and Stein fit a model to each clock, and this file leaves one of them too few samples at which it is read with two others. AT1, AT2 and EWFA form a scale from it."],
  [/^the KLTS interval needs/, "The KLTS interval needs a file without missing samples. The hat and the Groslambert covariance separate the clocks of this one."],
  [/^the frequency record is missing/, "Frequency columns with missing samples cannot be read."],
];
const FEW = [/^hat needs three or more simulated/, "Simulated clocks have a perfect reference, which counts as no clock. Simulate three or more, on the Data or the Synthesis page."];
const inWords = (message) => GAPS.find(([re]) => re.test(message))?.[1].concat(" The gap filling of the Data page does not apply here.")
  ?? (FEW[0].test(message) ? FEW[1] : message);

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
    a.push("--value-cols", clocks.map((k) => k.value_col).join(","), "--time-col", c.time_col, "--header", c.header, "--scaling", c.scaling);
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

function request() {
  if (state.ens.view === "hat") return hatRequest();
  const E = state.ens, clocks = fileClocks();
  const [a, file] = clockArgs("timescale", clocks);
  a.push("--algorithm", E.algorithm, ...(KALMAN.has(E.algorithm) ? ["--pivot", pivotOf(clocks.length || 4)] : []), ...(CAPPED.has(E.algorithm) ? ["--wmax", E.wmax] : []),
    ...(E.algorithm === "at2" && E.steps ? ["--steps", E.steps] : []), ...(FLICKER.has(E.algorithm) && E.flicker ? ["--flicker"] : []), "--format", "json", "--series", 2000);
  return { clocks, body: { args: a.map(String), file } };
}

const views = () => `<div class="seg ens-view" role="group" aria-label="Analysis">${[["scale", "Time scale"], ["hat", "N-cornered hat"]].map(([id, label]) =>
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
    <div class="rule"></div>
    <p class="muted small">${!clocks.length
      ? "No file with several clocks is loaded and none simulated, so four simulated clocks are used. Load a file with one column per clock, or simulate several clocks, on the Data page."
      : !file ? `${clocks.length} clocks simulated on the Data page. Their reference is perfect and counts as no clock.`
      : pairs ? `The three columns of ${h(clocks[0].source.name)} as the differences A−B, B−C and C−A of clocks A, B and C, each measured on its own.`
      : `${clocks.length} clocks from ${h(clocks[0].source.name)}, each read against the same reference${corners > clocks.length ? ", which is one more clock" : ""}.`}</p>
    <p class="muted small">${ABOUT_HAT[method]} Grid and confidence level are those of the Stability page.</p>`;
}

export function controls() {
  if (state.ens.view === "hat") return hatControls();
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
    ${FLICKER.has(E.algorithm) ? `<label class="opt"><span>Flicker states</span><select data-ctl="flicker">${[["", "Off"], ["on", "On"]].map(([id, label]) =>
      `<option value="${id}" ${E.flicker === (id === "on") ? "selected" : ""}>${label}</option>`).join("")}</select></label>` : ""}
    <div class="rule"></div>
    ${KALMAN.has(E.algorithm) && (clocks.length || 4) < 3 ? `<p class="small ens-warn">${ALGS.find(([id]) => id === E.algorithm)[1]} needs three or more clocks.</p>` : ""}
    <p class="muted small">${clocks.length
      ? (clocks[0].source.type === "synth"
        ? `${clocks.length} clocks simulated on the Data page, each against a perfect reference.`
        : `${clocks.length} clocks from ${h(clocks[0].source.name)}, each read against the same reference.`)
      : "No file with several clocks is loaded and none simulated, so four simulated clocks are used. Load a file with one column per clock, or simulate several clocks, on the Data page."}</p>
    <p class="muted small">${ABOUT[E.algorithm]}${FLICKER.has(E.algorithm) ? ABOUT_FLICKER : ""}</p>`;
}

export function onControl(e) {
  const t = e.target.closest("[data-ctl]");
  if (!t) return;
  const k = t.dataset.ctl;
  if (k === "view" && e.type === "click") update((s) => { s.ens.view = t.dataset.v; });
  if (e.type !== "change") return;
  if (k === "reference") update((s) => { s.ens.hat.reference = t.checked; });
  const v = t.value.trim();
  if (k === "method") update((s) => { s.ens.hat.method = v; });
  if (k === "hat-dev") update((s) => { s.ens.hat.dev = v; });
  if (k === "pairs") update((s) => { s.ens.hat.pairs = v === "on"; });
  if (k === "alg") update((s) => { s.ens.algorithm = v; });
  if (k === "steps") update((s) => { s.ens.steps = v; });
  if (k === "flicker") update((s) => { s.ens.flicker = v === "on"; });
  if (k === "pivot") update((s) => { s.ens.pivot = Number(v); });
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
      <div class="card plot-card ens-weights"><div class="card-row"><h2>Weights</h2></div><div class="plot" data-ref="weights"></div></div>
    </div>
    <div class="card side-card ens-table" data-ref="table"></div>`;
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

// Ask the core for the scale or the hat when it is missing or out of date.
function ensure() {
  const { body } = request();
  const key = JSON.stringify(body);
  if (state.ens.result?.key === key || state.ens.error?.key === key || inflight === key || !state.core.ok) return;
  inflightCtl?.abort();
  const ctl = inflightCtl = new AbortController();
  inflight = key;
  state.ens.error = null;
  ensemble(body, ctl.signal).then((doc) => {
    if (inflight !== key) return;
    inflight = null;
    update((s) => { s.ens.result = { key, ts: doc.timescale, hat: doc.hat }; });
  }).catch((err) => {
    if (err.name === "AbortError" || inflight !== key) return;
    inflight = null;
    update((s) => { s.ens.error = { key, message: inWords(err.message) }; s.ens.result = null; });
  });
}

const weight = (w) => (w === null ? "–" : num(w, 3));

// Without a scale the plots are hidden, and its table and the line that describes it go.
function showEmpty(text) {
  const e = root.querySelector("[data-ref=empty]");
  e.innerHTML = text ? `<p>${text}</p>` : "";
  e.hidden = !text;
  for (const p of root.querySelectorAll(".plot")) p.style.visibility = text ? "hidden" : "visible";
  if (text) for (const ref of ["table", "sub"]) root.querySelector(`[data-ref=${ref}]`).innerHTML = "";
}

// The hat: each clock's σ(τ) with its interval, over the pairs as measured.
function refreshHat(r, { clocks, pairs }) {
  const e = byId[r.results[0].deviation], res = r.results[0];
  const col = (c, k) => c.points.map((p) => p[k]);
  // The core names the file's clocks by column and its reference "reference".
  const names = r.clocks.map((n, i) => (pairs ? n : clocks[i]?.name ?? (n === "reference" ? clocks[0]?.against || "Reference" : n)));
  const color = (i) => (pairs || clocks[i] ? seriesColor(pairs ? i : clocks[i].slot) : clocks.length ? css("--ink-2") : seriesColor(i));
  const sigma = e.axis === "sigma_x" ? "σ<sub>x</sub>(τ) (s)" : "σ<sub>y</sub>(τ)";
  root.querySelector("[data-ref=title]").textContent = `${e.label} of each clock`;
  root.querySelector("[data-ref=sub]").textContent = `${HAT_METHODS.find(([id]) => id === (r.interval === "klts" ? "klts" : r.method))[1]} · ${r.clocks.length} clocks · ${r.n} samples · τ₀ = ${num(r.tau0, 4)} s · ${num(100 * r.confidence, 3)} % intervals`;
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
    ...res.clocks.map((c, i) => ({ x: tau, y: col(c, "dev"), name: names[i], mode: "lines+markers", line: { color: color(i), width: 2 }, marker: { color: color(i), size: 6 },
      error_y: bar(c), text: c.points.map(span), hovertemplate: `<b>${h(names[i])}</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}%{text}<extra></extra>` })),
    // KLTS gives an upper limit where the clock is not resolved.
    ...res.clocks.map((c, i) => ({ x: tau, y: c.points.map((p) => (p.dev === null ? p.ci_upper : null)), name: `${names[i]}, upper limit`, showlegend: false, mode: "markers",
      marker: { color: color(i), size: 9, symbol: "triangle-down-open" }, hovertemplate: `<b>${h(names[i])}</b>, not resolved<br>τ = %{x:.4g} s<br>σ < %{y:.3e}<extra></extra>` })),
  ], { ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis(sigma), showlegend: true, uirevision: `hat-${res.deviation}` }, CONFIG);
  Plotly.Plots.resize(root.querySelector("[data-ref=adev]"));
}

export function refresh() {
  if (!root) return;
  ensure();
  const req = request(), { clocks } = req;
  const hat = state.ens.view === "hat";
  const r = hat ? state.ens.result?.hat : state.ens.result?.ts;
  const table = root.querySelector("[data-ref=table]");
  if (state.ens.error) return showEmpty(h(state.ens.error.message));
  if (!r) return showEmpty(!state.core.ok ? "No analysis core. Start the explorer with serve.py." : hat ? "Separating the clocks…" : "Forming the time scale…");
  showEmpty("");
  if (hat) return refreshHat(r, req);
  root.querySelector("[data-ref=title]").textContent = "ADEV of the scale and of each clock against it";
  const names = r.clocks.map((n, i) => clocks[i]?.name ?? n);
  const color = (i) => seriesColor(clocks[i]?.slot ?? i);
  // Jones–Tryon has no weights: its plot goes and the ADEV plot takes the room.
  const weighted = r.weight.some((w) => w !== null);
  const wcard = root.querySelector(".ens-weights"), relaid = wcard.hidden === weighted;
  wcard.hidden = !weighted;
  const alg = ALGS.find(([id]) => id === r.algorithm)?.[1] ?? r.algorithm;
  root.querySelector("[data-ref=sub]").textContent = `${alg}${r.pivot ? ` at ${names[r.pivot - 1]}` : ""}${r.steps ? `, ${r.steps} steps` : ""}${r.flicker ? ", flicker states" : ""} · ${r.n} samples · τ₀ = ${num(r.tau0, 4)} s · against ${r.reference || "the file's reference"}`;
  // With flicker states on: the Markov states each clock's fitted model carries.
  const states = (i) => [[r.flicker.fm[i], "FM"], [r.flicker.pm[i], "PM"]].filter(([k]) => k).map(([k, kind]) => `${k} ${kind}`).join(" + ");
  const carried = r.flicker ? names.map((n, i) => [n, states(i)]).filter(([, s]) => s) : [];
  const flicker = !r.flicker ? "" : carried.length ? `Flicker states: ${carried.map(([n, s]) => `${h(n)} (${s})`).join(", ")}. ` : "No clock's fit has flicker that needs states. ";
  table.innerHTML = `<div class="card-row"><h2>Clocks</h2></div>
    <div class="tbl"><div class="tbl-row tbl-head"><span>Clock</span><span>Weight</span><span>Mean</span><span>RMS error</span></div>
    ${names.map((n, i) => `<div class="tbl-row tbl-body"><span class="ens-name"><span class="clock-dot" style="background:${color(i)}"></span>${h(n)}</span>
      <span>${weight(r.weight[i])}</span><span>${weight(r.weight_mean[i])}</span><span>${sci(r.error_rms[i], 2)} s</span></div>`).join("")}</div>
    <p class="muted small ens-note">${flicker}${weighted ? "" : "The Jones–Tryon scale has no per-clock weights. "}Mean weight and RMS prediction error cover the second half of the record.</p>`;
  if (!Plotly) return;
  const scale = css("--accent");
  Plotly.react(root.querySelector("[data-ref=adev]"), [
    { x: r.adev.tau, y: r.adev.scale, name: "Time scale", mode: "lines+markers", line: { color: scale, width: 3 }, marker: { color: scale, size: 6 },
      hovertemplate: "<b>Time scale</b><br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>" },
    ...r.adev.clocks.map((dev, i) => ({ x: r.adev.tau, y: dev, name: names[i], mode: "lines", line: { color: color(i), width: 1.5, dash: "dot" },
      hovertemplate: `<b>${h(names[i])}</b> − scale<br>τ = %{x:.4g} s<br>σ = %{y:.3e}<extra></extra>` })),
  ], { ...baseLayout(), xaxis: logAxis("Averaging time τ (s)"), yaxis: logAxis("σ<sub>y</sub>(τ)"), showlegend: true, uirevision: "ens" }, CONFIG);
  if (relaid) Plotly.Plots.resize(root.querySelector("[data-ref=adev]"));
  const s = r.series;
  Plotly.react(root.querySelector("[data-ref=weights]"), s && weighted ? s.weights.map((w, i) => ({
    x: s.t, y: w, name: names[i], mode: "lines", line: { color: color(i), width: 1.5 },
    hovertemplate: `<b>${h(names[i])}</b><br>t = %{x:.6g} s<br>w = %{y:.3f}<extra></extra>` })) : [],
  { ...baseLayout(), margin: { l: 62, r: 14, t: 8, b: 44 }, xaxis: linAxis("Time (s)"), yaxis: { ...linAxis("Weight"), range: [0, 1] }, showlegend: false, uirevision: "ens-w" }, CONFIG);
}

export function hint() {
  return state.ens.view === "hat" ? "Solid = each clock, separated · dotted = the pairs as measured · bars = the confidence interval"
    : "Solid = the scale against the reference · dotted = each clock against the scale";
}

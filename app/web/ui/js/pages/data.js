// Data page: load a file (with the import dialog) or synthesize a record,
// the record check, and the phase plot with outliers in orange, the
// frequency with a fitted drift model, or the histogram (with the statistics
// of the samples) or the autocorrelation of phase or frequency.
import { state, update, activeClock, multi, markStale } from "../state.js";
import { drift as fetchDrift, driftKey, driftMethod, distribution, distKey, distOf } from "../core.js";
import { noiseLabel, noiseTitle } from "../estimators.js";
import { loadPlotly, baseLayout, linAxis, CONFIG, css, seriesColor, timeUnit } from "../plot.js";
import { thousands, duration, sci, startText, escapeHtml as h } from "../format.js";

const FILE_ICON = `<svg width="22" height="26" viewBox="0 0 22 26" fill="none" stroke="currentColor" stroke-width="1.8" aria-hidden="true"><path d="M3 2h10l6 6v16H3z"/><path d="M13 2v6h6"/></svg>`;
// Fitted to the phase and removed before every analysis. On phase a line is a
// frequency offset and a parabola adds a linear frequency drift.
const DRIFTS = [
  ["none", "None", "Analyze the record as measured"],
  ["linear", "Linear", "Remove a frequency offset (a straight line in phase)"],
  ["quadratic", "Quadratic", "Remove a frequency offset and a linear frequency drift (a parabola in phase)"],
];

// Frequency drift models fitted on the "Frequency drift" view.
const MODELS = [
  ["linear", "Linear", "y = a + b·t"],
  ["quadratic", "Quadratic", "y = a + b·t + c·t²"],
  ["log", "Log", "y = a·ln(1 + b·t) + c, the stabilization of an oscillator"],
  ["diffusion", "Diffusion", "y = a + b·√(t + c)"],
];

// Estimators of the linear model (drift_fit's methods).
const METHODS = [
  ["frequency_fit", "Line through frequency", "Least-squares line through the frequency: least variance on white FM"],
  ["phase_fit", "Parabola through phase", "Least-squares parabola through the phase: least variance on white PM"],
  ["bisection", "Bisection (three-point)", "Mean frequency of the second half less that of the first, over half the record"],
  ["second_difference", "Mean second difference", "Mean of the second differences of the phase: least variance on random-walk FM"],
  ["four_point", "Four-point", "Greenhall's estimator from four values of the summed phase: near the least variance on every noise"],
];

// The plots of the record card.
const VIEWS = [
  ["phase", "Phase", "The phase record, with its outliers"],
  ["drift", "Frequency drift", "Fit a frequency drift model"],
  ["hist", "Histogram", "Histogram of the samples, with a normal curve and their statistics"],
  ["acf", "Autocorrelation", "Autocorrelation of the samples, with the noise type it points to"],
];

const MAX_SIM_CLOCKS = 8;
let root, plotEl, Plotly, ctx, relayoutBound = false;
const inflight = new Map();   // clock id (and view) → key of the request in flight

export function mount(el, c) {
  root = el; ctx = c;
  root.innerHTML = `
    <div class="col">
      <div class="toolbar">
        <div class="seg" role="group" aria-label="Source" data-ref="mode">
          <button type="button" data-mode="file">Load a file</button><button type="button" data-mode="synth">Synthesize</button></div>
        <span class="caption" data-ref="caption"></span>
      </div>
      <div class="data-top" data-ref="top"></div>
      <div class="card record-card">
        <div class="card-row">
          <div class="seg seg-sm" role="group" aria-label="Plot" data-ref="view">${VIEWS.map(([v, l, tip]) =>
            `<button type="button" data-view="${v}" title="${tip}">${l}</button>`).join("")}</div>
          <h2 data-ref="plot-title">Phase x(t)</h2><span class="muted small" data-ref="plot-hint"></span><span class="range-bar" data-ref="range"></span></div>
        <div class="drift-bar" data-ref="drift-bar" hidden></div>
        <div class="plot-row"><div class="plot" data-ref="plot"></div><div class="dist-aside" data-ref="aside" hidden></div><div class="empty" data-ref="empty" hidden></div></div>
      </div>
    </div>`;
  plotEl = root.querySelector("[data-ref=plot]");
  relayoutBound = false;
  root.addEventListener("click", onClick);
  root.addEventListener("change", onChange);
  // The drift view's toolbar takes height from the plot; follow its size.
  // offsetParent is null once the page has moved on (Plotly rejects a hidden plot).
  new ResizeObserver(() => { if (Plotly && plotEl.data && plotEl.offsetParent) Plotly.Plots.resize(plotEl); }).observe(plotEl);
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {});
}

function fileCard(c) {
  const isFile = c?.source.type === "file";
  const nameField = isFile
    ? `<input class="clock-name" data-field="name" value="${h(c.name)}" aria-label="Clock name" title="Rename this clock">`
    : `<b>${h(c?.name ?? "")}</b>`;
  return `<div class="card file-card">
    ${isFile ? `<div class="file-title">${FILE_ICON}<div style="flex:1;min-width:0;overflow-wrap:anywhere">${nameField}
        <span class="muted small">${h(c.source.name)} · ${h(c.sourceInfo)}</span></div></div>
      <div class="file-actions"><button class="btn" data-act="import">Import settings…</button><button class="btn" data-act="open">Open another file…</button>${addClockBtn()}</div>`
    : `<div class="file-actions"><button class="btn btn-primary" data-act="open">Open a file…</button>${addClockBtn()}</div>`}
    <div class="dropzone" data-ref="drop">Drop a data file anywhere in the window</div>
    ${c?.error ? `<p class="error small">${h(c.error)}</p>` : ""}
  </div>`;
}

const addClockBtn = () => state.clocks.some((c) => c.record) ? `<button class="btn" data-act="add-clock">+ Add clock</button>` : "";

function synthCard(c) {
  const isSynth = c?.source.type === "synth";
  const s = isSynth ? c.source : { n: 16384, seed: 42, tau0: 1 };
  return `<div class="card synth-card">
    <div><b>Synthetic record</b><div class="muted small">White plus random-walk frequency noise, made in the core. Same seed, same record.
      With several clocks, each is noisier than the last by a factor of two, all against one perfect reference; the Ensemble page forms its time scale from them.</div></div>
    <div class="synth-grid">
      <label><span>Clocks</span><input type="number" data-synth="clocks" min="1" max="${MAX_SIM_CLOCKS}" step="1" value="${s.clocks ?? 1}"></label>
      <label><span>Samples</span><input type="number" data-synth="n" min="64" max="10000000" step="1" value="${s.n}"></label>
      <label><span>τ<sub>0</sub> (s)</span><input type="number" data-synth="tau0" min="0" step="any" value="${s.tau0}"></label>
      <label><span>Seed</span><input type="number" data-synth="seed" min="0" step="1" value="${s.seed}"></label>
      ${isSynth ? `<label class="synth-name"><span>Name of this clock</span><input class="clock-name" data-field="name" value="${h(c.name)}" aria-label="Clock name"></label>` : ""}
    </div>
    <div class="file-actions"><button class="btn btn-primary" data-act="synth">${c?.source.type === "synth" ? "Make record" : "Use a synthetic record"}</button>${addClockBtn()}</div>
  </div>`;
}

// NaN samples from the core, which puts a reading missing from the time
// column in its slot as one.
function gapRow(k) {
  const tip = k.gaps
    ? `${thousands(k.gaps)} missing sample${k.gaps === 1 ? "" : "s"}. ADEV, MDEV, HDEV, MHDEV, TDEV, HTDEV, MTIE and TIE rms leave them out; the other estimators, the spectra and the dynamic maps need Fill gaps by reflection`
    : "No missing samples";
  return `<div class="kv" title="${tip}"><span>Gaps</span><span class="mono ${k.gaps ? "warnv" : ""}">${thousands(k.gaps)}</span></div>`;
}

function checkCard(c) {
  const k = c?.record?.check;
  // A record the core refused (a frequency record with gaps, unfilled) keeps the box that fills it.
  if (!k) return `<div class="card check-card"><h2>Record check</h2><p class="muted small">${c?.error ? "The core could not read this record." : c ? "Reading…" : "No record yet."}</p>${c?.error ? fillBox(c) : ""}</div>`;
  const drift = Number.isFinite(k.drift_per_day) ? `${sci(k.drift_per_day, 2)} /day` : "–";
  const out = c.remove_outliers
    ? `<span class="mono">removed · <button class="link-btn" data-act="keep-outliers">undo</button></span>`
    : k.outliers ? `<span class="mono warnv">${k.outliers} · <button class="link-btn" data-act="remove-outliers">remove</button></span>`
    : `<span class="mono">0</span>`;
  return `<div class="card check-card">
    <h2>Record check</h2>
    <div class="kv"><span>Samples</span><span class="mono">${thousands(k.samples)}</span></div>
    <div class="kv"><span>Span</span><span class="mono">${duration(k.span_s)}</span></div>
    ${k.start ? `<div class="kv" title="UTC time of the first sample, from the file's timestamps"><span>Start</span><span class="mono">${startText(k.start)}</span></div>` : ""}
    ${gapRow(k)}
    <div class="kv"><span>Points past 5σ</span>${out}</div>
    <div class="kv" title="Least-squares linear frequency drift of the record as analyzed"><span>${c.drift === "quadratic" ? "Drift left" : "Frequency drift"}</span><span class="mono">${drift}</span></div>
    <div class="rule"></div>
    ${driftBlock(c, k)}
    ${fillBox(c)}
  </div>`;
}

const fillBox = (c) => `<label class="check" title="Fill the missing samples and NaN readings from the data beside each gap, reflected (Howe and Schlossberger)">
      <input type="checkbox" data-field="fillgaps" ${c.fillgaps ? "checked" : ""}>Fill gaps by reflection</label>`;

function driftBlock(c, k) {
  const mode = c.drift ?? "none";
  const removed = [];
  if (Number.isFinite(k.removed_offset)) removed.push(`offset ${sci(k.removed_offset, 3)}`);
  if (Number.isFinite(k.removed_drift)) removed.push(`drift ${sci(k.removed_drift, 3)} /day`);
  return `<div class="drift-block">
    <div class="opt"><span title="Fitted to the phase by least squares and removed before every analysis">Drift removal</span>
      <div class="seg seg-sm" role="group" aria-label="Drift removal">${DRIFTS.map(([v, l, tip]) =>
        `<button type="button" data-drift="${v}" class="${v === mode ? "on" : ""}" title="${tip}">${l}</button>`).join("")}</div></div>
    <div class="muted small drift-note">${mode === "none" ? "Frequency offset and drift stay in the record." : removed.length ? `Removed ${removed.join(" · ")}` : "Fitting…"}</div>
  </div>`;
}

export function refresh() {
  if (!root) return;
  const c = activeClock();
  const mode = state.dataMode;
  root.querySelectorAll("[data-mode]").forEach((b) => b.classList.toggle("on", b.dataset.mode === mode));
  root.querySelector("[data-ref=caption]").innerHTML = mode === "file"
    ? "Phase or frequency records: text columns, Stable32, TimeLab, counter logs, BIPM UTC − UTC(k) and RINEX clock files"
    : "A test record with known noise, for trying estimators";
  const focus = document.activeElement?.dataset?.field === "name";
  if (!focus) root.querySelector("[data-ref=top]").innerHTML = (mode === "file" ? fileCard(c) : synthCard(c)) + checkCard(c);
  else root.querySelector(".check-card").outerHTML = checkCard(c);
  const view = state.dataView;
  root.querySelectorAll("[data-view]").forEach((b) => b.classList.toggle("on", b.dataset.view === view));
  const removed = c?.drift && c.drift !== "none" ? `, ${c.drift} drift removed` : "";
  const of = c ? (distOf(c) === "frequency" ? "frequency y" : "phase x") : "";
  const title = view === "drift" ? `Frequency y(t)${removed}`
    : view === "hist" ? `Histogram of ${of}${removed}`
    : view === "acf" ? `Autocorrelation of ${of}${removed}`
    : removed ? `Phase residual x(t)${removed}` : "Phase x(t)";
  root.querySelector("[data-ref=plot-title]").textContent = multi() && c ? `${title} · ${c.name}` : title;
  root.querySelector("[data-ref=plot-hint]").textContent = view === "drift" ? "Fit in crimson · drag to zoom"
    : view === "hist" ? (state.dist.normal ? "Normal curve in crimson" : "")
    : view === "acf" ? (c?.range ? "Of the analysis range" : "")
    : c?.range ? "" : state.dataDrag === "range" ? "Outliers in orange · drag across to set the range" : "Outliers in orange · drag to zoom";
  root.querySelector("[data-ref=range]").innerHTML = view === "phase" ? rangeBar(c) : "";
  const bar = root.querySelector("[data-ref=drift-bar]");
  const dist = view === "hist" || view === "acf";
  bar.hidden = view !== "drift" && !(dist && c?.record);
  if (dist && c?.record) bar.innerHTML = distBar(c, view);
  const aside = root.querySelector("[data-ref=aside]");
  aside.innerHTML = dist ? distAside(c, view) : "";
  aside.hidden = !aside.innerHTML;
  if (view === "drift") {
    ensureFit(c);
    bar.innerHTML = driftBar(c);
    drawDrift(c);
  } else if (view === "hist" || view === "acf") {
    ensureDist(c, view);
    view === "hist" ? drawHist(c) : drawAcf(c);
  } else {
    drawRecord(c);
  }
}

// Ask the core for the drift fit when it is missing or out of date.
function ensureFit(c) {
  if (!c?.record) return;
  const key = driftKey(c);
  if (c.driftFit?.key === key || inflight.get(c.id) === key) return;
  inflight.set(c.id, key);
  fetchDrift(c).then((out) => {
    if (inflight.get(c.id) !== key) return;
    inflight.delete(c.id);
    update(() => { c.driftFit = { key, ...out.drift }; });
  }).catch((err) => {
    if (inflight.get(c.id) !== key) return;
    inflight.delete(c.id);
    update(() => { c.driftFit = { key, error: err.message }; });
  });
}

const fitReady = (c) => c?.record && c.driftFit?.key === driftKey(c) && !c.driftFit.error;

function driftBar(c) {
  const f = fitReady(c) ? c.driftFit : null;
  const err = c?.driftFit?.key === driftKey(c) ? c.driftFit.error : null;
  const names = ["a", "b", "c"];
  // A method without a standard error gives null.
  const pm = (v, digits) => (Number.isFinite(v) ? ` ± ${sci(v, digits)}` : "");
  const coef = f ? f.coef.map((v, i) => `${names[i]} = ${sci(v, 4)}${pm(f.stderr[i], 2)}`).join(" · ") : "";
  const linear = state.driftModel === "linear";
  const method = driftMethod();
  const result = err ? `<span class="error">${h(err)}</span>`
    : f ? `<span title="${h(coef)}"><b class="mono">${sci(f.drift_per_day, 3)}${pm(f.drift_per_day_err, 2)} /day</b> at the end of the record
        <span class="muted">· rms residual ${sci(f.rms, 3)}</span></span>`
    : `<span class="muted">Fitting…</span>`;
  return `<div class="seg seg-sm" role="group" aria-label="Drift model">${MODELS.map(([v, l, tip]) =>
      `<button type="button" data-model="${v}" class="${v === state.driftModel ? "on" : ""}" title="${tip}">${l}</button>`).join("")}</div>
    <label class="tool-field" title="${linear ? "How the linear drift is estimated" : "The estimators apply to the linear model"}">Estimator <select data-ref="method" ${linear ? "" : "disabled"}>${METHODS.map(([v, l, tip]) =>
      `<option value="${v}" title="${tip}" ${v === method ? "selected" : ""}>${l}</option>`).join("")}</select></label>
    <label class="check"><input type="checkbox" data-ref="resid" ${state.driftResid ? "checked" : ""}>Residuals</label>
    <span class="drift-result small">${result}</span>`;
}

function drawDrift(c) {
  const empty = root.querySelector("[data-ref=empty]");
  const f = fitReady(c) ? c.driftFit : null;
  const s = f?.series;
  if (!s || !s.t.length) {
    empty.hidden = false;
    empty.innerHTML = `<p>${c?.record ? (c.driftFit?.error ? "No fit for this record." : "Fitting the drift model…") : "Open a file or drop one anywhere in the window."}</p>`;
    plotEl.style.visibility = "hidden";
    return;
  }
  empty.hidden = true;
  plotEl.style.visibility = "visible";
  if (!Plotly) return;
  const resid = state.driftResid;
  const y = resid ? s.y.map((v, i) => (Number.isFinite(v) ? v - s.fit[i] : null)) : s.y.map((v) => (Number.isFinite(v) ? v : null));
  const traces = [{
    x: s.t, y, mode: "lines", name: resid ? "residual" : "frequency",
    line: { color: multi() ? seriesColor(c.slot) : css("--s1"), width: 1.2 },
    hovertemplate: `t = %{x:.6g} s<br>${resid ? "y − fit" : "y"} = %{y:.5g}<extra></extra>`,
  }];
  if (!resid) {
    traces.push({
      x: s.t, y: s.fit, mode: "lines", name: `${state.driftModel} fit`,
      line: { color: css("--accent"), width: 2 },
      hovertemplate: `t = %{x:.6g} s<br>fit = %{y:.5g}<extra></extra>`,
    });
  }
  const yaxis = { ...linAxis(resid ? "y − fit" : "Fractional frequency y"), exponentformat: "e" };
  Plotly.react(plotEl, traces, { ...baseLayout(), margin: { l: 70, r: 18, t: 10, b: 46 }, xaxis: linAxis("Time t (s)"), yaxis, showlegend: false, uirevision: `${c.id}-drift` }, CONFIG);
}

// Ask the core for the statistics with the histogram or the autocorrelation
// when they are missing or out of date.
function ensureDist(c, what) {
  if (!c?.record) return;
  const key = distKey(c, what), slot = `${c.id}:${what}`;
  if (c[what]?.key === key || inflight.get(slot) === key) return;
  inflight.set(slot, key);
  const done = (fields) => {
    if (inflight.get(slot) !== key) return;
    inflight.delete(slot);
    update(() => { c[what] = { key, ...fields }; });
  };
  distribution(c, what).then(done, (err) => done({ error: err.message }));
}

const distReady = (c, what) => c?.record && c[what]?.key === distKey(c, what) && !c[what].error ? c[what] : null;

// Beside the plot: the statistics that Stable32's Statistics function
// shows, or the readings of the autocorrelation.
const kv = ([k, v, tip]) => `<div class="kv"${tip ? ` title="${h(tip)}"` : ""}><span>${k}</span><span class="mono">${v}</span></div>`;

function distAside(c, what) {
  const err = c?.[what]?.key === distKey(c, what) ? c[what].error : null;
  const got = distReady(c, what);
  return err ? `<p class="error small">${h(err)}</p>` : got ? distReadings(got, what) : "";
}

function distReadings(got, what) {
  if (what === "hist") {
    const st = got.statistics;
    return `<div aria-label="Statistics">${[["Points", thousands(st.n), "Samples measured"], ["Gaps", thousands(st.gaps), "Missing samples"],
      ["Min", sci(st.min, 4)], ["Max", sci(st.max, 4)], ["Mean", sci(st.mean, 4)], ["Median", sci(st.median, 4)],
      ["Std dev", sci(st.std, 4), "Standard deviation of the samples, with n − 1 in the denominator"],
      ...(st.quantity === "phase" ? [["Std dev (freq.)", sci(st.std_y, 4), "Standard deviation of the fractional frequency, the Std Dev of Stable32"]] : []),
      ["σ<sub>y</sub>(τ<sub>0</sub>)", sci(st.adev, 4), "Allan deviation at the sample interval"],
      ["Noise", noiseLabel(st.noise_type), `${noiseTitle(st.noise_type)}, identified at the sample interval`]].map(kv).join("")}</div>`;
  }
  const a = got.acf;
  return `<div aria-label="Autocorrelation">${[["τ", duration(a.tau), "m·τ₀"], ["Points", thousands(a.n), "Points in the series"],
    ["r<sub>1</sub>", a.r.length > 1 ? a.r[1].toFixed(3) : "–", "Autocorrelation at lag 1"],
    ["α", Number.isFinite(a.alpha) ? a.alpha.toFixed(2) : "–", "Lag-1 estimate of the power-law exponent (Riley and Greenhall); it needs a record of 30·m samples"],
    ["Noise", noiseLabel(a.noise_type), `${noiseTitle(a.noise_type)}, identified at τ`]].map(kv).join("")}</div>`;
}

function distBar(c, what) {
  const of = distOf(c), d = state.dist;
  const seg = `<div class="seg seg-sm" role="group" aria-label="Samples">${[["phase", "Phase"], ["frequency", "Frequency"]].map(([v, l]) =>
    `<button type="button" data-of="${v}" class="${v === of ? "on" : ""}">${l}</button>`).join("")}</div>`;
  const field = (k, label, tip, extra) => `<label class="opt" title="${tip}"><span>${label}</span>
    <input type="number" class="num-sm" data-dist="${k}" min="1" step="1" ${extra} value="${d[k] ?? ""}"></label>`;
  const controls = what === "hist"
    ? `${field("bins", "Bins", "Number of bins; empty for Scott's rule", `max="1000" placeholder="auto"`)}
      <label class="check"><input type="checkbox" data-dist="normal" ${d.normal ? "checked" : ""}>Normal curve</label>`
    : `${field("m", "Averaging factor m", "τ = m·τ₀", "")}
      ${field("lags", "Lags", "Largest lag; empty for half the series, at most 1000", `placeholder="auto"`)}`;
  return `${seg}${controls}`;
}

function showEmpty(msg) {
  const empty = root.querySelector("[data-ref=empty]");
  empty.hidden = !msg;
  if (msg) empty.innerHTML = `<p>${msg}</p>`;
  plotEl.style.visibility = msg ? "hidden" : "visible";
}

function distEmpty(c, what) {
  if (!c?.record) return "Open a file or drop one anywhere in the window.";
  if (distReady(c, what)) return "";
  return c[what]?.key === distKey(c, what) ? `No ${what === "hist" ? "histogram" : "autocorrelation"} for this record.` : "Computing…";
}

// Bars of the histogram and the normal curve of the samples' mean and
// standard deviation. Values far from zero next to their spread are shown
// less their mean, which the axis title gives. Phase is in a time unit.
function drawHist(c) {
  const msg = distEmpty(c, "hist");
  showEmpty(msg);
  if (msg || !Plotly) return;
  const hg = c.hist.histogram, freq = distOf(c) === "frequency";
  const shift = Math.abs(hg.mean) > 10 * hg.std ? hg.mean : 0;
  const [f, u] = freq ? [1, ""] : timeUnit(Math.max(Math.abs(hg.edges[0] - shift), Math.abs(hg.edges.at(-1) - shift)) || 1);
  const v = (x) => (x - shift) / f;
  const centers = hg.counts.map((_, k) => v((hg.edges[k] + hg.edges[k + 1]) / 2));
  const w = hg.edges.length > 1 ? (hg.edges[1] - hg.edges[0]) : 0;
  const traces = [{
    type: "bar", x: centers, y: hg.counts, width: w / f || undefined, name: "count",
    marker: { color: multi() ? seriesColor(c.slot) : css("--s1"), line: { width: 0 } },
    hovertemplate: `%{x:.5g}${u ? ` ${u}` : ""}<br>%{y} samples<extra></extra>`,
  }];
  if (state.dist.normal && hg.std > 0 && w > 0) {
    const lo = Math.min(hg.edges[0], hg.mean - 4 * hg.std), hi = Math.max(hg.edges.at(-1), hg.mean + 4 * hg.std);
    const xs = Array.from({ length: 241 }, (_, i) => lo + (hi - lo) * i / 240);
    const k = hg.n * w / (hg.std * Math.sqrt(2 * Math.PI));
    traces.push({
      x: xs.map(v), y: xs.map((x) => k * Math.exp(-0.5 * ((x - hg.mean) / hg.std) ** 2)), mode: "lines", name: "normal",
      line: { color: css("--accent"), width: 2 }, hovertemplate: "normal: %{y:.4g}<extra></extra>",
    });
  }
  const name = freq ? "Fractional frequency y" : "Phase x";
  const xaxis = { ...linAxis(shift ? `${name} − ${sci(shift, 4)}${freq ? "" : " s"}` : name), exponentformat: "e", ...(u ? { ticksuffix: ` ${u}` } : {}) };
  Plotly.react(plotEl, traces, {
    ...baseLayout(), margin: { l: 70, r: 18, t: 10, b: 46 }, bargap: 0.04, xaxis, yaxis: linAxis("Samples per bin"),
    uirevision: `${c.id}-hist`,
  }, CONFIG);
}

function drawAcf(c) {
  const msg = distEmpty(c, "acf");
  showEmpty(msg);
  if (msg || !Plotly) return;
  const a = c.acf.acf;
  Plotly.react(plotEl, [{
    x: a.lag, y: a.r, mode: a.lag.length > 200 ? "lines" : "lines+markers", name: "acf",
    line: { color: multi() ? seriesColor(c.slot) : css("--s1"), width: 1.4 }, marker: { size: 4 },
    customdata: a.lag.map((k) => k * a.tau),
    hovertemplate: "lag %{x} (%{customdata:.4g} s)<br>r = %{y:.4f}<extra></extra>",
  }], {
    ...baseLayout(), margin: { l: 70, r: 18, t: 10, b: 46 },
    xaxis: linAxis(`Lag k (τ = ${duration(a.tau)})`), yaxis: { ...linAxis("r(k)"), rangemode: "tozero", zeroline: true, zerolinecolor: css("--line-strong") },
    uirevision: `${c.id}-acf`,
  }, CONFIG);
}

// The analysis range: drag across the plot in "Select range" mode. The check
// above and this plot keep the whole record; the runs use the range.
function rangeBar(c) {
  if (!c?.record) return "";
  const r = c.range;
  return `${r ? `<span class="small">Analyzing <span class="mono">${duration(r[0])} – ${duration(r[1])}</span> · <button class="link-btn" data-act="clear-range">whole record</button></span>` : ""}
    <span class="seg seg-sm" role="group" aria-label="Drag on the plot to">${[["zoom", "Zoom", "Drag to zoom"], ["range", "Select range", "Drag across the plot to analyze only that part of the record"]].map(([v, l, tip]) =>
      `<button type="button" data-drag="${v}" class="${state.dataDrag === v ? "on" : ""}" title="${tip}">${l}</button>`).join("")}</span>`;
}

function setRange(c, r) {
  update((s) => { c.range = r; c.result = null; c.dynamic = null; markStale(s, "record"); });
}

function drawRecord(c) {
  const empty = root.querySelector("[data-ref=empty]");
  const s = c?.record?.series;
  if (!s || !s.t.length) {
    empty.hidden = false;
    empty.innerHTML = `<p>${c?.error ? h(c.error) : c ? "Reading the record…" : "Open a file or drop one anywhere in the window."}</p>`;
    plotEl.style.visibility = "hidden";
    return;
  }
  empty.hidden = true;
  plotEl.style.visibility = "visible";
  if (!Plotly) return;
  const k = c.record.check;
  const finite = s.x.filter(Number.isFinite);
  const maxAbs = Math.max(...finite.map(Math.abs), ...(k?.outlier_x ?? []).map(Math.abs));
  const [f, u] = timeUnit(maxAbs || 1);
  const traces = [{
    x: s.t, y: s.x.map((v) => (Number.isFinite(v) ? v / f : null)), mode: "lines", line: { color: multi() ? seriesColor(c.slot) : css("--s1"), width: 1.4 },
    hovertemplate: `t = %{x:.6g} s<br>x = %{y:.5g} ${u}<extra></extra>`, name: "phase",
  }];
  if (k?.outlier_t?.length && !c.remove_outliers) {
    traces.push({
      x: k.outlier_t, y: k.outlier_x.map((v) => v / f), mode: "markers", name: "outliers",
      marker: { color: css("--outlier") || css("--s2"), size: 8 },
      hovertemplate: `outlier<br>t = %{x:.6g} s<br>x = %{y:.5g} ${u}<extra></extra>`,
    });
  }
  const picking = state.dataDrag === "range";
  const yaxis = { ...linAxis(""), ticksuffix: ` ${u}`, fixedrange: picking };
  const shapes = c.range ? [{ type: "rect", xref: "x", yref: "paper", x0: c.range[0], x1: c.range[1], y0: 0, y1: 1,
    fillcolor: css("--accent-soft") || "rgba(158,27,50,0.12)", line: { width: 0 }, layer: "below" }] : [];
  Plotly.react(plotEl, traces, {
    ...baseLayout(), margin: { l: 70, r: 18, t: 10, b: 46 }, xaxis: linAxis("Time t (s)"), yaxis, shapes,
    dragmode: "zoom", uirevision: `${c.id}-${picking}`,
  }, { ...CONFIG, showTips: false });
  if (!relayoutBound) {
    relayoutBound = true;
    // In "Select range" mode a horizontal drag is a zoom we turn into the range.
    plotEl.on("plotly_relayout", (ev) => {
      const k = activeClock();
      const a = ev["xaxis.range[0]"], b = ev["xaxis.range[1]"];
      if (state.dataDrag !== "range" || state.dataView === "drift" || !k?.record || a === undefined || b === undefined) return;
      const end = (k.record.n - 1) * k.record.tau0;
      const t0 = Math.max(0, Math.min(a, b)), t1 = Math.min(end, Math.max(a, b));
      const snap = (t) => Math.round(t / k.record.tau0) * k.record.tau0;
      Plotly.relayout(plotEl, { "xaxis.autorange": true });
      if (snap(t1) - snap(t0) >= 16 * k.record.tau0) setRange(k, [snap(t0), snap(t1)]);
      else ctx.status("warn", "The range needs at least 16 samples.");
    });
  }
}

function onClick(e) {
  const v = e.target.closest("[data-view]")?.dataset.view;
  if (v) { update((s) => { s.dataView = v; }); return; }
  const of = e.target.closest("[data-of]")?.dataset.of;
  if (of) { update((s) => { s.dist.of = of; }); return; }
  const md = e.target.closest("[data-model]")?.dataset.model;
  if (md) { update((s) => { s.driftModel = md; }); return; }
  const dr = e.target.closest("[data-drift]")?.dataset.drift;
  if (dr) { const c = activeClock(); if (c && c.drift !== dr) ctx.setRecord(c, { drift: dr }); return; }
  const dg = e.target.closest("[data-drag]")?.dataset.drag;
  if (dg) { update((s) => { s.dataDrag = dg; }); return; }
  const m = e.target.closest("[data-mode]")?.dataset.mode;
  if (m) { update((s) => { s.dataMode = m; }); return; }
  const act = e.target.closest("[data-act]")?.dataset.act;
  const c = activeClock();
  if (act === "open") ctx.openFile(c);
  else if (act === "import") ctx.reimport(c);
  else if (act === "add-clock") ctx.addClock();
  else if (act === "remove-outliers") ctx.setRecord(c, { remove_outliers: true });
  else if (act === "keep-outliers") ctx.setRecord(c, { remove_outliers: false });
  else if (act === "clear-range") setRange(c, null);
  else if (act === "synth") {
    const v = (k) => Number(root.querySelector(`[data-synth=${k}]`).value);
    const n = Math.round(v("n")), tau0 = v("tau0"), seed = Math.round(v("seed")), clocks = Math.round(v("clocks"));
    if (!(n >= 64) || !(tau0 > 0) || !(seed >= 0)) { ctx.status("error", "Samples must be at least 64 and τ₀ above zero."); return; }
    if (!(clocks >= 1 && clocks <= MAX_SIM_CLOCKS)) { ctx.status("error", `Simulate between 1 and ${MAX_SIM_CLOCKS} clocks.`); return; }
    ctx.useSynth(c, { n, tau0, seed, clocks });
  }
}

function onChange(e) {
  if (e.target.dataset.ref === "resid") { update((s) => { s.driftResid = e.target.checked; }); return; }
  if (e.target.dataset.ref === "method") { update((s) => { s.driftMethod = e.target.value; }); return; }
  const dk = e.target.dataset.dist;
  if (dk === "normal") { update((s) => { s.dist.normal = e.target.checked; }); return; }
  if (dk) {
    // Bins and lags: empty for the default; m: 1 or more.
    const n = Math.round(Number(e.target.value));
    const ok = e.target.value === "" ? dk !== "m" : n >= 1 && (dk !== "bins" || n <= 1000);
    if (!ok) { ctx.status("error", dk === "bins" ? "Bins must be between 1 and 1000." : "Enter a whole number of 1 or more."); refresh(); return; }
    update((s) => { s.dist[dk] = e.target.value === "" ? null : n; });
    return;
  }
  const f = e.target.dataset.field;
  const c = activeClock();
  if (!c || !f) return;
  if (f === "name") { const v = e.target.value.trim() || c.name; update(() => { c.name = v; c.typedName = true; }); }
  else if (f === "fillgaps") ctx.setRecord(c, { fillgaps: e.target.checked });
}

export function hint() {
  return state.dataDrag === "range"
    ? "Drag across the plot to choose the part of the record to analyze"
    : "Drop a file anywhere in the window · drag on the plot to zoom · double-click to reset";
}

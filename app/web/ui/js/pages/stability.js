// Stability page: σ(τ) plot and, on the right, the per-τ table (one clock) or
// "Compare at τ" (several clocks, overlaid). Mounted once per visit;
// refresh() redraws from state.
import { state, update, visibleClocks, multi, markStale } from "../state.js";
import { byId, label, gridName, CONFIDENCES, noiseLabel, noiseTitle, noiseSlope, BIAS_CORRECTED, biasEntry } from "../estimators.js";
import { loadPlotly, baseLayout, logAxis, CONFIG, seriesColor, DASHES, hexAlpha, css } from "../plot.js";
import { sci, num, escapeHtml as h } from "../format.js";

// ITU-T wander masks the core checks (sigmatau --mask), by its names.
const MASKS = [
  ["G811", "G.811 PRC"], ["G813_opt1", "G.813 SEC opt 1"], ["G813_opt2", "G.813 SEC opt 2"],
  ["G8262_opt1", "G.8262 EEC opt 1"], ["G8262_opt2", "G.8262 EEC opt 2"],
  ["G823_PRC", "G.823 PRC interface"], ["G823_SSU", "G.823 SSU interface"],
  ["G823_SEC", "G.823 SEC interface"], ["G823_PDH", "G.823 PDH interface"],
];

const AXES = { sigma_y: "Frequency σ<sub>y</sub>(τ)", sigma_x: "Time σ<sub>x</sub>(τ)" };
const Y_TITLE = { sigma_y: "σ<sub>y</sub>(τ)", sigma_x: "σ<sub>x</sub>(τ) (s)" };

let root, plotEl, Plotly, clickBound = false;

export function mount(el, ctx) {
  root = el;
  root.innerHTML = `
    <div class="col">
      <div class="toolbar">
        <div class="seg" role="group" aria-label="Quantity" data-ref="axis">${Object.entries(AXES).map(([k, v]) =>
          `<button type="button" data-axis="${k}">${v}</button>`).join("")}</div>
        <span class="caption" data-ref="caption"></span>
        <span class="pin-list" data-ref="pins"></span>
        <span class="toggles">
          <label class="check" title="Noise type identified at each τ (lag-1 autocorrelation), with the power-law slope it implies"><input type="checkbox" data-tog="noiseId">Noise ID</label>
          <label class="check" title="Draw the confidence intervals as a shaded band instead of error bars"><input type="checkbox" data-tog="band">CI band</label>
          <label class="check" title="The power-law noise model fitted to MHDEV: its curve for each estimator, and h_α with the ADEV it gives at round τ"><input type="checkbox" data-tog="fit">Fit</label>
          <select data-ref="mask" aria-label="Telecom mask" title="ITU-T wander mask drawn over MTIE and TDEV, with pass or fail">
            <option value="">No mask</option>${MASKS.map(([v, l]) => `<option value="${v}">${l}</option>`).join("")}</select>
        </span>
      </div>
      <div class="card plot-card"><div class="plot" data-ref="plot"></div><div class="empty" data-ref="empty" hidden></div></div>
    </div>
    <div class="card side-card" data-ref="table"></div>`;
  plotEl = root.querySelector("[data-ref=plot]");
  clickBound = false;
  root.addEventListener("click", (e) => {
    const ax = e.target.closest("[data-axis]")?.dataset.axis;
    if (ax) update((s) => { s.axis = ax; });
    if (e.target.closest("[data-act=load]")) ctx.goData();
    if (e.target.closest("[data-act=run]")) ctx.run();
    if (e.target.closest("[data-act=pin]")) pinRun();
    const act = e.target.closest("[data-fit]")?.dataset.fit;
    if (act) editModel(act, ctx);
    const un = e.target.closest("[data-unpin]")?.dataset.unpin;
    if (un) update((s) => { s.pins = s.pins.filter((p) => p.id !== Number(un)); });
  });
  root.addEventListener("change", (e) => {
    const tog = e.target.dataset.tog;
    if (tog) { update((s) => { s[tog] = e.target.checked; }); return; }
    const r = e.target.dataset.ref;
    if (r === "table-dev") update((s) => { s.tableDev = e.target.value; });
    if (r === "cmp-tau") update((s) => { s.compareTau = Number(e.target.value); });
    if (r === "mask") update((s) => { s.mask = e.target.value; });
  });
  loadPlotly().then((P) => { Plotly = P; refresh(); }).catch(() => {
    showEmpty("Plotly could not be loaded. Run <code>python serve.py --fetch-plotly</code> once while online.");
  });
}

const devsOnAxis = () => state.devs.filter((d) => byId[d].axis === state.axis);
const resultOf = (clock, dev) => clock.result?.results.find((r) => r.dev === dev) ?? null;
const sameTau = (a, b) => Math.abs(a - b) <= 1e-9 * Math.max(Math.abs(a), Math.abs(b));

function showEmpty(html) {
  const e = root.querySelector("[data-ref=empty]");
  e.innerHTML = html;
  e.hidden = !html;
  plotEl.style.visibility = html ? "hidden" : "visible";
}

export function refresh() {
  if (!root) return;
  root.querySelectorAll("[data-axis]").forEach((b) => b.classList.toggle("on", b.dataset.axis === state.axis));
  const conf = state.ci ? `${(CONFIDENCES.find(([v]) => v === state.confidence) ?? [0, state.confidence])[1]} confidence · ` : "";
  const biased = state.devs.some((d) => BIAS_CORRECTED.has(d));
  const bias = biased ? ` · ${biasEntry(state.bias)[2]}` : "";
  root.querySelector("[data-ref=caption]").textContent = `${conf}${gridName(state.grid).toLowerCase()} τ grid${bias}${state.noise === "auto" ? "" : ` · noise ${state.noise === "fit" ? "from the fit" : noiseLabel(state.noise)}`}`;
  const maskSel = root.querySelector("[data-ref=mask]");
  maskSel.value = state.mask;
  maskSel.hidden = state.axis !== "sigma_x";
  // Noise ID and the band need the intervals: noise identification runs with them.
  root.querySelectorAll("[data-tog]").forEach((el) => {
    const lab = el.closest("label");
    lab.dataset.tip ??= lab.title;
    el.checked = !!state[el.dataset.tog];
    const off = !state.ci && el.dataset.tog !== "fit";
    el.disabled = off;
    lab.title = off ? "Turn on confidence intervals in the panel; noise identification runs with them" : lab.dataset.tip;
  });

  const clocks = visibleClocks();
  const devs = devsOnAxis();
  const hasRecord = state.clocks.some((c) => c.record);
  const hasResult = clocks.some((c) => c.result);
  const other = state.axis === "sigma_y" ? "sigma_x" : "sigma_y";

  if (!hasRecord) showEmpty(`<h2>No record yet</h2><p>Load a phase or frequency file, or synthesize one, on the Data page.</p><button class="btn" data-act="load">Go to Data</button>`);
  else if (!hasResult) showEmpty(`<h2>Not analyzed yet</h2><p>Pick estimators in the panel and run the analysis.</p><button class="btn btn-primary" data-act="run">Run analysis</button>`);
  else if (!devs.length) showEmpty(`<h2>No ${state.axis === "sigma_y" ? "σ<sub>y</sub>" : "σ<sub>x</sub>"} estimator selected</h2><p>${state.axis === "sigma_x" ? "Add TDEV, MTIE or another time estimator" : "Add ADEV or another frequency estimator"}, or switch to ${AXES[other]}.</p>`);
  else if (!clocks.length) showEmpty(`<h2>All clocks are hidden</h2><p>Show a clock with ◉ in the panel.</p>`);
  else showEmpty("");

  renderPins();
  renderTable(clocks, devs);
  if (Plotly && hasResult && devs.length && clocks.length) draw(clocks, devs);
}

// Pinned runs: earlier results kept for comparison, drawn thin and dashed
// under the current ones.
const PIN_DASH = ["dot", "dashdot", "longdash", "longdashdot"];
function pinTraces(devs) {
  const out = [];
  state.pins.forEach((p, k) => {
    for (const c of p.clocks) for (const r of c.results) {
      if (!devs.includes(r.dev)) continue;
      const di = state.devs.indexOf(r.dev);
      const color = hexAlpha(seriesColor(multi() ? c.slot : Math.max(di, 0)), 0.75);
      const name = `${p.label} · ${multi() ? `${c.name} · ` : ""}${label(r.dev)}`;
      out.push({
        x: r.tau, y: r.sigma, name, mode: "lines", line: { color, width: 1.3, dash: PIN_DASH[k % PIN_DASH.length] },
        hovertemplate: `<b>${h(name)}</b><br>τ = %{x:.4g} s<br>σ = %{y:.4e}<extra></extra>`,
      });
    }
  });
  return out;
}

// What a run was, from the core's own record of it (the settings may have
// changed since): shown on its pin.
function runSummary(c) {
  const a = c.result.doc.analysis ?? {};
  const parts = [(a.deviations ?? []).map(label).join(", "), a.tau_mode === "m" ? "explicit m" : `${gridName(a.tau_mode ?? "").toLowerCase()} τ`];
  const [r0, r1] = a.range ?? [null, null];
  if (r0 !== null || r1 !== null) parts.push(`range ${r0 === null ? "start" : `${num(r0, 4)} s`} – ${r1 === null ? "end" : `${num(r1, 4)} s`}`);
  if (a.remove_drift && a.remove_drift !== "none") parts.push(`${a.remove_drift} drift removed`);
  if (a.tau_min !== null && a.tau_min !== undefined || a.tau_max !== null && a.tau_max !== undefined) parts.push(`τ ${a.tau_min ?? "…"}–${a.tau_max ?? "…"} s`);
  return parts.join(" · ");
}

export function pinRun() {
  const clocks = visibleClocks().filter((c) => c.result);
  if (!clocks.length) return false;
  update((s) => {
    const n = (s.pinCount = (s.pinCount ?? 0) + 1);
    s.pins.push({
      id: n, label: `Pin ${n}`, title: clocks.map((c) => `${c.name}: ${runSummary(c)}`).join("\n"),
      clocks: clocks.map((c) => ({ name: c.name, slot: c.slot, results: c.result.results.map((r) => ({ dev: r.dev, tau: r.tau, sigma: r.sigma })) })),
    });
  });
  return true;
}

function renderPins() {
  const el = root.querySelector("[data-ref=pins]");
  const has = visibleClocks().some((c) => c.result);
  el.innerHTML = `<button class="btn btn-sm" data-act="pin" ${has ? "" : "disabled"} title="Keep these results on the plot to compare with the next run">Pin run</button>${state.pins.map((p) =>
    `<span class="chip" title="${h(p.title)}">${h(p.label)}<button class="link-btn" data-unpin="${p.id}" aria-label="Remove ${h(p.label)}">×</button></span>`).join("")}`;
}

function traces(clocks, devs) {
  const several = multi();
  const out = pinTraces(devs);
  for (const c of clocks) {
    devs.forEach((d) => {
      const r = resultOf(c, d);
      if (!r) return;
      const di = state.devs.indexOf(d);
      const color = several ? seriesColor(c.slot) : seriesColor(di);
      const dash = several ? DASHES[di % DASHES.length] : "solid";
      const ci = state.ci && r.ci_lo.length === r.tau.length;
      const name = several ? `${c.name} · ${label(d)}` : label(d);
      if (ci && state.band) {
        out.push({ x: r.tau, y: r.ci_hi, mode: "lines", line: { width: 0 }, hoverinfo: "skip", showlegend: false });
        out.push({ x: r.tau, y: r.ci_lo, mode: "lines", line: { width: 0 }, fill: "tonexty", fillcolor: hexAlpha(color, css("--ci-alpha")), hoverinfo: "skip", showlegend: false });
      }
      out.push({
        x: r.tau, y: r.sigma, name, mode: several ? "lines" : "lines+markers",
        line: { color, width: 2, dash }, marker: { color, size: 5 },
        error_y: ci && !state.band ? {
          type: "data", symmetric: false, color, thickness: 1, width: 0,
          array: r.ci_hi.map((v, i) => (v ?? NaN) - r.sigma[i]),
          arrayminus: r.sigma.map((v, i) => v - (r.ci_lo[i] ?? NaN)),
        } : undefined,
        customdata: r.tau.map((_, i) => [r.edf[i], noiseLabel(r.noise[i] ?? "")]),
        hovertemplate: `<b>${h(name)}</b><br>τ = %{x:.4g} s<br>σ = %{y:.4e}${r.edf.length ? "<br>EDF %{customdata[0]:.4g} · %{customdata[1]}" : ""}<extra></extra>`,
      });
      // While the model is built by hand: the MHDEV its components leave.
      const res = state.fit && state.fitEdit && !several && d === "mhdev" && c.result.powerlaw?.residual;
      if (res) out.push({
        x: res.tau, y: res.dev.map((v) => v || null), name: "MHDEV left by the model", mode: "lines+markers", showlegend: false,
        line: { color: css("--muted"), width: 1.2, dash: "dash" }, marker: { color: css("--muted"), size: 4 },
        hovertemplate: `<b>MHDEV left by the model</b><br>τ = %{x:.4g} s<br>σ = %{y:.4e}<extra></extra>`,
      });
      // The ITU-T mask's limit for MTIE or TDEV, once per estimator.
      const mask = maskOf(c, d);
      if (mask && c === clocks[0]) {
        const k = mask.tau.map((_, i) => i).filter((i) => mask[d][i] !== null);
        out.unshift({
        x: k.map((i) => mask.tau[i]), y: k.map((i) => mask[d][i]), name: `${maskName()} ${label(d)} mask`, mode: "lines",
        line: { color: css("--accent"), width: 1.6, dash: d === "mtie" ? "dash" : "dashdot" },
        hovertemplate: `<b>${h(mask.title)}</b><br>${label(d)} limit<br>τ = %{x:.4g} s<br>%{y:.4e} s<extra></extra>`,
        });
      }
      // The power-law model's curve, dotted under the measurement.
      const fit = state.fit && c.result.powerlaw?.curves?.find((k) => k.deviation === d);
      if (fit) out.splice(out.length - 1, 0, {
        x: fit.tau, y: fit.dev, name: `${name} model`, mode: "lines", showlegend: false,
        line: { color, width: 1.6, dash: "dot" },
        hovertemplate: `<b>${h(name)}, noise model</b><br>τ = %{x:.4g} s<br>σ = %{y:.4e}<extra></extra>`,
      });
    });
  }
  return out;
}

const maskName = () => (MASKS.find(([v]) => v === state.mask) ?? ["", ""])[1];
// The chosen mask from a clock's run, when it limits this estimator.
function maskOf(clock, dev) {
  if (!state.mask || (dev !== "mtie" && dev !== "tdev")) return null;
  return clock.result?.masks?.find((m) => m.name === state.mask) ?? null;
}

// Pass or fail of each MTIE / TDEV result against the chosen mask.
function maskBox(c, devs) {
  const m = state.mask && c.result?.masks?.find((k) => k.name === state.mask);
  const checks = m ? m.checks.filter((k) => devs.includes(k.deviation)) : [];
  if (!checks.length) return "";
  return `<div class="model-box"><h3>${h(m.title)}</h3><div class="model-grid">${checks.map((k) => {
    const ok = k.pass, q = k.margin;
    const text = q === null ? "no τ in the mask's range" : ok ? `passes, ${num(q, 3)}× under` : `fails, ${num(1 / q, 3)}× over`;
    return `<span>${label(k.deviation)}</span><span></span><span class="${ok ? "" : "warnv"}">${text}</span>`;
  }).join("")}</div></div>`;
}

// The σ(τ) figure for one axis: traces, layout and the noise overlay. The
// report export draws the same figure off screen.
export function figure(clocks, devs, axis = state.axis) {
  const layout = {
    ...baseLayout(),
    xaxis: logAxis("Averaging time τ (s)"),
    yaxis: logAxis(Y_TITLE[axis]),
    showlegend: !multi(),
    shapes: [], annotations: [],
    uirevision: `${axis}`,
  };
  const data = traces(clocks, devs);
  if (state.noiseId && state.ci && axis === "sigma_y") noiseOverlay(layout, data, clocks, devs);
  return { data, layout };
}

function draw(clocks, devs) {
  const { data, layout } = figure(clocks, devs);
  if (multi() && state.compareTau) {
    const accent = css("--accent");
    layout.shapes.push({ type: "line", xref: "x", yref: "paper", x0: state.compareTau, x1: state.compareTau, y0: 0, y1: 1, line: { color: accent, width: 1.2, dash: "dash" } });
    layout.annotations.push({ xref: "x", yref: "paper", x: Math.log10(state.compareTau), y: 1, xanchor: "left", yanchor: "top", showarrow: false, text: `τ = ${num(state.compareTau, 4)} s`, font: { color: accent, size: 11.5 } });
  }
  Plotly.react(plotEl, data, layout, CONFIG);
  if (!clickBound) {
    clickBound = true;
    plotEl.on("plotly_click", (ev) => {
      if (!multi()) return;
      const x = ev.points?.[0]?.x;
      const taus = compareTaus(visibleClocks(), devsOnAxis());
      if (!Number.isFinite(x) || !taus.length) return;
      const best = taus.reduce((a, b) => (Math.abs(Math.log(b / x)) < Math.abs(Math.log(a / x)) ? b : a));
      update((s) => { s.compareTau = best; });
    });
  }
}

// ── Noise identification overlay ──
// The strip along the top shows, per clock, the noise type identified at each
// τ, as runs of equal type. With one clock, dashed guides show the slope each
// run implies for the estimator the strip comes from.

const NOISE_VAR = { WHPM: "--noise-wpm", FLPM: "--noise-fpm", WHFM: "--noise-wfm", FLFM: "--noise-ffm", RWFM: "--noise-rwfm", FWFM: "--noise-fwfm", RRFM: "--noise-rrfm" };
const STRIP_PX = 18, STRIP_GAP = 3;

// The first selected σ_y estimator with noise types: the strip's source.
function noiseSource(clock, devs) {
  for (const d of devs) {
    const r = resultOf(clock, d);
    if (r && r.noise.some((v) => v)) return r;
  }
  return null;
}

// Runs of equal noise type, with edges halfway (in log τ) between points.
export function noiseRuns(r) {
  const runs = [];
  const n = r.tau.length;
  const edge = (i) => (i <= 0 ? r.tau[0] / Math.SQRT2 : i >= n ? r.tau[n - 1] * Math.SQRT2 : Math.sqrt(r.tau[i - 1] * r.tau[i]));
  for (let i = 0; i < n; i++) {
    const v = r.noise[i];
    if (!v) continue;
    const last = runs[runs.length - 1];
    if (last && last.noise === v && last.end === i) { last.end = i + 1; continue; }
    runs.push({ noise: v, start: i, end: i + 1 });
  }
  for (const u of runs) { u.lo = edge(u.start); u.hi = edge(u.end); }
  return runs;
}

function noiseOverlay(layout, tr, clocks, devs) {
  const several = clocks.length > 1;
  const withNoise = clocks.map((c) => [c, noiseSource(c, devs)]).filter(([, r]) => r);
  if (!withNoise.length) return;
  const ink = css("--ink-2");
  // Strips sit above the plot area in the Clocks list's order (first on
  // top), sized in pixels.
  withNoise.forEach(([c, r], k) => {
    const p0 = 4 + (withNoise.length - 1 - k) * (STRIP_PX + STRIP_GAP), p1 = p0 + STRIP_PX;
    for (const u of noiseRuns(r)) {
      layout.shapes.push({
        type: "rect", xref: "x", yref: "paper", ysizemode: "pixel", yanchor: 1, x0: u.lo, x1: u.hi, y0: p0, y1: p1,
        fillcolor: css(NOISE_VAR[u.noise] ?? "--table-head"), line: { width: several ? 1.5 : 0, color: several ? seriesColor(c.slot) : "rgba(0,0,0,0)" },
      });
      layout.annotations.push({
        xref: "x", yref: "paper", x: Math.log10(Math.sqrt(u.lo * u.hi)), y: 1, yanchor: "middle", yshift: (p0 + p1) / 2, showarrow: false,
        text: noiseLabel(u.noise), font: { size: 10.5, color: ink }, hovertext: `${several ? `${c.name}: ` : ""}${noiseTitle(u.noise)}`,
      });
      if (!several) {
        // Faint run boundaries through the plot.
        if (u.start > 0) layout.shapes.push({ type: "line", xref: "x", yref: "paper", x0: u.lo, x1: u.lo, y0: 0, y1: 1, layer: "below", line: { color: css("--line-strong"), width: 1, dash: "dot" } });
        const slope = noiseSlope(r.dev, u.noise);
        const pts = r.tau.slice(u.start, u.end).map((t, i) => [t, r.sigma[u.start + i]]).filter(([, v]) => v > 0);
        if (slope === null || !pts.length) continue;
        // Through the run's center of mass in log-log, across the run.
        const lx = pts.reduce((a, [t]) => a + Math.log10(t), 0) / pts.length;
        const ly = pts.reduce((a, [, v]) => a + Math.log10(v), 0) / pts.length;
        const at = (t) => 10 ** (ly + slope * (Math.log10(t) - lx)) * 1.6;
        const a = pts.length > 1 ? pts[0][0] : u.lo, b = pts.length > 1 ? pts[pts.length - 1][0] : u.hi;
        tr.push({
          x: [a, b], y: [at(a), at(b)], mode: "lines", hoverinfo: "text", showlegend: false,
          text: `${noiseTitle(u.noise)}: ${label(r.dev)} ∝ τ^${slope}`,
          line: { color: css("--muted"), width: 1.2, dash: "dash" },
        });
      }
    }
  });
  // Room for the strips above the curves.
  layout.margin = { ...layout.margin, t: 10 + withNoise.length * (STRIP_PX + STRIP_GAP) };
}

function compareTaus(clocks, devs) {
  const all = [];
  for (const c of clocks) for (const d of devs) for (const t of resultOf(c, d)?.tau ?? []) {
    if (!all.some((u) => sameTau(u, t))) all.push(t);
  }
  return all.sort((a, b) => a - b);
}

function renderTable(clocks, devs) {
  const box = root.querySelector("[data-ref=table]");
  const withResult = clocks.filter((c) => c.result);
  if (!withResult.length || !devs.length) { box.hidden = true; return; }
  box.hidden = false;
  if (multi()) return renderCompare(box, withResult, devs);

  const c = withResult[0];
  if (!devs.includes(state.tableDev)) state.tableDev = devs[0];
  const r = resultOf(c, state.tableDev);
  const rows = (r?.tau ?? []).map((t, i) => {
    const tip = r.ci_lo.length ? `${sci(r.ci_lo[i], 4)} … ${sci(r.ci_hi[i], 4)} (${(state.confidence * 100).toFixed(1)} %), n_eff ${r.neff[i] ?? "–"}` : `n_eff ${r.neff[i] ?? "–"}`;
    return `<div class="tbl-row tbl-body" title="${tip}"><span>${num(t, 5)}</span><span>${sci(r.sigma[i])}</span><span>${r.edf.length ? num(r.edf[i], 4) : "–"}</span>${r.noise[i] ? `<span class="badge" title="${noiseTitle(r.noise[i])}">${noiseLabel(r.noise[i])}</span>` : "<span></span>"}</div>`;
  }).join("");
  box.className = "card side-card per-tau";
  box.innerHTML = `
    <div class="card-head"><h2>Per-τ results</h2>
      <select data-ref="table-dev" aria-label="Estimator">${devs.map((d) => `<option value="${d}" ${d === state.tableDev ? "selected" : ""}>${label(d)}</option>`).join("")}</select></div>
    <div class="tbl"><div class="tbl-row tbl-head"><span>τ (s)</span><span>σ</span><span>EDF</span><span>noise</span></div>${rows || `<div class="tbl-row muted">No result for ${label(state.tableDev)}.</div>`}</div>
    ${state.axis === "sigma_x" ? maskBox(c, devs) : ""}
    ${state.fit ? modelBox(c) : ""}
    <div class="card-foot">Hover a row for its confidence interval. Every value is in the TSV and JSON exports.</div>`;
}

// The power-law model on MHDEV: the h_α it found and the ADEV it gives at
// round τ, the numbers a datasheet quotes. Edit builds it by hand instead,
// one noise type over a τ span at a time, each fitted to what the ones
// before it leave.
const NOISES = ["WHPM", "FLPM", "WHFM", "FLFM", "RWFM", "FWFM", "RRFM"];

function modelBox(c) {
  const p = c.result?.powerlaw;
  if (!p) return "";
  const hand = c.components.length > 0;
  const head = `<h3><span>Noise model <span class="muted">· S<sub>y</sub>(f) = Σ h<sub>α</sub> f<sup>α</sup>, ${hand ? "built by hand" : "fitted"} on MHDEV</span></span>
    <button type="button" class="link-btn" data-fit="toggle">${state.fitEdit ? "Done" : "Edit"}</button></h3>`;
  const editor = state.fitEdit ? modelEditor(c) : "";
  if (p.error) return `<div class="model-box">${head}<p class="muted">No noise model: ${h(p.error)}</p>${editor}</div>`;
  const comps = p.alphas.map((a, i) => [a, p.h[i], p.noise_types[i]]).filter(([, v]) => v > 0);
  const at = p.sigma_at;
  return `<div class="model-box">${head}
    <div class="model-grid">${comps.map(([a, v, n]) => `<span class="badge" title="${noiseTitle(n)}">${noiseLabel(n)}</span><span>h<sub>${a}</sub></span><span>${sci(v, 3)}</span>`).join("")}</div>
    ${editor}
    ${at.tau.length ? `<h3>ADEV from the model</h3>
    <div class="model-grid">${at.tau.map((t, i) => `<span>τ = ${t >= 86400 ? `${num(t / 86400, 3)} d` : `${num(t, 4)} s`}</span><span></span><span>${sci(at.adev[i], 3)}</span>`).join("")}</div>` : ""}
  </div>`;
}

function modelEditor(c) {
  const span = (k) => `${k.from ?? "start"} – ${k.to ?? "end"} s`;
  const used = new Set(c.components.map((k) => k.noise));
  return `<div class="model-edit">
    ${c.components.length ? `<ol class="model-steps">${c.components.map((k) => `<li><span class="badge">${noiseLabel(k.noise)}</span> ${span(k)}</li>`).join("")}</ol>`
      : `<p class="muted">Pick the noise that dominates a stretch of MHDEV and fit its level there. The dashed line is what is left for the next one.</p>`}
    <div class="model-add">
      <select data-ref="fit-noise" aria-label="Noise type">${NOISES.filter((n) => !used.has(n)).map((n) => `<option value="${n}" title="${noiseTitle(n)}">${noiseLabel(n)}</option>`).join("")}</select>
      <input data-ref="fit-from" inputmode="decimal" placeholder="from τ" aria-label="From τ in seconds">
      <input data-ref="fit-to" inputmode="decimal" placeholder="to τ" aria-label="To τ in seconds"><span class="muted">s</span>
      <button type="button" class="btn btn-sm" data-fit="add">Add</button>
    </div>
    ${c.components.length ? `<div class="model-add"><button type="button" class="btn btn-sm" data-fit="undo">Undo</button>
      <button type="button" class="link-btn" data-fit="auto">Back to the automatic fit</button></div>` : ""}
  </div>`;
}

function editModel(act, ctx) {
  const c = visibleClocks().find((k) => k.result);
  if (!c) return;
  if (act === "toggle") return update((s) => { s.fitEdit = !s.fitEdit; });
  const value = (ref) => {
    const v = root.querySelector(`[data-ref=${ref}]`)?.value.trim();
    return v && Number.isFinite(Number(v)) ? Number(v) : null;
  };
  if (act === "add") {
    const noise = root.querySelector("[data-ref=fit-noise]")?.value;
    if (!noise) return;
    update((s) => { c.components = [...c.components, { noise, from: value("fit-from"), to: value("fit-to") }]; markStale(s, "inputs"); });
  }
  if (act === "undo") update((s) => { c.components = c.components.slice(0, -1); markStale(s, "inputs"); });
  if (act === "auto") update((s) => { c.components = []; markStale(s, "inputs"); });
  ctx.run();
}

function renderCompare(box, clocks, devs) {
  const taus = compareTaus(clocks, devs);
  if (!taus.length) { box.hidden = true; return; }
  if (!taus.some((t) => state.compareTau && sameTau(t, state.compareTau))) {
    state.compareTau = taus[Math.floor((taus.length - 1) / 2)];
  }
  const cols = devs.slice(0, 3);
  const valueAt = (c, d) => {
    const r = resultOf(c, d);
    const i = r ? r.tau.findIndex((t) => sameTau(t, state.compareTau)) : -1;
    return i >= 0 ? r.sigma[i] : null;
  };
  const best = Object.fromEntries(cols.map((d) => {
    const vals = clocks.map((c) => valueAt(c, d)).filter((v) => v !== null);
    return [d, vals.length > 1 ? Math.min(...vals) : null];
  }));
  box.className = "card side-card cmp";
  box.style.setProperty("--ncols", cols.length);
  box.innerHTML = `
    <div class="card-head"><h2>Compare at τ</h2>
      <select data-ref="cmp-tau" aria-label="τ">${taus.map((t) => `<option value="${t}" ${sameTau(t, state.compareTau) ? "selected" : ""}>${num(t, 4)} s</option>`).join("")}</select></div>
    <div class="tbl">
      <div class="tbl-row tbl-head"><span>Clock</span>${cols.map((d) => `<span>${label(d)}</span>`).join("")}</div>
      ${clocks.map((c) => `<div class="tbl-row"><span class="name"><span class="clock-dot" style="background:${seriesColor(c.slot)};width:9px;height:9px"></span>${h(c.name)}</span>${cols.map((d) => {
        const v = valueAt(c, d);
        return `<span class="${v !== null && v === best[d] ? "best" : ""}">${v === null ? "–" : sci(v, 3)}</span>`;
      }).join("")}</div>`).join("")}
    </div>
    <div class="card-foot" style="margin-top:0;border-top:0">Bold marks the most stable clock at this τ. Click the plot to move τ.${devs.length > 3 ? ` Showing the first three ${state.axis === "sigma_y" ? "σ<sub>y</sub>" : "σ<sub>x</sub>"} estimators.` : ""}</div>`;
}

export function hint() {
  if (multi()) {
    const devs = devsOnAxis();
    const styles = ["Solid", "dashed", "dotted", "dash-dot", "long dash", "long dash-dot"];
    return `${devs.map((d) => `${styles[state.devs.indexOf(d) % 6]} = ${label(d)}`).join(" · ")} · ◉ hides a clock · click the plot to set τ`;
  }
  return "Hover a point for its value · drag to zoom · double-click to reset";
}

// Stability page: σ(τ) plot and, on the right, the per-τ table (one clock) or
// "Compare at τ" (several clocks, overlaid). Mounted once per visit;
// refresh() redraws from state.
import { state, update, visibleClocks, multi } from "../state.js";
import { byId, label, gridName, CONFIDENCES, noiseLabel, noiseTitle, noiseSlope, BIAS_CORRECTED } from "../estimators.js";
import { loadPlotly, baseLayout, logAxis, CONFIG, seriesColor, DASHES, hexAlpha, css } from "../plot.js";
import { sci, num, escapeHtml as h } from "../format.js";

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
        <span class="toggles">
          <label class="check" title="Noise type identified at each τ (lag-1 autocorrelation), with the power-law slope it implies"><input type="checkbox" data-tog="noiseId">Noise ID</label>
          <label class="check" title="Draw the confidence intervals as a shaded band instead of error bars"><input type="checkbox" data-tog="band">CI band</label>
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
  });
  root.addEventListener("change", (e) => {
    const tog = e.target.dataset.tog;
    if (tog) { update((s) => { s[tog] = e.target.checked; }); return; }
    const r = e.target.dataset.ref;
    if (r === "table-dev") update((s) => { s.tableDev = e.target.value; });
    if (r === "cmp-tau") update((s) => { s.compareTau = Number(e.target.value); });
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
  const bias = biased ? (state.bias ? " · bias-corrected" : " · not bias-corrected") : "";
  root.querySelector("[data-ref=caption]").textContent = `${conf}${gridName(state.grid).toLowerCase()} τ grid${bias}`;
  // Both need the intervals: noise identification runs with them.
  root.querySelectorAll("[data-tog]").forEach((el) => {
    const lab = el.closest("label");
    lab.dataset.tip ??= lab.title;
    el.checked = !!state[el.dataset.tog];
    el.disabled = !state.ci;
    lab.title = state.ci ? lab.dataset.tip : "Turn on confidence intervals in the panel; noise identification runs with them";
  });

  const clocks = visibleClocks();
  const devs = devsOnAxis();
  const hasRecord = state.clocks.some((c) => c.record);
  const hasResult = clocks.some((c) => c.result);
  const other = state.axis === "sigma_y" ? "sigma_x" : "sigma_y";

  if (!hasRecord) showEmpty(`<h2>No record yet</h2><p>Load a phase or frequency file, or synthesize one, on the Data page.</p><button class="btn" data-act="load">Go to Data</button>`);
  else if (!hasResult) showEmpty(`<h2>Not analysed yet</h2><p>Pick estimators in the panel and run the analysis.</p><button class="btn btn-primary" data-act="run">Run analysis</button>`);
  else if (!devs.length) showEmpty(`<h2>No ${state.axis === "sigma_y" ? "σ<sub>y</sub>" : "σ<sub>x</sub>"} estimator selected</h2><p>${state.axis === "sigma_x" ? "Add TDEV, MTIE or another time estimator" : "Add ADEV or another frequency estimator"}, or switch to ${AXES[other]}.</p>`);
  else if (!clocks.length) showEmpty(`<h2>All clocks are hidden</h2><p>Show a clock with ◉ in the panel.</p>`);
  else showEmpty("");

  renderTable(clocks, devs);
  if (Plotly && hasResult && devs.length && clocks.length) draw(clocks, devs);
}

function traces(clocks, devs) {
  const several = multi();
  const out = [];
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
        out.push({ x: r.tau, y: r.ci_lo, mode: "lines", line: { width: 0 }, fill: "tonexty", fillcolor: hexAlpha(color, 0.14), hoverinfo: "skip", showlegend: false });
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
    });
  }
  return out;
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

const NOISE_VAR = { WHPM: "--noise-wpm", FLPM: "--noise-fpm", WHFM: "--noise-wfm", FLFM: "--noise-ffm", RWFM: "--noise-rwfm" };
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
        // Through the run's centre of mass in log-log, across the run.
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
    <div class="card-foot">Hover a row for its confidence interval. Every value is in the TSV and JSON exports.</div>`;
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

// Plot element and table rows for the export menu.
export function plotElement() { return plotEl; }
export function hint() {
  if (multi()) {
    const devs = devsOnAxis();
    const styles = ["Solid", "dashed", "dotted", "dash-dot", "long dash", "long dash-dot"];
    return `${devs.map((d) => `${styles[state.devs.indexOf(d) % 6]} = ${label(d)}`).join(" · ")} · ◉ hides a clock · click the plot to set τ`;
  }
  return "Hover a point for its value · drag to zoom · double-click to reset";
}
